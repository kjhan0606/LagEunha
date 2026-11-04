static int gfs_kdk_adiabatic_volume_on(void)
{
    const char *s=getenv("GFS_KDK_ADIABATIC_VOLUME");
    return s && strcmp(s,"1")==0;
}

/* Exact fixed-mass ideal-gas PdV update: U V^(Gamma-1) is invariant. */
static int gfs_kdk_adiabatic_volume_update(double U0,double V0,double V1,
                                           double Gamma,double *U1)
{
    if(!U1 || !(U0>0.0) || !(V0>0.0) || !(V1>0.0) || !(Gamma>1.0) ||
       !isfinite(U0) || !isfinite(V0) || !isfinite(V1) || !isfinite(Gamma))
        return 0;
    const double ratio=V0/V1;
    const double value=U0*pow(ratio,Gamma-1.0);
    if(!(value>0.0) || !isfinite(value) || !isfinite((float)value) ||
       !((float)value>0.0f))
        return 0;
    *U1=value;
    return 1;
}

typedef struct {
    double x, y, vx, vy, U, V, mass;
    double ax, ay, Udot;
    unsigned period;
    long start_tick;
    int valid;
} GfsKdkLtsBlock2State;

/* Symmetric KDK subcycling of the external Kepler potential only.
 * The hydro RHS remains held at its block value during each drift. */
static int gfs_kdk_lts_kepler_drift(double x, double y, double vx, double vy,
        double elapsed, double mass, int wall_events, double eta, double cx, double cy, double eps,
        double accx_ext, double accy_ext,
        double *xout, double *yout, double *vxout, double *vyout,
        GfsKdkFixedPolyWallStats *wall_stats)
{
    double t = 0.0;
    int n = 0;
    if (wall_stats) memset(wall_stats, 0, sizeof(*wall_stats));
    if (!(elapsed >= 0.0) || !isfinite(elapsed) ||
        !(eta > 0.0) || !isfinite(eta) ||
        !(mass > 0.0) || !isfinite(mass) ||
        !(eps >= 0.0) || !isfinite(eps) ||
        !isfinite(x) || !isfinite(y) || !isfinite(vx) || !isfinite(vy) ||
        !xout || !yout || !vxout || !vyout || !wall_stats)
        return -1;
    while (t < elapsed) {
        if (n >= 4096) return -10;
        double dx = x - cx, dy = y - cy;
        double r2 = dx*dx + dy*dy + eps*eps;
        double hmax = eta * pow(r2, 0.75);
        double remaining = elapsed - t;
        double h = fmin(hmax, remaining);
        if (!(h > 0.0) || !isfinite(h) || !(t + h > t)) return -11;

        postype ax0 = (postype)accx_ext, ay0 = (postype)accy_ext;
        kepler_accel((postype)x, (postype)y, &ax0, &ay0);
        double vhalf_x = vx + 0.5*h*(double)ax0;
        double vhalf_y = vy + 0.5*h*(double)ay0;
        double xn, yn;
        if (wall_events) {
            GfsKdkFixedPolyWallStats one = {0};
            xn = x;
            yn = y;
            if (!gfs_kdk_fixedpoly_wall_drift(&xn,&yn,&vhalf_x,&vhalf_y,
                    h,mass,&one)) {
                wall_stats->failure_code=one.failure_code;
                wall_stats->failure_side=one.failure_side;
                wall_stats->failure_remaining=one.failure_remaining;
                wall_stats->failure_gap=one.failure_gap;
                return -(100+one.failure_code);
            }
            wall_stats->impacts += one.impacts;
            wall_stats->wall_jx += one.wall_jx;
            wall_stats->wall_jy += one.wall_jy;
            wall_stats->wall_torque += one.wall_torque;
            wall_stats->kinetic_roundoff += one.kinetic_roundoff;
        } else {
            xn = x + h*vhalf_x;
            yn = y + h*vhalf_y;
        }
        postype ax1 = (postype)accx_ext, ay1 = (postype)accy_ext;
        kepler_accel((postype)xn, (postype)yn, &ax1, &ay1);
        vx = vhalf_x + 0.5*h*(double)ax1;
        vy = vhalf_y + 0.5*h*(double)ay1;
        x = xn;
        y = yn;
        t += h;
        ++n;
    }
    if (!isfinite(x) || !isfinite(y) || !isfinite(vx) || !isfinite(vy))
        return -12;
    *xout=x; *yout=y; *vxout=vx; *vyout=vy;
    return n;
}
/*
 * Two-level correctness prototype. All cells are predicted to each fine
 * event and the existing CUDA Voronoi operator is evaluated globally. Each
 * cell receives a KDK velocity correction only at its block endpoint.
 * This intentionally does not claim active-only geometry or a speedup.
 */
static int gfs_kdk_lts_block2(
        SimParameters *simpar,
        treevorostressrk4particletype *sbp,
        postype dt_base,
        postype dt_ceiling,
        postype xmin, postype ymin, postype xmax, postype ymax,
        postype OrderOfAccuracy, postype Courant, postype Gamma,
        void (*paddingAllTreeParticles)(SimParameters *, postype),
        Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
        treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
        void (*mkLinkedList2D)(SimParameters *, postype, postype, postype,
                               postype, postype,
                               void (*)(SimParameters *, postype)),
        void (*postStage)(SimParameters *),
        int (*targetBP)(treevorork4particletype *, postype, postype),
        postype Lx, postype Ly,
        postype *dt_macro_out,
        double *update_seconds, double *force_seconds, double *post_seconds)
{
    const int np0 = (int)VORO_NP(simpar);
    size_t max_id = 0;
    int n_target = 0;
    int n_period1 = 0, n_period2 = 0;
    int i, tick;
    double base = (double)dt_base;
    double macro, fine;
    double eta = 0.0, cx = 2.0, cy = 2.0, eps = 0.0;
    double subcycle_eta = 0.0;
    int subcycle_active = 0;
    const int wall_mode = gfs_kdk_fixedpoly_wall_events();
    const int adiabatic_volume = gfs_kdk_adiabatic_volume_on();
    const char *ic = getenv("EUNHA_IC");
    const char *s;
    GfsKdkLtsBlock2State *state = NULL;
    unsigned char *seen = NULL;

    if (!dt_macro_out || !(base > 0.0) || !isfinite(base) ||
        np0 <= 0 || !postStage || !targetBP)
        return -1;

    /* A caller's exact-stop/diagnostic ceiling already equal to the global
     * base step leaves no macro interval to subdivide; retain ordinary KDK. */
    if (dt_ceiling > 0.0 &&
        (double)dt_ceiling <= base * (1.0 + 1.0e-12))
        return 0;
    if (dt_ceiling > 0.0 &&
        (double)dt_ceiling < 2.0 * base * (1.0 - 1.0e-12))
        return 0;

    if (ic && strcmp(ic, "kepler") == 0) {
        s = getenv("EUNHA_KEPLER_DTETA");
        if (s && s[0]) eta = atof(s);
        s = getenv("EUNHA_KEPLER_CX");
        if (s && s[0]) cx = atof(s);
        s = getenv("EUNHA_KEPLER_CY");
        if (s && s[0]) cy = atof(s);
        s = getenv("EUNHA_KEPLER_EPS");
        if (s && s[0]) eps = atof(s);
        if (!(eta >= 0.0) || !isfinite(eta) || !isfinite(cx) ||
            !isfinite(cy) || !(eps >= 0.0) || !isfinite(eps))
            return -1;
    }
    s = getenv("GFS_KDK_KEPLER_SUBCYCLE_ETA");
    if (s && s[0]) {
        subcycle_eta = atof(s);
        if (!(subcycle_eta > 0.0) || !isfinite(subcycle_eta) ||
            !ic || strcmp(ic, "kepler") != 0)
            return -1;
        subcycle_active = 1;
        postype warm_x=0, warm_y=0;
        kepler_accel((postype)0, (postype)0, &warm_x, &warm_y);
    }
    if (wall_mode && !subcycle_active) {
        fprintf(stderr,"[KDK_LTS_BLOCK2] KDK wall events require central-force KDK subcycling\n");
        return -1;
    }

    for (i = 0; i < np0; ++i) {
        treevorork4particletype *p =
            (treevorork4particletype *)(sbp + i);
        if (!targetBP(p, Lx, Ly)) continue;
        size_t id = (size_t)PINDX(p);
        if (id > max_id) max_id = id;
        ++n_target;
    }
    if (n_target <= 0 || max_id > (size_t)np0 * 16u + 1024u ||
        max_id == (size_t)-1)
        return -1;

    size_t nslot = max_id + 1u;
    state = (GfsKdkLtsBlock2State *)calloc(nslot, sizeof(*state));
    seen = (unsigned char *)calloc(nslot, sizeof(*seen));
    if (!state || !seen) {
        free(state);
        free(seen);
        return -1;
    }

    for (i = 0; i < np0; ++i) {
        treevorork4particletype *p =
            (treevorork4particletype *)(sbp + i);
        if (!targetBP(p, Lx, Ly)) continue;
        size_t id = (size_t)PINDX(p);
        if (id >= nslot || seen[id]) {
            fprintf(stderr, "[KDK_LTS_BLOCK2] duplicate/out-of-range particle id=%zu\n", id);
            MPI_Abort(MPI_COMM(simpar), 96);
        }
        seen[id] = 1;

        double candidate = (double)sbp[i].dt;
        if (eta > 0.0) {
            double dx = (double)sbp[i].x - cx;
            double dy = (double)sbp[i].y - cy;
            double r2 = dx*dx + dy*dy + eps*eps;
            double dtg = eta * pow(r2, 0.75);
            if (!(dtg > 0.0) || !isfinite(dtg)) {
                fprintf(stderr, "[KDK_LTS_BLOCK2] invalid local Kepler dt id=%zu\n", id);
                MPI_Abort(MPI_COMM(simpar), 96);
            }
            if (dtg < candidate) candidate = dtg;
        }
        if (!isfinite(candidate) || !(candidate > 0.0) ||
            candidate < base * (1.0 - 1.0e-6)) {
            fprintf(stderr,
                    "[KDK_LTS_BLOCK2] local dt violates synchronized base id=%zu dt=%.17g base=%.17g\n",
                    id, candidate, base);
            MPI_Abort(MPI_COMM(simpar), 96);
        }

        GfsKdkLtsBlock2State *q = state + id;
        q->x = (double)sbp[i].x;
        q->y = (double)sbp[i].y;
        q->vx = (double)sbp[i].vx;
        q->vy = (double)sbp[i].vy;
        q->U = sbp[i].weight_thermal.U;
        q->V = sbp[i].weight_thermal.V;
        q->mass = (double)sbp[i].mass;
        q->Udot = sbp[i].weight_thermal.Udot;
        postype gx = 0, gy = 0;
        if (!subcycle_active)
            kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
        q->ax = (double)sbp[i].ax +
            (subcycle_active ? 0.0 : (double)GAS_ACCX(simpar) + (double)gx);
        q->ay = (double)sbp[i].ay +
            (subcycle_active ? 0.0 : (double)GAS_ACCY(simpar) + (double)gy);
        q->period = candidate >= 2.0 * base ? 2u : 1u;
        q->start_tick = 0;
        q->valid = 1;
        if (q->period == 1u) ++n_period1;
        else ++n_period2;
        if (!(q->mass > 0.0) || !(q->U > 0.0) || !(q->V > 0.0) ||
            !isfinite(q->mass) || !isfinite(q->U) || !isfinite(q->V) ||
            !isfinite(q->Udot) || !isfinite(q->ax) || !isfinite(q->ay)) {
            fprintf(stderr, "[KDK_LTS_BLOCK2] invalid initial state id=%zu\n", id);
            MPI_Abort(MPI_COMM(simpar), 96);
        }
    }

    /* No asynchronous work exists when all cells select the base step. */
    if (n_period2 == 0) {
        free(state);
        free(seen);
        return 0;
    }

    macro = 2.0 * base;
    if (dt_ceiling > 0.0 && macro > (double)dt_ceiling)
        macro = (double)dt_ceiling;
    fine = 0.5 * macro;
    if (!(macro > 0.0) || !(fine > 0.0) ||
        !isfinite(macro) || !isfinite(fine)) {
        free(state);
        free(seen);
        return -1;
    }

    const size_t entry_bytes=(size_t)np0*sizeof(*sbp);
    treevorostressrk4particletype *entry_backup=
        (treevorostressrk4particletype*)malloc(entry_bytes);
    if(!entry_backup){free(state);free(seen);return -1;}
    memcpy(entry_backup,sbp,entry_bytes);
    const double saved_floor=gfs_stage_floor_l;
    const int saved_floor_n=gfs_stage_floor_n_l;
    const long long saved_pairs=gfs_pair_face_ends;
    const double saved_eaud[5]={gfs_eaud_hyd,gfs_eaud_hyd_abs,gfs_eaud_wg,
                                gfs_eaud_eh0,gfs_eaud_ep0};
    unsigned long long macro_wall_hits=0;
    double macro_wall_jx=0.0,macro_wall_jy=0.0;
    double macro_wall_torque=0.0,macro_wall_dk=0.0;

    if (MYID(simpar) == 0) {
        double active_fraction =
            ((double)n_period1 * 0.5 + (double)n_period2) / (double)n_target;
        fprintf(stderr,
                "[KDK_LTS_BLOCK2] cells=%d bin0=%d bin1=%d dt_base=%.9e dt_macro=%.9e fine_events=2 active_cell_fraction=%.6e full_geometry=all_cells\n",
                n_target, n_period1, n_period2, fine, macro, active_fraction);
        fflush(stderr);
    }

    for (tick = 1; tick <= 2; ++tick) {
        int bad = 0, nsub_max = 0, prediction_logged = 0;
        int predictor_recoverable = 0, predictor_fatal = 0;
        unsigned long long wall_hits = 0;
        double wall_jx = 0.0, wall_jy = 0.0;
        double wall_torque = 0.0, wall_dk = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:bad,nsub_max,predictor_recoverable,predictor_fatal) reduction(+:wall_hits,wall_jx,wall_jy,wall_torque,wall_dk)
#endif
        for (i = 0; i < np0; ++i) {
            treevorork4particletype *p =
                (treevorork4particletype *)(sbp + i);
            if (!targetBP(p, Lx, Ly)) {
                sbp[i].vx = sbp[i].vy = 0;
                continue;
            }
            size_t id = (size_t)PINDX(p);
            if (id >= nslot || !state[id].valid ||
                tick < state[id].start_tick) {
                bad = 1;
                predictor_fatal = 1;
#ifdef _OPENMP
#pragma omp critical(gfs_kdk_lts_bad_cell)
#endif
                if (!prediction_logged) {
                    fprintf(stderr,"[KDK_LTS_BAD_PREDICTION_STATE] tick=%d id=%zu reason=id_or_start_tick nslot=%zu valid=%d start_tick=%ld\n",
                        tick,id,nslot,(id<nslot)?state[id].valid:0,
                        (id<nslot)?state[id].start_tick:-1L);
                    prediction_logged=1;
                }
                continue;
            }
            const GfsKdkLtsBlock2State *q = state + id;
            const double elapsed =
                (double)(tick - q->start_tick) * fine;
            double x, y, vx, vy;
            const double U = adiabatic_volume ? q->U : q->U + q->Udot * elapsed;
            GfsKdkFixedPolyWallStats wall_piece = {0};
            if (subcycle_active) {
                const double hblock = (double)q->period * fine;
                const double vx_start = q->vx + 0.5*hblock*q->ax;
                const double vy_start = q->vy + 0.5*hblock*q->ay;
                int nsub = gfs_kdk_lts_kepler_drift(q->x, q->y,
                    vx_start, vy_start, elapsed, q->mass, wall_mode,
                    subcycle_eta, cx, cy, eps,
                    (double)GAS_ACCX(simpar), (double)GAS_ACCY(simpar),
                    &x, &y, &vx, &vy, &wall_piece);
                if (nsub < 0) {
                    bad = 1;
                    predictor_fatal = 1;
#ifdef _OPENMP
#pragma omp critical(gfs_kdk_lts_bad_cell)
#endif
                    if (!prediction_logged) {
                        fprintf(stderr,"[KDK_LTS_BAD_PREDICTION] tick=%d id=%zu drift_status=%d elapsed=%.17g period=%u x=%.17g y=%.17g vx=%.17g vy=%.17g wall_failure=%d side=%d remaining=%.17g gap=%.17g\n",
                            tick,id,nsub,elapsed,q->period,q->x,q->y,q->vx,q->vy,
                            wall_piece.failure_code,wall_piece.failure_side,
                            wall_piece.failure_remaining,wall_piece.failure_gap);
                        prediction_logged=1;
                    }
                    continue;
                }
                if (nsub > nsub_max) nsub_max = nsub;
            } else {
                const double half = 0.5 * elapsed * elapsed;
                x = q->x + q->vx * elapsed + q->ax * half;
                y = q->y + q->vy * elapsed + q->ay * half;
                vx = q->vx + q->ax * elapsed;
                vy = q->vy + q->ay * elapsed;
            }
            int prediction_flags =
                (!isfinite(x) ? 1 : 0) | (!isfinite(y) ? 2 : 0) |
                (!isfinite(vx) ? 4 : 0) | (!isfinite(vy) ? 8 : 0) |
                (!(U > 0.0) ? 16 : 0) | (!isfinite(U) ? 32 : 0) |
                ((float)U <= 0.0f || !isfinite((float)U) ? 64 : 0) |
                ((double)sbp[i].mass != q->mass ? 128 : 0);
            if (prediction_flags) {
                bad = 1;
                if (prediction_flags == 16) predictor_recoverable = 1;
                else predictor_fatal = 1;
#ifdef _OPENMP
#pragma omp critical(gfs_kdk_lts_bad_cell)
#endif
                if (!prediction_logged) {
                    fprintf(stderr,"[KDK_LTS_BAD_PREDICTION_STATE] tick=%d id=%zu flags=%d elapsed=%.17g U0=%.17g Udot=%.17g Upred=%.17g mass_now=%.17g mass_start=%.17g\n",
                        tick,id,prediction_flags,elapsed,q->U,q->Udot,U,
                        (double)sbp[i].mass,q->mass);
                    prediction_logged=1;
                }
                continue;
            }
            sbp[i].x = (postype)x;
            sbp[i].y = (postype)y;
            sbp[i].vx = (postype)vx;
            sbp[i].vy = (postype)vy;
            sbp[i].weight_thermal.U = U;
            sbp[i].ie = (float)U;
            if (wall_mode && ((unsigned)tick % q->period == 0u)) {
                wall_hits += wall_piece.impacts;
                wall_jx += wall_piece.wall_jx;
                wall_jy += wall_piece.wall_jy;
                wall_torque += wall_piece.wall_torque;
                wall_dk += wall_piece.kinetic_roundoff;
            }
        }
        if (bad) {
            if (predictor_recoverable && !predictor_fatal) {
                if ((int)VORO_NP(simpar) != np0) {
                    fprintf(stderr,"[KDK_LTS_FALLBACK] cannot restore after particle-count change\n");
                    MPI_Abort(MPI_COMM(simpar),96);
                }
                memcpy(VORORK4_TBP(simpar),entry_backup,entry_bytes);
                gfs_stage_floor_l=saved_floor;
                gfs_stage_floor_n_l=saved_floor_n;
                gfs_pair_face_ends=saved_pairs;
                gfs_eaud_hyd=saved_eaud[0];gfs_eaud_hyd_abs=saved_eaud[1];
                gfs_eaud_wg=saved_eaud[2];gfs_eaud_eh0=saved_eaud[3];
                gfs_eaud_ep0=saved_eaud[4];
                fprintf(stderr,"[KDK_LTS_FALLBACK] reason=nonpositive_predictor_U tick=%d dt_macro=%.17g retry=global_KDK dt=%.17g state_and_ledgers_restored=1\n",
                    tick,macro,base);
                free(entry_backup);free(state);free(seen);
                return 0;
            }
            fprintf(stderr, "[KDK_LTS_BLOCK2] invalid prediction or mass changed at tick=%d\n", tick);
            MPI_Abort(MPI_COMM(simpar), 96);
        }
        if (subcycle_active && MYID(simpar) == 0) {
            fprintf(stderr,
                "[KDK_LTS_KEPLER_SUBCYCLE] tick=%d eta=%.6g max_substeps_per_cell=%d\n",
                tick, subcycle_eta, nsub_max);
        }

        double t0 = MPI_Wtime();
        postStage(simpar);
        if (post_seconds) *post_seconds += MPI_Wtime() - t0;
        sbp = (treevorostressrk4particletype *)VORORK4_TBP(simpar);
        if ((int)VORO_NP(simpar) != np0) {
            fprintf(stderr, "[KDK_LTS_BLOCK2] particle count changed in one-rank event\n");
            MPI_Abort(MPI_COMM(simpar), 96);
        }

        for (i = 0; i < np0; ++i) {
            sbp[i].w2 = sbp[i].rk4.w2backup;
            if (GAS_AVMODE(simpar) >= 1)
                sbp[i].stress.vsig_max = 0;
        }

        t0 = MPI_Wtime();
        updateDenW2Pressure2DBlend(simpar, xmin, ymin, xmax, ymax,
                Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
                mkLinkedList2D, fine);
        if (update_seconds) *update_seconds += MPI_Wtime() - t0;
        if (adiabatic_volume) {
            sbp=(treevorostressrk4particletype*)VORORK4_TBP(simpar);
            int adiabat_bad=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:adiabat_bad)
#endif
            for(i=0;i<np0;i++){
                treevorork4particletype *p=(treevorork4particletype*)(sbp+i);
                if(!targetBP(p,Lx,Ly)) continue;
                size_t id=(size_t)PINDX(p);
                if(id>=nslot || !state[id].valid){adiabat_bad=1;continue;}
                double Unew=0.0;
                if(!gfs_kdk_adiabatic_volume_update(state[id].U,state[id].V,
                        sbp[i].weight_thermal.V,(double)Gamma,&Unew)){
                    adiabat_bad=1;continue;
                }
                sbp[i].weight_thermal.U=Unew;
                sbp[i].ie=(float)Unew;
            }
            if(adiabat_bad){
                fprintf(stderr,"[KDK_LTS_ADIABATIC_VOLUME] invalid volume or thermal state at tick=%d\n",tick);
                MPI_Abort(MPI_COMM(simpar),96);
            }
        }

        t0 = MPI_Wtime();
#ifdef USE_CUDA
        if (GAS_GPU_ENABLED(simpar))
            (void)getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
                    OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
                    find2DNeighborBP, find2DCellBP, mkLinkedList2D);
        else
#endif
        {
            fprintf(stderr, "[KDK_LTS_BLOCK2] CUDA force path is required\n");
            MPI_Abort(MPI_COMM(simpar), 96);
        }
        if (force_seconds) *force_seconds += MPI_Wtime() - t0;
        sbp = (treevorostressrk4particletype *)VORORK4_TBP(simpar);

        bad = 0;
        int recoverable_bad = 0, fatal_bad = 0;
        int bad_logged = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:bad,recoverable_bad,fatal_bad)
#endif
        for (i = 0; i < np0; ++i) {
            treevorork4particletype *p =
                (treevorork4particletype *)(sbp + i);
            if (!targetBP(p, Lx, Ly)) continue;
            size_t id = (size_t)PINDX(p);
            if (id >= nslot || !state[id].valid ||
                (double)sbp[i].mass != state[id].mass) {
                bad = 1;
                fatal_bad = 1;
#ifdef _OPENMP
#pragma omp critical(gfs_kdk_lts_bad_cell)
#endif
                if (!bad_logged) {
                    fprintf(stderr,"[KDK_LTS_BAD_ENDPOINT] tick=%d id=%zu reason=id_state_or_mass nslot=%zu valid=%d mass_now=%.17g mass_start=%.17g\n",
                        tick,id,nslot,(id<nslot)?state[id].valid:0,
                        (double)sbp[i].mass,(id<nslot)?state[id].mass:NAN);
                    bad_logged = 1;
                }
                continue;
            }
            GfsKdkLtsBlock2State *q = state + id;
            if ((unsigned)tick % q->period != 0u) continue;

            const double h = (double)q->period * fine;
            double ax_end, ay_end, vx_end, vy_end;
            if (subcycle_active) {
                ax_end = (double)sbp[i].ax;
                ay_end = (double)sbp[i].ay;
                vx_end = (double)sbp[i].vx + 0.5*h*ax_end;
                vy_end = (double)sbp[i].vy + 0.5*h*ay_end;
            } else {
                postype gx = 0, gy = 0;
                kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
                ax_end = (double)sbp[i].ax +
                    (double)GAS_ACCX(simpar) + (double)gx;
                ay_end = (double)sbp[i].ay +
                    (double)GAS_ACCY(simpar) + (double)gy;
                vx_end = q->vx + 0.5*h*(q->ax + ax_end);
                vy_end = q->vy + 0.5*h*(q->ay + ay_end);
            }
            const double Udot_end = sbp[i].weight_thermal.Udot;
            const double U_end = adiabatic_volume ? sbp[i].weight_thermal.U :
                q->U + 0.5 * h * (q->Udot + Udot_end);
            const double V_end = sbp[i].weight_thermal.V;
            if (!(U_end > 0.0) || !isfinite(U_end) || !(V_end>0.0) ||
                !isfinite(V_end) || !isfinite(vx_end) || !isfinite(vy_end) ||
                !isfinite(ax_end) || !isfinite(ay_end) ||
                !isfinite(Udot_end) || !isfinite((float)U_end)) {
                bad = 1;
                int flags = (!(U_end > 0.0) ? 1 : 0) |
                    (!isfinite(U_end) ? 2 : 0) |
                    (!isfinite(vx_end) || !isfinite(vy_end) ? 4 : 0) |
                    (!isfinite(ax_end) || !isfinite(ay_end) ? 8 : 0) |
                    (!isfinite(Udot_end) ? 16 : 0) |
                    (!isfinite((float)U_end) ? 32 : 0);
                if(flags==1) recoverable_bad=1;
                else fatal_bad=1;
#ifdef _OPENMP
#pragma omp critical(gfs_kdk_lts_bad_cell)
#endif
                if (!bad_logged) {
                    fprintf(stderr,"[KDK_LTS_BAD_ENDPOINT] tick=%d id=%zu flags=%d period=%u h=%.17g U0=%.17g Udot0=%.17g Udot1=%.17g U1=%.17g v=(%.17g,%.17g) a=(%.17g,%.17g)\n",
                        tick,id,flags,q->period,h,q->U,q->Udot,Udot_end,U_end,
                        vx_end,vy_end,ax_end,ay_end);
                    bad_logged = 1;
                }
                continue;
            }

            sbp[i].vx = (postype)vx_end;
            sbp[i].vy = (postype)vy_end;
            sbp[i].weight_thermal.U = U_end;
            sbp[i].ie = (float)U_end;
            if (adiabatic_volume) {
                q->U=U_end;
                q->V=V_end;
                q->Udot=Udot_end;
            }
            if (tick < 2) {
                q->x = (double)sbp[i].x;
                q->y = (double)sbp[i].y;
                q->vx = vx_end;
                q->vy = vy_end;
                q->U = U_end;
                q->ax = subcycle_active ? (double)sbp[i].ax : ax_end;
                q->ay = subcycle_active ? (double)sbp[i].ay : ay_end;
                q->Udot = Udot_end;
                q->start_tick = tick;
            }
        }
        if (bad) {
            if(recoverable_bad && !fatal_bad){
                if((int)VORO_NP(simpar)!=np0){
                    fprintf(stderr,"[KDK_LTS_FALLBACK] cannot restore after particle-count change\n");
                    MPI_Abort(MPI_COMM(simpar),96);
                }
                memcpy(VORORK4_TBP(simpar),entry_backup,entry_bytes);
                gfs_stage_floor_l=saved_floor;
                gfs_stage_floor_n_l=saved_floor_n;
                gfs_pair_face_ends=saved_pairs;
                gfs_eaud_hyd=saved_eaud[0];gfs_eaud_hyd_abs=saved_eaud[1];
                gfs_eaud_wg=saved_eaud[2];gfs_eaud_eh0=saved_eaud[3];
                gfs_eaud_ep0=saved_eaud[4];
                fprintf(stderr,"[KDK_LTS_FALLBACK] reason=nonpositive_corrector_U tick=%d dt_macro=%.17g retry=global_KDK dt=%.17g state_and_ledgers_restored=1\n",
                    tick,macro,base);
                free(entry_backup);free(state);free(seen);
                return 0;
            }
            fprintf(stderr, "[KDK_LTS_BLOCK2] invalid endpoint update or mass changed at tick=%d\n", tick);
            MPI_Abort(MPI_COMM(simpar), 96);
        }
        macro_wall_hits+=wall_hits;
        macro_wall_jx+=wall_jx;macro_wall_jy+=wall_jy;
        macro_wall_torque+=wall_torque;macro_wall_dk+=wall_dk;
    }

    if (macro_wall_hits && MYID(simpar) == 0) {
        fprintf(stderr,
            "[KDK_LTS_WALL_LEDGER] step=%d t_end=%.17g hits=%llu wall_J=(%.17g,%.17g) wall_torque=%.17g kinetic_roundoff=%.17g\n",
            gfs_rk4_context_step,gfs_rk4_context_time+macro,macro_wall_hits,
            macro_wall_jx,macro_wall_jy,macro_wall_torque,macro_wall_dk);
        fflush(stderr);
    }
    *dt_macro_out = (postype)macro;
    free(entry_backup);free(state);free(seen);
    return 1;
}
