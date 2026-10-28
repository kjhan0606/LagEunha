"""Attach an opt-in experimental fixed-mass weighted RK4 mode.

The input must be a fresh, byte-verified accepted bounded near32 source
worktree. This patch adds the validated GPU Jacobi rate and the global
pressure-work-neutral projection at each RK4 stage. It is deliberately
restricted to one GPU rank and no padding for its first native gate.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil


DOCS = Path(__file__).resolve().parent


def replace_once(source, old, new, label):
    count = source.count(old)
    if count != 1:
        raise ValueError(f"{label}: expected one anchor, got {count}: {old[:80]!r}")
    return source.replace(old, new)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("worktree", type=Path)
    args = parser.parse_args()
    root = args.worktree.resolve()
    if not (root/".kepler_bounded_overlay.json").exists():
        raise ValueError("Accepted bounded source overlay is missing")
    marker = root/".kepler_energy_neutral_mode.json"
    if marker.exists():
        raise FileExistsError("Weighted source patch already applied")

    paths = {
        name: root/"Exam"/name for name in
        ("exam.c", "exam_gpu.cu", "exam_gpu_extract.c", "exam_gpu.h", "Makefile")
    }
    original = {name: path.read_text() for name,path in paths.items()}
    updated = dict(original)

    c = updated["exam.c"]
    old = "static int gfs_w2lin_active = 0;\nint gfs_w2lin_is_active(void){ return gfs_w2lin_active; }\n"
    new = old + (
        "/* Experimental GPU volume-rate weights; no change when disabled. */\n"
        "static int gfs_volw2_active = 0;\n"
        "int gfs_volw2_is_active(void){ return gfs_volw2_active; }\n"
        "static int gfs_volw2_on(void){\n"
        "\tstatic int on = -1;\n"
        "\tif(on < 0){ const char *s=getenv(\"GFS_GPU_VOLW2\"); on=(s && atoi(s)>0) ? 1 : 0; }\n"
        "\treturn on;\n}\n")
    c = replace_once(c,old,new,"exam.c mode flag")
    old = "static void dyn_face_push(dyn_face_t *buf, size_t n){\n"
    helper = r'''/* Rates are stored by stable particle ID, including across RK4 stage motion. */
int gfs_volw2_stage_current(void){ return dyn_stage_cur; }
void gfs_volw2_store_stage(SimParameters *simpar, int np, const double *rate){
    treevorostressrk4particletype *bp=(treevorostressrk4particletype*)VORORK4_TBP(simpar);
    size_t maxid=0;
    if(!gfs_volw2_active || dyn_stage_cur<0 || dyn_stage_cur>3 || np!=VORO_NP(simpar))
        MPI_Abort(MPI_COMM(simpar),97);
    for(int i=0;i<np;i++){
        size_t id=(size_t)PINDX((treevorork4particletype*)(bp+i));
        if(id==MAX_INDEX || !isfinite(rate[i])) MPI_Abort(MPI_COMM(simpar),97);
        if(id>maxid) maxid=id;
    }
    dyn_ensure(maxid);
    for(int i=0;i<np;i++){
        size_t id=(size_t)PINDX((treevorork4particletype*)(bp+i));
        dyn_k[dyn_stage_cur][id]=rate[i];
        bp[i].w2hydro=(postype)rate[i];
    }
}
'''
    c = replace_once(c,old,helper+old,"exam.c stage store")
    c = replace_once(c,
        "if(GAS_AVMODE(simpar) == 5 && GAS_Kappa(simpar) == 0){\n",
        "if(GAS_AVMODE(simpar) == 5 && GAS_Kappa(simpar) == 0 && !gfs_volw2_on()){\n",
        "exam.c restart weight guard")
    old = "\tint w2lin = gfs_w2lin_on() && GAS_Kappa(simpar) > 0 && !sedov_lagvol_on();\n"
    mode = r'''	int volw2 = gfs_volw2_on();
#ifdef USE_CUDA
	if(volw2 && !(GAS_GPU_ENABLED(simpar) && NID(simpar)==1 &&
		av_mode==5 && GAS_Kappa(simpar)==0 && !sedov_lagvol_on() &&
		!gfs_w2lin_on() && !gfs_dynw2_on() &&
		GAS_FCENTROID(simpar)==0 && GAS_XSPHEPS(simpar)==0 &&
		GAS_VISCOSITY(simpar)==0)){
		fprintf(stderr,"[GPU_VOLW2] requires one CUDA rank, av_mode=5, kappa=0, no steering or physical viscosity\n");
		MPI_Abort(MPI_COMM(simpar),97);
	}
#else
	if(volw2){ fprintf(stderr,"[GPU_VOLW2] CUDA required\n"); MPI_Abort(MPI_COMM(simpar),97); }
#endif
	{
		static int banner=0;
		if(volw2 && !banner++ && MYID(simpar)==0)
			fprintf(stderr,"[GPU_VOLW2] experimental pressure-work-neutral Jacobi16 RK4 weights enabled\n");
	}
'''
    c = replace_once(c,old,mode+old,"exam.c mode gate")
    c = replace_once(c,"\tgfs_w2lin_active = w2lin;\n",
        "\tgfs_w2lin_active = w2lin || volw2;\n\tgfs_volw2_active = volw2;\n",
        "exam.c stage face-speed flag")
    c = replace_once(c,
        "&& NID(simpar) == 1 && !w2lin && !sedov_lagvol_on();",
        "&& NID(simpar) == 1 && !w2lin && !volw2 && !sedov_lagvol_on();",
        "exam.c dynamic mode exclusion")
    if c.count("if(dynw2) gfs_dynw2_set(") != 3:
        raise ValueError("Expected three dynamic-weight RK4 stage setters")
    c = c.replace("if(dynw2) gfs_dynw2_set(",
                  "if(dynw2 || volw2) gfs_dynw2_set(")
    c = replace_once(c,
        "/* GFS_DYN_W2: W_{n+1} = W_n + Dt (k1 + 2k2 + 2k3 + k4)/6. */\n\tif(dynw2){",
        "/* GFS_DYN_W2 / GPU_VOLW2: RK4 combination of stage rates. */\n\tif(dynw2 || volw2){",
        "exam.c final weighted RK4 combination")
    c = replace_once(c,
        "\tgfs_w2lin_active = 0;\n\tgfs_dynw2_active = 0;\n",
        "\tgfs_w2lin_active = 0;\n\tgfs_dynw2_active = 0;\n\tgfs_volw2_active = 0;\n",
        "exam.c mode reset")
    updated["exam.c"] = c

    cu = updated["exam_gpu.cu"]
    old = "if (jbp_is_ghost || ibp_w2 != 0.0f || pw2[j] != 0.0f) {"
    new = ("if (jbp_is_ghost || ((ibp_w2 != 0.0f || pw2[j] != 0.0f) && "
           "!(w2lin && laguerre_rot && geom_fv))) {")
    cu = replace_once(cu,old,new,"GPU weighted variational gate")
    cu = replace_once(cu,
        'extern "C" int gpu_varpress_w0_enabled(void) {',
        'extern "C" int gfs_volw2_is_active(void);\n'
        'extern "C" int gpu_varpress_w0_enabled(void) {',
        "GPU weighted runtime declaration")
    cu = replace_once(cu,
        '        params->w2lin || params->laguerre_rot) {',
        '        ((params->w2lin || params->laguerre_rot) &&\n'
        '         !(gfs_volw2_is_active() && params->w2lin &&\n'
        '           params->laguerre_rot && params->geom_fv))) {',
        "W0 variational guard permits only the gated weighted mode")
    cu = replace_once(cu,
        '        fprintf(stderr, "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active\\n");',
        '        fprintf(stderr, gfs_volw2_is_active() ?\n'
        '            "[VARPRESS_W0] gated weighted fixed-mass variational pressure active\\n" :\n'
        '            "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active\\n");',
        "weighted pressure banner")
    updated["exam_gpu.cu"] = cu

    header = updated["exam_gpu.h"]
    old = "void gpu_download_results(GPUContext *ctx, int n_particles);\n"
    new = old + (
        "/* Experimental fixed-mass power-weight rate, global work projection. */\n"
        "int gpu_probe_energy_neutral_rate(GPUContext *ctx,\n"
        "    const unsigned char *host_active, int n_particles, int n_total,\n"
        "    int iterations, double omega, double *host_rate,\n"
        "    double *host_gradient, double *host_scalar,\n"
        "    unsigned long long *host_unsupported);\n"
        "int gpu_probe_upload_weight_rate(GPUContext *ctx,\n"
        "    const double *host_rate, int n_particles, int n_total);\n")
    updated["exam_gpu.h"] = replace_once(header,old,new,"GPU header prototypes")

    make = updated["Makefile"]
    make = replace_once(make,
        "EXOBJ      += exam_gpu.o exam_gpu_extract.o\n",
        "EXOBJ      += exam_gpu.o exam_gpu_extract.o gfs_energy_neutral_probe.o\n",
        "GPU object list")
    old = ("exam_gpu_extract.o: exam_gpu_extract.c exam_gpu.h gfs_thermal64.h gfs_pair.h "
           "gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h "
           "gfs_enriched_euler.h gfs_enriched_coupling.h\n"
           "\t$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@\n")
    new = old + ("\ngfs_energy_neutral_probe.o: gfs_energy_neutral_probe.cu "
                 "kepler_gpu_jacobi_probe_20261027.cu exam_gpu.h\n"
                 "\t$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@\n")
    updated["Makefile"] = replace_once(make,old,new,"GPU probe object rule")

    extract = updated["exam_gpu_extract.c"]
    old = "int gfs_w2lin_is_active(void);\n"
    new = old + ("int gfs_volw2_is_active(void);\n"
                 "int gfs_volw2_stage_current(void);\n"
                 "void gfs_volw2_store_stage(SimParameters *, int, const double *);\n")
    extract = replace_once(extract,old,new,"GPU extract stage externs")
    old = "static int gpu_call_count = 0;\n"
    helper = r'''/* Opt-in weighted RK4 stage: stable-ID mask, work-neutral Wdot, no padding. */
static void gpu_volw2_stage(SimParameters *simpar, int nbp, int n_total)
{
    if(!gfs_volw2_is_active()) return;
    static unsigned char *by_id=0;
    static size_t mask_count=0;
    if(NID(simpar)!=1 || n_total!=nbp || MYID(simpar)!=0)
        MPI_Abort(MPI_COMM(simpar),97);
    if(!by_id){
        const char *path=getenv("GFS_GPU_VOLW2_MASK");
        if(!path || !path[0]) MPI_Abort(MPI_COMM(simpar),97);
        FILE *input=fopen(path,"rb");
        if(!input || fseek(input,0,SEEK_END)!=0) MPI_Abort(MPI_COMM(simpar),97);
        long count=ftell(input);
        if(count<=0 || count>100000000 || fseek(input,0,SEEK_SET)!=0)
            MPI_Abort(MPI_COMM(simpar),97);
        mask_count=(size_t)count;
        by_id=(unsigned char*)malloc(mask_count);
        if(!by_id || fread(by_id,1,mask_count,input)!=mask_count)
            MPI_Abort(MPI_COMM(simpar),97);
        fclose(input);
    }
    unsigned char *active=(unsigned char*)calloc((size_t)n_total,1);
    double *rate=(double*)malloc((size_t)nbp*sizeof(double));
    double *gradient=(double*)malloc((size_t)nbp*sizeof(double));
    if(!active || !rate || !gradient) MPI_Abort(MPI_COMM(simpar),97);
    const size_t stride=TVORORK4_DDINFO(simpar)[0].n_size;
    char *base=(char*)VORORK4_TBP(simpar);
    int nactive=0;
    for(int i=0;i<nbp;++i){
        treevorork4particletype *bp=(treevorork4particletype*)(base+(size_t)i*stride);
        unsigned long long id=((unsigned long long)PINDX(bp)) & ((1ULL<<56)-1);
        if(id>=mask_count) MPI_Abort(MPI_COMM(simpar),97);
        active[i]=by_id[id] ? 1 : 0;
        nactive+=active[i];
    }
    if(nactive<=0) MPI_Abort(MPI_COMM(simpar),97);
    unsigned long long unsupported=0;
    double alpha=0.;
    const double started=MPI_Wtime();
    int status=gpu_probe_energy_neutral_rate(&g_gpu_ctx,active,nbp,n_total,
        16,.8,rate,gradient,&alpha,&unsupported);
    if(status || unsupported){
        fprintf(stderr,"[GPU_VOLW2] stage=%d rate status=%d unsupported=%llu\n",
            gfs_volw2_stage_current(),status,unsupported);
        MPI_Abort(MPI_COMM(simpar),97);
    }
    gfs_volw2_store_stage(simpar,nbp,rate);
    status=gpu_probe_upload_weight_rate(&g_gpu_ctx,rate,nbp,n_total);
    if(status) MPI_Abort(MPI_COMM(simpar),97);
    long double work=0.;
    for(int i=0;i<nbp;++i){
        g_gpu_ctx.h_parts.w2old[i]=(float)rate[i];
        if(active[i]) work+=(long double)gradient[i]*rate[i];
    }
    fprintf(stderr,"[GPU_VOLW2] stage=%d active=%d faces=%d alpha=%.9g work=%.3Le elapsed=%.6f s\n",
        gfs_volw2_stage_current(),nactive,g_gpu_ctx.d_faces.n_faces_total,
        alpha,work,MPI_Wtime()-started);
    free(gradient); free(rate); free(active);
}

'''
    extract = replace_once(extract,old,helper+old,"GPU stage call helper")
    old = "    double t4_kernel = MPI_Wtime();\n    /* --- Launch force kernel --- */\n"
    extract = replace_once(extract,old,
        "    gpu_volw2_stage(simpar, nbp, n_total);\n"+old,
        "GPU force stage attachment")
    old = "params.laguerre_rot = (lr && lr[0] && atoi(lr) > 0) ? 1 : 0;"
    if extract.count(old) != 2:
        raise ValueError("Expected two GPU/validation centroid-speed assignments")
    extract = extract.replace(old,
        "params.laguerre_rot = gfs_volw2_is_active() ? 1 : (lr && lr[0] && atoi(lr) > 0) ? 1 : 0;")
    old = "params.geom_fv = (gv && gv[0] && atoi(gv) > 0) ? 1 : 0;"
    if extract.count(old) != 2:
        raise ValueError("Expected two GPU/validation geometric-work assignments")
    extract = extract.replace(old,
        "params.geom_fv = gfs_volw2_is_active() ? 1 : (gv && gv[0] && atoi(gv) > 0) ? 1 : 0;")
    updated["exam_gpu_extract.c"] = extract

    source_cu = DOCS/"kepler_gpu_energy_neutral_probe_20261027.cu"
    base_cu = DOCS/"kepler_gpu_jacobi_probe_20261027.cu"
    targets = (root/"Exam/gfs_energy_neutral_probe.cu",
               root/"Exam/kepler_gpu_jacobi_probe_20261027.cu")
    if any(path.exists() for path in targets):
        raise FileExistsError("New CUDA sources already exist in worktree")
    for name,path in paths.items():
        path.write_text(updated[name])
    shutil.copy2(source_cu,targets[0])
    shutil.copy2(base_cu,targets[1])
    receipt = {
        "worktree": str(root),
        "mode": "experimental one-rank fixed-mass pressure-work-neutral weighted RK4",
        "original_sha256": {name:hashlib.sha256(original[name].encode()).hexdigest()
                            for name in original},
        "patched_sha256": {name:digest(path) for name,path in paths.items()},
        "new_cuda_sha256": {str(path.relative_to(root)):digest(path) for path in targets},
        "scope": "Source attachment only. Weighted hydro is not accepted until full build and one-step native conservation and geometry gates pass.",
    }
    with marker.open("x") as stream:
        json.dump(receipt,stream,indent=2)
        stream.write("\n")
    print(json.dumps(receipt,indent=2))


if __name__ == "__main__":
    main()
