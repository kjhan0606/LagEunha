/* Same RK4 equations, one float cast from accepted U for each stage.
 * E_inv_xx is a step-local double anchor only under the narrow av5 contract.
 * It travels with the particle through migration and is reinitialized at
 * every accepted step boundary, including native restart. No entropy repair.
 */
static int gfs_rk4_direct_u_on(SimParameters *simpar){
    const char *s=getenv("GFS_RK4_DIRECT_U");
    if(!s || atoi(s)<=0)return 0;
    if(gfs_thermal64_requested()){
        fprintf(stderr,"[RK4_DIRECT_U] conflicts with authoritative thermal64 payload\n");
        MPI_Abort(MPI_COMM(simpar),99);return 0;
    }
    gfs_es_config();
#ifdef USE_CUDA
    if(GAS_GPU_ENABLED(simpar) && gpu_varpress_w0_enabled() &&
       getenv("EUNHA_IC") && !strcmp(getenv("EUNHA_IC"),"kepler") &&
       getenv("GFS_GPU_RAW_FACES") && atoi(getenv("GFS_GPU_RAW_FACES"))==1 &&
       getenv("GFS_GPU_FINAL_VOLUME") && atoi(getenv("GFS_GPU_FINAL_VOLUME"))==1 &&
       GAS_AVMODE(simpar)==5 && GAS_Kappa(simpar)==0 && GAS_ENTROPY_MODE(simpar)==0 &&
       GAS_XSPHEPS(simpar)==0 && GAS_FCENTROID(simpar)==0 &&
       sedov_phase1_on() && gfs_dual_energy_eta()==0 && !gfs_es_on && !gfs_half_on &&
       !gfs_dynw2_on() && !gfs_w2lin_on() && !sedov_lagvol_on() && !gfs_vaud_on() &&
       TVORORK4_DDINFO(simpar)[0].n_size==sizeof(treevorostressrk4particletype))return 1;
#endif
    fprintf(stderr,"[RK4_DIRECT_U] requires native smooth W0/av5 with no conflicting entropy/volume/weight/stress representation\n");
    MPI_Abort(MPI_COMM(simpar),99);return 0;
}
