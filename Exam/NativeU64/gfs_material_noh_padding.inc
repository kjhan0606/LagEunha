/* Physical exterior samples from immutable IC grid coordinates. The cells
 * being evolved remain the original N² material cells. Exterior samples are
 * prescribed boundary data, never counted in owned mass/checkpoints. */
static void gfs_material_noh_padding(SimParameters*s,postype unused_width){
 (void)unused_width;GfsNohExterior b;
 if(!gfs_material_noh_configuration(s,&b)||
    TVORORK4_DDINFO(s)[0].n_size!=sizeof(treevorostressrk4particletype)){
  fprintf(stderr,"[NOH_NATIVE_BOUNDARY] invalid configuration or missing stage time\n");MPI_Abort(MPI_COMM(s),99);return;
 }
 int N=NX(s),n=N*N,npad=3*n;double h=Xmax_HydroExam(s)/N,V=h*h;
 /* Match the native padding callback ownership contract: the caller frees
  * prior padding after the RHS/endpoint. The initial padding slot is not an
  * owned allocation and must not be inspected or freed here. */
 VORORK4_TBPP(s)=NULL;VORO_NPAD(s)=0;
 treevorostressrk4particletype*p=my_malloc((size_t)npad*sizeof(*p));memset(p,0,(size_t)npad*sizeof(*p));
 int k=0;
 for(int y=-N/2;y<3*N/2;y++)for(int x=-N/2;x<3*N/2;x++){
  if(x>=0&&x<N&&y>=0&&y<N)continue;
  double x0[3]={(x+.5)*h,(y+.5)*h,0},pos[3];GfsFvPrimitive q;
  if(!gfs_noh_exterior_position(&b,x0,pos)||!gfs_noh_exterior_state(&b,pos,&q)){
   fprintf(stderr,"[NOH_NATIVE_BOUNDARY] exterior left upstream regime at t=%.17g\n",b.time);MPI_Abort(MPI_COMM(s),99);return;
  }
  treevorostressrk4particletype*a=p+k;
  a->x=pos[0];a->y=pos[1];a->vx=q.u[0];a->vy=q.u[1];a->mass=q.rho*V;
  a->u4if.indx=n+k;SET_FLAG(a,VOROflag);SET_FLAG(a,BoundaryGhostflag);
  a->weight=(GfsWeightState){0,0,1,GFS_WEIGHT_STATE_TAG};
  a->weight_thermal=(GfsWeightThermal){q.p*V/(b.gamma-1),0,V};
  a->volume=V;a->ie=a->weight_thermal.U;a->den=q.rho;a->pressure=q.p;
  a->csound=sqrt(b.gamma*q.p/q.rho);a->w2=a->w2old=h*h;
  a->stress.K=q.p/pow(q.rho,b.gamma);
  a->stress.E_inv_xx=a->weight_thermal.U;a->stress.E_inv_xy=V;
  a->stress.E_inv_yx=0;a->stress.E_inv_yy=2026102603.;k++;
 }
 if(k!=npad)MPI_Abort(MPI_COMM(s),99);
 VORORK4_TBPP(s)=(void*)p;VORO_NPAD(s)=npad;
}
static void mkLinkedList2D_material_noh(SimParameters*s,postype cell,postype xmin,
 postype ymin,postype xmax,postype ymax,void(*unused)(SimParameters*,postype)){
 (void)unused;
 mkLinkedList2D_oExam(s,cell,xmin,ymin,xmax,ymax,gfs_material_noh_padding);
}
