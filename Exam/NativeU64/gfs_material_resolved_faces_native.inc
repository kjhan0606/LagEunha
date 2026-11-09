/* One resolved-face contract for pressure work and volume Jacobians.
 * This is the existing pressure adapter's construction, without new cutoffs. */
static int gfs_material_resolved_faces2d(GPUContext *ctx,int n,int nt,
 double lx,double ly,const int *alias,const unsigned char *image_reflection,
 GfsTraceDirected2D *edge,unsigned char *walls,int *out_count){
 ParticleSoA *p=&ctx->h_parts;FaceCSR *f=&ctx->h_faces;int count=0;
 for(int i=0;i<n;i++)for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
  int j=f->neighbor_idx[k];double area=hypot(f->c2x[k]-f->c1x[k],f->c2y[k]-f->c1y[k]);
  if(j<0||j>=nt){return 0;}
  double canonical_dx=p->x[j]-p->x[i],canonical_dy=p->y[j]-p->y[i];
  if(!image_reflection||!image_reflection[j]){
   canonical_dx=material_canonical_delta(p->x[i],p->x[alias[j]],p->x[j],lx);
   canonical_dy=image_reflection?p->y[alias[j]]-p->y[i]:material_canonical_delta(p->y[i],p->y[alias[j]],p->y[j],ly);
  }
  /* Direct FV must retain resolved short faces: a fixed absolute cutoff
   * removes finite volume on refined meshes. Use the CSR's existing
   * coordinate-roundoff scale, retaining all resolved face contributions. */
#ifdef GFS_MATERIAL_DIRECT_FV
  double area_cut=256*DBL_EPSILON*hypot(canonical_dx,canonical_dy);
#else
  double area_cut=1e-12;
#endif
  if(area<area_cut)continue;
  int wall=image_reflection?image_reflection[j]:0;
  if((wall&&alias[j]!=i)||(!wall&&(f->is_ghost[k]||alias[j]==i))){return 0;}
  double dx=canonical_dx,dy=canonical_dy;
  double cx=.5*(f->c1x[k]+f->c2x[k]),cy=.5*(f->c1y[k]+f->c2y[k]);
  double d=hypot(dx,dy),plane=.5*(dx*dx+dy*dy+p->weight_W[i]-p->weight_W[alias[j]]);
  if(!(d>0)||fabs(dx*cx+dy*cy-plane)>1e-11*d){fprintf(stderr,"[MATERIAL_TRACE_FAIL] phase=plane i=%d j=%d error=%.17g distance=%.17g\n",i,j,dx*cx+dy*cy-plane,d);return 0;}
  /* A shared face must not survive on only one side of the same cutoff.
   * Search only the roundoff-width band, keeping the existing geometry
   * tolerance and cutoff. A genuinely unmatched face still fails CSR. */
  if(!wall&&area<=area_cut+256*DBL_EPSILON*d){
   double tol=256*DBL_EPSILON*d;int omit=0;
   for(int h=f->face_offset[alias[j]];h<f->face_offset[alias[j]+1];h++){
    int b=f->neighbor_idx[h];if(b<0||b>=nt||alias[b]!=i)continue;
    double bx=material_canonical_delta(p->x[alias[j]],p->x[alias[b]],p->x[b],lx);
    double by=image_reflection?p->y[alias[b]]-p->y[alias[j]]:material_canonical_delta(p->y[alias[j]],p->y[alias[b]],p->y[b],ly);
    if(fabs(bx+dx)>tol||fabs(by+dy)>tol)continue;
    double other=hypot(f->c2x[h]-f->c1x[h],f->c2y[h]-f->c1y[h]);
    double bcx=.5*(f->c1x[h]+f->c2x[h]),bcy=.5*(f->c1y[h]+f->c2y[h]);
    if(other<area_cut&&fabs(other-area)<=tol&&fabs(bcx-(cx-dx))<=tol&&fabs(bcy-(cy-dy))<=tol)omit=1;
   }
   if(omit)continue;
  }
  if(walls)walls[count]=(unsigned char)wall;
  edge[count++]=(GfsTraceDirected2D){i,alias[j],dx,dy,cx,cy,area};
 }
 *out_count=count;return 1;
}
