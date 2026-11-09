/* FaceCSR/ParticleSoA bridge shared by both retained native variants. */
int gfs_material_import_facecsr(const FaceCSR *faces,const ParticleSoA *parts,
 struct GfsMaterialCsr2D *mesh){
 if(!faces||!parts||!mesh||faces->n_particles<1||parts->n_total<faces->n_particles)return 0;
 int nf=faces->n_faces_total,nc=faces->n_particles;
 if(nf!=faces->face_offset[nc])return 0;
 int *boundary=(int*)malloc((size_t)nf*sizeof(int));if(!boundary)return 0;
 for(int f=0;f<nf;f++){
  boundary[f]=faces->is_ghost[f]||faces->neighbor_idx[f]<0;
  /* Exterior rank/periodic images require explicit vertex ownership. */
  if(!boundary[f]&&faces->neighbor_idx[f]>=nc){free(boundary);return 0;}
 }
 int ok=gfs_material_csr_import(nc,faces->face_offset,parts->x,parts->y,
     faces->c1x,faces->c1y,faces->c2x,faces->c2y,boundary,mesh);
 free(boundary);return ok;
}
