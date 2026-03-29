#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "color.h"
#include "makemap.h"
typedef HydroTreeLinkedCell CellType;

int makemap(SimParameters *simpar, int icount, int nbin, char *outfile, char *saofile){
    postype cellsize = KP_GridSize(simpar);
    postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
    postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
    int nximg = NX(simpar)/nbin;
    int nyimg = NY(simpar)/nbin;
    float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
    float *img = (float*)malloc(sizeof(float)*nximg*nyimg);
    float *ndist = (float*)malloc(sizeof(float)*nximg*nyimg);
    int i,j;
    int ii,jj;
    for(i=0;i<nximg*nyimg;i++) map[i] = 0;
    for(i=0;i<nximg*nyimg;i++) ndist[i] = 1.e20;
    postype pixsize = Lx/nximg;
    postype xmin,ymin,xmax,ymax;
    xmin = KP_XMIN(simpar);
    ymin = KP_YMIN(simpar);
    xmax = KP_XMAX(simpar);
    ymax = KP_YMAX(simpar);
    int mx,my;
    mx = ceil((xmax-xmin)/cellsize);
    my = ceil((ymax-ymin)/cellsize);
    CellType *cells = (CellType*)malloc(sizeof(CellType)*mx*my);
    for(i=0;i<mx*my;i++){
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
    treevoroparticletype *bp = VORO_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        int index = ix+mx*iy;
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }

    for(j=0;j<nyimg;j++){
        postype yp = (j+0.5)*pixsize; // img pixel position
        if(yp < ymin || yp >= ymax) continue;
        for(i=0;i<nximg;i++){
            postype xp = (i+0.5)*pixsize; // img pixel position
            if(xp < xmin || xp >= xmax) continue;
            postype nearden = 0;
            postype idist = 1.e20;
            int jy,ix;
            ix = (xp-xmin)/cellsize;
            jy = (yp-ymin)/cellsize;
            for(jj=jy-1;jj<=jy+1;jj++){
                if(jj<0 || jj>=my) continue;
                for(ii=ix-1;ii<=ix+1;ii++){
                    if(ii<0 || ii>=mx) continue;
                    size_t ipixel = ii + mx*jj;
                    struct linkedlisttype *tmp = cells[ipixel].link;
                    while(tmp){
                        treevoroparticletype *tt = (treevoroparticletype*)tmp;
                        if(tt->x >= xmin && tt->x < xmax && tt->y >=ymin && tt->y < ymax){
                            postype distx = fabs(tt->x - xp);
                            postype disty = fabs(tt->y - yp);
                            postype dist2 = distx*distx + disty*disty;
                            if(dist2 < idist) {
                                idist = dist2;
                                nearden = tt->den;
                            }
                        }
                        tmp = tmp->next;
                    }

                }
            }
            map[i+nximg*j] = nearden;
            ndist[i+nximg*j] = idist;

        }
    }
    int id;
    if(MYID(simpar)==0){
        float *idist = (float*)malloc(sizeof(float)*nximg*nyimg);
        MPI_Status status;

        for(id=1;id<NID(simpar);id++){
            MPI_Recv(idist, nximg*nyimg,MPI_FLOAT,id,id, MPI_COMM_WORLD, &status);
            MPI_Recv(img, nximg*nyimg,MPI_FLOAT,id,id, MPI_COMM_WORLD, &status);
            for(j=0;j<nximg*nyimg;j++){
                if(idist[j] < ndist[j]){
                    map[j] = img[j];
                    ndist[j] = idist[j];
                }
            }
        }
        free(idist);
    }
    else {
        MPI_Send(ndist, nximg*nyimg,MPI_FLOAT,0,MYID(simpar), MPI_COMM_WORLD);
        MPI_Send(map, nximg*nyimg,MPI_FLOAT,0,MYID(simpar), MPI_COMM_WORLD);
    }

    if(MYID(simpar)==0){
        float *img2 = (float*)malloc(sizeof(float)*nximg*nyimg);;
        for(j=0;j<nyimg;j++){
            for(i=0;i<nximg;i++){
                img2[i+nximg*(nyimg-j-1)] = map[i+nximg*j];
            }
        }
        colorizeit(7,img2, nximg,nyimg,saofile, outfile, 0., 2.);
        free(img2);
    }
    free(cells);
    free(map); free(img);
    free(ndist);
    return 0;
}

int maketscmap(SimParameters *simpar, int icount, int nbin, char *outfile, char *saofile){
    int nx = NX(simpar);
    int ny = NY(simpar);
    int np = VORO_NP(simpar);
    int nximg, nyimg;
    nximg = nx/nbin;
    nyimg = ny/nbin;
    treevoroparticletype *bp = VORO_TBP(simpar);
    float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
    float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
    float cellsize = Lx/nximg;
    int i;
    float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
    float *tmap = (float*)malloc(sizeof(float)*nximg*nyimg);
    for(i=0;i<nximg*nyimg;i++) map[i] = 0;
    for(i=0;i<nximg*nyimg;i++) tmap[i] = 0;
    postype mscale = (nximg/Lx)*(nyimg/Ly);

    for(i=0;i<np;i++){
        postype pmass = mscale*bp[i].mass;
        postype xp = bp[i].x/cellsize;
        postype yp = bp[i].y/cellsize;
        int nearx = rint(xp);
        int neary = rint(yp);
        int ic = (nearx + nximg) %nximg;
        int jc = (neary + nyimg) %nyimg;
        postype xmin = xp -nearx;
        postype ymin = yp -neary;
        int icc = (nximg+nearx + (int)(copysign(1., xmin))) %nximg;
        int jcc = (nyimg+neary + (int)(copysign(1., ymin))) %nyimg;
        int iccc = (nximg+nearx - (int)(copysign(1., xmin))) %nximg;
        int jccc = (nyimg+neary - (int)(copysign(1., ymin))) %nyimg;
        float xd1 = fabs(xmin);
        float yd1 = fabs(ymin);
        float wx1 = (0.75-xd1*xd1)*pmass;
        float wy1 = (0.75-yd1*yd1);
        float wx3 = 0.5*pmass * (0.25+xd1*(xd1-1.));
        float wx2 = wx3 + pmass * xd1;
        float wy3 = 0.5 * (0.25+yd1*(yd1-1.));
        float wy2 = wy3 + yd1;

        map[ic+nximg*jc] += wx1*wy1;
        map[icc+nximg*jc] += wx2*wy1;
        map[iccc+nximg*jc] += wx3*wy1;
        map[ic+nximg*jcc] += wx1*wy2;
        map[icc+nximg*jcc] += wx2*wy2;
        map[iccc+nximg*jcc] += wx3*wy2;
        map[ic+nximg*jccc] += wx1*wy3;
        map[icc+nximg*jccc] += wx2*wy3;
        map[iccc+nximg*jccc] += wx3*wy3;

    }
    MPI_Reduce(map, tmap, nximg*nyimg, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(MYID(simpar)==0){
        colorizeit(7,tmap, nximg,nyimg,saofile, outfile, 0., 2.);
    }
    free(map);
    free(tmap);

    return 0;
}

