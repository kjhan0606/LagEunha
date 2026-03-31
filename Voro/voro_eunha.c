#include <stdio.h>
#include "eunha.h"
#include "voro.h"
#include "voro_eunha.h"


postype findSDen3D(Voro3D_Vertex *vertex, int ip, Voro3D_point *neighwork,
		int *intwork, int mp, Voro3D_GasParticle *bp){
	int i;
	for(i=0;i<mp;i++) intwork[i] = 0;
	for(i=0;i<ip;i++){
		intwork[vertex[i].related[0]] = 1;
		intwork[vertex[i].related[1]] = 1;
		intwork[vertex[i].related[2]] = 1;
	}
	postype tweight = 0;
	postype hsml=0.1;
	postype den = getWeight3D(0.,0.,0., hsml);
	for(i=0;i<mp;i++){
		if(intwork[i] ==1){
			int jndx = neighwork[i].indx;
			postype weight = getWeight3D(neighwork[i].x, neighwork[i].y,neighwork[i].z, hsml);
			den += bp[jndx].mass/bp[jndx].volume * weight;
			tweight += weight;
		}
	}
	den = den/tweight;
	return den;
}

postype getAvgPressure2D(int indx, Voro2D_Corner *voro, Voro2dGasParticle *bp, 
		Voro2D_point *neigh){
	int jndx = neigh[voro->upperrelated].indx;
	int j1id = neigh[voro->upperlink->upperrelated].indx;
	int j2id = neigh[voro->lowerrelated].indx;
	postype res = (5./12.)*(bp[jndx].pressure + bp[indx].pressure)
		+ (1./12.)*(bp[j1id].pressure + bp[j2id].pressure);
	return res;
}
void get2dAreaAvgNeighorPressure(
        treevorork4particletype *ibp,
        Voro2D_Corner *vorocorner,
        Voro2D_point *neigh,
        treevorork4particletype *bp
        ){
    int mp=1;
    ibp->volume = Area2DPolygon(vorocorner, mp); // mp is obsolete here.
    postype avgNeighboringPressure = 0;
    postype Length = 0;
    Voro2D_Corner *tmp = vorocorner;
    do {
        if(tmp->upperrelated >= 0){
            postype dx = tmp->x-(tmp->upperlink)->x;
            postype dy = tmp->y-(tmp->upperlink)->y;
            postype ds = sqrt(dx*dx+dy*dy);
            Length += ds;
            treevorork4particletype *jbp = neigh[tmp->upperrelated].bp;
            avgNeighboringPressure += ds*jbp->pressure;
        }
		tmp = tmp->upperlink;
    } while(vorocorner != tmp);
    if(Length > 0)
        ibp->avgNeighboringPressure = avgNeighboringPressure/Length;
    else
        ibp->avgNeighboringPressure = ibp->pressure;
}


Voro2D_point get2dUpqradRk4(treevorork4particletype *p, treevorork4particletype *q,
		postype dtold){
	Voro2D_point uradpq,upq, er;
	upq.x = q->vx - p->vx;
	upq.y = q->vy - p->vy;
	postype absupq = sqrt(upq.x*upq.x + upq.y*upq.y);
	er.x = q->x - p->x;
	er.y = q->y - p->y;
	postype dpq2 = (er.x*er.x + er.y*er.y);
	postype dpq = sqrt(dpq2);
	er.x /= dpq;
	er.y /= dpq;
	postype wp2,wq2,wpold2,wqold2;
	wp2 = p->w2;
	wq2 = q->w2;
	wpold2 = p->w2old;
	wqold2 = q->w2old;
	postype fact1 = 0.5*(1+(wp2-wq2)/dpq2);
	postype dwpdt = (sqrt(wp2)-sqrt(wpold2))/dtold;
	postype dwqdt = (sqrt(wq2)-sqrt(wqold2))/dtold;
	postype vpw, vqw;
	if(dwpdt >0) vpw = (p->csound < dwpdt ? p->csound: dwpdt);
	else vpw = (-p->csound > dwpdt ? p->csound: dwpdt);
	if(dwqdt >0) vqw = (q->csound < dwqdt ? q->csound: dwqdt);
	else vqw = (-q->csound > dwqdt ? q->csound: dwqdt);

	postype fact2 = (sqrt(wp2)*vpw - sqrt(wq2)*vqw)/dpq;
	fact2 = fact2 - (wp2-wq2)/dpq2 * absupq; 
	uradpq.x = fact1*upq.x + fact2*er.x;
	uradpq.y = fact1*upq.y + fact2*er.y;
	return uradpq;
}

Voro2D_point get2dUpqrad(treevoroparticletype *p, treevoroparticletype *q,
		postype dtold){
	Voro2D_point uradpq,upq, nr;
	upq.x = q->vx - p->vx;
	upq.y = q->vy - p->vy;
	postype absupq = sqrt(upq.x*upq.x + upq.y*upq.y);
	nr.x = q->x - p->x;
	nr.y = q->y - p->y;
	postype dpq2 = (nr.x*nr.x + nr.y*nr.y);
	postype dpq = sqrt(dpq2);
	nr.x /= dpq;
	nr.y /= dpq;
	postype wp2,wq2,wpold2,wqold2;
	wp2 = p->w2;
	wq2 = q->w2;
	wpold2 = p->w2old;
	wqold2 = q->w2old;
	postype fact1 = 0.5*(1+(wp2-wq2)/dpq2);
    postype dwpdt = (sqrt(wp2)-sqrt(wpold2))/dtold;
    postype dwqdt = (sqrt(wq2)-sqrt(wqold2))/dtold;
    postype vpw = (p->csound < dwpdt ? p->csound: dwpdt);
    postype vqw = (q->csound < dwqdt ? q->csound: dwqdt);
    postype fact2 = (sqrt(wp2)*vpw - sqrt(wq2)*vqw)/(dpq+dpq);
    fact2 = fact2 - (wp2-wq2)/dpq2 * absupq;
	uradpq.x = fact1*upq.x + fact2*nr.x;
	uradpq.y = fact1*upq.y + fact2*nr.y;
	return uradpq;
}

Voro3D_point get3dUpqradRk4(treevorork4particletype *p, treevorork4particletype *q,
		postype dtold){
	Voro3D_point uradpq,upq, nr;
	upq.x = q->vx - p->vx;
	upq.y = q->vy - p->vy;
	upq.z = q->vz - p->vz;
	nr.x = q->x - p->x;
	nr.y = q->y - p->y;
	nr.z = q->z - p->z;
	postype dpq2 = (nr.x*nr.x + nr.y*nr.y + nr.z*nr.z);
	postype dpq = sqrt(dpq2);
	nr.x /= dpq;
	nr.y /= dpq;
	postype wp2,wq2,wpold2,wqold2;
	wp2 = p->w2;
	wq2 = q->w2;
	wpold2 = p->w2old;
	wqold2 = q->w2old;
	postype fact1 = 0.5*(1-(wp2-wq2)/dpq2);
	postype fact2 = 0.5/(dpq*dtold)* (wp2-wpold2 - (wq2-wqold2));
	uradpq.x = fact1*upq.x + fact2*nr.x;
	uradpq.y = fact1*upq.y + fact2*nr.y;
	uradpq.z = fact1*upq.z + fact2*nr.z;
	return uradpq;
}
Voro3D_point get3dUpqrad(treevoroparticletype *p, treevoroparticletype *q,
		postype dtold){
	Voro3D_point uradpq,upq, nr;
	upq.x = q->vx - p->vx;
	upq.y = q->vy - p->vy;
	upq.z = q->vz - p->vz;
	nr.x = q->x - p->x;
	nr.y = q->y - p->y;
	nr.z = q->z - p->z;
	postype dpq2 = (nr.x*nr.x + nr.y*nr.y + nr.z*nr.z);
	postype dpq = sqrt(dpq2);
	nr.x /= dpq;
	nr.y /= dpq;
	postype wp2,wq2,wpold2,wqold2;
	wp2 = p->w2;
	wq2 = q->w2;
	wpold2 = p->w2old;
	wqold2 = q->w2old;
	postype fact1 = 0.5*(1.-(wp2-wq2)/dpq2);
	postype fact2 = 0.5/(dpq*dtold)* (wp2-wpold2 - (wq2-wqold2));
	uradpq.x = fact1*upq.x + fact2*nr.x;
	uradpq.y = fact1*upq.y + fact2*nr.y;
	uradpq.z = fact1*upq.z + fact2*nr.z;
	return uradpq;
}

