#include "eunha.h"
#include "voro.h"

void voro3D_EulerRot(Voro3D_point *rin, int np, ptype theta, ptype phi, ptype psi){
	Voro3D_point res,*r;
	int i;
	r = rin;
	for(i=0;i<np;i++){
		res.x = cos(theta)*cos(psi)*r->x + 
			(-cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi))*r->y + 
			( sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi))*r->z;
		res.y = cos(theta)*sin(psi)*r->x + 
			( cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi))*r->y + 
			(-sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi))*r->z;
		res.z = -sin(theta)        *r->x + 
			( sin(phi)*sin(theta)                             )*r->y + 
			( cos(phi)*cos(theta)                             )*r->z;
		rin[i] = res;
		r ++;
	}
}


Voro3D_point voro3D_norm(Voro3D_point *rin, int np){
	Voro3D_point a, b,res;
	b = Vec3DSub(rin+1, rin  );
	a = Vec3DSub(rin+2, rin+1);
	res = Vec3DCross(&b,&a);
	ptype amp = sqrt(Vec3DInnP(&res,&res));
	res.x /= amp; res.y /= amp; res.z /= amp;
	return res;
}
Voro2D_point voro2D_norm(Voro2D_point *rin){
	Voro3D_point a,b,c;
	Voro2D_point res;
	b.x =  rin->y;
	b.y = -rin->x;
	b.z = 0;
	a.x =  rin->x;
	a.y =  rin->y;
	a.z = 0;
	c = Vec3DCross(&b, &a);
	if(c.z<0) {
		b.x = - b.x;
		b.y = - b.y;
	}
	ptype amp =  sqrt(Vec3DInnP(&b,&b));
	res.x = b.x/amp; res.y = b.y/amp;
	return res;
}

/*
const Voro3D_Vertex p3dguess[8] = {
	{1, {-5 -1 -3}, {0x4,0x1,0x3}, 0,{0, 0, 0}, -1, -1, -1}, 
	{1, {-5 -3 -2}, {0x7,0x2,(0x0)}, 0,{0, 0, 0}, -1, -1, 1}, 
	{1, {-5 -2 -4}, {0x6,0x3,0x1}, 0,{0, 0, 0}, -1, 1, 1}, 
	{1, {-5 -4 -1}, {0x5,(0x0),0x2}, 0,{0, 0, 0}, -1, 1, -1}, 
	{1, {-6 -3 -1}, {(0x0),0x5,0x7}, 0,{0, 0, 0}, 1, -1, -1}, 
	{1, {-6 -1 -4}, {0x3,0x6,0x4}, 0,{0, 0, 0}, 1, 1, -1}, 
	{1, {-6 -4 -2}, {0x2,0x7,0x5}, 0,{0, 0, 0}, 1, 1, 1}, 
	{1, {-6 -2 -3}, {0x1,0x4,0x6}, 0,{0, 0, 0}, 1, -1, 1},

};

Voro2D_Corner Voro2D_FindInterSec(Voro2D_point *ctrCenter, Voro2D_point *v, Voro2D_Corner *a, Voro2D_Corner *b){
	Voro2D_point umv,bma,upv,amupv;

	upv.x = ctrCenter->x + v->x;
	upv.y = ctrCenter->y + v->y;

	umv.x = ctrCenter->x - v->x;
	umv.y = ctrCenter->y - v->y;

	bma.x = b->x - a->x;
	bma.y = b->y - a->y;

	amupv.x = a->x - 0.5*upv.x;
	amupv.y = a->y - 0.5*upv.y;

	ptype t = (umv.x*amupv.x + umv.y*amupv.y)/(bma.x*umv.x+bma.y*umv.y);

	Voro2D_Corner res;
	res.x = a->x + t*(b->x-a->x);
	res.y = a->y + t*(b->y-a->y);
	return res;
}
*/




void InactivateLowerAndUpper(Voro2D_Corner *lower, Voro2D_Corner *upper){
	lower = lower->lowerlink;
	while(lower != upper){
		lower->status = Inactive;
		lower = lower->lowerlink;
	}
}


int Voro3D_Trim_Vertex(Voro3D_Vertex *voroVertex, int mp, ptype *maxdist2, int *Old2New){
	int newmp=0;
	int i,j,k;
	ptype mdist2=-1.e20, dist2;
	for(i=0;i<mp;i++){
		if(voroVertex[i].status == Active){
			dist2 = Vec3DInnP(voroVertex+i, voroVertex+i);
			if(dist2>mdist2) mdist2 = dist2;
			Old2New[i] = newmp;
			newmp ++;
		}
	}
	newmp = 0;
	for(i=0;i<mp;i++){
		if(voroVertex[i].status == Active){
			voroVertex[newmp] = voroVertex[i];
			voroVertex[newmp].indx = newmp;
			voroVertex[newmp].link[0] = Old2New[voroVertex[i].link[0] - voroVertex] + voroVertex;
			voroVertex[newmp].link[1] = Old2New[voroVertex[i].link[1] - voroVertex] + voroVertex;
			voroVertex[newmp].link[2] = Old2New[voroVertex[i].link[2] - voroVertex] + voroVertex;
			newmp ++;
		}
	}
	*maxdist2 = mdist2;
	return newmp;
}

int Voro2D_Trim_Corner(Voro2D_Corner *voroCorner, int mp, ptype *maxdist2){
	int Old2New[mp];
	int newmp=0;
	int i,j,k;
	ptype mdist2=-1.e20, dist2;
	for(i=0;i<mp;i++){
		if(voroCorner[i].status == Active){
			dist2 = Vec2DInnP(voroCorner+i, voroCorner+i);
			if(dist2>mdist2) mdist2 = dist2;
			Old2New[i] = newmp;
			newmp ++;
		}
	}
	newmp = 0;
	for(i=0;i<mp;i++){
		if(voroCorner[i].status == Active){
			voroCorner[newmp] = voroCorner[i];
			voroCorner[newmp].upperlink = Old2New[voroCorner[i].upperlink - voroCorner] + voroCorner;
			voroCorner[newmp].lowerlink = Old2New[voroCorner[i].lowerlink - voroCorner] + voroCorner;
			newmp ++;
		}
	}
	*maxdist2 = mdist2;
	return newmp;
}

/*
 *
 *
 *
 * */

/*
Voro3D_point Vec3DSub(Voro3D_point *a, Voro3D_point *b){
	Voro3D_point res;
	res.x = a->x - b->x;
	res.y = a->y - b->y;
	res.z = a->z - b->z;
	return res;
}

Voro2D_point Vec2DSub(Voro2D_point *a, Voro2D_point *b){
	Voro2D_point res;
	res.x = a->x - b->x;
	res.y = a->y - b->y;
	return res;
}
*/

/* Basically we assume r1 is the origin in the coordinate */
Voro2D_Corner Voro2D_FindIntersection(Voro2D_Corner *a, Voro2D_Corner *b, Voro2D_point *neigh, ptype *tt){
	Voro2D_Corner res;
	Voro2D_point ba = Vec2DSub(b,a);
	Voro2D_point *pba = &ba;

	ptype mod = Vec2DInnP(pba, neigh);
	ptype den = 0.5*Vec2DInnP(neigh,neigh) - Vec2DInnP(a, neigh);
	{
		ptype t = den/mod;
		/* to avoid the FP exceptions caused by the round-off errors */
		t = MAX(0., t); 
		t = MIN(1., t); 
		if(t>=0 && t<=1) {
			*tt = t;
			return Vec2DMulAdd(t,pba,a);
		}
		else {
			fprintf(stderr,"Error found in t: out of bound: t=%g\n", t);
			exit(9);
		}
	}
}
Voro3D_Vertex Voro3D_FindIntersection(Voro3D_Vertex *a, Voro3D_Vertex *b, Voro3D_point *neigh){
	Voro3D_point ba = Vec3DSub(b,a);
	Voro3D_point *pba = &ba;
	ptype mod = Vec3DInnP(&ba, neigh);
	ptype den = 0.5*Vec3DInnP(neigh,neigh) - Vec3DInnP(a, neigh);
	if(mod ==0) {
		fprintf(stderr,"Error found : mod= 0\n");
		/*
		exit(9);
		*/
		ptype t = 0.5L;
		return Vec3DMulAdd(t,pba,a);
	}
	else {
		ptype t = den/mod;
		t = MIN(t, 1.L-EPS);
		t = MAX(t, EPS);
		if(t>0 && t<1) {
			return Vec3DMulAdd(t,pba,a);
		}
		else {
			fprintf(stderr,"Error found in t: out of bound: t=%g\n", t);
			exit(9);
		}
	}
}

int Voro2D_CutOrStay(Voro2D_point *a, Voro2D_Corner *neigh){
	ptype dist = Vec2DInnP(neigh,neigh);
	Voro2D_point b = Vec2DSub(neigh,a);
	ptype dist2 = Vec2DInnP(&b,&b);
	if(dist > dist2) return Outside;
	else if (dist == dist2) return Onborder;
	else return Inside;
}
int Voro3D_CutOrStay(Voro3D_point *a, Voro3D_Vertex *neigh, ptype *offset){
	ptype d = Vec3DInnP(neigh,neigh);
	Voro3D_point b = Vec3DSub(neigh,a);
	ptype d2 = Vec3DInnP(&b,&b);
	ptype ratio  = d/d2;
	if(ratio > 1 + EPS2_InOut) {
		*offset = d-d2;
		return Outside;
	}
	else {
		*offset = -1;
		return Inside;
	}
}


int Voro2D_InitializeCorner(Voro2D_Corner *Corner,int np, ptype boxsize){
	int i;
	ptype halfbox = 0.5*boxsize;
	for(i=0;i<4;i++){
		Corner[i].x = boxsize*(((i+1)%4)/2) - halfbox;
		Corner[i].y = boxsize*(i/2) - halfbox;
		Corner[i].status = Active;
		Corner[i].lowerlink = Corner+ (i-1+4)%4;
		Corner[i].upperlink = Corner+ (i+1+4)%4;
		Corner[i].neighlowerlink = NULL;
		Corner[i].neighupperlink = NULL;
		Corner[i].upperrelated = -(i + 1);
		Corner[i].lowerrelated = -( (i+3)%4 + 1) ;
	}

	/*
	Corner[0].x = -halfbox; Corner[0].y = -halfbox;
	Corner[1].x =  halfbox; Corner[1].y = -halfbox;
	Corner[2].x =  halfbox; Corner[2].y =  halfbox;
	Corner[3].x = -halfbox; Corner[3].y =  halfbox;
	
	for(i=0;i<4;i++){
		Corner[i].status = Active;
		Corner[i].lowerlink = Corner+ (i-1+4)%4;
		Corner[i].upperlink = Corner+ (i+1+4)%4;
		Corner[i].related = None;
	}
	*/

	for(i=4;i<np;i++) Corner[i].status = Inactive;
	return 4;
}

int Voro3D_InitializeVertex(Voro3D_Vertex *Vertex,int np, ptype boxsize){
	int i,k;
	ptype halfbox = 0.5*boxsize;
	for(i=0;i<8;i++){
		int j = i%4;
		int k = i/4;
		Vertex[i].x = boxsize*(((j+1)%4)/2) - halfbox;
		Vertex[i].y = boxsize*(j/2) - halfbox;
		Vertex[i].z = boxsize*(i/4) - halfbox;
		Vertex[i].status = Active;
		Vertex[i].link[0] = Vertex + (j-1+4)%4 + k*4;
		if(i<4){
			Vertex[i].related[0] = -(j+1);
			Vertex[i].related[1] = -((j+3)%4+1);
			Vertex[i].related[2] = -5;

			Vertex[i].link[1] = Vertex + (j+1+4)%4 + k*4;
			Vertex[i].link[2] = Vertex + (i + 4)%8;
		}else {
			Vertex[i].related[0] = -(j+1);
			Vertex[i].related[1] = -6;
			Vertex[i].related[2] = -((j+3)%4+1);

			Vertex[i].link[2] = Vertex + (j+1+4)%4 + k*4;
			Vertex[i].link[1] = Vertex + (i + 4)%8;
		}
		for(k=0;k<3;k++) {
//			Vertex[i].related[k] = None;
			Vertex[i].facecount[k] = 0;
		}
	}
	for(i=8;i<np;i++) Vertex[i].status = Inactive;
	return 8;
}


Voro2D_Corner *FindingUpperLimit( Voro2D_point *neigh, Voro2D_Corner *a){
	a = a->lowerlink;
	while(Voro2D_CutOrStay(neigh,a) ==Outside){
		a = a->lowerlink;
	}
	return a;
}
Voro2D_Corner *FindingLowerLimit(Voro2D_point *neigh, Voro2D_Corner *a){
	a = a->upperlink;
	while(Voro2D_CutOrStay(neigh,a) ==Outside){
		a = a->upperlink;
	}
	return a;
}

Voro3D_Vertex *Finding3DBoundaries(Voro3D_point *neigh, Voro3D_Vertex *a){
	Voro3D_Vertex *check = a;
	return check;
}

int Voro3D_CreateOneVertex(Voro3D_Vertex *voroVertex, Voro3D_Vertex *now,Voro3D_Vertex *target, int ip, 
		Voro3D_point *neighbor, int pp){
	int k;
	voroVertex[ip] = Voro3D_FindIntersection(now, target, neighbor+pp);
	voroVertex[ip].link[0] = target;
    voroVertex[ip].link[1] = voroVertex + ip + 1;
    voroVertex[ip].link[2] = voroVertex + ip - 1;
    voroVertex[ip].related[0] = pp;
    voroVertex[ip].status = Active;
	for(k=0;k<3;k++) voroVertex[ip].facecount[k] = 0;
    /* This is to update the disconnected inner vertex to link to the new vertex. */
    for(k=0;k<3;k++) if(target->link[k] == now) { 
		target->link[k] = voroVertex+ip;
		voroVertex[ip].related[2] = target->related[(k+1)%3];
		voroVertex[ip].related[1] = target->related[(k-1+3)%3];
	}
	ip ++;
	return ip;
}

int CreateNewVertices(Voro3D_point *neighbor, int pp, Voro3D_Vertex2 *start, 
		Voro3D_Vertex *voroVertex, int ip){
	ptype offset;
	int ibreak=1;
	int i,j,k;
	int jp = ip; /* backup for the starting ip */
	Voro3D_Vertex *now = start->vertex;
	Voro3D_Vertex *target = start->link;
	int ilink = start->ilink;

	ip = Voro3D_CreateOneVertex(voroVertex, now, target, ip, neighbor, pp);

	ilink ++;

	do {
		for(i=0;i<3;i++){
			j = (ilink+i)%3;
			target = now->link[j];

			if(now == start->vertex && target == start->link) {
				/* to close the new loop of vertices */
				voroVertex[ip-1].link[1] = voroVertex+jp  ;
				voroVertex[jp  ].link[2] = voroVertex+ip-1;
				ibreak = 0;
				break;
			}

			int inoutflag = Voro3D_CutOrStay(neighbor+pp,target,&offset);
			switch(inoutflag){
				case Outside: 
					for(k=0;k<3;k++){ 
						if(now == target->link[k]){ 
							now = target; 
							ilink = (k+1)%3; 
//							target = now->link[ilink];
							break; 
						} 
					} 
					break;
				default:
					ip = Voro3D_CreateOneVertex(voroVertex, now, target, ip, neighbor, pp);
			}
			if(inoutflag == Outside) break;
		}
	} while(ibreak);

	return ip;
}
int CreateNewCorner(Voro2D_Corner *lower, Voro2D_Corner *upper, Voro2D_Corner *voroCorner, int ip, 
		Voro2D_point *neighbor, int pp){
	int ilow,iupper;

	Voro2D_Corner *nextlower,*nextupper;

	nextlower = lower->lowerlink;
	nextupper = upper->upperlink;

	ptype tl, tu;

	voroCorner[ip] = Voro2D_FindIntersection(lower, nextlower, neighbor, &tl);
	voroCorner[ip].upperlink = lower;
	voroCorner[ip].lowerlink = voroCorner + (ip+1);
	voroCorner[ip].neighupperlink = lower->neighlowerlink;
	voroCorner[ip].neighlowerlink = neighbor;
	voroCorner[ip].status = Active;
	voroCorner[ip].lowerrelated = pp;
	voroCorner[ip].upperrelated = nextlower->upperrelated;
	lower->lowerlink = voroCorner+ip;
	ip++;

	voroCorner[ip] = Voro2D_FindIntersection(upper, nextupper, neighbor, &tu);
	voroCorner[ip].upperlink = voroCorner + (ip-1);
	voroCorner[ip].lowerlink = upper;
	voroCorner[ip].neighlowerlink = upper->neighupperlink;
	voroCorner[ip].neighupperlink = neighbor;
	voroCorner[ip].status = Active;
	voroCorner[ip].upperrelated = pp;
	voroCorner[ip].lowerrelated = nextupper->lowerrelated;
	upper->upperlink = voroCorner+ip;
	ip++;


	Voro2D_Corner *old= upper;
	Voro2D_Corner *new= upper->upperlink;
	do {
		ptype dist = sqrt(neighbor->x*neighbor->x + neighbor->y*neighbor->y)*0.5;
		Voro2D_Corner line  = Vec2DSubCorner(old,new);
		ptype linelength = sqrt(Vec2DInnP(&line,&line));
		if(linelength < EPS*dist )
		{
			if(old == lower->lowerlink){
				Voro2D_Corner *upperalive = lower->upperlink;
				old->upperlink = upperalive;
				old->upperrelated = new->upperrelated;
				old->neighlowerlink = new->neighlowerlink;
				upperalive->lowerlink = old;
				new->status = Inactive;
			}
			else {
				Voro2D_Corner *lastalive = old->lowerlink;
				new->lowerlink = lastalive;
				new->lowerrelated = old->lowerrelated;
				new->neighlowerlink = old->neighlowerlink;
				lastalive->upperlink = new;
				old->status = Inactive;
			}
		}
		old = new;
		new = old->upperlink;
	} while(old != lower);

	return ip;
}

int OLD_CreateNewCorner(Voro2D_Corner *lower, Voro2D_Corner *upper, 
		Voro2D_Corner *voroCorner, int ip, Voro2D_point *neighbor, int pp){
	int ilow,iupper;

	Voro2D_Corner *nextlower,*nextupper;

	nextlower = lower->lowerlink;
	nextupper = upper->upperlink;

	ptype tl, tu;
	Voro2D_Corner cl,cu;
	cl = Voro2D_FindIntersection(lower, nextlower, neighbor, &tl);
	cu = Voro2D_FindIntersection(upper, nextupper, neighbor, &tu);

	if(tl !=0 && tu !=0){
		voroCorner[ip] = cl;
		voroCorner[ip].upperlink = lower;
		voroCorner[ip].lowerlink = voroCorner + (ip+1);
		voroCorner[ip].neighupperlink = lower->neighlowerlink;
		voroCorner[ip].neighlowerlink = neighbor;
		voroCorner[ip].status = Active;
		voroCorner[ip].lowerrelated = pp;
		voroCorner[ip].upperrelated = nextlower->upperrelated;
		lower->lowerlink = voroCorner+ip;
		ip++;
		voroCorner[ip] = cu;
		voroCorner[ip].upperlink = voroCorner + (ip-1);
		voroCorner[ip].lowerlink = upper;
		voroCorner[ip].neighlowerlink = upper->neighupperlink;
		voroCorner[ip].neighupperlink = neighbor;
		voroCorner[ip].status = Active;
		voroCorner[ip].upperrelated = pp;
		voroCorner[ip].lowerrelated = nextupper->lowerrelated;
		upper->upperlink = voroCorner+ip;
		
		ip++;
	}
	else if(tl==0 && tu !=0){
		lower->lowerlink = voroCorner + ip;
		lower->neighlowerlink = neighbor;
		lower->lowerrelated = pp;
        upper->upperlink = voroCorner+ip;
		voroCorner[ip] = cu;
        voroCorner[ip].upperlink = lower;
        voroCorner[ip].lowerlink = upper;
        voroCorner[ip].neighlowerlink = upper->neighupperlink;
        voroCorner[ip].neighupperlink = neighbor;
        voroCorner[ip].status = Active;
        voroCorner[ip].upperrelated = pp;
        voroCorner[ip].lowerrelated = nextupper->lowerrelated;

        ip++;
	}
	else if(tl !=0 && tu ==0){
		voroCorner[ip] = cl;
        voroCorner[ip].upperlink = lower;
        voroCorner[ip].lowerlink = upper;
        voroCorner[ip].neighupperlink = lower->neighlowerlink;
        voroCorner[ip].neighlowerlink = neighbor;
        voroCorner[ip].status = Active;
        voroCorner[ip].lowerrelated = pp;
        voroCorner[ip].upperrelated = nextlower->upperrelated;
        lower->lowerlink = voroCorner+ip;
		upper->upperlink = voroCorner + ip;
		upper->neighupperlink = neighbor;
		upper->upperrelated = pp;

        ip++;
	}
	else if(tl ==0 && tu ==0){
		lower->lowerlink = upper;
		lower->neighlowerlink = neighbor;
		lower->lowerrelated = pp;
		upper->upperlink = lower;
		upper->neighupperlink = neighbor;
		upper->upperrelated = pp;
	}

	return ip;
}

Voro3D_point getPolygonCentroid(Voro3D_Vertex *start, Voro3D_Vertex *next){
	Voro3D_point centroid;
	centroid.x = centroid.y = centroid.z = 0;
	Voro3D_Vertex *nnext,*now = start;
	do {
		int i = 0;
		while(next->link[i]!=now) {
            i = (i+1)%3;
        }
        i = (i-1+3)%3;
        nnext = next->link[i];

        Voro3D_Vertex *a,*b,*c;
		a = start;
		b = next;
		c = nnext;
        Voro3D_point r1,r2;
        r1.x = b->x - a->x;
        r1.y = b->y - a->y;
        r1.z = b->z - a->z;
        r2.x = c->x - a->x;
        r2.y = c->y - a->y;
        r2.z = c->z - a->z;

		Voro3D_point dS = Vec3DCross(&r1,&r2);
		Voro3D_point ab = AddPoint3D(a,b);
		Voro3D_point bc = AddPoint3D(b,c);
		Voro3D_point ca = AddPoint3D(c,a);
		centroid.x += dS.x * (ab.x*ab.x + bc.x*bc.x + ca.x*ca.x);
		centroid.y += dS.y * (ab.y*ab.y + bc.y*bc.y + ca.y*ca.y);
		centroid.z += dS.z * (ab.z*ab.z + bc.z*bc.z + ca.z*ca.z);

		next->considered[i] = Yes;

		now = next;
		next = nnext;
	} while(now != start);
	return centroid;
}

Voro3D_point findPolyhedronCentroid(Voro3D_Vertex *vertex, int np){
	Voro3D_point centroid,res;
	int i,j,k;
	centroid.x = centroid.y = centroid.z = 0;
	for(i=0;i<np;i++) for(j=0;j<3;j++) vertex[i].considered[j] = No;
	for(i=0;i<np;i++){
		Voro3D_Vertex *start = vertex+i;
		for(j=0;j<3;j++){
			if(start->considered[j] == No){
				Voro3D_Vertex *next = start->link[j];
				res = getPolygonCentroid(start, next);
				centroid.x += res.x;
				centroid.y += res.y;
				centroid.z += res.z;

				start->considered[j] = Yes;
			}
		}
	}
	centroid.x /= 48.;
	centroid.y /= 48.;
	centroid.z /= 48.;
	return centroid;
}


typedef struct delaunay2D {
	Voro2D_point *delaunay[3];
} delaunay2D;

void findNewCorner(Voro2D_Corner *old, Voro2D_point *neighbors, Voro2D_point *ctrcenter){
	Voro2D_point center = *ctrcenter;
	center.x = center.y = 0;
	center.csound = ctrcenter->csound;
	Voro2D_point new;
	int imax=0;
	ptype maxcsound = center.csound;
	delaunay2D triangle;
	if(old->upperrelated <0 || old->lowerrelated <0) return;
	triangle.delaunay[0] = &center;
	triangle.delaunay[1] = neighbors + old->upperrelated;
	triangle.delaunay[2] = neighbors + old->lowerrelated;
	{
		ptype csound;
		csound = triangle.delaunay[1]->csound;
		if(csound > maxcsound) {
			maxcsound = csound;
			imax = 1;
		}
		csound = triangle.delaunay[2]->csound;
		if(csound > maxcsound) {
			maxcsound = csound;
			imax = 2;
		}
	}
	int j;
	j = (imax+1)%3;
	Voro2D_point r1, r2;
	ptype d1, d2, t;

	ptype wi, wmax;

	wmax = pow(triangle.delaunay[imax]->csound,0.25L);
	wi = pow(triangle.delaunay[j]->csound,0.25L);
	t = wmax/(wi+wmax);
	t = MIN(t, 0.7);
	findLine2D(triangle.delaunay[imax], triangle.delaunay[j], t,  &r1, &d1);
	j = (imax+2)%3;
	wi = pow(triangle.delaunay[j]->csound,0.25L);
	t = wmax/(wi+wmax);
	t = MIN(t, 0.7);
	findLine2D(triangle.delaunay[imax], triangle.delaunay[j], t, &r2, &d2);

	ptype slope1, slope2;
	slope1 = -r1.x/r1.y;
	slope2 = -r2.x/r2.y;
	d1 = d1 / r1.y;
	d2 = d2 / r2.y;
	new = calculateIntersection2D(slope1, d1, slope2, d2);
	old->x = new.x;
	old->y = new.y;
}

/* finding normal intersection * input: r1, r2 (two vectors) & output: r (normal vector), d(residual) */
void findLine2D(Voro2D_point *r1, Voro2D_point *r2, ptype t, Voro2D_point *r, ptype *d){
	Voro2D_point a = Vec2DSub(r2, r1);
	r->x = a.x; r->y = a.y;
	*d = dotProduct2D(r, r1) + t * dotProduct2D(r,r);
}
Voro2D_point calculateIntersection2D(ptype slope1, ptype c1, ptype slope2, ptype c2) {
    Voro2D_point intersection;
    intersection.x = (c2 - c1) / (slope1 - slope2);
    intersection.y = slope1 * intersection.x + c1;
    return intersection;
}

typedef struct delaunay3D {
    Voro2D_point *delaunay[4];
} delaunay3D;

void findNewVertex(Voro3D_Vertex *old, Voro3D_point *neighbors, Voro3D_point *ctrcenter){
    Voro3D_point center = *ctrcenter;
    center.x = center.y = center.z = 0;
    center.csound = ctrcenter->csound;
    Voro3D_point new;
    int imax=0;
    ptype maxcsound = center.csound;
    delaunay3D triangle;
    if(old->related[0] <0 || old->related[1] <0 || old->related[2]<0) return;
    triangle.delaunay[0] = &center;
    triangle.delaunay[1] = neighbors + old->related[0];
    triangle.delaunay[2] = neighbors + old->related[1];
    triangle.delaunay[3] = neighbors + old->related[2];
    int i;
	for(i=0;i<4;i++){
        ptype csound;
        csound = triangle.delaunay[i]->csound;
        if(csound > maxcsound) {
            maxcsound = csound;
            imax = 1;
        }
    }
    int j;
    Voro3D_point r1, r2, r3;
    ptype d1, d2, d3, t;

    ptype wi, wmax;
    wmax = pow(triangle.delaunay[imax]->csound,1./4.);

    j = (imax+1)%4;
    wi = pow(triangle.delaunay[j]->csound,1./4.);
    t = wmax/(wi+wmax);
    t = MIN(t, 0.7);
    findPlane3D(triangle.delaunay[imax], triangle.delaunay[j], t,  &r1, &d1);

    j = (imax+2)%4;
    wi = pow(triangle.delaunay[j]->csound,1./4.);
    t = wmax/(wi+wmax);
    t = MIN(t, 0.7);
    findPlane3D(triangle.delaunay[imax], triangle.delaunay[j], t, &r2, &d2);

    j = (imax+3)%4;
    wi = pow(triangle.delaunay[j]->csound,1./4.);
    t = wmax/(wi+wmax);
    t = MIN(t, 0.7);
    findPlane3D(triangle.delaunay[imax], triangle.delaunay[j], t, &r3, &d3);

    new = calculateIntersection3D(&r1,d1, &r2, d2, &r3, d3);
    old->x = new.x;
    old->y = new.y;
    old->z = new.z;
}

/* finding normal intersection * input: r1, r2 (two vectors) & output: r (normal vector), d(residual) */
void findPlane3D(Voro3D_point *r1, Voro3D_point *r2, ptype t, Voro3D_point *r, ptype *d){
	Voro3D_point a = Vec3DSub(r2, r1);
	r->x = a.x; r->y = a.y; r->z = a.z;
	*d = dotProduct3D(r, r1) + t * dotProduct3D(r,r);
}

Voro3D_point calculateIntersection3D(
		Voro3D_point *r1, ptype d1, 
		Voro3D_point *r2, ptype d2,
		Voro3D_point *r3, ptype d3
		) { 
	Voro3D_point intersection; 
	ptype determinant = r1->x * (r2->y * r3->z - r3->y * r2->z) 
		- r1->y * (r2->x * r3->z - r3->x * r2->z) 
		+ r2->z * (r2->x * r3->y - r3->x * r2->y);
    if (determinant == 0) {
        printf("Planes are parallel, no intersection point.\n");
        intersection.x = NAN;
        intersection.y = NAN;
        intersection.z = NAN;
    } else {
        intersection.x = (
				d1 * (r2->y * r3->z - r3->y * r2->z) 
				- r1->y * (d2 * r3->z - d3 * r2->z) 
				+ r2->z * (d2 * r3->y - d3 * r2->y)
				) / determinant;
        intersection.y = (
				r1->x * (d2 * r3->z - d3 * r2->z) 
				- d1 * (r2->x * r3->z - r3->x * r2->z) 
				+ r2->z * (r2->x * d3 - r3->x * d2)
				) / determinant;
        intersection.z = (
				r1->x * (r2->y * d3 - r3->y * d2) 
				- r1->y * (r2->x * d3 - r3->x * d2) 
				+ d1 * (r2->x * r3->y - r3->x * r2->y)
				) / determinant;
    }

    return intersection;
}

Voro2D_Corner Centroid2DPolygon( Voro2D_Corner *c, int np){
	Voro2D_Corner res;
	ptype A=0;
	res.x = 0;
	res.y = 0;
	Voro2D_Corner *tmp = c;
	Voro2D_Corner *tmp2;

	do {
		tmp2 = tmp->upperlink;
		A += tmp->x*tmp->y - tmp2->x*tmp2->y;
		res.x += (tmp->x+tmp2->x)*(tmp->x*tmp2->y - tmp2->y*tmp->y);
		res.y += (tmp->y+tmp2->y)*(tmp->x*tmp2->y - tmp2->y*tmp->y);
	} while(tmp != c);
	A = A*0.5;
	res.x = res.x/(6*A);
	res.y = res.y/(6*A);
	return res;
}
ptype Voro3D_pyramid_volume(Voro3D_Vertex *start, Voro3D_Vertex *next){
	int i;
	Voro3D_point dS;
	ptype volume=0;
	Voro3D_Vertex *nnext,*now;
	now = start;

	do {
		i=0;
		while(next->link[i]!=now) {
			i = (i+1)%3;
		}
		i = (i-1+3)%3;
		nnext = next->link[i];

		Voro3D_point r1,r2;
		r1.x = next->x-start->x; 
		r1.y = next->y-start->y; 
		r1.z = next->z-start->z;
		r2.x = nnext->x - next->x;
		r2.y = nnext->y - next->y;
		r2.z = nnext->z - next->z;
		dS= Vec3DCross(&r1,&r2);
		volume += start->x*dS.x;
		volume += start->y*dS.y;
		volume += start->z*dS.z;
		next->considered[i] = Yes;

		now = next;
		next = nnext;
	} while(now != start);
	volume = (1./6.)*volume;
	return volume;
}
ptype Voro3D_Volume_Polyhedron(Voro3D_Vertex *vorovertex,int nvertex){
	int i,j;
	Voro3D_Vertex *start;
	ptype volume = 0;
	int istart = 0;
	int inow = istart;
	for(i=0;i<nvertex;i++) 
		vorovertex[i].considered[0] = vorovertex[i].considered[1] = vorovertex[i].considered[2] = No;
	for(i=0;i<nvertex;i++){
		start = vorovertex+i;
		for(j=0;j<3;j++){
			if(start->considered[j] == No){
				Voro3D_Vertex *next = start->link[j];
				volume += Voro3D_pyramid_volume(start, next);
				start->considered[j] = Yes;
			}
		}
	}
	return volume;
}
Voro3D_point Voro3D_norm_polygon(Voro3D_Vertex *start, Voro3D_Vertex *next){
    int i;
    Voro3D_point norm,dnorm;
    ptype volume=0;
    Voro3D_Vertex *nnext,*now;
	now = start;

	norm.x = norm.y = norm.z = 0;

    do {
        i=0;
        while(next->link[i]!=now) {
            i = (i+1)%3;
        }
        i = (i-1+3)%3;
        nnext = next->link[i];


        Voro3D_point r1,r2;
        r1.x = next->x-start->x;
        r1.y = next->y-start->y;
        r1.z = next->z-start->z;
        r2.x = nnext->x - next->x;
        r2.y = nnext->y - next->y;
        r2.z = nnext->z - next->z;
        dnorm= Vec3DCross(&r1,&r2);
        norm.x += dnorm.x*0.5L;
        norm.y += dnorm.y*0.5L;
        norm.z += dnorm.z*0.5L;
        next->considered[i] = Yes;
        now = next;
        next = nnext;
    } while(now != start);
    return norm;
}


Voro3D_point Voro3D_denGrad_Polyhedron(Voro3D_Vertex *vorovertex,int nvertex, 
		int indx, Voro3D_HessianParticle *bp, Voro3D_point *neighwork){
	Voro3D_point res;
	res.x = res.y = res.z = 0;
    int i,j;
    Voro3D_Vertex *start;
    ptype volume = 0;
    int istart = 0;
    int inow = istart;
    for(i=0;i<nvertex;i++)
        vorovertex[i].considered[0] = vorovertex[i].considered[1] = vorovertex[i].considered[2] = No;
    for(i=0;i<nvertex;i++){
        start = vorovertex+i;
        for(j=0;j<3;j++){
            if(start->considered[j] == No){
				int iface = (start->related[j]+2)%3;
				int jndx = neighwork[iface].indx;
                Voro3D_Vertex *link = start->link[j];
                Voro3D_point norm = Voro3D_norm_polygon(start, link);
				ptype avgDen;
				if(iface>=0){
					avgDen = 0.5*(bp[indx].den + bp[jndx].den);
				}
				else {
					avgDen = bp[indx].den;
				}
				res.x += norm.x * avgDen;
				res.y += norm.y * avgDen;
				res.z += norm.z * avgDen;
                start->considered[j] = Yes;
            }
        }
    }
	res.x = res.x/bp[indx].volume;
	res.y = res.y/bp[indx].volume;
	res.z = res.z/bp[indx].volume;
    return res;
}
Voro3D_tensor Voro3D_denHessian_Polyhedron(Voro3D_Vertex *vorovertex,int nvertex, 
		int indx, Voro3D_HessianParticle *bp, Voro3D_point *neighwork){
	Voro3D_tensor res;
	res.xx = res.xy = res.xz = res.yx = res.yy = res.yz = res.zx = res.zy = res.zz = 0;
    int i,j;
    Voro3D_Vertex *start;
    ptype volume = 0;
    int istart = 0;
    int inow = istart;
    for(i=0;i<nvertex;i++)
        vorovertex[i].considered[0] = vorovertex[i].considered[1] = vorovertex[i].considered[2] = No;
    for(i=0;i<nvertex;i++){
        start = vorovertex+i;
        for(j=0;j<3;j++){
            if(start->considered[j] == No){
				int iface = (start->related[j]+2)%3;
				int jndx = neighwork[iface].indx;
                Voro3D_Vertex *link = start->link[j];
                Voro3D_point norm = Voro3D_norm_polygon(start,  link);
				Voro3D_point avgDen;
				if(iface>=0){
					avgDen.x = 0.5*(bp[indx].denGrad.x + bp[jndx].denGrad.x);
					avgDen.y = 0.5*(bp[indx].denGrad.y + bp[jndx].denGrad.y);
					avgDen.z = 0.5*(bp[indx].denGrad.z + bp[jndx].denGrad.z);
				}
				else {
					avgDen.x = bp[indx].denGrad.x;
					avgDen.y = bp[indx].denGrad.y;
					avgDen.z = bp[indx].denGrad.z;
				}
				res.xx += norm.x * avgDen.x;
				res.xy += norm.x * avgDen.y;
				res.xz += norm.x * avgDen.z;
				res.yy += norm.y * avgDen.y;
				res.yz += norm.y * avgDen.z;
				res.zz += norm.z * avgDen.z;
                start->considered[j] = Yes;
            }
        }
    }
	res.xx = res.xx/bp[indx].volume;
	res.xy = res.xy/bp[indx].volume;
	res.xz = res.xz/bp[indx].volume;
	res.yy = res.yy/bp[indx].volume;
	res.yz = res.yz/bp[indx].volume;
	res.zz = res.zz/bp[indx].volume;
    return res;
}


ptype Area2DPolygon(Voro2D_Corner *c, int np){
	int i;
	ptype parea = 0;
	ptype marea = 0;
	ptype area = 0;

	Voro2D_Corner *tmp = c;
	do {
		parea += tmp->x * (tmp->upperlink)->y;
		marea += tmp->y * (tmp->upperlink)->x;
		tmp = tmp->upperlink;
	} while(tmp != c);

	area = 0.5*fabs(parea-marea);
	return area;
}

int sort2Ddist(const void *a, const void *b){
	Voro2D_point *aa,*bb;
	aa = (Voro2D_point*)a;
	bb = (Voro2D_point*)b;
	if(aa->dist2>bb->dist2) return 1;
	else if(aa->dist2<bb->dist2) return -1;
	else return 0;
}

int Voro2D_FindVC(Voro2D_point *ctrCenter, Voro2D_point *inputneighbors, Voro2D_point *neighbors,int np, 
		Voro2D_Corner *voroCorner, int mp, ptype boxsize){
	ptype maxdist2=boxsize/2*sqrt(3.);
	int i,j;
	int ip;

	j =0;

	for(i=0;i<np;i++){
		if(inputneighbors[i].indx != ctrCenter->indx){ 
			neighbors[j] = Vec2DSub(inputneighbors+i , ctrCenter);
			neighbors[j].csound = inputneighbors[i].csound;
			neighbors[j].pressure = inputneighbors[i].pressure;
			neighbors[j].dist2 = (Vec2DInnP(neighbors+j, neighbors+j));
			neighbors[j].bp = inputneighbors[i].bp;
			j++;
		}
	}
	np = j;

	qsort(neighbors, np, sizeof(Voro2D_point), sort2Ddist);

	ip = Voro2D_InitializeCorner(voroCorner,mp,boxsize);
	for(i=0;i<np;i++){
		if(neighbors[i].dist2>4*maxdist2) break;
		j = 0;
		while(j<ip){
			if(voroCorner[j].status==Active && Voro2D_CutOrStay(neighbors+i,voroCorner+j) == Outside){
				/* Warning!!! Here the lower and upper are the boundary of the activated corners */
				Voro2D_Corner *lower = voroCorner+j;
				lower = FindingLowerLimit(neighbors+i, lower);
				Voro2D_Corner *upper = voroCorner+j;
				upper = FindingUpperLimit(neighbors+i,upper);
				InactivateLowerAndUpper(lower,upper); 
				ip = CreateNewCorner(lower,upper, voroCorner,ip,neighbors+i,i);
				ip = Voro2D_Trim_Corner(voroCorner,ip, &maxdist2);
				j = ip;
				break;
			}
			else { 
				j++;
			}
		}
	}
	return ip;
}

int InactivateOutsideVertex(Voro3D_point *neighbors, Voro3D_Vertex *voroVertex, int np){
	int i,icount=0;

	ptype offset;
	for(i=0;i<np;i++){
		if(Voro3D_CutOrStay(neighbors,voroVertex+i,&offset) == Outside){
			voroVertex[i].status = Inactive;
			icount ++;
		}
	}
	return icount;
}

Voro3D_Vertex2 FindNearestBoundaryVertex(Voro3D_point *neighbor, Voro3D_Vertex *vertex){
	Voro3D_Vertex2 res;
	while(1){
		ptype minoffset = 1.e27;
		ptype offset;
		int i,j;
		for(i=0;i<3;i++){
			Voro3D_Vertex *link;
			link = vertex->link[i];
			if(Voro3D_CutOrStay(neighbor, link, &offset) != Outside){
				res.vertex = vertex;
				res.link = link;
				res.ilink = i;
				return res;
			}
			if(offset < minoffset) {
				minoffset = offset;
				j = i;
			}
		}
		vertex = vertex->link[j];
	}

}

int sort3Ddist(const void *a, const void *b){
	Voro3D_point *aa,*bb;
	aa = (Voro3D_point*)a;
	bb = (Voro3D_point*)b;
	if(aa->dist2>bb->dist2) return 1;
	else if(aa->dist2<bb->dist2) return -1;
	else return 0;
}

int Voro3D_FindVC(Voro3D_point *ctrCenter, Voro3D_point *inputneighbors, Voro3D_point *neighbors, int np, 
		Voro3D_Vertex *voroVertex, int mp, ptype boxsize, int ishrink, int *intwork){
	ptype maxdist2=(boxsize/2)*(boxsize/2)*3;
	int i,j;
	int ip;

	j =0;

	for(i=0;i<np;i++){
		if(inputneighbors[i].indx != ctrCenter->indx){ 
			neighbors[j] = Vec3DSub(inputneighbors+i,ctrCenter);
			neighbors[j].csound = inputneighbors[i].csound;
			neighbors[j].pressure = inputneighbors[i].pressure;
			neighbors[j].dist2 = (Vec3DInnP(neighbors+j,neighbors+j));
			j++;
		}
	}
	np = j;
	qsort(neighbors, np, sizeof(Voro3D_point), sort3Ddist);

	ip = Voro3D_InitializeVertex(voroVertex,mp,boxsize);


	for(i=0;i<np;i++){
		if(neighbors[i].dist2>4*maxdist2) break;
		j = 0;
		while(j<ip){
			ptype offset;
			if(voroVertex[j].status==Active && Voro3D_CutOrStay(neighbors+i,voroVertex+j,&offset) == Outside){
				Voro3D_Vertex2 start;
				start = FindNearestBoundaryVertex(neighbors+i, voroVertex+j);
				int icount = InactivateOutsideVertex(neighbors+i, voroVertex,ip);
				ip = CreateNewVertices(neighbors,i,&start,voroVertex, ip);
				ip = Voro3D_Trim_Vertex(voroVertex,ip, &maxdist2, intwork);
				j = ip;
				break;
			}
			else { 
				j++;
			}
		}
	}
	/*
	if(ishrink){
		ptype ww = 1./ctrCenter->csound;
		for(i=0;i<ip;i++){
			ptype xc,yc,zc;
			xc = yc = zc = 0;
			ptype wi,wj,wk;
			wi = sqrt(1./neighbors[voroVertex[i].related[0]].csound);
			wj = sqrt(1./neighbors[voroVertex[i].related[1]].csound);
			wk = sqrt(1./neighbors[voroVertex[i].related[2]].csound);
			ww += wi + wj + wk;
			xc += wi * neighbors[voroVertex[i].related[0]].x;
			xc += wj * neighbors[voroVertex[i].related[1]].x;
			xc += wk * neighbors[voroVertex[i].related[2]].x;
			yc += wi * neighbors[voroVertex[i].related[0]].y;
			yc += wj * neighbors[voroVertex[i].related[1]].y;
			yc += wk * neighbors[voroVertex[i].related[2]].y;
			zc += wi * neighbors[voroVertex[i].related[0]].z;
			zc += wj * neighbors[voroVertex[i].related[1]].z;
			zc += wk * neighbors[voroVertex[i].related[2]].z;
			voroVertex[i].x = xc/ww*4.L;
			voroVertex[i].y = yc/ww*4.L;
			voroVertex[i].z = zc/ww*4.L;
		}
	}
	*/


	return ip;
}

int Voro3D_FaceExtract(Voro3D_Vertex *Vertex, int mp,int fid, Voro3D_point *outcorner){
	int i,j,k;
	int np = 0;
	for(i=0;i<mp;i++){
		for(j=0;j<3;j++){
			if(Vertex[i].related[j] == fid){
				Voro3D_Vertex *onow = Vertex+i;
				Voro3D_Vertex *now = onow;
				k = (j+1)%3;
				Voro3D_Vertex *otarget = (Vertex+i)->link[k];
				Voro3D_Vertex *target = otarget;
				DumpTo3DPoint(outcorner+np, now); np++;
				do{
					for(k=0;k<3;k++){
						if(target->related[k] == fid){
							now = target;
							target = target->link[ (k+1)%3];
							break;
						}
					}
					DumpTo3DPoint(outcorner+np, now); np++;
				} while(now != onow || target != otarget);
				return np;
			}
		}
	}
	return np;
}
int Voro2D_LineExtract(Voro2D_Corner *corner, int mp,Voro2D_point *outpolygon){
	int i,j,k;
	int np = 0;
	Voro2D_Corner *ptr;
	ptr = corner;
	for(i=0;i<mp;i++){
		outpolygon[i].x = ptr->x;
		outpolygon[i].y = ptr->y;
		outpolygon[i].indx = ptr->upperrelated;
		ptr = ptr->upperlink;
	}
	return np;
}


ptype getWeight3D(ptype x, ptype y, ptype z, ptype hsml){
	ptype r = sqrt(x*x+y*y+z*z);
	ptype rh = r/hsml;
	if(rh<1.){
	}
	else if(rh<2.){
	}
	else return 0.;
	return 0;
}

ptype findSDen3D(Voro3D_Vertex *vertex, int ip, Voro3D_point *neighwork,
		int *intwork, int mp, Voro3D_GasParticle *bp){
	int i,j;
	for(i=0;i<mp;i++) intwork[i] = 0;
	for(i=0;i<ip;i++){
		intwork[vertex[i].related[0]] = 1;
		intwork[vertex[i].related[1]] = 1;
		intwork[vertex[i].related[2]] = 1;
	}
	ptype tweight = 0;
	ptype hsml;
	ptype den = getWeight3D(0.,0.,0., hsml);
	for(i=0;i<mp;i++){
		if(intwork[i] ==1){
			int jndx = neighwork[i].indx;
			ptype weight = getWeight3D(neighwork[i].x, neighwork[i].y,neighwork[i].z, hsml);
			den += bp[jndx].mass/bp[jndx].volume * weight;
			tweight += weight;
		}
	}
	den = den/tweight;
	return den;
}

ptype getAvgPressure2D(int indx, Voro2D_Corner *voro, Voro2dGasParticle *bp, 
		Voro2D_point *neigh){
	int jndx = neigh[voro->upperrelated].indx;
	int j1id = neigh[voro->upperlink->upperrelated].indx;
	int j2id = neigh[voro->lowerrelated].indx;
	ptype res = (5./12.)*(bp[jndx].pressure + bp[indx].pressure)
		+ (1./12.)*(bp[j1id].pressure + bp[j2id].pressure);
	return res;
}
ptype getAvgPressureOnSurface3D(ptype Pi, ptype Pj, Voro3D_point *neigh, Voro3D_Vertex *start, Voro3D_Vertex *next){
	int i;
    ptype totSij=0;
    ptype totPij=0;
    Voro3D_Vertex *now;
    now = start;

	Voro3D_point com = getPolygonCentroid(start, next);
	Voro3D_point *centroid = &com;
	ptype Pij = 0.5*(Pi+Pj);

	do {
//        Voro3D_Vertex *nnext = getNextPolygonVertex(now,next); 
		i=0;
        while(next->link[i]!=now) {
            i = (i+1)%3;
        }
        i = (i-1+3)%3;
        Voro3D_Vertex *nnext = next->link[i];

		Voro3D_point r1,r2;
        r1.x = now->x-centroid->x;
        r1.y = now->y-centroid->y;
        r1.z = now->z-centroid->z;
        r2.x = next->x - centroid->x;
        r2.y = next->y - centroid->y;
        r2.z = next->z - centroid->z;
		
		Voro3D_point dnorm= Vec3DCross(&r1,&r2);
		ptype dSij = 0.5*sqrt(dnorm.x*dnorm.x + dnorm.y*dnorm.y + dnorm.z*dnorm.z);
		ptype Pj1 = 0.25L*(Pi+neigh[now->related[0]].pressure + neigh[now->related[1]].pressure + neigh[now->related[2]].pressure);
		ptype Pj2 = 0.25L*(Pi+neigh[next->related[0]].pressure + neigh[next->related[1]].pressure + neigh[next->related[2]].pressure);
		ptype Pj12 = 0.3333333L*(Pij + Pj1 + Pj2);
		totSij += dSij;
		totPij += Pj12*dSij;
		


		next->considered[i] = Yes;
		now = next; next = nnext;

	} while(now != start);
	/* (1/3, 1/6, 1/12) */
	ptype avPressure = totPij/totSij;

	return avPressure;
}
