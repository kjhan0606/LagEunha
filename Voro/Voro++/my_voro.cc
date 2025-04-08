// Author   : Juhan Kim at KIAS
// Email    : kjhan0606@gmail.com
// Date     : Febrary 27th 2019

#include "voro++.hh"
using namespace voro;

using std::vector;

extern "C" int myvoro (int np, int *id,double *x, double *y, double *z, 
	int *c_iface, double *c_area, double *c_norm) 
{ 
	int j;
	voronoicell_neighbor v; 
	
	v.init(-1,1,-1,1,-1,1); 
	for(j=0;j<np;j++) { 
		v.nplane(x[j],y[j],z[j],id[j]); 
	} 

	vector<int> iface;
	vector<double> area;
	vector<double> norm;

	v.neighbors(iface);
	v.face_areas(area);
	v.normals(norm);
	unsigned int i;
	for(i=0;i<iface.size();i++){
		c_iface[i] = iface[i];
		c_area[i] = area[i];
		c_norm[3*i] = norm[3*i];
		c_norm[3*i+1] = norm[3*i+1];
		c_norm[3*i+2] = norm[3*i+2];
	}
	/*
    v.draw_pov(0,0,0,"single_cell_v.pov");
    v.draw_pov_mesh(0,0,0,"single_cell_m.pov");
	*/

	return iface.size();
}

extern "C" int view_myvoro (int np, int *id,double *x, double *y, double *z, 
	int *c_iface, double *c_area, double *c_norm) 
{ 
	int j;
	voronoicell_neighbor v; 
	
	v.init(-1,1,-1,1,-1,1); 
	for(j=0;j<np;j++) { 
		v.nplane(x[j],y[j],z[j],id[j]); 
	} 

    v.draw_pov(0,0,0,"single_cell_v.pov");
    v.draw_pov_mesh(0,0,0,"single_cell_m.pov");

	return 0;
}
