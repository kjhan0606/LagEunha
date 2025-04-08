#define OBSERVER_ES {\
	void LightConeES(treepmparticletype *, long, int, float, float, float, int);\
	if(nowTsubdiv == 0){\
		LightConeES(simpar.dm.u.tbp,simpar.dm.np,simpar.indtinfo.tnumcount,simpar.anow,  pre_da, damin,1);\
	}\
	else if(nowTsubdiv != (1<<maxTsubpower) ){\
		LightConeES(simpar.dm.u.tbp,simpar.dm.np,simpar.indtinfo.tnumcount,simpar.anow,  pre_da, damin,0);\
	}\
}
#define OBSERVER_S {\
	void LightConeS(treepmparticletype *, long, int, float, float, float, float );\
	if(nowTsubdiv != (1<<maxTsubpower) && nowTsubdiv != 0){\
		LightConeS(simpar.dm.u.tbp,simpar.dm.np,simpar.indtinfo.tnumcount,simpar.amax,simpar.anow,pre_da,damin);\
	}\
}
