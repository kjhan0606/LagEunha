#ifndef FV_ORIGINAL_CASES_H
#define FV_ORIGINAL_CASES_H
/* Original nine 1D problem definitions. Open boundaries use an enlarged
 * moving material interval; accuracy is evaluated on the physical interval.
 * This fixed-mass 1D driver has no global weight solve. */
typedef struct {
    double lo,hi,cut[2],gamma,end,rho[3],u[3],p[3];
    int regions,lower,upper,sinusoidal;
} FvOriginalCase;
static int original_case(int problem,FvOriginalCase *s){
    *s=(FvOriginalCase){0,1,{.5,1},1.4,0,{1,1,1},{0,0,0},{1,1,1},2,GFS_FV_MOVING_TRANSMIT,GFS_FV_MOVING_TRANSMIT,0};
    switch(problem){
    case BLAST:s->regions=3;s->cut[0]=.1;s->cut[1]=.9;s->p[0]=1000;s->p[1]=.01;s->p[2]=100;s->lower=s->upper=GFS_FV_REFLECT;s->end=.038;break;
    case COLLISION:s->cut[0]=.4;s->rho[0]=5.99924;s->rho[1]=5.99242;s->u[0]=19.5975;s->u[1]=-6.19633;s->p[0]=460.894;s->p[1]=46.095;s->end=.035;break;
    case SHU_OSHER:s->lo=-5;s->hi=5;s->cut[0]=-4;s->rho[0]=3.857143;s->u[0]=2.629369;s->p[0]=10.33333;s->end=1.8;s->sinusoidal=1;break;
    case NOH:s->gamma=5./3;s->regions=1;s->u[0]=-1;s->p[0]=1e-6;s->lower=GFS_FV_REFLECT;s->end=.6;break;
    case LAX:s->rho[0]=.445;s->rho[1]=.5;s->u[0]=.698;s->p[0]=3.528;s->p[1]=.571;s->end=.13;break;
    case DOUBLE_FAN:s->u[0]=-2;s->u[1]=2;s->p[0]=s->p[1]=.4;s->end=.15;break;
    case CONTACT_OPEN:s->rho[1]=.125;s->end=.25;break;
    case MOVING_CONTACT_OPEN:s->rho[1]=.125;s->u[0]=s->u[1]=1;s->end=.25;break;
    default:return 0;
    }return 1;
}

#endif
