#ifndef MATERIAL_CURVED_BOUNDS_H
#define MATERIAL_CURVED_BOUNDS_H
#include <limits>
/* Outward-rounded long-double intervals for the analytic P3 map represented
 * by the stored double coefficients. Input coefficients are not re-fitted.
 * sqrt endpoint rounding is widened; no fast-math is permitted. Native
 * double evaluation roundoff is a separate explicit envelope, not a claim
 * that a rounded polynomial is the exact underlying geometric oracle.
 */
namespace curved_bounds {
using R=long double;
static R down(R x){return nextafterl(x,-std::numeric_limits<R>::infinity());}
static R up(R x){return nextafterl(x,std::numeric_limits<R>::infinity());}
struct I {R lo,hi;I(R x=0):lo(x),hi(x){}I(R a,R b):lo(a),hi(b){need(a<=b,"ordered interval");}};
static bool exact_zero(I a){return a.lo==0 && a.hi==0;}
static I operator+(I a,I b){if(exact_zero(a))return b;if(exact_zero(b))return a;if(a.lo==a.hi && b.lo==b.hi && a.lo==-b.lo)return I(0);return {down(a.lo+b.lo),up(a.hi+b.hi)};}
static I operator-(I a){return {-a.hi,-a.lo};}
static I operator-(I a,I b){return a+(-b);}
static I operator*(I a,I b){if(exact_zero(a) || exact_zero(b))return I(0);R v[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};return {down(*std::min_element(v,v+4)),up(*std::max_element(v,v+4))};}
static I inverse(I a){need(a.lo>0 || a.hi<0,"interval denominator separated from zero");return {down(1/a.hi),up(1/a.lo)};}
static I operator/(I a,I b){return a*inverse(b);}
static I square(I a){if(exact_zero(a))return I(0);if(a.lo<=0 && a.hi>=0)return {0,up(std::max(a.lo*a.lo,a.hi*a.hi))};R p=a.lo*a.lo,q=a.hi*a.hi;return {down(std::min(p,q)),up(std::max(p,q))};}
static I root(I a){need(a.lo>=0,"nonnegative sqrt interval");if(exact_zero(a))return I(0);R l=sqrtl(a.lo),h=sqrtl(a.hi);for(int i=0;i<4;i++){l=std::max(0.L,down(l));h=up(h);}return {l,h};}
static R absmax(I a){return std::max(fabsl(a.lo),fabsl(a.hi));}
static I widen(I a,R e){return {down(a.lo-e),up(a.hi+e)};}
struct D {I v,x,y;D(I a=I(0),I b=I(0),I c=I(0)):v(a),x(b),y(c){}};
static D operator+(D a,D b){return {a.v+b.v,a.x+b.x,a.y+b.y};}
static D operator-(D a,D b){return {a.v-b.v,a.x-b.x,a.y-b.y};}
static D operator-(D a){return D()-a;}
static D operator*(D a,D b){return {a.v*b.v,a.x*b.v+a.v*b.x,a.y*b.v+a.v*b.y};}
static D operator/(D a,D b){I inv=inverse(b.v);return {a.v*inv,(a.x-a.v*inv*b.x)*inv,(a.y-a.v*inv*b.y)*inv};}
static D sqrt_dual(D a){I s=root(a.v),factor=inverse(I(2)*s);return {s,a.x*factor,a.y*factor};}
using Tri=std::array<std::array<R,2>,3>;
#define MATERIAL_CURVED_VELOCITY_GRADIENT_BOUNDS
struct Enclosure {I displacement[2],velocity[2],gradient[4],velocity_reference_gradient[4],J,s;R lipschitz,native_position_envelope,native_gradient_envelope,native_velocity_gradient_envelope;};
static Enclosure enclose(const curved::Record &r,const Tri &tri,R eps){
    I range[2];for(int k=0;k<2;k++){R l=tri[0][k],h=l;for(int a=1;a<3;a++){l=std::min(l,tri[a][k]);h=std::max(h,tri[a][k]);}range[k]=I(l,h);}
    D X(range[0],I(1),I(0)),Y(range[1],I(0),I(1)),dx=X-D(I(r.triangle.geometry.x[0][0])),dy=Y-D(I(r.triangle.geometry.x[0][1]));
    D bary[3];bary[1]=dx*D(I(r.triangle.geometry.grad[1][0]))+dy*D(I(r.triangle.geometry.grad[1][1]));bary[2]=dx*D(I(r.triangle.geometry.grad[2][0]))+dy*D(I(r.triangle.geometry.grad[2][1]));bary[0]=D(I(1))-bary[1]-bary[2];
    D N[10];for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=D(I(1));for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=N[a]*(bary[k]*D(I(3))-D(I(j)))/D(I(index[k]-j));}
    I ax(r.triangle.geometry.x[0][0]),ay(r.triangle.geometry.x[0][1]);I za=square(ax)+square(ay)+square(I(eps));D L=(dx*D(I(2)*ax)+dy*D(I(2)*ay))/D(za),Q=(dx*dx+dy*dy)/D(za);
    D A3=D(I(1))-D(I(.75L))*L-D(I(.75L))*Q+D(I(.65625L))*L*L+D(I(1.3125L))*L*Q-D(I(.6015625L))*L*L*L;
    D A2=D(I(1))-D(I(.75L))*L-D(I(.75L))*Q+D(I(.65625L))*L*L;
    I rz=root(root(za)),oa=inverse(rz*rz*rz);D z=X*X+Y*Y+D(square(I(eps)));
    /* X*X interval alone may include a tiny negative lower endpoint at zero.
     * Supply the exact nonnegative square enclosure for the radial value. */
    z.v=square(range[0])+square(range[1])+square(I(eps));D rz0=sqrt_dual(sqrt_dual(z)),omega=D(I(1))/(rz0*rz0*rz0);
    D q[2];R q_roundoff_scale[2]={};D scalar_s=L+Q;R abs_s=absmax(scalar_s.v);
    if(abs_s<.75L){
        /* Exact same initial remainder, evaluated without subtracting two
         * broad enclosures of almost identical functions. The finite series
         * is enclosed with an explicit infinite-tail value/derivative bound. */
        D term=D(I(.56396484375L))*scalar_s*scalar_s*scalar_s*scalar_s,tail=term;
        for(int n=4;n<24;n++){term=term*D(I(-.75L-n)/I(n+1))*scalar_s;tail=tail+term;}
        R t=absmax(term.v),value_tail=up(t*abs_s/(1-abs_s)),derivative_tail=up(t*(25/(1-abs_s)+abs_s/((1-abs_s)*(1-abs_s))));
        tail.v=widen(tail.v,value_tail);tail.x=widen(tail.x,up(derivative_tail*absmax(scalar_s.x)));tail.y=widen(tail.y,up(derivative_tail*absmax(scalar_s.y)));
        D A=D(I(.65625L))*Q*Q-D(I(.6015625L))*(D(I(3))*L*L*Q+D(I(3))*L*Q*Q+Q*Q*Q)+tail;
        D B=D(I(.65625L))*(D(I(2))*L*Q+Q*Q)-D(I(.6015625L))*scalar_s*scalar_s*scalar_s+tail;
        q[0]=D(oa)*((-D(ay))*A+(-dy)*B);q[1]=D(oa)*(D(ax)*A+dx*B);
        R lm=absmax(L.v),qm=absmax(Q.v),am=.65625L*qm*qm+.6015625L*(3*lm*lm*qm+3*lm*qm*qm+qm*qm*qm)+absmax(tail.v),bm=.65625L*(2*lm*qm+qm*qm)+.6015625L*abs_s*abs_s*abs_s+absmax(tail.v);
        q_roundoff_scale[0]=absmax(oa)*(absmax(ay)*am+absmax(dy.v)*bm);q_roundoff_scale[1]=absmax(oa)*(absmax(ax)*am+absmax(dx.v)*bm);
    }else{
        q[0]=(-Y)*omega-D(oa)*((-D(ay))*A3+(-dy)*A2);q[1]=X*omega-D(oa)*(D(ax)*A3+dx*A2);
    }
    Enclosure out={};out.s=(L+Q).v;D disp[2];R basis_magnitude=1,basis_gradient_magnitude=1;
    for(int k=0;k<2;k++){
        I coefficient[10];for(int a=0;a<10;a++){
            int index[3];gfs_material_p3_index(a,index);I node(0);for(int j=0;j<3;j++)node=node+I(index[j])*I(r.triangle.geometry.x[j][k])/I(3);
            coefficient[a]=I(r.map.x[k][a])-node+I(r.origin[k]);
        }
        /* Remove the complete affine field before interval evaluation of the
         * higher polynomial modes. Affine reproduction is an algebraic
         * identity of the complete P3 basis, not an approximation or clipping. */
        disp[k]=D(coefficient[0])+bary[1]*D(coefficient[1]-coefficient[0])+bary[2]*D(coefficient[2]-coefficient[0]);
        for(int a=3;a<10;a++){
            int index[3];gfs_material_p3_index(a,index);I affine(0);for(int j=0;j<3;j++)affine=affine+I(index[j])*coefficient[j]/I(3);disp[k]=disp[k]+D(coefficient[a]-affine)*N[a];
        }
        for(int j=0;j<2;j++){
            D residual=q[j]-D(I(r.triangle.remainder[0][j]));for(int a=1;a<10;a++)residual=residual-D(I(r.triangle.remainder[a][j])-I(r.triangle.remainder[0][j]))*N[a];
            I coeff=I(r.map.x[k][10+j])/I(r.scale[j]);disp[k]=disp[k]+D(coeff)*residual;
            basis_magnitude+=absmax(coeff)*(absmax(q[j].v)+1);basis_gradient_magnitude+=absmax(coeff)*(absmax(q[j].x)+absmax(q[j].y)+1);
        }
        /* Linear reference interpolant differs from literal identity only by
         * the stored geometry-gradient arithmetic, retained explicitly. */
        I base[2]={I(0),I(0)};for(int d=0;d<2;d++)for(int j=1;j<3;j++)base[d]=base[d]+(I(r.triangle.geometry.x[j][k])-I(r.triangle.geometry.x[0][k]))*I(r.triangle.geometry.grad[j][d]);
        out.displacement[k]=disp[k].v;out.gradient[2*k]=disp[k].x+base[0]-I(k==0);out.gradient[2*k+1]=disp[k].y+base[1]-I(k==1);
        I vc[3]={I(r.map.v[k][0]),I(r.map.v[k][1]),I(r.map.v[k][2])};D velocity=D(vc[0])+bary[1]*D(vc[1]-vc[0])+bary[2]*D(vc[2]-vc[0]);
        for(int a=3;a<10;a++){int index[3];gfs_material_p3_index(a,index);I affine(0);for(int j=0;j<3;j++)affine=affine+I(index[j])*vc[j]/I(3);velocity=velocity+D(I(r.map.v[k][a])-affine)*N[a];}
        for(int j=0;j<2;j++){D residual=q[j]-D(I(r.triangle.remainder[0][j]));for(int a=1;a<10;a++)residual=residual-D(I(r.triangle.remainder[a][j])-I(r.triangle.remainder[0][j]))*N[a];velocity=velocity+D(I(r.map.v[k][10+j])/I(r.scale[j]))*residual;}
        out.velocity[k]=widen(velocity.v,4096.L*DBL_EPSILON*(1+absmax(velocity.v)));
        R operation_scale=1;for(int a=0;a<10;a++)operation_scale+=fabsl(r.map.v[k][a])*(absmax(N[a].x)+absmax(N[a].y));
        for(int j=0;j<2;j++){R residual_scale=absmax(q[j].x)+absmax(q[j].y);for(int a=1;a<10;a++)residual_scale+=fabsl((R)r.triangle.remainder[a][j]-r.triangle.remainder[0][j])*(absmax(N[a].x)+absmax(N[a].y));operation_scale+=fabsl((R)r.map.v[k][10+j]/r.scale[j])*residual_scale;}
        R velocity_rounding=16384.L*DBL_EPSILON*operation_scale;out.native_velocity_gradient_envelope=std::max(out.native_velocity_gradient_envelope,velocity_rounding);out.velocity_reference_gradient[2*k]=widen(velocity.x,velocity_rounding);out.velocity_reference_gradient[2*k+1]=widen(velocity.y,velocity_rounding);
    }
    /* Explicit conservative engineering envelope for double input-helper and
     * basis evaluation; distinct from long-double interval rounding. The
     * exact analytic map enclosure above already includes current coefficient
     * cancellation. This envelope is recorded and must stay far below margin. */
    out.native_position_envelope=4096.L*DBL_EPSILON*basis_magnitude;
    out.native_gradient_envelope=16384.L*DBL_EPSILON*basis_gradient_magnitude;
    R sum=0;for(int d=0;d<4;d++){out.gradient[d]=widen(out.gradient[d],out.native_gradient_envelope);sum+=square(I(absmax(out.gradient[d]))).hi;}out.lipschitz=up(sqrtl(sum));
    if(abs_s<.75L){
        /* Absolute operation-error envelope for the native position path.
         * The anchor is added once; counting it in every basis operation
         * would overestimate errors by orders of magnitude. Unit roundoff
         * is eps/2. Native pow is assumed within8 ulps; the4096-operation
         * stable-remainder allowance and the explicit tail cover that scale.
         * This is an engineering floating-point model, independently checked
         * against native values; it is not a formal libm proof. */
        R u=DBL_EPSILON/2;auto gamma=[&](int n){return up(n*u/(1-n*u));};R dl[3]={},dxmax=absmax(dx.v),dymax=absmax(dy.v);
        for(int j=1;j<3;j++)dl[j]=gamma(5)*(fabsl(r.triangle.geometry.grad[j][0])*dxmax+fabsl(r.triangle.geometry.grad[j][1])*dymax);dl[0]=dl[1]+dl[2]+gamma(2)*3;
        R shape_bound[10],shape_error[10];for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);if(a<3){shape_bound[a]=1.000000001L;shape_error[a]=6*dl[a]+gamma(16)*8;}else if(a==9){shape_bound[a]=1.000000001L;shape_error[a]=7*(dl[0]+dl[1]+dl[2])+gamma(8)*2;}else{int first=-1,second=-1;for(int j=0;j<3;j++)if(index[j]){if(first<0)first=j;else second=j;}shape_bound[a]=1.100000001L;shape_error[a]=32*(dl[first]+dl[second])+gamma(16)*16;}}
        R residual_scale[2],residual_error[2];for(int j=0;j<2;j++){residual_scale[j]=absmax(q[j].v)+fabsl(r.triangle.remainder[0][j]);residual_error[j]=gamma(4096)*q_roundoff_scale[j];for(int a=1;a<10;a++){R delta=fabsl((R)r.triangle.remainder[a][j]-r.triangle.remainder[0][j]);residual_scale[j]+=delta*shape_bound[a];residual_error[j]+=delta*shape_error[a];}residual_error[j]+=gamma(32)*residual_scale[j];}
        R component[2]={};for(int k=0;k<2;k++){R absolute_sum=0,error=0;for(int a=0;a<10;a++){R c=fabsl(r.map.x[k][a]);absolute_sum+=c*shape_bound[a];error+=c*shape_error[a];}for(int j=0;j<2;j++){R c=fabsl((R)r.map.x[k][10+j]/r.scale[j]);absolute_sum+=c*residual_scale[j];error+=c*(residual_error[j]+gamma(2)*residual_scale[j]);}error+=gamma(24)*absolute_sum;component[k]=up(error+gamma(1)*(fabsl(r.origin[k])+absolute_sum+error));}
        R input_radius=hypotl(absmax(range[0]),absmax(range[1]));out.native_position_envelope=up(hypotl(component[0],component[1])+(1+out.lipschitz)*gamma(1)*input_radius);
    }
    for(int k=0;k<2;k++)out.displacement[k]=widen(out.displacement[k],out.native_position_envelope);
    out.J=(I(1)+out.gradient[0])*(I(1)+out.gradient[3])-out.gradient[1]*out.gradient[2];return out;
}
static I divergence(const Enclosure &b){need(b.J.lo>0,"positive J for velocity-divergence enclosure");I F00=I(1)+b.gradient[0],F01=b.gradient[1],F10=b.gradient[2],F11=I(1)+b.gradient[3];return (b.velocity_reference_gradient[0]*F11-b.velocity_reference_gradient[1]*F10-b.velocity_reference_gradient[2]*F01+b.velocity_reference_gradient[3]*F00)/b.J;}
static std::array<Tri,4> split(const Tri &t){std::array<R,2> a,b,c;for(int k=0;k<2;k++){a[k]=(t[0][k]+t[1][k])/2;b[k]=(t[1][k]+t[2][k])/2;c[k]=(t[2][k]+t[0][k])/2;}return {Tri{t[0],a,c},Tri{a,t[1],b},Tri{c,b,t[2]},Tri{a,b,c}};}
struct Certificate {R min_J=std::numeric_limits<R>::infinity(),max_J=-std::numeric_limits<R>::infinity(),max_lipschitz=0,max_displacement=0,native_position_envelope=0,native_gradient_envelope=0,max_abs_s=0;uint64_t panels=0,leaves=0;int max_depth=0;bool complete=true;};
static void certify(const curved::Record &r,const Tri &tri,R eps,int depth,int cap,R max_lipschitz,Certificate &out){
    auto e=enclose(r,tri,eps);out.panels++;bool valid=e.J.lo>0 && e.lipschitz<max_lipschitz && absmax(e.s)<.75L;
    if(!valid && depth<cap){for(const auto &s:split(tri))certify(r,s,eps,depth+1,cap,max_lipschitz,out);return;}
    out.complete=out.complete && valid;out.leaves++;out.max_depth=std::max(out.max_depth,depth);out.min_J=std::min(out.min_J,e.J.lo);out.max_J=std::max(out.max_J,e.J.hi);out.max_lipschitz=std::max(out.max_lipschitz,e.lipschitz);out.max_abs_s=std::max(out.max_abs_s,absmax(e.s));
    out.max_displacement=std::max(out.max_displacement,up(sqrtl(square(I(absmax(e.displacement[0]))).hi+square(I(absmax(e.displacement[1]))).hi)));out.native_position_envelope=std::max(out.native_position_envelope,e.native_position_envelope);out.native_gradient_envelope=std::max(out.native_gradient_envelope,e.native_gradient_envelope);
}
static Enclosure positive_enclosure(const curved::Record &r,const Tri &tri,R eps){
    auto out=enclose(r,tri,eps);if(out.J.lo>0)return out;
    /* A wide Cartesian box may overestimate even a J=1 quadratic shear.
     * Tighten only the Jacobian range using the existing bounded source
     * triangle certificate. Keep the original conservative field/evaluation
     * bounds; no state, coefficient or manufactured deformation is changed. */
    Certificate c;certify(r,tri,eps,0,6,.24L,c);need(c.complete && c.min_J>0,"positive map after bounded source-triangle interval refinement");out.J=I(c.min_J,c.max_J);return out;
}
}
#endif
