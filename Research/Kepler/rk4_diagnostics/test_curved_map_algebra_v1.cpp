/* Actual curved packet algebra for future remap. No overlap/weight solve. */
#include <cfloat>
#include <algorithm>
#include <tuple>
#define GFS_ENRICH_TEST_LIBRARY
#include "test_enriched_euler_v1.cpp"
#include "material_curved_packet.h"
using R=curved::R;
static bool contains(const curved::Record &r,double X,double Y){
    const auto &g=r.triangle.geometry;double x=X-g.x[0][0],y=Y-g.x[0][1],l1=g.grad[1][0]*x+g.grad[1][1]*y,l2=g.grad[2][0]*x+g.grad[2][1]*y,l0=1-l1-l2,tol=64*DBL_EPSILON*(1+fabs(l1)+fabs(l2));return l0>=-tol && l1>=-tol && l2>=-tol;
}
int main(int argc,char **argv){
    need(argc==4,"curved packet, original native polygon directory, new output");auto packet=curved::read(argv[1]);std::set<uint64_t> ids;for(const auto &r:packet.record)ids.insert(r.cell);auto original=load(argv[2],ids);
    std::map<uint64_t,std::array<double,2>> sites;std::string sitepath=std::string(argv[3])+"/advected_sites.csv";FILE *sf=fopen(sitepath.c_str(),"w");need(sf,"new actual advected generator sites");fprintf(sf,"id,reference_x,reference_y,mapped_x,mapped_y,vx,vy\n");
    for(auto id:ids){double X=original.at(id).gx-12,Y=original.at(id).gy-12;bool found=false;for(const auto &r:packet.record)if(r.cell==id && contains(r,X,Y)){
        auto p=curved::native_point(r,X,Y,packet.parameters[3]);sites[id]={p.position[0],p.position[1]};fprintf(sf,"%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",(unsigned long long)id,X,Y,p.position[0],p.position[1],p.velocity[0],p.velocity[1]);found=true;break;
    }need(found,"every original generator mapped by actual fluid map");}fclose(sf);
    R position_error=0,velocity_error=0,radial_error=0,tangent_error=0,Jacobian_error=0,plane_error=0,plane_derivative_error=0,max_curve=0,max_curve_relative=0;
    uint64_t queries=0,planes=0;double q[4],w[4];gauss(4,q,w);
    for(const auto &r:packet.record){
        const auto &g=r.triangle.geometry;R phase=atan2l((R)g.x[0][1]+g.x[1][1]+g.x[2][1],(R)g.x[0][0]+g.x[1][0]+g.x[2][0]);R length=sqrtl(2*(R)g.area);
        for(int i=0;i<4;i++)for(int j=0;j<4;j++){
            double l1=q[i],l2=(1-q[i])*q[j],l0=1-l1-l2,X=l0*g.x[0][0]+l1*g.x[1][0]+l2*g.x[2][0],Y=l0*g.x[0][1]+l1*g.x[1][1]+l2*g.x[2][1];R radius=hypotl(X,Y);if(radius==0)continue;R theta=atan2l(Y,X),delta=atan2l(sinl(theta-phase),cosl(theta-phase)),u=tanl(delta/2);auto ray=curved::ray(r,radius,phase,packet.parameters[3]);auto native=curved::native_point(r,X,Y,packet.parameters[3]);std::array<R,3> x[2];
            for(int k=0;k<2;k++){
                x[k]=curved::evaluate(ray.position[k],u);auto v=curved::evaluate(ray.velocity[k],u);position_error=std::max(position_error,fabsl(x[k][0]-native.position[k])/length);velocity_error=std::max(velocity_error,fabsl(v[0]-native.velocity[k])/std::max(1.L,fabsl(native.velocity[k])));
                R nr=(native.F[2*k]*X+native.F[2*k+1]*Y)/radius,nt=-native.F[2*k]*Y+native.F[2*k+1]*X;radial_error=std::max(radial_error,fabsl(x[k][1]-nr));tangent_error=std::max(tangent_error,fabsl(x[k][2]-nt)/radius);
            }
            R J=(x[0][1]*x[1][2]-x[0][2]*x[1][1])/radius;Jacobian_error=std::max(Jacobian_error,fabsl(J-native.J)/native.J);auto own=sites.at(r.cell);int wi=std::distance(sites.begin(),sites.find(r.cell));
            for(const auto &entry:sites)if(entry.first!=r.cell){int wj=std::distance(sites.begin(),sites.find(entry.first));R ax=(R)entry.second[0]-own[0],ay=(R)entry.second[1]-own[1],b0=.5L*((R)entry.second[0]*entry.second[0]+(R)entry.second[1]*entry.second[1]-(R)own[0]*own[0]-(R)own[1]*own[1]),norm=hypotl(ax,ay);
                for(int weights=0;weights<2;weights++){
                    R b=b0+(weights?.5e-7L*(sinl(.7L*wi)-sinl(.7L*wj)):0);auto polynomial=curved::power_plane(ray,ax,ay,b);auto value=curved::evaluate(polynomial,u);R exact=ax*native.position[0]+ay*native.position[1]-b;plane_error=std::max(plane_error,fabsl(value[0]-exact)/(norm*length));R er=(ax*native.F[0]+ay*native.F[2])*X/radius+(ax*native.F[1]+ay*native.F[3])*Y/radius,et=-(ax*native.F[0]+ay*native.F[2])*Y+(ax*native.F[1]+ay*native.F[3])*X;
                    plane_derivative_error=std::max({plane_derivative_error,fabsl(value[1]-er)/norm,fabsl(value[2]-et)/(norm*radius)});planes++;
                }
            }queries++;
        }
        for(int e=0;e<3;e++){
            int k=(e+1)%3;auto a=curved::native_point(r,g.x[e][0],g.x[e][1],packet.parameters[3]),b=curved::native_point(r,g.x[k][0],g.x[k][1],packet.parameters[3]),mid=curved::native_point(r,.5*(g.x[e][0]+g.x[k][0]),.5*(g.x[e][1]+g.x[k][1]),packet.parameters[3]);R gap=hypotl((R)mid.position[0]-.5L*(a.position[0]+b.position[0]),(R)mid.position[1]-.5L*(a.position[1]+b.position[1]));max_curve=std::max(max_curve,gap);max_curve_relative=std::max(max_curve_relative,gap/length);
        }
    }
    printf("{\"scope\":\"actual nonuniform P3 curved-map packet and arbitrary current coefficients; degree6 power-plane preimage algebra and radial/angular derivatives only; no root isolation, overlap integration, mass/weight solve or re-expression claim\",\"triangles\":%zu,\"physical_cells\":%llu,\"queries\":%llu,\"power_plane_comparisons\":%llu,\"polynomial_degree_cap\":6,\"position_difference_over_triangle_length\":%.17Lg,\"velocity_relative_difference\":%.17Lg,\"radial_derivative_difference\":%.17Lg,\"angular_derivative_difference_over_radius\":%.17Lg,\"Jacobian_relative_difference\":%.17Lg,\"power_plane_difference_over_normal_length\":%.17Lg,\"power_plane_derivative_normalized_difference\":%.17Lg,\"actual_edge_midpoint_chord_gap_max\":%.17Lg,\"actual_edge_chord_gap_over_triangle_length_max\":%.17Lg}\n",packet.record.size(),(unsigned long long)packet.cells,(unsigned long long)queries,(unsigned long long)planes,position_error,velocity_error,radial_error,tangent_error,Jacobian_error,plane_error,plane_derivative_error,max_curve,max_curve_relative);
    need(std::max({position_error,velocity_error,radial_error,tangent_error,Jacobian_error,plane_error,plane_derivative_error})<1.e-9L,"actual curved-map rational algebra gate");need(max_curve_relative>1.e-10L,"actual curved geometry, not affine stand-in");return 0;
}
