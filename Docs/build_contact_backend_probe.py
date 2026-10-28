"""Extract existing CPU/CUDA contact-state routines without changing their arithmetic."""
from pathlib import Path
import re,json,hashlib
repo=Path('/tmp/lageunha-push');out=Path('/tmp/contact_backend_probe_20261018');out.mkdir(exist_ok=False)
def extract(path,name,qualifier):
 s=path.read_text();clean=re.sub(r'/\*.*?\*/|//[^\n]*','',s,flags=re.S);m=re.search(r'void\s+'+name+r'\s*\(',clean);assert m,name;start=m.start();b=clean.index('{',start);level=1;j=b+1
 while level:
  level+=(clean[j]=='{')-(clean[j]=='}');j+=1
 return qualifier+' '+clean[start:j]+'\n'
cpu=repo/'Exam/exam.c';gpu=repo/'Exam/exam_gpu.cu'
preamble='#include <cmath>\n#include <cstdio>\n#include <cstdlib>\n#include <cstring>\n#include <omp.h>\n#include <cuda_runtime.h>\ntypedef double postype;\n#define HLLC_VACUUM_PRATIO 100.0\n#define HLLC_VACUUM_PMIN_GUARD .01\nstatic int sedov_phase1_on(){return 0;}\n'
s=preamble
for n in ('hll_star_state','hllc_face_2d','hllc_face_2d_rest_frame'):s+=extract(cpu,n,'static inline')
for n in ('dev_hll_star_state','dev_hllc_face_2d','dev_hllc_face_2d_rest_frame'):s+=extract(gpu,n,'__device__')
s+='''
#define CK(x) do{cudaError_t e=(x);if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\\n",cudaGetErrorString(e));return 2;}}while(0)
__global__ void run(const double *a,double *o,int N){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=N)return;const double *p=a+7*i;dev_hllc_face_2d_rest_frame(p[0],p[2],p[3],p[4],p[1],p[2],p[3],p[5],p[6],5./3.,0,o+2*i,o+2*i+1);}
int main(){const int N=4096;double *a=new double[7*N],*one=new double[2*N],*four=new double[2*N],*gpu=new double[2*N];
for(int i=0;i<N;i++){double *p=a+7*i;p[0]=(i%2)?.01:1.;p[1]=(i%3)?.3:1.;p[2]=(i%5)?1e-5:1.;p[3]=(i-2048)*.013;p[4]=sqrt((5./3.)*p[2]/p[0]);p[5]=sqrt((5./3.)*p[2]/p[1]);p[6]=p[3]+((i%7)-3)*.017;}
for(int pass=0;pass<2;pass++){omp_set_num_threads(pass?4:1);double *out=pass?four:one;
#pragma omp parallel for
for(int i=0;i<N;i++){double *p=a+7*i;hllc_face_2d_rest_frame(p[0],p[2],p[3],p[4],p[1],p[2],p[3],p[5],p[6],5./3.,out+2*i,out+2*i+1);}}
double *da,*do_;CK(cudaMalloc(&da,7*N*sizeof(double)));CK(cudaMalloc(&do_,2*N*sizeof(double)));CK(cudaMemcpy(da,a,7*N*sizeof(double),cudaMemcpyHostToDevice));run<<<(N+127)/128,128>>>(da,do_,N);CK(cudaGetLastError());CK(cudaDeviceSynchronize());CK(cudaMemcpy(gpu,do_,2*N*sizeof(double),cudaMemcpyDeviceToHost));
double ep=0,ev=0,backend=0;int finite=1;for(int i=0;i<N;i++){ep=fmax(ep,fabs(one[2*i]-a[7*i+2]));ev=fmax(ev,fabs(one[2*i+1]-a[7*i+3]));for(int j=0;j<2;j++){finite=finite&&std::isfinite(one[2*i+j])&&std::isfinite(gpu[2*i+j]);backend=fmax(backend,fabs(one[2*i+j]-gpu[2*i+j]));}}
printf("{\\"cases\\":%d,\\"finite\\":%d,\\"CPU_OMP_bitwise\\":%s,\\"CPU_CUDA_bitwise\\":%s,\\"backend_max_abs\\":%.17g,\\"pressure_reference_error\\":%.17g,\\"normal_velocity_reference_error\\":%.17g}\\n",N,finite,memcmp(one,four,2*N*sizeof(double))?"false":"true",memcmp(one,gpu,2*N*sizeof(double))?"false":"true",backend,ep,ev);
return !(finite&&ep<1e-12&&ev<1e-12&&backend<1e-12&&memcmp(one,four,2*N*sizeof(double))==0);}
'''
(out/'probe.cu').write_text(s)
(out/'provenance.json').write_text(json.dumps({'scope':__doc__,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (cpu,gpu)},'phase1':0,'type':'postype double','limits':['Extracted helper routines only, not full production kernel.','Exact contact star pressure/normal velocity only; no full ALE flux or MPI.']},indent=2)+'\n')
print(out)
