from pathlib import Path
s=Path('/tmp/contact_backend_probe_20261018/probe.cu').read_text().split('#define CK')[0];s=s.replace('static int sedov_phase1_on(){return 0;}','static thread_local int diagnostic_phase1=0;\nstatic int sedov_phase1_on(){return diagnostic_phase1;}')
s+='''
#define CK(x) do{cudaError_t e=(x);if(e!=cudaSuccess){fprintf(stderr,"%s\\n",cudaGetErrorString(e));return 2;}}while(0)
__global__ void kernel(const double*a,double*o,int N){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=N)return;const double*p=a+11*i;double*z=o+6*i;
dev_hllc_face_2d_rest_frame(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],(int)p[10],z,z+1);
dev_hllc_face_2d_rest_frame(p[0],p[1],p[2]+7.,p[3],p[4],p[5],p[6]+7.,p[7],p[8]+7.,p[9],(int)p[10],z+2,z+3);
dev_hllc_face_2d_rest_frame(p[4],p[5],-p[6],p[7],p[0],p[1],-p[2],p[3],-p[8],p[9],(int)p[10],z+4,z+5);}
int main(){FILE*f=fopen("inputs.bin","rb");if(!f)return 3;fseek(f,0,SEEK_END);long size=ftell(f);rewind(f);if(size%(11*sizeof(double)))return 4;int N=size/(11*sizeof(double));double*a=new double[11*N];if(fread(a,sizeof(double),11*N,f)!=(size_t)(11*N))return 5;fclose(f);double*one=new double[6*N],*four=new double[6*N],*gpu=new double[6*N];int teams[2]={0,0};omp_set_dynamic(0);
for(int pass=0;pass<2;pass++){omp_set_num_threads(pass?4:1);double*out=pass?four:one;
#pragma omp parallel
{
#pragma omp single
teams[pass]=omp_get_num_threads();
#pragma omp for
for(int i=0;i<N;i++){const double*p=a+11*i;double*z=out+6*i;diagnostic_phase1=(int)p[10];
hllc_face_2d_rest_frame(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],z,z+1);
hllc_face_2d_rest_frame(p[0],p[1],p[2]+7.,p[3],p[4],p[5],p[6]+7.,p[7],p[8]+7.,p[9],z+2,z+3);
hllc_face_2d_rest_frame(p[4],p[5],-p[6],p[7],p[0],p[1],-p[2],p[3],-p[8],p[9],z+4,z+5);}}}
double*da,*dout;CK(cudaMalloc(&da,size));CK(cudaMalloc(&dout,6*N*sizeof(double)));CK(cudaMemcpy(da,a,size,cudaMemcpyHostToDevice));kernel<<<(N+127)/128,128>>>(da,dout,N);CK(cudaGetLastError());CK(cudaDeviceSynchronize());CK(cudaMemcpy(gpu,dout,6*N*sizeof(double),cudaMemcpyDeviceToHost));double diff=0;int finite=1;for(int i=0;i<6*N;i++){diff=fmax(diff,fabs(one[i]-gpu[i]));finite=finite&&std::isfinite(one[i])&&std::isfinite(gpu[i]);}f=fopen("cpu_outputs.bin","wb");fwrite(one,sizeof(double),6*N,f);fclose(f);
printf("{\\"N\\":%d,\\"teams\\":[%d,%d],\\"finite\\":%d,\\"cpu_omp_bitwise\\":%s,\\"cpu_cuda_bitwise\\":%s,\\"backend_max_abs\\":%.17g}\\n",N,teams[0],teams[1],finite,memcmp(one,four,6*N*sizeof(double))?"false":"true",memcmp(one,gpu,6*N*sizeof(double))?"false":"true",diff);return !(finite&&teams[0]==1&&teams[1]==4&&diff<1e-12);}
'''
p=Path('/tmp/existing_interface_probe_20261018.cu');p.write_text(s)
