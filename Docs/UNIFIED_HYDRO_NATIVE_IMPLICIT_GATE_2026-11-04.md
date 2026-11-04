# 공통 implicit 점성의 실제 native 결합 및 유한 열화 검사

전체 목표는 진행 중이다. 1D contact 원인 규명과 공통 대수 검증을
완료했고, 공통 연산자로 실제 native 2D 진화에 처음 연결했다.
**일반 2D 전체 시험, 3D 진화, GIZMO 대비 실행시간 2배 기준은 미완료다.**
기존 20pi 무점성 기준선과 새 점성 후보는 서로 다른 실행이다.

## 코드와 재현 위치

- 개발 checkout /gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1
- 실험 /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1
- 동결 결합 실행 파일 native_coupled_build_v5/unified_native_coupled.exe
- SHA256 547dd5622a42782949d85e87a14f5c732b3bf59985a306fbd66b52db322dd6cf
- 공통 graph 연산자 Exam/gfs_viscosity_implicit.{h,c}
- native adapter Exam/NativeU64/gfs_unified_viscosity.inc
- native KDK Exam/NativeU64/gfs_unified_kdk.inc

기존 OpenMP+MPI+CUDA 본체와 696-byte restart ABI를 사용했다.
기하는 기존 GPU Voronoi, 압력은 공통 CUDA 함수, 점성은 공통 C/OpenMP
PCG이다. Python은 실행 조직과 결과 분석만 수행한다.
기존 검증 archive에 새 object를 부분 재링크했으며 전체 clean build는 아니다.
관련 build.sh도 새 공통 C object를 포함하도록 갱신했다.

## 동일 소산식의 차원 공통 시간 갱신

중앙 pair 계수 D>=0을 half-step 시작에 고정하고, 종점에서 접근하는
pair만 활성화하는 BE 최소화를 푼다. 활성집합별 SPD 선형계는
(M+hL) delta_v = -h L v0 이다. 압력 kick을 점성 열에 넣지 않는다.
열 배분은 pair의 h*D*g(v1)^2와 셀의 0.5*m*|delta_v|^2를 합한다.
벽은 유체에 전체 소산 열을 주고 외부 충격량과 토크를 분리한다.
BE 점성은 시간 1차이며 대칭 half-step 배치만으로 2차가 되지 않는다.

공통 C graph 구현의 1D 결과와 기존 tridiagonal 결과의 최대 속도 차이는
7.43e-13이다. 1D/2D/3D graph 24개 사례에서 에너지, 운동량,
각운동량과 양의 열을 검사했다. 큰 Galilean boost의 반올림 오차도
별도 기록했다. graph 대수 통과는 3D 유체 진화 통과가 아니다.

1D 실제 점성 단계 9개 검사에서 에너지와 경계 반작용을 포함한 운동량
잔차가 반올림 범위이고 질량과 위치가 바뀌지 않았다.
Sod/blast/Noh의 explicit/implicit 각각 dt, dt/2, dt/4, dt/8로
총 24개 전체 지정 종점 계산이 완료됐다. Noh implicit 전체 에너지 오차는
8.98e-4 -> 2.54e-4 -> 6.63e-5 -> 1.62e-5로 감소했다.
같은 반이산 해로 수렴하는지 실제 공간 구간 적분으로도 비교했다.
이는 충격 해가 시간 2차라는 주장이 아니다.

## 실제 native 2D 결합 범위

P/2 -> D/2 -> 연화 중력 KDK subcycle + drift/단열 체적 갱신
-> D/2 -> P/2 순서를 사용했다. W=0, 고정질량, 물리점성 0,
추가 복원력·floor·조향은 없다. 원래 ID로 drift 전후 체적을 대응시켜
migration이 입자 순서를 바꾸어도 단열 갱신의 셀 대응을 유지했다.
현재 entry gate는 저장된 one-rank 고정 polygon Kepler 상태다.
periodic alias, 일반 IC/경계, MPI 분할과 3D adapter는 미완료다.
압력 벽 반작용의 전체 토크 장부도 후속 검증 대상이다.

초기 기하에는 반올림 크기의 점 접촉 선분이 한쪽에만 남거나,
같은 이웃의 면이 긴 선분과 극소 선분으로 나뉜 경우가 있었다.
같은 이웃의 선분 면적을 합하고 양쪽 면적 평균으로 하나의 중앙 pair를
만들었다. 해상도 있는 비상호 면은 거부한다. clipping 반올림 경계는
절대 좌표의 뺄셈 오차와 셀 크기를 함께 사용한다.
실제 t20pi 최대 면적 차이는 1.03e-14이며 기록한 반올림 경계 이내다.
해당 변경은 숨겨진 pair force나 입자 이동을 추가하지 않는다.

- t0 -> .002, 90,352셀, 2 accepted steps, dErel=2.17e-9.
- 기존 20pi 상태 -> 20pi+1e-4, 2 accepted steps, dErel=-1.01e-11.
- 후자의 첫 dt는 7.84e-5로 explicit diffusion 제한 8.95e-8보다 약 876배 크다.
- 두 실행 모두 ID와 질량이 정확히 유지됐고 W=pi=0, U/V>0이다.

기하와 전송 비용을 포함한 production 성능 통과는 아니다.
짧은 선형계 solve 시간이나 큰 dt 비율을 GIZMO 가속비로 부르지 않는다.

## 유한 열화와 압력 반응 gate

Astra의 단계 평가에 따라, 기존 .002 결과부터 이어서 .01/.03/.1까지
계산했다. 독립 dt/2 비교는 t0에서 시작했다. 두 branch 모두 같은
동결 source와 계수를 사용했고 총 6개 최종 스냅샷만 추가했다.
원래 ID의 core R0<.5, disk .5<=R0<=2, ambient R0>2를 계속 추적했다.
반사 event는 없었고 기존 입력과 출력은 보존했다.

| 측정 | 기본 dt | dt/2 |
|---|---:|---:|
| core U(t=.1), 초기 U=1.17522e-6 | 4.06636e-4 | 3.55262e-4 |
| disk 초기 체적에서 환산한 엔트로피 열 증가 | 1.18310e-7 | 1.18297e-7 |
| disk 초기 U | 1.78048e-5 | 1.78048e-5 |
| paper 현재 annulus 속도 L1(t=.1) | 2.20355e-7 | 2.20240e-7 |
| 원래 disk 전체 속도 L1(t=.1) | 1.73155e-7 | 1.72952e-7 |
| t0 기준 전체 dErel(t=.1) | 1.12567e-5 | 2.98342e-7 |

실제 D substep 열 장부의 disk 열률은 전 구간 약 1.183e-6으로
거의 일정하고 dt/2에서도 같다. core 열률은 .01-.03에서
약 .00264/.00263, .03-.1에서 .00525/.00446으로 다시 증가한다.
따라서 일시적인 초기 잡음 감쇠로 결론내리지 않는다.
core 열량과 전체 에너지 오차는 아직 시간수렴하지 않았다.
작은 disk 궤도 오차만으로 장기 정확도나 전체 밀도장 통과를 주장하지 않는다.

모든 단열 drift의 메모리 내 U*V^(gamma-1) 상대 잔차는 최대 4.44e-16이다.
소산 열은 비음수이고 점성 에너지 잔차는 약 1e-18이다.
처음 발견한 출력 엔트로피의 약 -2.4e-8 잔차는 초기 파일 저장 V와
첫 실제 GPU V의 최대 상대차 5.716e-8에서 왔다.
실제 초기 V로 비교하면 최소 엔트로피 비율은 1-3e-16이다.
GasInfo.gamma의 동일 float 값을 CPU/CUDA가 double로 승격한다.
실제 Gamma 불일치는 없었으며 Gamma 보정은 추가하지 않았다.

## 저장 상태만 이용한 소산 원인 분리

원래 q=max(reconstruction branch, compression branch)에서
compression-only 열률과 나머지 추가 열률을 분리했다.
두 branch를 독립적으로 더하지 않아 중복 집계를 피했다.

초기 core 열률의 약 95.9%는 reconstruction 추가 경로다.
초기 disk 열률의 약 77.2%는 fitted divergence 압축 경로다.
따라서 reconstruction 경로만 제거하면 해결된다고 결론내릴 수 없다.

동일 m,V,U,기하를 유지한 진단용 속도/gradient 대조를 수행했다.
정확 gradient는 원인 분리용 oracle이며 생산 코드에 넣지 않았다.

| t0 disk 열률 | 값 |
|---|---:|
| 실제 속도 + 기존 LSQ | 1.18679e-6 |
| 같은 위치의 IC 원운동 속도 + 기존 LSQ | 1.18679e-6 |
| IC 원운동 속도 + 정확 gradient | 3.80376e-10 |

.03과 .1 상태에서도 정확 gradient의 disk 열률은 약 3.9e-10이다.
초기 core 열률은 .003581에서 정확 gradient로 .0003744로 줄지만
곡률에 의한 남은 소산이 있다. LSQ condition bound는 core/disk에서
약 2.5로 작아 illconditioning으로 설명되지 않는다.

가중치 A/d를 주는 면적 LSQ도 **frozen RHS 진단에서만** 비교했다.
disk 열률은 1.1868e-6 -> 7.5702e-7, core는 약 11% 감소에 그친다.
이 대안을 생산 fit에 추가하거나 해결책으로 승격하지 않았다.

## 일반 유동의 공간 척도 검사

실제 프로젝트 Voro가 생성한 정사각 주기 셀과 같은 공유 연산자로
psi=sin(2pi*x)*sin(4pi*y)/(4pi)의 divergence-free 유동을 검사했다.
P=1e-6, rho=1이며 새 중력이나 Kepler 가정은 없다.

| 한 축 셀 수 | 기존 LSQ 열률 | 정확 gradient 진단 열률 |
|---:|---:|---:|
| 16 | .0519788 | .0106539 |
| 32 | .00532563 | .000694327 |
| 64 | .000507812 | .0000440514 |

기존 LSQ 소산은 해상도 두 배에서 약 9.76/10.49배 감소한다.
반면 격자 교대 모드는 H*t_cross/K=8.00258로 해상도에 무관하다.
매끄러운 곡률 유동과 교대 모드를 공간 척도로 구별하는 증거다.
해당 시험은 고정장 RHS 검사이며 주기 유체 시간진화 검증은 아니다.

## 다음 우선순위

1. 같은 실제 native 본체에서 일반 IC, periodic 이미지/소유 셀 대응,
   경계 및 KDK를 연결해 이동 contact -> affine shear -> 평면 shock을 실행한다.
2. 실제 면과 V에 맞는 보존량 IC를 사용하고 회전된 평면도 비교한다.
   1D 해와의 비교를 차원 공통 기법 검증으로 사용한다.
3. 동일 scheme으로 KH/RT/Gresho/Sedov/Noh/cylinder/Kepler 전체 2D를 재수행하고
   실제 3D 진화 adapter와 평면 삽입 검증을 완료한다.
4. 물리 정확도 판정 뒤 중복 geometry/graph 재사용, 그 뒤 필요할 때
   implicit GPU 이전을 한 가지씩 평가한다. 동등한 GIZMO 정확도·IC·자원·출력
   조건의 전체 처리량으로 2배 실행시간 기준을 판정한다.

현재 후보의 장기 Kepler와 속도 승격은 보류한다. 새 detector, 고차 fit,
Kepler 전용 교정과 여러 closure를 동시에 추가하지 않는다.
모든 파일은 /gpfs에 두며 /scratch는 사용하지 않는다.
