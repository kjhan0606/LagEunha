# 차원 공통 유체역학 진행 기록

현재 목표는 1D contact 원인 규명, 1D/2D/3D 기법 통일, 같은 기법으로
2D 시험 재수행, 동등 조건 GIZMO 실행시간의 2배 이하 달성이다.
**목표는 진행 중이다. 공통 기법의 기본 1D/2D/3D 실제 진화와 축 shock
횡속도 수정은 검증했다. 전체 2D suite, 새 방법의 장기 Kepler와 동등
조건 GIZMO 속도 기준은 남아 있다.**

## 작업 위치 및 보존한 기준선

- 현재 개발 checkout
  /gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1
- branch codex/native20pi-contact-audit-20261104
- 개발 실행 및 분석
  /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1
- 별도 초기 인과 감사
  /gpfs/kjhan/LagForce/Contact1D_CausalAudit_20261104_v1
- 기존 Kepler 20pi 소스와 원고 snapshot은 a1b5247로 GitHub에 푸시했다.
  Overleaf main은 6cb0e1a이고 원격 ref까지 확인했다.
  GitHub SSH의 서버 오류/인증 오류 후 credential helper를 이용한 HTTPS
  push가 성공했다. force push와 기존 dirty checkout 변경은 하지 않았다.
- 이후 공통 solver 개발 파일은 위 checkout의 추가 변경이다.
  기존 20pi 실행을 새 공통 solver 결과로 재분류하지 않는다.

## 확정한 원인과 결과

1. 기존 1D contact는 등질량·불균등 간격 초기화 뒤 실제 Voronoi
   체적으로 rho를 바꾸면서 e를 유지했다. 첫 RHS 전에 P=1 평형이
   약 0.379/1.303으로 깨진다. 33개 C 대조 계산으로 확인했다.
2. 공통 W0 연산자는 실제 V, rho=m/V, P=(gamma-1)U/V,
   B=dV/dx, F=B^T P, Udot=-P Bv를 쓴다.
   실제 프로젝트 Voro로 1D/2D/3D 기하 미분을 검사했다.
   최대 유한차분 오차는 4.03e-11, 전체 힘과 에너지 RHS 잔차는
   약 1e-16이다. 이는 연산자 검사이며 다차원 진화 검증이 아니다.
3. 같은 단일 선형대수 velocity fit을 1D/2D/3D에 제공한다.
   affine fit 최대 오차는 2.22e-16이고 rank 부족은 명시적으로 거부한다.
4. 동일 공유 헤더의 CPU 및 실제 Ada CUDA 검사가 통과했다.
   device_audit.json은 세 차원의 fit, 단열 갱신, pair 열화 및
   affine shear 소산을 검사한다. 이 검사도 완전한 유체 진화가 아니다.

## 충격 closure 선택 과정

첫 acoustic drag 및 유한 속도차 기반 AV 후보는 contact를 보존해도
Sod/blast에서 생성점 간격이 0으로 닫혔다. 단순히 계수를 계속 올리거나
복원력을 추가하지 않았다. 현재 계열의 실패 자료는 아래에 보존했다.

- first_evolution_results.json
- monaghan_candidate_results.json
- compression_switch_candidate_results.json

Sod의 accepted t=0.014068225066433274 상태를 저장하고 같은 RHS를
KDK/RK4, dt/dt2/dt4로 이어 계산했다. 모두 t=0.0168730517로 수렴하는
붕괴를 보였다. 체적은 양수이고 접근속도는 약 0.183738이다.
RK4 에너지 잔차가 2.3e-13이어도 붕괴한다.
따라서 해당 후보의 실패는 KDK 분할 오차로 설명되지 않는다.
integrator_causal_comparison에 checkpoint와 pair별 힘 분해가 있다.

다음 후보는 h*delta_v/d라는 실제 변형률을 쓰는 하나의 선형/이차
압축 점성이다. alpha=1, beta=2를 모든 문제에 동일하게 사용한다.
힘은 중앙 pair force이고 운동에너지 소산을 양의 열화와 짝지었다.
정확 affine shear에서 재구성 잔차가 0이 되도록 공통 velocity fit을
사용한다. 압축/회전 switch는 물리 속도 fit에서 계산한다.
Voronoi 체적 변화율 자체를 물리 divergence로 부르지 않는다.

재구성은 실제 접근속도에 곱하는 [0,1]의 점성 계수로 제한한다.
실제 접근속도가 0으로 갈 때 유한 힘이 남는 constitutive discontinuity와
q/approach의 발산을 피한다. uncapped_fit_source에는 해당 실패를 남겼다.
추가 pair restoring force, 에너지 floor, 입자 조향은 사용하지 않는다.

## 1D 진화 상태

초기화는 실제 균일 셀에 m, momentum, total energy를 적분한다.
각 셀 질량은 진화 중 고정이다. 유입 경계에서 삽입한 새 셀은 별도다.
global KDK에서 고정 위치 kick은 midpoint RK2를 사용하고,
drift의 열에너지는 정확 단열 체적 관계로 갱신한다.
조기 종료는 실패 status와 비영 exit code로 보고한다.

N=200에서 정지 contact, 이동 contact, Sod, blast, Noh, Lax,
double rarefaction, collision, Shu--Osher 모두 전체 지정 시각에 도달했다.
열 floor나 생성점 교차는 없다. bounded_strain_fit_results.json에 있다.
기존 minmod fit 대조군의 N=400/800도 전체 시각에 도달했다.
새 공통 fit N=200과 minmod 대조군은 최종 출력이 동일했다.

| N | Sod L1rho | blast L1rho | Noh L1rho |
|---:|---:|---:|---:|
| 200 | 0.00643925 | 0.0562251 | 0.0157278 |
| 400 | 0.00335522 | 0.0327945 | 0.00787103 |
| 800 | 0.00170778 | 0.0156568 | 0.00394206 |

contact N=200의 L1rho는 약 4.62e-16이고 max|P-1|=2.22e-14다.
이동 contact의 max|P-1|는 약 3.36e-12다.
에너지 잔차는 shock에서도 검사하며 N=200 Sod 약 1.37e-5,
blast 약 1.53e-5, Noh 약 3.22e-5다.
완료 및 수렴은 속도 목표 통과와 구분한다.

## Shu--Osher 기준해 오류

기존 ref.dat의 실제 좌표 범위는 [-0.3185,4.9908]이고 초기 영역은
[-5,5]다. 기존 기준해는 유입 입자를 넣지 않아 upstream driver가
고갈됐고 leftmost P는 9.6185까지 낮아졌다. 새 N=3200 기준 계산은
같은 위치의 upstream P=10.33333을 유지한다.
따라서 기존 ref 대비 약 1.46의 밀도 norm이 정체되는 현상을
새 solver의 비수렴으로 해석할 수 없다.

동일 새 모델 N=3200 해와 공통 구간 [-4.9,4.9]에서 실제 구간 적분
L1rho를 비교하면 N=200/400/800/1600에서
0.5956/0.3339/0.1670/0.06646으로 감소한다.
shu_reference_audit.json에 계산 방법과 범위가 있다.
이는 self-convergence이며 독립 Euler/GIZMO 기준해 검증은 남아 있다.

## 속도와 남은 구현

불필요한 RHS 재평가를 제거해 한 step의 RHS를 7회에서 4회로 줄였다.
9종 N=200 최종 snapshot은 최적화 전후 byte-identical이다.
cached_rhs_results.json에 검증했다. 실행시간은 동시 실행 영향을 받으므로
이 짧은 시험을 GIZMO 대비 속도 달성으로 보고하지 않는다.

새 공통 CUDA adapter를 Exam/NativeU64/gfs_unified_gpu.inc에 추가했다.
현재 최초 연결 범위는 fully owned one-rank W0 및 고정 벽 face의 RHS다.
production global KDK, periodic image synchronization, 일반 boundary,
IC projection 및 3D force adapter는 후속 구현/검증 대상이다.
기존 Kepler 전용 KDK guard를 삭제해서 충격 검증을 대신하지 않는다.

다음으로 native CUDA에서 같은 초기상태의 RHS·소산·기하 장부를 검사하고
일반 global KDK와 경계/IC 경로에 연결한다. 그 뒤 같은 모델로
2D contact/KH/RT/Gresho/Sedov/Noh/cylinder/Kepler를 다시 실행한다.
과거 3.339배 GIZMO 실행시간은 2배 목표 미달이며 새 matched benchmark는
아직 수행하지 않았다. /scratch는 사용하지 않는다.


## 2026-11-04 추가 결과와 공개 범위

새 native CUDA 모듈을 별도 compilation unit으로 빌드하고 기존
OpenMP+MPI+CUDA 실행 파일에 부분 재링크했다. 기존 696-byte
restart ABI는 유지했다. 전체 clean rebuild나 다중 랭크 검증은 아니다.
저장된 t=0 및 t=20pi 상태에서 RHS만 평가했고 accepted evolution step은
0이다. 해당 결과를 새 방법의 Kepler 장기 진화로 집계하지 않는다.

새 explicit 점성을 기존 t=20pi 상태에 평가하면 저밀도 외곽의 가까운
쌍 때문에 최소 dt가 8.9503e-8로 작아진다. 따라서 같은 방법의 장기
explicit 계산은 제출하지 않았다. 열화율과 실제 유한 시간의 열 증가를
구분하며, 가열 문제를 해결했다고 보고하지 않는다.

같은 중앙 점성 계수의 frozen-geometry BE 선형계만 별도로 풀었다.
t=20pi, h=4.4515e-5에서 PCG 22회, CPU 8 threads solve 약 0.0342초,
유한 열 증가 8.4096e-7, 에너지 identity 잔차 2.24e-14다.
이는 그래프 구축·전송·기하 갱신을 제외한 대수 시험이다.
일부 쌍은 갱신 끝에서 팽창 상태가 되어 frozen 선형계가 인력을 준다.
따라서 해당 probe는 압축 전용 점성의 완성 구현이 아니다.

후속 1D 후보 unified_implicit은 종점 접근 부호를 active-set으로 풀고
P/2 -> D/2 -> drift+adiabat -> D/2 -> P/2로 갱신한다.
열은 h*v1^T L v1 + 0.5*delta_v^T M delta_v로 계산하며
압력 kick의 속도 변화는 점성 열에 포함하지 않는다.
BE 점성이 활성화된 구간에서는 시간 1차이고, 대칭 배치만으로
전체 시간 2차가 된다고 주장하지 않는다. 무점성 KDK 극한은 별도다.

N=200의 9종 모두 지정 종점에 도달했고 floor, retry, crossing은 0이다.
contact L1rho=1.53e-14, max|P-1|=4.04e-13이다.
Sod L1rho=0.00643825, blast L1rho=0.0604734, Noh L1rho=0.0189555다.
Noh 에너지 잔차 약 8.98e-4는 explicit 후보보다 크므로
시간간격 수렴과 경계 일 장부를 추가 검증해야 한다.
implicit_1d_v1_results.json의 Shu norm은 알려진 불완전한 과거 ref를
사용하므로 정확도 판정에 사용하지 않는다.
새 implicit 후보의 공간/시간 수렴, 공통 다차원 진화,
GIZMO 대비 2배 시간 기준은 아직 미검증이다.

공개 snapshot은 검증 완료한 이전 Kepler 기준선과 실험 중인 새 후보를
분리한다. 원고의 기존 결과를 새 공통 solver의 결과로 대체하지 않는다.


## Native 결합 이후의 최신 상태

UNIFIED_HYDRO_NATIVE_IMPLICIT_GATE_2026-11-04.md에 실제 2D 결합,
1D 시간수렴, 2D .1 유한 열화 gate, 원인 분리와 다음 순서를 기록했다.
그 문서가 위 RHS-only 개발 상태를 갱신한다. 전체 목표는 미완료다.

## 주기 2D 진화 이후의 최신 상태

UNIFIED_HYDRO_PERIODIC_2D_2026-11-04.md에 기존 native 프로그램의
주기 이미지/소유 셀 연결, 보존량 IC 적분, 실제 접촉·전단·Sod 진화를
기록했다. 접촉 여섯 경우와 수정된 전단 세 해상도는 지정 시각에
도달했다. 전단에서 발견한 active-set의 반올림 cycling을 보정변수의
속도차 평가로 수정했다. 새 물리 switch나 epsilon cutoff는 추가하지
않았다. 공통 대수 검사와 1D 9종 회귀는 통과했고 1D 최종장은 동일하다.

축 Sod N=64/128의 최대 횡속도는 0.07087/0.04285다. 수렴 추세가 있어도
정확한 평면 대칭 붕괴 원인은 미해결이다. 범용 shock 통과로 간주하지
않는다. 다음은 최초 횡방향 RHS 진단과 작은 3D 기하·실제 진화 연결,
이후 나머지 2D suite다. 새 결과를 원고의 기존 방법 결과로 바꾸지 않는다.

## 면 기울기와 3D 연결 이후의 최신 상태

UNIFIED_HYDRO_FACE_FIT_3D_2026-11-04.md가 위 주기 2D 단계의 미해결
항목을 갱신한다. 소멸 면의 유한 LSQ 비중이 축 Sod 대칭을 깨는 원인을
확인하고 공통 A/r^(D-1) 면 기하 비중으로 수정했다. 1D 9종은 동일한
최종장을 유지했고 새 실제 2D 14개 시험이 종료 시각에 도달했다.

기존 Voro3D/OpenMP에 같은 압력·fit·점성·BE·KDK를 연결했다. 3D의
유한 clipping 허용치로 생긴 reciprocal 면 불일치도 실제 실패 좌표의
재생으로 분리했다. 명시적 strict double entry와 반올림 범위 안의
endpoint 계약을 사용하며 jitter·floor·유한 면적 cutoff는 없다.
실제 3D 15개 시험과 추가 시간간격 2개 시험이 완료됐다. 32^3 Sod의
에너지 오차는 dt/dt2/dt4에서 8.76e-4/2.84e-4/7.87e-5로 감소했다.

다음은 Gresho/KH, 2D 중복 기하 생성 제거와 회귀·profiling이다.
공통 수식의 기본 차원 연결을 전체 물리 suite 또는 속도 목표 통과로
확대하지 않는다. 최종 목표는 계속 진행 중이다.

## 기하 재사용 및 Gresho/KH 단계

UNIFIED_HYDRO_ROTATING_SHEAR_2026-11-04.md에 후속 구현과 결과를 기록했다.
세 실제 native A/B 시험에서 최종 물리장은 동일하고 중복 기하 생성을
제거했다. Gresho 64²/128²는 t=1, KH 64×128/128×256은 t=2에 도달했다.
KH 256×512는 공개 시점에도 진행 중이며 종점 완료로 집계하지 않는다.
Gresho 시간간격 축소는 에너지 오차를 줄이지만 속도 오차와 소산은 남는다.
동결 음향 Hessian 양의 부분 검사는 전체 기하 Hessian을 포함하지 않으며
새 시간간격 조건을 채택한 것이 아니다. 새 closure의 Kepler 중심 가열도
미해결이다. 원고는 기존 검증 기준선을 유지하며 전체 목표는 진행 중이다.

## 2026-11-05: KH 종점 및 전체 압력 곡률 검사

UNIFIED_HYDRO_ACOUSTIC_CFL_2026-11-05.md가 위 단계의 진행 상태를 갱신한다.
KH 256×512도 t=2에 도달해 세 해상도 종점과 셀별 고정질량을 확인했다.
전체 등엔트로피 압력 Hessian 방향미분에서 기존 CFL이 빠른 압력 mode를
놓치는 상태를 확인했다. 차원 공통 block 상한의 대수 검사는 통과했고,
2D 선택 옵션으로 Gresho t=1까지 비교했다. 에너지 오차는 약 5e-6으로
줄지만 속도 오차와 약 1/3의 운동에너지 손실은 남는다. 스텝 수도 늘어
기본값으로 승격하지 않았다. 동일 기하의 정확 속도 대조는 이미 발생한
비원운동이 소산을 크게 증폭함을 보인다. 다음은 체적 일·압력 잡음·공간
힘 오차의 분리이며, 전체 물리 및 GIZMO 속도 목표는 미완료다.

## 2026-11-05: 물질 체적 수송과 국소 가중치 구속 후보

UNIFIED_HYDRO_MATERIAL_VOLUME_2026-11-05.md가 위 원인 분리의 결과다.
실제기하의 K:G는 비압축성 affine 전단에서도 국소 체적 변화를 만든다.
선형 압력·회전·전역 에너지 검사의 통과와 구분해야 한다. 정확 원운동
geometry-only 대조에서도 실용 해상도의 밀도 왜곡을 확인했다.
새 C=B+LR 후보는 작은 1D/2D/3D 대수·가변 strain·비영 가중치 기하
미분 검사를 통과했지만, 첫 native CPU 계수 풀이가 제한된 추가비용
안에서 수렴하지 못했다. 추가 OpenMP 교차검사에서도 ACTIVE는 반복당
비용을 줄이지만 256회 안의 수렴 조건을 만족하지 못했다. 현재 CPU
후보 적용은 보류한다. 기존 native를 같은 t=0.5 체크포인트에서 .51까지
비교하면 ACTIVE는 phase 시간 기준 1.486배 빠르고 최종 물리장은 완전히
동일하다. 생산 물리 본체는 변경하지 않았다. 물리 정확도·시간적분·전체
suite·GIZMO 속도 목표는 미완료다.

## 2026-11-05: 계수 잔차와 공간 연산자 응답의 분리

UNIFIED_HYDRO_COEFFICIENT_ERROR_2026-11-05.md가 다음 진단 묶음의 결과다.
같은 cap 256, tol 1e-10, 공통 시간 예산 안에서 두 계수를 한 번 저장했다.
둘 다 미수렴 상태를 유지한다. 독립 기준해와의 차이를 부피 일·힘·면 변위로
전파하고 2D/3D 동일 잔차 모드 대조 14개를 확인했다. tree 기반 상계는
너무 느슨하므로 새 허용오차를 임의로 선택하지 않았다. 더 중요한 점은 정확한
계수로도 기존 상태의 속도 편차에 대한 체적 응답이 1.681→6.556, 압력 편차에
대한 힘밀도 응답이 4.923→21.890으로 커진다는 점이다. 고정 상태의 반응이므로
불안정성 증명은 아니지만, 계수 solver를 더 비싸게 만드는 것보다 공간 연산자
자체의 반응을 먼저 평가해야 한다. 생산 물리와 원고 결과는 유지한다.

## 2026-11-05: 장파장 일관성과 변형 기하의 실용 stiffness 분리

UNIFIED_HYDRO_MODE_RESPONSE_2026-11-05.md가 다음 mode 검사 묶음의 결과다.
1D에서 C=B 환원은 정확하다. 불규칙 2D/3D 주기 symbol의 장파장 음속은
맞고, 고정 파장의 횡방향 가짜 압축은 반복패턴 refinement에서 감소한다.
기존 native 기하에 정확한 파장·진폭의 섭동을 넣어도 매끄러운 부피 응답은
개선된다. 그러나 압력 힘은 악화하며, 변형된 128² 기하의 양의 압력 block
최고 주파수는 B의 약 49.56배다. 독립 face 재조립과 제한된 power 대조로
확인했다. 이 값은 성장률이나 전체 KDK 안정한계가 아니지만 속도 목표에
중대한 경고다. 후보 계수 최적화와 새 장기 진화는 진행하지 않는다.
자기일관·단순성·물리 정확도·전체 2D suite·GIZMO 처리량 목표는 미완료다.

## 2026-11-05: 균일압력 평형 대조와 현 보정 후보 종료

UNIFIED_HYDRO_EQUILIBRIUM_RETIREMENT_2026-11-05.md가 마지막 최소 대조다.
같은 실제 기하와 고정질량에 v=0, 균일 P0를 두면 배경 기하 미분항이
소거되어 양의 acoustic block이 정확한 무점성 압력 선형화가 된다.
B의 엄밀한 행렬 norm 부등식과 기존 C Rayleigh 하한을 평가하면 주파수
비의 보수적 수치 하한도 약 36.02다. interval 인증과는 구분한다.
새 eigen 계산·계수 풀이·PDE 없이 이 차이를 확인했다. Astra의 종료 조건에
따라 현 psi_i G_i 구성은 생산후보에서 제외했다. 소스와 결과는 보존한다.
전체 C=B+LR 계열의 불가능성으로 확대하지 않는다. 새 후보는 정상 부피의
변형 기하에서도 acoustic norm이 폭증하지 않는 구조와 낮은 비용을 함께
확보해야 한다. 원래 물리 정확도·차원 공통화·2D 재검증·속도 목표는 유지한다.

## 2026-11-05: midpoint divergence의 압력 gate와 위치 표현 결함

UNIFIED_HYDRO_MIDPOINT_PRESSURE_2026-11-05.md가 후속 묶음 결과다.
과거 재구성 식의 (Gi-Gj)r/4 항을 제거한 D는 1D W=0에서 B로 환원하고
1D/2D/3D의 부피·adjoint·cut ledger를 통과했다. 그러나 같은 불규칙 unit
pattern을 두 배 정밀하게 반복해도 셀 평균 선형 압력의 힘밀도 오차가 줄지
않았다. 모든 6개 차원/축 검사가 불통과하여 Bloch·native stiffness·weight
풀이·후보 진화로 진행하지 않았다. 1D에서는 기존 B에도 같은 평균압력 결함이
있다. 실제 셀 중심과 생성점의 차이로 오차식을 유도하고 독립 C 및 기존
native pressure helper로 확인했다. 정지 contact의 초기화 원인과는 별개다.
다음은 보존 상태·열역학·생성점 이동의 표현을 함께 정하는 단계다. 생산
코드와 원고 기준선은 유지하며 전체 목표는 미완료다.
