# Astra — 비용 실패 뒤의 단일 compact-space 결정

2026-10-22. SHARED_CONTEXT 전체와 실제 paired step contract, selected cost, 네 stage quadrature difference, endpoint, v3 restart를 읽었다. MATERIAL_P3_RESULTS, MATERIAL_SHARED_RESULTS와 geometric-pressure 결과를 제한적으로 대조했다. 원본 LagEunha HEAD8e25dfe 및 연구 ROOT 분리 계약을 유지한다. 계산·폴링·소스 변경 없이 이 새 문서만 작성한다.

### 동시 계획 파일에 대한 실행 범위 명시

운전자가 별도 `ASTRA_CONSERVATIVE_COMPACT_REZONING_PLAN_20261022.md`의 존재를 알렸다. 이 파일은 본 평가자가 이번 turn에 작성한 문서가 아니며 작성 주체를 확인하지 않았다. 이를 삭제하거나 덮어쓰지 않는다. **이번 요청에 대한 선택·실행 계약은 이 문서 한 개다.** step5를 straight compact P3/Bernstein rho/e로 옮기는 nonlinear per-cell p/L/K projection/rezoning은 이번 평가에서 승인하지 않았다. 그 대안과 이 문서의 fresh compact branch를 병행 launch하지 않는다. 이전 거대공간의 long continuation 결정은 본 문서의 비용 실패 종료로 대체된다.

## 1. 완료 판정: 정확도 채택, 현 고정공간 경로 종료

accepted step5는 t=.00011263377801805954이며 실제 low/high 네 stage와 endpoint, v3 restart가 통과했다. stage1의 `high_k1_force_and_solve_reused=true`와 step-level `cached_initial_k1_reused=false`는 서로 다른 캐시 항목이다. 이를 high k1 재계산으로 오독하지 않는다. 에너지/U -1.409e-12, 각운동/2K 1.263e-17, thermal/U -2.267e-18과 원래 mass 보존을 채택한다.

그러나 selected-only130.547초/step은 paired reference 시간을 제외해도 기존225초 대비1.72배다. 원래 node dt bound를 사용한 낙관적 t1 비용은4.368일이다. 따라서 기존10배/24시간 비용 관문에 실패했다. **10.8M scalar DOF 공간의 추가 최적화·guard 개발·장시간 실행을 여기서 끝낸다.** 마지막 checkpoint, positive rules, GPU DD mass와 모든 증거는 보존한다. 이 실패는 양성 cubature의 정확도 실패나 물리 모델의 불가능성 증명이 아니다.

## 2. 다음 하나의 후보: core-protected nested P3 + 전역 초기 residual 두 개

기존 원본 CUDA/MPI/OpenMP 실행 구조와 geometric-pressure/entropy/hybrid 계약을 유지하면서 **trial space 자체를 줄이는 한 후보**를 선택한다. 생성점·물리 source cell·연속 IC·중력·eps·P0를 바꾸지 않는다.

- 71,080개 source cell과 canonical reference geometry는 그대로 둔다.
- 원래 fan triangle의 h/2 내부 분할은 **초기 source polygon이 r<.5와 만나는 cell에서만 유지**한다. 밖에서는 저장된 원래 parent fan의 완전 P3로 제한한다. 이 판정은 실제 polygon과 circle의 교차/거리로 고정하며 pressure 결과를 보고 이동시키지 않는다. 신규 radial IC refinement가 아니다.
- coarse/fine 내부 경계는 coarse P3 trace를 fine trace에 정확한 interpolation constraint로 주입한다. hanging node를 독립 운동시키거나 좌표를 평균하지 않는다. 제한은 x와 v에 동일하다.
- 각 실제 connected component에 두 scalar 함수

      R_k(X)=v_initial,k(X)-I_compact,P3(v_initial,k)(X),  k=x,y

  만 추가한다. 각각의 x/v 계수는 자유로운 동적 DOF다. 시간마다 초기 궤도를 대입하거나 미래 원운동 함수를 추가하지 않는다. 대규모 edge-residual/local-residual 자유도는 제거한다. shared canonical interpolation으로 R_k는 연속이고 안정한 `q-I3q` 평가를 사용한다.

과거 coarse P3 max chi=.005376이 .002 관문에 실패한 사실을 알고 있다. 따라서 **전 영역 uniform coarse P3 복귀는 선택하지 않았다.** core의 기존 polynomial h/2 해상도를 보호한다. 그렇더라도 두 전역 enrichment가 초기 v를 정확히 표현한다는 사실만으로 이후 acceleration/pressure curvature나 감긴 map 정확도를 보장하지 않는다. 이 한 후보의 실제 full initial 관문이 실패하면 core 범위/order를 다시 조정하는 탐색을 하지 않고 종료한다.

구조적 규모 목표는 대략1.5–2.0M scalar DOF이며 상한2.5M이다. 이는 아직 실제 count가 아니다. 저장 parent/child와 interface constraints로 정확한 owned/ghost/DOF count를 **구축 전에** 산출한다. parent 계보가 없으면 저장 정확 reference 연결에서만 복원하며 새 Voronoi 계산을 하지 않는다. parent 관계가 증명되지 않은 triangle을 임의 병합하지 않는다. 상한을 못 맞추면 full 실행에 들어가지 않는다.

## 3. 이전 계산을 재사용하는 정확한 선형대수

compact coefficients c와 기존 fine coefficients q 사이의 주입 T를 구성한다:

    q=T c,    N_compact=N_fine T,
    M_compact=T^T M_fine T,
    F_compact(q)=T^T F_fine(q).

polynomial restriction은 저장 reference의 정확 interpolation 관계다. 전역 residual column은 기존에 정확히 표현한 초기 v의 coefficient에서 compact interpolant의 fine coefficient를 뺀 것이다. 기존 normalizations를 무시한 단순 좌표 차감은 안 된다. coefficient 관계와 analytic canonical 함수 관계의 rounding 오차를 기록한다.

M_compact는 저장된 conditioned binary80 M의 congruence로 만든다. 초기 force도 저장 초기 F 또는 그 authoritative element blocks를 T^T로 축약한다. 저장 acceleration을 단순 평균해 새 acceleration이라고 하지 않는다. 두 global residual column 때문에 생기는 결합은 sparse polynomial block + rank2 border로 표현하고 작은 Schur complement를 사용한다. 전역 dense matrix를 만들지 않는다. 원래 binary80 true residual1e-16, DD GPU local arithmetic을 유지한다.

이 재사용은 **같은 x/v 함수에서의 weak operator restriction**에 대해 정확하다. 새 공간에서 독립적으로 진화한 다음 시각의 force를 과거 fine force에서 꺼낼 수 있다는 뜻은 아니다. 새로운 K/M metric은 새 공간에 대한 M_compact이며 오래된 M과 bitwise 같다는 주장도 하지 않는다.

기존 positive subset은 fine-space mass/derivative-product constraints를 지켰으므로 주입 T 아래 compact mass/volume 제약의 기준 자료로 재사용한다. coarse parent의 child subset을 합친 뒤, 동일한 positive moment elimination으로 compact constraints에 맞춰 **한 번만** 다시 줄일 수 있다. 새 sampling/차수/precision 탐색은 없다. 원래20억 점 또는5시간 전처리를 그대로 재생하지 않는다. 원점 rho/logistic/eps1e-6 구조를 삭제하지 않는다.

M_compact, determinant products, constant measure, 생성 thermal function을 보존하는 실제 residual과 cold work/torque contraction을 확인한다. positivity가 unresolved이거나 원래 허용치를 못 맞추는 parent는 실패다. weight 음수 clamp, rank 삭제, tolerance 변경은 없다. 실제 총 점수 목표는 기존 selected rule의1/4 이하이며 이 조건도 전체 집계로 확인한다. 점 수만 줄고 연산비가 줄지 않을 수 있으므로 최종 판정은 native 시간이다.

## 4. Step5를 압축해 이어간다고 주장하지 않는다

이미 진화한 fine-space step5의 일반 x/v는 compact subspace에 속하지 않는다. 따라서 이를 평균/LS projection한 뒤 잃은 K를 U로 넣거나 별도 spin reservoir에 보관하는 restart를 만들지 않는다. 기존 accepted state는 영구 보존한다.

새 후보는 **동일 연속 IC에서 t=0으로 시작하는 별도 공간 이산화 branch**다. 성공한 동일 알고리즘 구간을 단순 재생하는 것이 아니라, old space에 없는 cost viability를 얻기 위해 자유도를 바꾸는 필요한 비교다. 초기 위치는 같은 reference map, 초기 velocity는 위 두 residual로 원래 연속장을 표현한다. physical U 및 원장 연속장 E를 유지한다. point-generator legacy K와 연속 적분 K 차이를 열로 보상하지 않는다.

같은 source grouping, 원래 quadrature 정의와 generating jhat0/physical U0 함수를 유지할 수 있는 범위는 명확하다. 이는 새 물리 cell로 remap하지 않았기 때문이다. 생성 measure를 실제 새 양의 rule에서 평가하고 원래 생성 장부와 discrepancy를 남긴다. old step5의 current J를 새 jhat0로 삼지 않는다. compact initial x가 원래 reference map과 같고 determinant-product/constant-measure constraints가 유지될 때만 기존 생성 entropy 함수가 동일한 물리 초기 상태를 정의한다. 이를 실제 초기 U/P/volume로 확인하며 차이를 reset으로 없애지 않는다.

## 5. 한 실행 묶음: 초기 restriction부터 실제 연속 trajectory까지

### 고정 사전 관문

하나의 compact topology/space만 구성한다. 실제 DOF, selected points, M/preconditioner/halo bytes, GPU working set을 먼저 기록한다. 저장 M/F의 restriction으로 초기 acceleration을 구하고, 이미 보존한 원래 reference 함수/점 및 기존 witness에서 초기 pressure curvature와 velocity/energy 일치를 확인한다. 새 primitive별 토이 시험이나 전체 initial mass 적분은 하지 않는다.

기존 max chi<=.002 및 모든 셀 평균<=.001, arbitrary-constant pressure nullspace, 원래 physical U/K/M/torque 관문을 그대로 적용한다. 새 공간의 local reaction/virtual-work는 그 공간에 맞게 조립해야 하며 이전 공간의 reaction을 그대로 복사하지 않는다. 실제 reference points 이외 연속 supremum을 증명했다고 표시하지 않는다.

이 후보가 공간 정확도 또는 비용 규모 관문에 실패하면 원인 source와 실제 count를 보존하고 묶음을 끝낸다. core반경, h, degree, enrichment 수를 바꿔서 작은 시험을 반복하지 않는다.

### 실제 full-disk 실행

통과한 같은 native mode에서 fresh compact branch를 시작한다. 첫 step은 이어지는 trajectory의 일부다. 첫 full RHS의 CPU/DD spot 및 전 셀 work/torque, true residual, 초기/accepted restart를 확인하고 그 구간을 재실행하지 않는다. 이후 전체 고차 참조를 매 step 동시 실행하지 않는다.

actual native minimum-dt source/branch를 저장한다. compact 공간이라고 원점 회전 제한이 사라지지 않는다. 사전 비용은 현재 sustainable h≈.000345838을 **개선 없이** 사용하고, 실제 더 작은 h가 나오면 그 값으로 갱신한다. factor7을 이유로 h를 키우거나 low point 최대값만으로 CFL을 완화하지 않는다. 기존 timestep 법칙과 종료 조건을 유지한다.

기준 진행 시각은 .1/.2/1/6/8/10이다. .1/.2는 옛 고갈과 첫 가열을 확인하는 중간 출력이며 종료 목표가 아니다. 같은 Eulerian radial bins의 density/P/phase/velocity, U/K 및 압축 Jensen 분해, 실제 resolved L/I/Krot, 총 E/L와 외압 work, map quality를 남긴다. 고정 entropy로 entropy residual이 작다는 것만으로 물리 정확도를 선언하지 않는다.

기존보다 느슨한 물리 gate는 허용하지 않는다. 새 branch의 공개된 continuation engineering ceiling으로 r<8 각 기존 영역의 density/P 상대 RMS1%, velocity 상대 RMS1e-3, phase RMS.2도 이내를 사용하고, 기존 더 엄격한 기준이 있으면 그것이 우선한다. 이는 paper validation 기준이나 오류증명이 아니라 대규모 잘못된 경로를 조기에 끝내기 위한 상한이다. 실제 오류가 초기·quadrature 오차보다 큰 경우 그 증가를 반드시 기록한다. outer material boundary 연구 차이도 별도 표기한다.

### 실제 비용과 고정 wall 한도

실제 DOF/점수/메모리 및 첫 정상 연속 step들의 비용으로 t1뿐 아니라 t6/t10을 보고한다. setup·reference-only 시간과 steady time을 분리하고 같은 device에서 평가한다. 낙관적인 큰 h만으로 통과시키지 않는다.

- 첫 batch 최대24시간. 처음 실측에서 t1>24시간 또는 t10>7일이면 장기 실행으로 승격하지 않는다. 이 경우 마지막 accepted state에서 비용 실패로 종료한다.
- 통과하면 같은 trajectory를24시간 단위 durable restart로 이어 최대 누적7일 또는 t10까지 진행할 수 있다. 실제 dt/비용이 악화되어 t10 예산을 넘거나 물리/mesh 관문이 깨지면 종료한다. 작은 구간마다 Astra 승인 요청은 필요 없다.
- t6/8/10에 도달해도 그 자체가 해결 판정은 아니다. 밀도/압력·궤도·E/L·경계 의존성을 함께 평가한다. 논문 최장 약1000공전은 별도 정규화와 장기 수명 문제이며7일 budget으로 해결됐다고 표시하지 않는다.

## 7. Material 수명과 W0: 이번 축소가 해결하지 않는 것

두 global 초기 residual은 초기 함수 보존 장치다. 시간에 따라 감기는 원운동 map을 정확하게 저장하는 장치가 아니다. core cond(F)가 t8에10^5를 넘는 정확 전단 증거를 그대로 유지한다. 전체 enriched continuous-map guard의 개념적 계약은 이전 shear-lifetime 문서와 같지만, 종료한10.8M 공간에 대한 새 구현을 하지 않고 이 compact space가 초기/비용 관문을 통과했을 때만 실제 필요한 guard를 연결한다. source/vertex/외부 boundary의 unresolved와 실제 fold를 구별한다.

실제 감김 때문에 이 branch가 t6 전에 막히면 마지막 state에서 중단한다. dt/order scan은 공간 수명 해결책이 아니다. 그때 필요한 것은 새 reference로 물리 full field·mass·U·resolved K/L를 보존해 옮기고 새 geometric thermal functional과 연결하는 실제 재표현이다. 현재 compacting은 **같은 reference topology의 trial-space 제한**이므로 그 remap 문제를 해결한 것이 아니다.

W0에서 개별 material mass는 고정되지만 current material 경계가 Voronoi라는 뜻은 아니다. weight-linked mass relaxation/동적 Laguerre/Voronoi 환원은 아직 미완료다. 원래 원본 loop/GPU/MPI 사용과 원래 geometric model 동일성은 구별한다. 궁극적 사용자 요구를 충족했다고 선언하려면 이 간극도 닫아야 한다.

## 8. 원래 t8 CUDA baseline과 회전 증거의 사용

기존 t8/t10 snapshot은 원래 모델의 실패 지점·가열 분포·ring 정확도와 실행 비용을 비교하는 실제 기준으로 보존한다. 그것을 새로운 material coefficients로 변환하면서 이미 발생한 열을 제거하거나, spin 에너지를 만들어 E를 재분류하는 restart는 이번에 하지 않는다. 원래 생성점 구조와 native 실행 환경은 재사용하지만 legacy pressure/work 혼합 및 압력비100 branch를 수정 없이 다시 장기 실행하는 것은 새로운 해결 후보가 아니다.

이번 compact 후보의 rotation은 동일 velocity field의 antisymmetric gradient와 mass-centred internal angular momentum이며 실제 stress-work와 함께 진화한다. independent kinetic tensor/reservoir 또는 damping을 추가하지 않는다. shock entropy는 별도 미완료로 남으며 smooth 단열식을 충격 해법으로 포장하지 않는다.

## 9. 최종 실행 결정

현5step trajectory는 정확도 PASS/비용 FAIL로 종료·보존한다. 다음 단 하나의 고정 compact 후보를 위 규모·initial curvature·energy/thermal·실제 비용 관문으로 판단하고, 통과하면 **실제 연속 full-disk t1→6/8/10**을 진행한다. 실패 시 그 한 묶음을 종료하며 또 다른 작은 공간/parameter 변형을 자동 생성하지 않는다. 코딩·루틴 실행에 사용자 재승인을 묻지 않는다. 결과 보고는 compact space의 실제 규모, 초기/물리/비용 판정, 도달한 physical time과 종료 이유를 한 번에 묶는다.
