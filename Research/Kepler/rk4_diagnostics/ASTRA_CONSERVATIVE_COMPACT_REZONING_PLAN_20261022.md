# Astra — 비용 실패 이후 compact P3 보존 재격자 한 묶음

2026-10-22. SHARED_CONTEXT 전체, 원본 HEAD8e25dfe, ASTRA_LONG_EVOLUTION_COST_DECISION, 실제 step5/endpoint/restart/cost 결과를 확인했다. 계산·폴링·코드 수정 없이 작성한다. 운전자가 구현·실행한다.

## 1. 종료 판정

양의 rule 전처리와 실제 전환 step의 정확도·보존·restart 성공은 채택한다. accepted t=0.00011263377801805954와 모든 원래 payload를 보존한다. 같은 입력의 성공 점/구간을 반복하지 않는다.

그러나 선택 경로130.547초/step, 낙관적인 dt에서도 t1 남은4.3685일이며 향후 guard 비용조차 제외했다. 같은 Ada의 유리한 비교도3.387배에 그친다. 이는 H200 지속성능 비교가 아니지만 채택한 두 비용 대안 모두 실패했다는 결론에는 충분하다. **동일10.8M DOF의 추가 최적화·차수/허용치/dt/rank 탐색 및 장기 재개를 종료한다.** 기존 전처리 시간을 회수하려고 비현실적인 장기 작업을 제출하지 않는다.

## 2. 다음 후보 하나: compact P3 material rezoning

다음 한 묶음은 **현재 accepted 물리장을 새 power-cell의 straight reference fan에 보존 재표현하고, 초기해 전용 edge/bubble residual이 없는 shared P3 공간으로 규모를 줄이는 실제 native transaction**이다. 이후에는 같은 pressure/work compatible material evolution을 사용한다. 시간에 따라 누적되는 material 전단은 동일 보존 rezoning으로 reference를 갱신하는 경로를 택한다.

이것은 mean compression이나 단순 계수 삭제가 아니다. 저차 공간의 근사 오차가 생기는 새로운 numerical representation이다. exact IC를 다시 대입하거나 현재 속도를 초기 원궤도로 바꾸지 않는다. 보존을 위해 생긴 부족한 kinetic energy를 열로 보내지도 않는다. 이 후보가 실현가능하다는 보장은 없으며, 아래 동시 제약을 실제 저장 상태에서 만족할 때만 전환한다.

선택 공간은 고정한다:

- 각 실제 target convex polygon을 한 번 centroid fan으로 나눈 straight triangle mesh. 기존 모든 작은 reference triangle과 h/2 subdivision을 그대로 승계하지 않는다.
- 위치/속도: conforming P3, canonical vertex/edge/triangle-interior DOF. **IC 전용 residual basis 없음**. 처음 위치는 이 새 reference의 identity map이다.
- reference density와 reference-volume physical internal-energy density: 각 fan triangle의 비음수 Bernstein P3 계수로 표현한다. 속도와 달리 density/thermal은 cell 경계에서 불연속일 수 있다. 양수 상태에 맞지 않는 작은 음수 계수 clamp는 없다.
- 초기 정확해의 cusp/전단을 삭제하지 않는다. 그 원장을 현재 accepted field로부터 얼마나 근사했는지를 별도 오차 장부로 판정한다. 원장 profile 변경이 아니다.

공간 규모는 실제 생성한 DOF/삼각형 수로 산출한다. P3 mass의 polynomial 부분은 rho P3 × N P3 × N P3이므로 degree9의 양의 고정 적분으로 정확하게 조립할 수 있다. 기존 nonpolynomial 초기해 residual과 64×28 적분을 production 공간에 남겨 놓고 compact라고 부르지 않는다. 새 evolved heat/gravity는 여전히 비다항식이므로 별도 실제 적분오차 검증이 필요하다. polynomial exactness가 그 검증을 대신하지 않는다.

## 3. 실제 물리 전달 계약

source는 step5의 full field, fixed M, 생성 entropy와 현재 geometric-pressure context다. target별 실제 overlap에서 다음 적분을 각각 얻는다:

    m, q=∫rho x, Q=∫rho xxᵀ,
    p=∫rho v, L=∫rho(x cross v),
    K=1/2∫rho|v|², U=∫e, Phi=∫rho Phi_gravity.

기존12 moments와 겹치는 항은 같은 canonical packet을 source-outgoing/target-incoming 양쪽에 사용한다. 이 moment 목록만으로 old velocity를 재구성하지 않는다. 실제 old full function을 projection RHS와 오차 적분에도 사용한다. 이번 geometry는 현재 step5의 새 geometry이므로 새로 필요한 overlap 계산은 수행하되 과거 finiteCONTROL의 완료 packet을 현재 상태라고 재사용하지 않는다.

### 밀도와 열

positive Bernstein density를 현재 source rho에 대한 constrained approximation으로 구한다. target의 m/q/Q 및 Phi를 원래 conservation tolerance 안에서 유지한다. physical thermal field도 실제 source e에서 positive constrained approximation하며 target U를 별도로 보존한다. 두 선형 양성 제약문제가 infeasible이면 target 실패다. density/thermal을 유리한 평균값으로 리셋하지 않는다.

새 reference에서 jhat_new,0와 선택한 실제 e_new,0를 사용해 entropy 생성자료를 **한 번** 정의한다. 이것은 현재 U를 초기 P0로 되돌리는 heat reset이 아니다. 기존/새 e와 pressure의 공간적 차이 및 ∫rho s의 차이를 numerical representation error로 기록한다. 새 delta가 달라졌다는 이유만으로 old entropy 함수를 그대로 넣어 U/P를 바꾸지 않는다. 전환 후 생성자료는 다시 immutable이다.

### 속도·회전·kinetic energy

새 density의 consistent mass metric Mnew에 대해 가장 가까운 shared P3 속도를 구하되, 각 target의 p/L/K를 동시에 제약한다:

    minimize 1/2 ∫rho_new |v_new-v_old|²
    subject to p_new=p_old, L_new=L_old, K_new=K_old.

rho가 바뀌므로 목적함수와 보존 제약에 사용하는 measure를 명시하고, 원래 measure에서도 velocity error를 독립 집계한다. 원래 U와 Phi도 보존하므로 E=K+U+Phi 장부를 별도로 닫을 수 있다.

먼저 선형 p/L 제약하의 minimum K를 계산한다. K_old<Kmin이면 공간/밀도/제약이 실현불가능하다. K_old≥Kmin만으로 충분하다고 선언하지 않는다. shared DOF로 연결된 여러 셀의 K 제약과 field-error 조건도 동시에 풀어야 한다. sparse constrained solve를 사용하고 global dense Schur는 만들지 않는다.

K 부족을 새 고주파 bubble 한 개에 넣거나 독립 scalar reservoir로 저장하지 않는다. 기존 compact 공간 안에서 old field와 가장 가까운 해만 후보이며, 방향이 정해지지 않는 nullspace에 임의 회전을 넣지 않는다. 해를 찾더라도 허용 field 오차를 넘으면 실패다. 이런 제약 투영은 physical dissipation이 아니라 수치 재표현이며, 정확 보존만으로 물리적으로 무해하다고 주장하지 않는다.

resolved intrinsic angular momentum은 L−x_cm cross p, kinetic energy는 실제 shared velocity와 Mnew로 계산한다. 별도 spin state나 damping이 없다. 같은 법선 pressure traction의 비중심 torque와 compatible work를 새 공간에도 유지한다.

## 4. 가중치·질량과 material 수명

generator는 현재 source map에서 운반된 원래 generator 위치를 사용한다. 분석적 궤도나 mesh quality force로 이동시키지 않는다. 필요시 승인된 screened finite law

    m_i(W) + mu_i W_i = m_i,reference

를 실제 old-field overlap mass로 푼다. mu와 기준 질량은 immutable manifest에 고정한다. target mass를 projection이 임의로 선택하지 않는다. 전체 질량과 외부 domain coverage, canonical plane/event와 weight gauge도 함께 닫아야 한다.

W=0일 때 m_i=m_i,reference는 **이 방정식이 충족되는 경우의 조건부 환원**이다. 일반 moving Voronoi overlap에서 W=0만 넣으면 개별 질량이 고정된다는 주장은 금지한다. 방정식의 geometry/positivity가 실패하면 μ 변경이나 entropy 조정으로 숨기지 않는다.

새 reference를 identity로 다시 만드는 것은 좌표 재설정이고 물리 generator steering이 아니다. x/v/rho/e의 실제 현재장을 위 조건으로 전달해야 한다. 이 재격자가 통과하면 누적 material shear를 영구히 같은 고차 map으로 담아야 하는 문제를 피할 수 있다. 그러나 반복 재격자의 오차·비용·weight feasibility는 아직 미증명이다. 이번 한 번의 성공으로 무한 수명이라고 하지 않는다.

## 5. 실행할 한 묶음과 중단선

이 묶음에는 작은 guard 선행실험을 따로 두지 않는다. 현재 전체 accepted state 한 개를 대상으로 다음 transaction을 구성한다.

1. **Native full-state import와 후보 생성.** 저장 step5를 직접 읽고 SHA를 확인한다. 현재 source map의 단사성/coverage가 실제 전달에 충분히 입증되지 않으면 조립을 중단한다. sampled J만으로 이 전제를 통과시키지 않는다. 이 검사는 전환의 전제이며 새로운 일반 proof 탐색 프로젝트로 확대하지 않는다.
2. **한 번의 전체 compact feasibility/전달.** 실제 target registry, weight-linked mass, positive rho/e fit, shared velocity 보존 투영을 같은 묶음으로 조립한다. 고정 공간·정해진 solver와 유한 반복 예산만 사용한다. failed target만 상세 witness를 남기고 전체 transaction은 미commit한다. 실패하면 자동 P상승/IC smoothing/추가 residual basis/열 전환은 없다.
3. **오차 예산 분리.** 기존 conserved-moment/thermal 전달 tau는 그대로다. 그것과 별도로 velocity/density/pressure field 근사오차와 entropy discrepancy를 평가한다. 기존 T=10의 normalized 물리오차1e-3 예산을 근사공간 변경 때문에 새로 리셋하지 않는다. quadrature, 시간, 반복 remap의 누적 장부와 분리하여 남은 예산을 사전 배정하고, 보존 오차가 작다는 이유로 field 오차를 면제하지 않는다. 이전 같은-function interval bound를 다른 함수의 근사오차 증명으로 전용하지 않는다.
4. **통과한 경우만 다음 실제 RK4 step.** 원래 accepted time에서 compact state를 commit하고 새 smooth compatible operator로 다음 한 step을 실행한다. 과거 구간은 replay하지 않는다. 고차 원장 경로의 실제 force/pressure 응답과 비교하는 데 필요한 같은 입력 평가를 이 전환 관문 안에서 한 번 수행한다. 과거 초기/곡률/atlas 계산은 반복하지 않는다. positivity, local E/L/work, physical pressure accuracy와 restart를 검사한다.
5. **비용 관문.** 실제 compact step+예상 필수 remap/guard의 amortized 비용으로 t1≤24시간을 요구한다. remap을 공짜로 놓은 예측은 쓰지 않는다. 새 mesh reset이 timestep을 키운다는 가정도 실제 제한 계산 전에는 쓰지 않는다. 비용 또는 정확도 실패 시 compact 상태를 production으로 승격하지 않는다.

대규모 overlap quadrature를 시작하기 전에 source coverage와 compact 공간의 선형 제약/필요 DOF 및 memory, projection solve 예산을 산출한다. 임의 target 수를 늘려 toy loop를 만들지 않는다. 이번 전처리·전달·조건부 한 step은 wall24시간의 한 번 bounded native 배치로 하며, 예산 끝에는 원본 step5 또는 이미 유효하게 채택된 새 state를 보존한다. 부족한 시간 때문에 미완료 stage를 accepted로 표시하지 않는다.

## 6. 성공·실패 후의 단일 경로

보존 투영이 불가능하거나 cold field 오차가 남은 예산을 넘으면 **compact P3 가설 기각**으로 종료한다. 그때 원래10.8M 공간으로 장기 복귀하거나 차수 메뉴를 돌리지 않는다. 실패 witness는 어떤 conserved constraints/양성/필요 kinetic content/field 오차가 공간 축소를 막는지를 명시해야 한다. 이는 실제 물리해가 틀렸다는 증거와 다르다.

정확도와 전체 비용을 통과하면 같은 native CUDA+MPI/OpenMP 경로를 t1까지 이어가며, remap trigger는 현재 deformation/거리/condition의 사전 고정 조건과 실제 전달가능성에 따른다. 열 상태를 refit하는 자동 repair는 없다. 이후 t6/8/10에 대해 실제 rho/P, orbit/phase, E/L/physical U와 remap 누적오차·질량/weight·비용을 평가한다. 현재 Pext 경계 차이, smooth entropy mode의 shock 미완료, 원래 periodic 환원은 계속 명시한다.

GIZMO 최장 목표는 `KEPLER_GIZMO_COMPLETION_CRITERIA_20261022.md`에 따른다. code time600만으로 약1000공전을 선언하지 않으며 GM/eps/reference radius와 T(r)을 고정해 code time과 orbit count를 함께 기록해야 한다. 새 연구 경계로 원래 논문과 직접 동일하다는 주장은 불가하다.

## 7. 운전자에게 넘기는 결정

다음 구현은 **실제 step5 → compact P3 보존 rezoning → 조건부 다음 실제 native step** 하나다. GPU/DD/MPI 기반, source-context 및 transfer 회계의 기존 성공 코드는 재사용하되 오래된 함수·packet을 현재 상태로 속이지 않는다. 추가 미세 평가 호출 없이 전체 feasibility/정확도/비용 결과 또는 구체적 구조적 실패를 한 번에 보고한다. 목표는 여전히 t6/8/10와 최장 공전이며, compact 전달의 PASS 자체가 완료 조건은 아니다.
