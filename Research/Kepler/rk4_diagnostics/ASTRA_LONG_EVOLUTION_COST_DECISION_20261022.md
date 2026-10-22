# Astra — 장기 진화 비용에 따른 전략 결정

2026-10-22. SHARED_CONTEXT 전체, RESULT_FULL_EVOLUTION, handover의 최신 절, 기존 계획과 실제 evolved/cache/mass host 및 CUDA 연산을 읽었다. 원본 identity는 운전자가 재확인한 LagEunha HEAD8e25dfe다. 계산·폴링·원본 수정 없이 작성한다. 현재 Ada 작업에는 개입하지 않는다.

## 1. 판정

실제 full RK4/restart와 DD mass continuation은 유효한 진전이다. immutable 감사로 분리한 archive-only wrapper 실패 때문에 성공한 구간을 재실행하지 않는다. 그러나 세 step, t≈5.8e-5는 역사적 t6–10의 해법 증거가 아니다.

현재 지속비용 약225초/step과 실제 dt bound는 장기 목표에 맞지 않는다. t1 약7.5일/t10 약75일이라는 추정은 현 bound가 유지된다는 낙관적인 조건까지 포함한다. 장기 material shear가 조건수와 dt를 악화시키면 더 느려진다. **이 경로를 그대로24시간 또는75일 달리는 계획은 철회한다.** 진행 중인 Ada cached step은 보존하고 운전자가 완료·감사한다. Ada stage1의 cache 일치와 시간은 채택 가능한 중간 증거이나 아직 accepted fourth step이나 H200 성능 증거는 아니다.

선택하는 다음 한 묶음은 **고정 state/DOF/M을 유지한 양의 reference cubature 압축과 device-resident DD mass iteration**이다. 둘은 독립 최적화 후보의 무한 탐색이 아니라, 현재 두 지배비용을 함께 제거하는 한 번의 production 비용 관문이다. 무조건 차수를 낮추거나 h를 늘리는 경로를 택하지 않는다.

## 2. 먼저 인정할 구조적 한계

10,816,441 scalar DOF와 source당 다수의 P3 residual basis는 71,080개 생성점이라는 표시보다 훨씬 큰 문제다. quadrature만 줄여도 이 DOF의 통신·mass iteration·변형 수명은 남는다. 따라서 이번 최적화가 1000공전을 실현한다고 약속할 수 없다. 공간 정확성만 계속 강화해 계산시간을 늘리는 방향은 중단한다.

현재 acoustic factor7은 polynomial P3에 기반한 추정이며, edge/bubble enrichment를 포함한 전체 consistent-M 연산자의 spectral CFL 정리가 아니다. `eta/gradv`, `eta/Omega`, altitude/sound 제한도 time accuracy와 전역 단사성을 함께 증명하지 않는다. 현 controller는 유지하되 더 큰 h를 허용하는 근거로 쓰지 않는다. 새 reduced rule이 기존보다 더 큰 dt bound를 내도 기존 bound와 비교해 상향 채택하지 않는다.

정확 원운동의 연속 material map도 시간이 갈수록 differential rotation에 의해 길어지고 굽어진다. 유한 basis 표현은 그 map을 영구히 정확히 담지 못한다. positive sampled J만으로 folding 부재를 선언할 수 없다. conservative repartition/recipient thermal context와 W0 동적 Voronoi 환원은 여전히 미완료다. 외압 material 경계와 원래 periodic 경계의 차이, smooth entropy 법칙과 shock entropy 생성의 차이도 유지한다.

## 3. 양의 cubature: 바꿀 것은 적분 rule, 유지할 것은 물리 함수

최신 감사된 accepted checkpoint 한 개를 anchor로 삼는다. 원래 basis, 모든 x/v 계수, source grouping, 고정 binary80 M, generating coefficient/jhat0/entropy 함수와 boundary 법칙을 보존한다. 원래 full64×28 rule을 oracle로 남긴다. 새 rule은 **원래 reference 점들의 부분집합과 비음수 weight**로 구성한다. 순서별 Gauss 시험이나 epsilon/차수 sweep을 하지 않는다.

각 source triangle의 원래 positive measure에 대해 deterministic scaled rank-revealing linear algebra와 positive moment elimination을 한 번 적용한다. constraints에는 다음을 넣는다.

1. reference measure 1, 그리고 고정 reference density를 곱한 basis products `rho*N_a*N_b`. 실제 rank를 사용하고 작은 singular direction을 임의 삭제하지 않는다. 이 조건은 고정 mass metric/K/L의 cubature 일치를 지킨다. authoritative M은 재조립하거나 교체하지 않는다.
2. `C_ab = N_a,X*N_b,Y - N_a,Y*N_b,X`의 독립 span. `J=sum_ab x_a*y_b*C_ab`이므로 이 적분을 보존하면 **모든 현재 계수**에서 QJ와 그 1·2차 configuration derivative가 원래 rule과 같다. 따라서 geometric delta와 adjoint volume operator를 quadrature 변경 때문에 바꾸지 않을 수 있다. constant reference measure도 같아야 한다.
3. anchor의 실제 heat integrand와 thermal gradient, 실제 softened potential과 gravity load의 필요한 독립 constraints. 이 부분은 모든 미래 상태에 대한 정확성 증명이 아니다. 현재 저장 상태를 바꾸지 않고 rule 전환의 에너지/force jump를 검사하기 위한 조건이다. 저장 aggregate가 없는 constraint만 원래 점에서 한 번 streaming 평가한다. 초기 mass/곡률/완료 transfer를 다시 돌리지 않는다.

weight positivity, affine-rank 관계와 각 보존 constraint를 원래 binary80 합과 비교한다. 원래 source의 작은 열 척도를 사용하며 kinetic 크기로 heat 오차를 숨기지 않는다. 원래 M 일치/volume/work/torque 관문은 유지한다. constraint residual을 state나 열에 더해 맞추지 않는다. 불가능하거나 압축률이 부족하면 그 source는 실패/원래 rule 유지로 명시한다. 자동 constraints 삭제나 tolerance 확대는 없다.

이 구성을 모든 source에 무작정 먼저 적용하지 않는다. 기존 full 비용/점 수 자료로 선택한 고정 대표 source 집합(원점 포함, 큰 density gradient, core, 두 contact, outer floor)을 대상으로 **한 번** rank/점 수/비용을 확인한다. 이는 새 trajectory가 아니라 full production rule을 구성할 수 있는지 확인하는 전처리다. 이를 통과한 동일 알고리즘만 전체 source에 적용하며, 미압축 source의 실제 비용까지 합산한다. 원점 density의1e-6 척도나 logistic 꼬리를 삭제하지 않는다.

### 열과 quadrature 전환의 의미

새 점에서 `U(X)=a(X)*jhat(X)^(1-gamma)`는 원래 함수다. point heat를 새 rule에 맞춰 재초기화하지 않는다. 적분 rule 변경으로 생기는 총 U/Phi 차이는 **numerical integral discrepancy**로 별도 기록하며 physical ledger를 조용히 zero-rebase하지 않는다. anchor 열/힘 constraints 및 원래 오차 기준을 못 맞추면 candidate rule을 commit하지 않는다.

이후 force는 새 positive thermal functional의 configuration derivative, heat work는 그 동일 adjoint를 사용한다. 기존 force만 더 적은 점으로 계산하고 원래 heat increment를 따로 쓰지 않는다. nonpolynomial thermal/gravity의 미래 quadrature 오차는 남는다. 이를 다음 실제 stage의 high-rule 대조로 판정한다. 상수압력/운동량/각운동 구조와 실제 pressure accuracy를 혼동하지 않는다.

### 유한 산술에서 constraint를 판정하는 방법

실제 저장 weight는 엄격한 양수여야 하며 제거된 점의 weight만 정확0이다. 작은 음수를0으로 clamp한 뒤 통과시키지 않는다. 각 보존 moment의 binary80 독립 합, 원래 target와의 signed residual 및 절대합에 따른 합산오차를 함께 저장한다. rank/positivity가 수치적으로 unresolved이면 실패로 남긴다. 원래 M을 새 합에 맞춰 바꾸지 않는다.

수학적 moment 보존과 유한 산술 일치를 구별한다. determinant constraint의 잔차를 `epsilon_ab`라 하면 실제 stage에서 `|delta QJ| <= sum_ab |x_a*y_b| epsilon_ab`이고, 같은 식을 계수에 미분해 G와 Hessian 오차를 전파할 수 있다. mass constraint 잔차 행렬 E_M도 실제 v/a로 `v^T E_M v`, `v^T E_M a`, 좌표와 속도의 angular contraction까지 전파한다. 이를 기존 point-vs-M, volume, work/torque 관문과 대조한다. 단순 각 entry의 작은 상대차 하나로 cold U 기준을 대신하지 않는다. 모든 계수에 대해 exact QJ라는 주장은 정확 moment 일치의 대수적 명제이며, 실제 저장 rule에는 이 잔차 장부가 붙는다.

기존 관문을 못 맞추면 그 rule로 실제 다음 stage를 채택하지 않는다. 단계 base는 그대로 보존하고 이번 묶음을 실패로 끝낸다. dt 축소·order 상승·rank 삭제로 통과시키지 않는다. 이 검사는 압축 전처리와 실제 한 전환 step 안에서 수행하며 별도 작은 certificate 루프를 만들지 않는다.

## 4. Mass: host 왕복을 inner iteration에서 제거

현재 `material_full_gpu_mass_host.h::operation`은 매 apply/precondition마다 전체 local vector를 binary80→DD로 변환·복사하고 host halo와 합을 수행한다. 54회 반복에서 GPU local kernel만 가속한 상태다. 다음 구현은 inner x/r/p/z와 local support reduction을 device DD에 유지한다.

- 원래 M/factor/partition의 정확 hi+lo 변환을 재사용한다. pure double mass 또는 tensor-core 저정밀 근사는 없다.
- DD AXPY/dot/reduction을 사용하고 MPI는 owned/ghost hi와 lo를 함께 교환한다. GPU-aware MPI가 없으면 ghost 경계만 pinned staging한다. 전체 vector를 매 iteration host로 이동시키지 않는다. ghost norm 중복, unordered double atomic 합산은 금지한다.
- outer correction과 최종 authoritative `b-Mx`는 기존 binary80 M으로 평가하고 원래1e-16 기준을 유지한다. DD inner와 binary80 outer를 동일 산술이라고 주장하지 않는다. 기존3 correction cap/양의 metric 검사와 zero-RHS 처리를 유지한다.
- x/y 두 RHS는 독립 PCG recurrence를 유지하면서 kernel/halo 패킹을 묶을 수 있다. 수렴 기준을 공유해 한 성분의 실패를 가리지 않는다.

첫 새 실제 RHS 한 개에서 이전 host-vector DD 경로와 새 device-vector 경로를 같은 b로 비교한다. true equation residual 및 기존 work/torque 기준으로 채택하고 이후 반복 비교하지 않는다. 실제 지배시간이 halo/true residual인지 kernel인지 분리한다. H200 한 장과 Ada의 결과를 섞지 않는다.

## 5. 한 번의 실제 전환 관문

현재 cached 실행이 완료하면 그 마지막 감사된 accepted state에서 시작한다. 아직 실행 중이면 변경하지 않는다. 원래 immutable entropy와 M을 가져오며 t0/성공 구간을 재생하지 않는다.

1. candidate rule을 구성하고 positivity·constraint·rank·총 point 수와 예상 비용을 보고한다. 다음단계에 가기 전에 전체 지배 point 처리량이 적어도4배 줄어드는지 확인한다. 이 수치는 물리 tolerance가 아니라 전처리 비용 중단 기준이다. 실패하면 rule/차수 탐색을 반복하지 않는다.
2. **실제 다음 RK4 step의 동일 stage 상태**에서 reduced RHS와 원래 full-rule RHS를 한 번씩 비교한다. 네 개 stage는 해당 step의 실제 stage이며 독립 trajectory sweep이 아니다. 오래된 initial/control을 다시 적분하지 않는다. 모든 셀의 U/Phi, pressure/volume rate, force와 mass-solved acceleration 차이를 보존한다. step 후보가 통과할 때만 commit한다. 비교 때문에 full-rule 원래 궤적을 따로 다시 진행하지 않는다.
3. 기존 algebraic/positivity/M/thermal/torque 관문 외에 새로운 quadrature truncation budget을 명시한다. 전체 계획 구간 T=10에서 새 적분오차가 단독으로 normalized pressure/density/velocity 오차1e-3을 소진하지 않도록, 해당 rate 차이의 dt 적분을 각 항목별 누적한다. 한 step 허용치를 dt와 무관한1e-3으로 주지 않는다. 기존 더 엄격한 기준이 있으면 그것이 우선이다. 이는 안정성 정리나 시간오차 보장이 아니라 **추가 quadrature 오차의 수치 예산**이다. 기존 해 자체의 물리 오차는 별도로 평가한다.
4. anchor 이후 geometric deformation/thermal contrast가 최초 검증 범위를 벗어나면 high-rule sentinel 평가가 필요하다. sentinel은 처음에는 실제 전환 step의 전 셀 대조이며, 이후 checkpoint의 보존/pressure/mesh 지표로 고정된 조건에서 수행한다. 실패 시 accepted state를 보존하고 멈춘다. rule을 조용히 다시 학습하거나 entropy/M을 refit하지 않는다.

DD mass와 cubature 둘 중 하나가 실패하면 통과한 성능 작업은 보존하되 전체 장기 승격은 보류한다. field representation/DOF 수를 줄이는 별도 재표현을 이번 묶음에서 몰래 추가하지 않는다.

## 6. 비용과 실행의 단호한 제한

새 다음 step의 실제 GPU/통신/solve/출력 시간을 사용한다. cache 첫 tile bitmatch는 전체 step 성능 대신 쓰지 않는다. 현재 dt bound가 미래에도 유지된다는 가정의 예측은 낙관적 하한으로 표시한다.

- 기존 대비 지속 step 비용10배 이상 감소를 이번 production 후보의 목표로 둔다. 또는 그것보다 작은 개선이라도 실측 h에서 t1까지24시간 내라는 기존 실행 관문을 만족해야 한다. 두 조건 중 어느 것도 못 맞추면24시간 blind continuation을 하지 않는다.
- 통과하면 동일 accepted trajectory를 t1 또는wall24시간까지 진행하고, .1/.2/1의 실제 Eulerian rho/P/궤도/E/L/회전 및 map condition을 기록한다. 작은 positive step 성공을 원래 열고갈 해결로 선언하지 않는다.
- t1 이후6/8/10으로 갈 때는 실제 변형에 따른 dt와 비용으로 다시 예측한다. 최장1000공전은 지금의 고정 material topology로 보장되지 않는다. current exterior BC가 원래 periodic와 다르다는 표시도 유지한다.

이번 알고리즘 고정 한 묶음이 비용/정확도 관문을 못 넘으면, **이 10.8M-DOF 공간을 유지한 채 최적화만 계속하는 전략은 종료**한다. 다음에는 실제 보존 고차 field 재표현/mesh 수명 문제를 포함한 공간 규모 변경이 필요하다. 그것은 새로운 수학 계약이며 mean compression, kinetic reservoir, heat repair로 대신할 수 없다. 현재 코드가 언젠가 빠를 것이라는 기대만으로 장기 작업을 제출하지 않는다.

현재 dt/step를 GIZMO 약1000공전까지 단순 연장하는 비용은 현실적이지 않다. 이번 압축과 resident mass는 우선 t1/t6의 실제 물리문제에 도달하기 위한 중간 조치다. 최종 장기 시간적분 효율, material mesh의 전단 수명 및 필요시 실제 보존 repartition은 남아 있다. 이번 가속의 통과를 그 문제들의 해결로 확대하지 않는다.

## 7. 종료 보고 범위

운전자는 재승인 질문 없이 위 전처리·native 구현·조건부 실제 전환 step을 수행한다. 보고는 압축 규칙의 실제 rank/점 수/오차, device mass의 실제 비용/true residual, 동일 step의 물리 차이와 누적 오차 예산, accepted continuation 여부를 한 번에 묶는다. 새 primitive별 평가를 반복하지 않는다. 기존 run/state/binary/결과는 보존한다. t6/8/10와 논문 최장 목표는 계속 미달이며, 이번 결정은 준비작업의 완성이 아니라 그 목표에 맞지 않는 계산비용을 실제로 기각·교체하는 선택이다.
