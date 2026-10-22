# Astra — material 전단 수명과 조건부 연속-map 검사

2026-10-22. SHARED_CONTEXT 전체, 최신 positive-rule connection 및 shear/compression 결과, exact shear/Jensen 유도, 현재 shared basis와 real-bound 구현을 제한적으로 읽었다. 프로젝트 identity는 해당 결과에 기록된 원본 LagEunha HEAD8e25dfe와 분리된 연구 ROOT를 따른다. 실행·규칙·checkpoint·runner를 변경하지 않았고 계산·폴링하지 않았다.

## 1. 이번 메모의 결정 범위

`ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md`가 계속 채택된 계획이다. 전처리 prefix의 압축률, resident DD mass의 성공, 저장 high k1은 각각의 증거이며 아직 없는 paired step의 정확도/비용 판정을 대신하지 않는다. 현재 실행을 그대로 두고 그 결과를 기다린다. 이 문서는 **paired accuracy/cost가 통과했을 경우** 필요한 guard와 장기 수명의 계약을 정한다. 새 실험, rule 변경, dt 완화, 추가 sweep의 지시가 아니다.

정확 map의 큰 condition number는 실제 numerical folding 관측이 아니다. 반대로 정확해가 J=1이라는 이유로 현재 유한 material 표현도 오래 단사일 것이라고 기대할 수는 없다. 초기 원점의 빠른 강체회전, r≈.054의 큰 전단, 원반 내부 경계 r=2는 별개의 제한이다.

## 2. A — guard의 세 수준과 채택 가능한 충분조건

### 2.1 진단과 증명을 분리한다

- 원래 quadrature-node J/jhat 양성, subset 양성, 원래/압축 rule의 force 일치는 **표본 진단**이다. 연속 source 내부 양성과 전역 단사성의 증명이 아니다.
- 전체 실제 enriched map의 domain enclosure로 얻은 J 하한, 공유 trace, boundary winding 검사는 **명시된 산술·함수 계약 아래의 충분조건 인증**이다. 불확실은 실제 folding과 다르게 기록한다.
- 교차의 interval existence 증명, 또는 오차를 넘는 실제 음의 J witness는 **실제 실패 증거**다. bound가0을 포함한다는 사실만으로 음의 J나 교차를 선언하지 않는다.

현재 `tight_real_bounds`의 `||Dchi-I||<1`은 초기 근처에는 쓸 수 있지만, 시간이 지나 회전한 map을 거부할 수 있다. 이를 전 시간의 필수 물리조건으로 삼지 않는다. 보강 없는 P3 determinant만 검사하는 것도 금지한다.

### 2.2 전체 basis를 포함한 local 계약

각 reference triangle에서 실제 저장 coefficient가 정의하는

    chi = P3 + sum(edge coefficient * s^4 r(t)) + local residual,
    t=lambda_B/s

전체 함수를 bound한다. 원래 radial remainder를 안정한 `q-I3q` 형태로 묶고 P3 derivative는 Bernstein control net으로 평가하는 기존 구조를 재사용할 수 있다. 단순히 두 큰 값을 따로 bound해 얻은 실패를 실제 mesh 실패로 해석하지 않는다.

apex s=0에서는 Cartesian complex analytic extension을 요구하지 않는다. 실제 real domain의 s,t blow-up과 s^4로 생기는 derivative decay를 이용한다. 이 함수는 apex에서 필요한 1차 미분을 갖지만 holomorphic이라는 뜻은 아니다. trace와 limit derivative를 포함해야 한다.

실용적인 첫 sufficient filter는 triangle별 비특이 기준행렬 A_T를 사용하여

    det(A_T)>0,
    sup_T || A_T^-1 Dchi(X) - I || < 1

을 인증하는 것이다. A_T는 실제 derivative의 수치 preconditioner이며 새 좌표/상태/궤도 fitting이 아니다. convex reference triangle에서 이 조건은 local bi-Lipschitz와 방향 보존을 준다. rigid rotation을 I와 비교해 부당하게 거부하는 문제를 피한다. 이 조건을 못 맞추면 곧바로 실제 folding이라고 하지 않는다. 직접 determinant interval 하한 및 제한된 reference proof subdivision을 사용할 수 있지만 물리 mesh/coefficients는 바꾸지 않는다.

동일 reference subdivision과 coefficient-dependent enclosure는 실제 stage의 범위를 덮어야 한다. 이미 얻은 bound는 새 coefficient 변화에 대한 엄격한 margin 전파가 있을 때만 재사용한다. 날짜/epoch만 바꿔 재사용하지 않는다. 정해진 finite work/memory 한도에서 미해결이면 `guard_unresolved`로 종료한다. 무제한 세분화·정밀도/차수 증가를 시작하지 않는다.

J>0의 연속 하한과 별도로 `jhat=J+delta_source>0`도 source delta의 계산오차를 포함해 확인한다. geometric map이 유효해도 corrected thermal measure가 양수라는 결론이 자동으로 나오지는 않는다.

### 2.3 전역 조건: 모든 triangle pair 검사 대신 planar degree

현재 reference는 고정된 disk와 위상동형인 square material 영역이다. 이 topology를 manifest로 확인한 경우 다음 충분조건 경로를 선택할 수 있다.

1. 각 triangle의 closure까지 방향 보존과 local 단사성, nonzero edge tangent를 위 방식으로 확인한다.
2. 공유 edge의 두 제한은 **같은 canonical 함수/DOF**이며 역방향 parameterization만 다르다는 구조적 identity를 확인한다. native cache 차이는 함수/forward-error 계약에 포함한다. 두 독립 근사 curve가 가깝다는 표본 결과로 exact C0를 대신하지 않는다.
3. shared edge의 양쪽 interior가 trace의 반대쪽으로 붙고, interior vertex의 cyclic fan이 한 번만 도는 local degree1을 확인한다. 원래 연결순서만 읽는 것으로는 현재 휘어진 fan 검사가 되지 않는다. 각 incident derivative wedge와 shared tangent를 이용하며 퇴화/unresolved는 남긴다.
4. 외부 canonical boundary가 단순한 oriented Jordan curve이고 winding degree1임을 확인한다. piecewise curved edge의 enclosure/BVH로 비인접 boundary 조각을 배제하고, 남은 후보만 interval curve intersection으로 판정한다. 인접 조각은 공통 endpoint 이외의 만남이 없어야 한다. 좌표를 평균해서 접합하지 않는다.

이때 piecewise map은 local homeomorphism이고 compact disk의 boundary degree가1이므로 내부의 여러 preimage를 허용하지 않는 전역 degree 논증을 적용할 수 있다. 이것은 내부 모든 source pair의 O(N²) 교차검색을 요구하지 않는 경로다. 다중 boundary/구멍/periodic torus에는 위 조건을 그대로 적용하지 않는다.

geometry proof가 실제 analytic canonical field에 관한 것인지, roundoff-enclosed native 평가에 관한 것인지 명시한다. 현 libm/engineering error premise를 썼다면 그 premise에 조건부인 인증이다. native 표본 일치만으로 무조건적인 실수 정리를 주장하지 않는다.

### 2.4 비용과 실행 시점

위는 paired step이 통과한 뒤 long continuation을 위한 구현 계약이다. 지금 별도 guard job을 제출하라는 뜻이 아니다. 기준 coefficient/shape enclosures와 canonical topology를 미리 재사용하고, stage별 변화 margin으로 빠르게 갱신하며 어려운 source/boundary만 자세히 처리하는 구성이 현실적이다. 전체 enriched bound의 full-disk 비용은 아직 측정되지 않았다. 값싼 GPU node 검사만으로 이를 대체할 수 없고, 비용이 현재 RHS보다 커도 괜찮다고 가정할 수도 없다.

accepted endpoint뿐 아니라 실제 RK stage를 검사해야 한다. endpoint 둘이 유효해도 사이 stage가 뒤집힐 수 있다. 시간 전체까지 유효하다는 주장은 별도의 시간구간 enclosure 없이는 하지 않는다. discrete stage guard와 연속시간 보증을 구별한다. 실패 report에는 source/edge/vertex ID, coefficient epoch, interval/margin, 실제-negative 또는 unresolved 구분과 비용을 남긴다. guard 미완료 표시를 paired 양성 결과만으로 지우지 않는다.

## 3. Spectral CFL의 남은 범위

현재 실제 dt와 원점 full-gradient 추정의1e-9 수준 일치는 제한 위치의 강한 단서다. 정확한 minimum source/branch를 저장하지 않았으므로 원인 확정은 아니다. 이후 이미 수행하는 reduction에서 최소 branch/source/point, grad·Omega·acoustic 세 후보와 h를 기록하면 새 RHS 없이 이 점을 명확히 할 수 있다.

P3의 factor7 acoustic 추정은 enrichment/consistent M/비선형 pressure tangent 전체의 spectral bound가 아니다. 회전 strain이0이라는 이유로 full-gradient 제한을 지우면 RK4가 여전히 적분해야 하는 빠른 회전/중력·위상 오차를 무시한다. 현재 dt 법칙을 바꾸지 않는다. 향후 증명이 필요하면 실제 전체 tangent와 M에 대한 보수적 spectrum 또는 energy stability 조건이 대상이며, 몇 개 Ritz 값·단일 frozen eigenvalue가 전체 nonlinear 안정성 증명은 아니다. 이번 메모는 그 계산을 새로 지시하지 않는다.

## 4. B — 고정 material 공간의 t6/8/10 현실성

정확 원운동에서 J=1이면서 t8 near-softening cond≈1.27e5인 것은 모순이 아니다. stretching과 compression은 다르다. 중심 r=0의 cond1도 그 인접 floor 영역 전체가 강체라는 뜻은 아니다. fixed basis의 해상도 요구는 cond 하나뿐 아니라 `t*Omega'(r)*Delta r` 및 더 높은 map derivative에 의해 커진다. 초기 radial velocity residual은 시간에 따라 감긴 정확 map 자체를 basis에 포함하지 않는다.

따라서 현재 상태/고정10.8M DOF로 t6/8/10을 정확하게 통과할 현실성은 **입증되지 않았고, 특히 central floor는 심각한 위험 구간**이다. 큰 cond만으로 불가능성이나 예상 folding 시각을 단정하지 않지만, dt나 quadrature 차수 조정이 공간 표현의 수명을 회복하지는 않는다. paired 비용/정확도 실패 시 이 공간을 무한 최적화하지 않는 기존 결정을 유지한다. 통과하더라도 압력/밀도·궤도 오차와 continuous-map margin으로 실제 수명을 판단해야 한다.

### 필요해질 수 있는 진짜 재표현 관문

이는 현재 새 방법의 채택이 아니라, representation 수명이 실제로 제한할 때 필요한 최소 조건이다.

- 단순 target labels 변경은 reference map의 감김, conditioning 또는 mass-space 규모를 줄이지 않는다. 실제 reset-reference/remesh는 물리장의 보존 재표현이어야 한다.
- 새 reference Y와 기존 X의 relabeling psi에 대해 같은 물리장을 `chi_new(Y)=chi_old(psi(Y))`, `v_new(Y)=v_old(psi(Y))`로 나타내야 한다. reference mass form은 Jacobian을 포함해 push-forward하고, 실제 thermal/entropy 분포도 같은 물질을 추적해야 한다. 단순 current-J entropy 재초기화가 아니다.
- p/L/총 E만 맞추고 velocity variance를 없애면 cold U를 오염시킬 수 있다. 기존 full field와 resolved rotation/strain/K를 정해진 오차 안에서 함께 보존해야 하며 누락 K를 별도 reservoir나 열로 넘기지 않는다. 새 공간에 담을 수 없으면 projection 자체가 실패다.
- 현재 `nu=w(J+delta_source)`는 quadrature와 source grouping에 의존한다. recipient별 새 delta를 적용하면서 A를 그대로 복사하는 것은 물리 U/P를 보존하지 않는다. source-context fragment를 유지하는 것과 새로운 well-conditioned thermodynamic space로 옮기는 것은 별개다. 새로운 nu/열 함수/force의 동일 변분·work 계약, positive U, same-state E/L/압력 및 restart equivalence가 필요하다.
- 새 고정 mass matrix가 필요한 경우 그것은 이전 M을 유지한 단순 가속이 아니다. 별도 epoch의 보존적 state 전환이며 원래 before/after discrepancy 장부가 있어야 한다. W0의 개별 고정질량/Voronoi 관계와 weight-linked mass relaxation도 그 새 계약에서 따로 닫아야 한다.

이 관문을 피하려고 floor 영역을 제외하거나 원궤도 강제, mesh steering, spin damping을 적용하지 않는다. 보존적 fragment 전달을 이미 완료했다는 사실도 recipient thermodynamics와 장기 비용의 자동 해결은 아니다. 현재 cost 묶음의 성공/실패가 아직 없으므로 이 재표현의 구현이나 새 실험을 지금 승인·지시하지 않는다.

## 5. C — Jensen 진단에 필요한 최소 실제 출력

고정 양의 rule마다 생성 measure를 명시한다:

    V0q=sum w*j0, alpha_q=w*j0/V0q,
    lambda_q=j/j0, m=sum alpha*lambda=Vq/V0q,
    f(lambda)=lambda^(1-gamma),
    B=f(m)-1,
    D=sum alpha*f(lambda)-f(m).

fresh constant-P 생성 law에서 U0q=P0*V0q/(gamma-1), U/U0q-1=B+D다. 정확산술에서 D>=0이다. physical compression energy의 분해이며 비가역 entropy 생산을 측정하는 식은 아니다. 정확 원형 정상해에는 lambda=1이므로, 이 문제에서 의미 있게 큰 B/D가 나타나면 궤도·map·수송 오차와 함께 추적해야 한다. D가 양수라는 사실만으로 numerical 원인 하나를 확정하지 않는다.

추가 point dump 없이, 이미 수행하는 native reduction에 다음을 셀별로 남기면 충분하다.

- 고정 생성자료: rule ID/SHA, gamma/P0, V0q/U0q, source/material-context ID, generating coeff/delta0 SHA. 원래 canonical V0와 quadrature V0q 차이를 별도 기록한다.
- 현재 U, Vq, canonical Vb, m, B, signed D와 각 합산/직렬화 오차 척도. 낮은 rule과 high rule은 각각 자기 생성 measure로 계산한다. 물리 U0를 덮거나 numerical transition ledger를 재기준화하지 않는다.
- `Var_alpha(lambda)`를 weighted centered reduction으로 계산한다. `mean(lambda²)-mean(lambda)²`만 쓰면 현재처럼1e-16 잔여에서 cancellation이 지배한다. 최소/최대 lambda 및 해당 source/point ID, 필요하면 centered third absolute moment를 함께 남긴다.
- J/jhat 및 generating J/jhat0의 extrema, delta/delta0와 volume identity 오차를 연결한다. 이 정보가 없으면 geometry 압축과 corrected-volume 항의 영향을 분리하기 어렵다. 가능하면 이미 같은 pass에서 계산하는 physical J/Jgen의 centered variance를 별도 항목으로 둔다.

안정적인 D 계산은 평균 m을 구한 뒤

    D=sum alpha*[ f(lambda)-f(m)-f'(m)*(lambda-m) ]

를 사용하는 방법이 있다. 기존 bulk/total 차감 결과도 signed 값으로 보존한다. 매우 작은 D를0으로 clamp하지 않는다. 두-pass reduction은 현재 dependent pass에 결합할 수 있으며 진단을 위해 새 force/mass solve가 필요하지 않다.

양의 lambda 범위가 [l,u]이면

    gamma*(gamma-1)/2 * u^(-gamma-1)*Var <= D
       <= gamma*(gamma-1)/2 * l^(-gamma-1)*Var

가 Taylor remainder로 성립한다. 이는 실제 sampled positive measure의 범위에 대한 관계다. 전체 연속장 extrema 인증으로 확대하지 않는다. m≈1일 때의5/9 계수는 gamma=5/3에만 해당하며 native에 기록된 실제 gamma를 써야 한다.

현재 saved U/V/U0만으로는 B와 signed residual을 계산할 수 있지만 lambda variance, extrema, tail 분포나 특정 subcell의 원인을 복원할 수 없다. 현재 잔여 약1e-16은 산술 잡음과 구분되지 않는다. 이를 현재 가짜 가열 검출이나 고갈 원인 해결로 보고하지 않는다. 같은 평균 U/V를 갖는 여러 lambda 분포가 가능하므로 Jensen gap을 variance의 정확한 역산값으로 쓰지도 않는다.

## 6. 조건부 후속 판단

1. 현재 paired step이 끝날 때까지 기존 계획과 실행을 유지한다. 결과 전 승격/기각하지 않는다.
2. paired 정확도·비용이 통과하면 guard_pending을 지우기 전에 위 전체 enriched-map/경계·fan 계약의 실제 범위와 비용을 닫아야 한다. 단순 추가 node 검사로 대체하지 않는다.
3. 실제 guard·physics가 제한하는 시점에는 마지막 accepted state를 보존한다. `actual_fold`, `guard_unresolved`, `spatial_accuracy_failure`, `cost_failure`를 분리한다. 같은 fixed 공간의 dt/order/precision sweep을 계속하지 않는다.
4. 장기 재표현이 필요하다는 판단과 특정 새 재표현 연산자가 이미 승인됐다는 주장은 다르다. 지금은 현재 cost 묶음의 실증 결과를 기다린다. t>6와 최장 GIZMO 목표, 경계·공전수 정규화 및 shock/W0 미완료 범위를 그대로 유지한다.

이 메모는 한 번의 전략·수학 보완이며 primitive별 추가 호출이나 현재 실행의 변경을 요구하지 않는다.
