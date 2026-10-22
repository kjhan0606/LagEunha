# Astra — compact topology, 유한 weight law와 실제 chart 결정

2026-10-22. SHARED_CONTEXT 전체, RESULT_COMPACT_POWER_DOMAIN_TOPOLOGY 및 실제 registry/power-curve/power-cell/domain 헤더를 정적으로 읽었다. 원본 HEAD8e25dfe와 accepted step5는 운전자의 확인된 identity를 계승한다. 실행·폴링·코드 수정 없이 작성한다. 작성 말미 운전자가 v13 다섯 합성 회귀와 component links/SHA의 완료 PASS를 전달했다. 실제 source 실행은 여전히 없다.

## 1. 종료 평가

현재 구현은 합성 geometry에서 검사한 **boundary graph 구성요소**로 채택 가능하다. 실제 source import, 전체 target coverage, positive chart, projection 또는 accepted step을 통과했다고 해석하지 않는다. 읽은 코드에서 명백한 방향 반전이나 unresolved event를 성공 처리하는 분기는 발견하지 않았다. 이는 모든 실제 입력의 무결성 증명은 아니다.

다음 구조는 타당하다:

- sorted site pair의 동일 plane, triple로 정의한 같은 power corner, canonical exterior curve/root ordinal을 연결한다.
- spatial exclusion은 현재 convex power polygon 전체와 node의 모든 generator에 대한 엄격한 하한이다. 이후 clipping이 polygon을 줄이므로 그 exclusion을 재사용할 수 있다. native high 값에 의한 탐색 순서는 허용하지만 실제 exclusion/plane은 canonical enclosure를 써야 한다.
- transverse root가 모두 격리되고 power-corner/endpoint/tangent/multiple 사건이 unresolved이면 거부한다. 그 조건하에서 연속 boundary 조각의 중간점 membership은 조각 전체의 membership을 결정할 수 있다. 구간 내부에 root가 없다는 completeness가 필수다.
- fixed-ray winding의 endpoint/tangent 불확실성을 실패시키는 정책은 보수적이다. ray 실패가 실제 domain 부재라는 의미는 아니다.

`close()`의 one-in/one-out cycle 검사만으로 geometric embedding을 증명하지 않는다. 실제 채택에는 source의 동일 epoch embedding, 모든 active plane/event의 completeness, 서로 다른 arc의 실제 교차가 사건 목록에 빠지지 않았음, oriented external cycle가 함께 필요하다. 같은 graph ID의 enclosure가 조금 다르면 합치는 tolerance 대신 동일 원함수/사건 정의와 enclosure 교집합으로 연결해야 한다. 현재 strict equality 때문에 거부되는 것은 false rejection일 수 있으나, 이를 좌표 평균으로 고쳐서는 안 된다.

`accepted_generator_anchors_bound` shortcut은 지금처럼 실제 global source certificate에 결합한 후만 사용한다. 그 certificate의 모든 rank coefficient identity와 registry hash를 transaction manifest에서 함께 묶는다. synthetic Registry 자체를 source embedding의 대용으로 쓰지 않는다.

실제 transaction의 `intersect_actual_domain` 입구는 shortcut 사용 여부와 무관하게 `bind_actual_source(g)` 완료를 강제해야 한다. low-level synthetic intersect가 unbound 입력을 허용하는 것 자체와 실제 physics 진입 정책을 구별한다. 운전자가 예정한 wrapper 강제는 필요한 보완이다.

v13의 high+low+tail/canonical_bound는 같은 canonical site를 더 정확하게 둘러싸는 표현으로 타당하다. tag가 없는 site의 uncertainty를 임의로 줄이지 않는다. canonical root 정의는 이 enclosure의 중심점이 아니라 원함수다.

## 2. 이번에 봉인할 finite mass law

제안한 선택을 채택한다:

    mref_i = 이전 accepted fixed-M source cell mass,
    Vsource_i = 같은 accepted source의 양의 canonical 물리 면적,
    mu_i = mref_i / Vsource_i,
    F_i(W) = m_i(W) + mu_i W_i - mref_i = 0.

mref/Vsource/실제 gamma와 epoch, 단위 및 산술 표현을 첫 manifest에 봉인한다. W의 단위는 길이²이고 mu의 단위는 이 2D 문제에서 질량/면적이다. W Newton 중에 mu나 mref를 갱신하지 않는다. W0 overlap mass를 mref로 재정의하지 않는다. Vsource는 새 target 면적이나 곡선의 chord 면적이 아니다.

q와 Ω* 및 old physical rho를 고정한 각 Newton 문제에서, i-j 내부 face의 계수를

    c_ij = ∫face rho_old / (2 |q_j-q_i|) dS

로 두면 dm_i/dW_i=Σ_j c_ij, dm_i/dW_j=-c_ij이고

    J_F = diag(mu) + L(c)

이다. pair를 한 번 적분해 두 행에 반대 offdiagonal로 넣는다. rho는 새 fitted target density가 아니라 실제 source field다. source 경계를 지나면 face 적분을 source 조각으로 나눈다. face가 density 불연속면과 정확히 겹쳐 고전적 미분이 정해지지 않는 경우 one-sided/semismooth 처리를 구현했다고 가정하지 말고 이번 유한 정책에서는 거부한다.

양의 mu 때문에 regular 상태의 행렬은 SPD이며 common-weight gauge도 이 law가 정한다. Newton 후 W에 임의 상수를 빼면 law가 달라진다. 완전 coverage와 정확한 총질량에서는 Σmu_i W_i=0이 따라온다. 실제 수치 장부에서는 Σmref−∫Ω*rho의 기존 허용 오차도 함께 기록하고 gauge로 그 차이를 숨기지 않는다.

고정 rho/domain의 convex power dual에 mu_i W_i²/2를 더하면 coercive strictly-convex 문제를 얻는다는 것은 유용하다. 그러나 실제 root에서 nonempty target/원래 cold budget/finite event 정책이 통과한다는 보장은 아니다. sparse Newton/CG와 고정 유한 line-search budget을 사용하고 geometry 또는 positivity 실패를 unresolved로 보존한다. 이 존재 논의를 과거 hard-material DAE의 장기 geometry 실패를 뒤집는 근거로 쓰지 않는다.

W=0에서 개별 mass 고정은 F(0)=0일 때의 조건부 명제다. 일반 Voronoi overlap에서 자동 성립한다고 하지 않는다.

## 3. 선택하는 chart 알고리즘 하나

제안한 **curved arc 1회 midpoint 분할 → chord polygon의 인증 ear triangulation → collar ear triangle의 centroid fan3 → patch당 최대 한 curved edge**를 실제 전체 transaction의 고정 chart 후보로 선택한다. midpoint는 실제 parameter 구간의 중점 정의이며 새로운 exact rounded vertex가 아니다.

순서는 다음과 같다.

1. 각 actual graph component의 oriented arc를 보존한다. curved arc마다 원래 암시적 두 parameter endpoint로 tmid=(ta+tb)/2를 정의한다. 새 node ID는 parent arc+midpoint 정의이며 power root와 혼동하지 않는다. Gamma(tmid)의 동일 함수를 enclosure로 평가한다. 역방향 arc도 같은 midpoint를 공유한다.
2. 이 node들을 연결한 chord polygon이 단순하고 방향이 확정되며 zero edge가 없음을 인증한다. 두-arc cap에서 midpoint가 필요하지만 그것만으로 positive polygon이 보장되지는 않는다. collinear/교차/unresolved이면 실패다. source의 곡선 영역을 chord 영역으로 바꾸는 단계가 아니라 **reference connectivity**를 만드는 단계다.
3. 고정 global-ID 순서로 ear를 선택한다. convex orientation, 열린 ear 내부의 다른 vertex 배제, diagonal의 polygon 내부 조건을 모든 endpoint enclosure에 대해 증명한다. 경계 equality를 epsilon으로 깨지 않는다. 인증 가능한 ear가 없으면 이번 후보는 종료한다.
4. curved boundary edge를 가진 ear triangle을 그 실제 chord triangle의 centroid로 세 분할한다. centroid도 원래 암시적 vertex의 대수적 평균 정의다. 각 새 patch는 최대 한 curved edge를 가진다. bulk는 기존 compact fan을 유지한다. 인접 patch의 straight interface는 같은 node/function ID를 사용한다.
5. 곡선 patch는 이전 부록의 동일 extension을 사용한다:

       G(s,t)=A+s[(1-t)B+tC]+s^4[Gamma(t)-((1-t)B+tC)].

   source field나 curve를 polynomial refit하지 않는다. 각 chart의 실제 J>0를 전체 reference 영역에서 인증한다. 단순 quadrature-node J>0는 부족하다.

## 4. 국소 J와 전체 coverage 사이의 필수 연결

ear triangulation이 chord polygon을 한 번 덮는 것과 곡선 patch가 실제 T_i를 한 번 덮는 것은 별도다. 여기서는 다음 하나의 연결을 사용한다: chart들이 공통 내부 edge에서 정확히 C0이고, 전체 외부 trace가 graph의 실제 단순 oriented boundary를 한 번 따르며, 각 chart가 orientation-preserving이고 유한 source/target component의 topology와 맞는다면 degree=1과 양의 local Jacobian으로 중복/누락을 배제한다. vertex/edge의 compatible local topology도 포함해야 한다. 이 전제를 확인할 수 없으면 coverage 미증명이다.

graph component가 여러 개면 component별로 처리한다. component 방향/포함관계가 hole을 나타내는데 현재 단순 ear 경로가 이를 지원하지 않으면 거부한다. `component.size()`만 보고 모든 cycle을 양의 filled region으로 적분하지 않는다. 서로 다른 target의 직선 경계는 동일 canonical plane/vertex/event와 반대 orientation을 가져야 한다. 외부 curved arc는 전체 target registry에서 정확히 한 번 분할 소유한다.

이렇게 degree/trace 연결을 완성하면 전 target의 모든 patch 쌍을 거대한 일반 교차 문제로 다시 만들 필요는 없다. 다만 전제가 unresolved인 부분을 global conservation 하나로 면제해서도 안 된다.

## 5. geometry modes와 열·회전 계약

현재 component별 고정 geometry carrier h1,h2를 실제 합성된 curved charts와 straight P3 reference의 차이로 만든다. 위치·속도·mass·force는 이전 부록의 동일 geometry-completed P3 공간을 쓴다. curved geometry만 유지하고 velocity에서 h를 빠뜨리면 rotation-closed 상태공간이 아니므로 거부한다.

전역 최대 두 scalar mode/component의 rank를 실제 양의 mass metric으로 확인한다. exact dependency와 수치 unresolved를 구별한다. 작은 고유값을 tolerance로 잘라서 rank를 만들지 않는다. geometry mode를 초기화 이후 고정시키거나 외부 위치를 고정하는 경로는 없다. 계수는 실제 방정식으로 진화한다.

곡선 collar는 degree9 polynomial mass exactness의 적용 대상이 아니다. 실제 J/geometry mode를 포함한 positive integration과 기존 오차 장부로 처리한다. interior bulk의 기존 degree9 코드는 재사용한다. physical/reference density tag, 생성 jhat0/entropy, 같은 boundary G/H·traction/work는 유지한다. 순간 재격자에 물리 경계 일을 추가하지 않는다. velocity 투영으로 변한 이후 boundary-work rate는 기존 representation budget으로 평가한다.

## 6. 다음 실행 범위와 단일 종료 보고

운전자는 v13 합성 묶음의 최종 결과를 보존한 뒤 별도 primitive 평가 호출 없이 전체 bounded transaction을 연결할 수 있다. actual source import/context → embedding → current sites/immutable mu → 실제 target domain/weight solve → 위 fixed chart/coverage/rank → overlap/positive rho,e/shared p,L,K projection → 조건부 실제 다음 RK4 step 순서다. 기존 통과 자료와 완성된 코드만 재사용하며 성공한 step5 또는 옛 전체 적분을 반복하지 않는다.

전체24시간 예산, 기존 cold/conserved tau, 누적 T=10의1e-3 representation·적분 예산 및 actual amortized t1≤24h 비용 관문을 유지한다. root/ear/J/rank/weight/projection 중 unresolved이면 accepted step5와 실패 witness를 보존하고 종료한다. midpoint 추가를 무한 세분화로 확장하거나 plane/차수/정밀도/허용치를 바꿔 통과시키지 않는다.

이번 종료 보고는 실제 full transaction의 채택 또는 구체적 구조적 실패 한 묶음이어야 한다. 합성 topology PASS나 chart 후보 완성으로 t6/8/10 및 최장 GIZMO를 대체하지 않는다. Re∞, no steering, physical U/E/L와 resolved rotation, 조건부 W0 환원 및 원래 native CUDA/MPI/OpenMP 경로를 유지한다.
