# Astra: registry 완료를 위한 다음 묶음

2026-10-21. 운전자가 구현·실행한다. SHARED_CONTEXT와 bridged-domain 장부를 읽었다. 본 평가는 저장 증거와 제안 수식에 대한 평가이며 새로운 수치검사·폴링은 하지 않았다.

## 판정과 범위

9216 arrangements, internal power-interface cancellation, 1152 source-edge chain 및536 paired edge의 저장 사건 장부 완료를 인정한다. 96 external fragment는 실제 source support 경계로 보존한다. 원래 native exit2와 별도 bridge의 유효 완료를 구별한 provenance도 유지한다.

다음 한 묶음은 **global source-image disjointness + 실제 전 target atlas/cost + 다른23 target normalizer manifest**다. 전체 quadrature·Newton·압축·시간진화는 아직 실행하지 않는다. target68526 출판 및 전체 원본 state/heat/함수는 불변이다.

## 1. 제안한 global pair certificate 채택 조건

384개 reference triangle의 unordered73536쌍을 빠짐없이 분류하고, 아래 충분조건을 고정 순서로 적용한다. 실패는 실제 겹침 판정이 아니라 pending이다.

### Strict image boxes

기존 whole-image enclosure가 한 축에서 엄격히 분리되면 제외한다. 단순 touching은 이 관문에서 제외하지 않는다. binary64 reference geometry, canonical map 및 enclosure epoch/hash를 묶는다.

### 공유 edge: convex union의 bi-Lipschitz

두 양의 면적 triangle이 정확히 한 reference edge를 공유하고 내부는 겹치지 않으며 union이 convex임을 exact predicate로 확인한다. union의 외곽에 collinear vertex가 있어 triangle 모양으로 축약되는 경우는 허용하지만 source 자체의 면적0은 허용하지 않는다.

같은 lab coordinate에서 chi=X+d이고 두 조각의 exact canonical trace가 일치하며 ||Dd||_infinity<=L=max(L_A,L_B)<1이면, union 안의 직선 경로에 대해

    ||chi(X)-chi(Y)||_infinity >= (1-L)||X-Y||_infinity

이다. segment가 shared edge를 통과해도 piecewise derivative bound와 C0가 있으므로 성립한다. 따라서 두 image의 교집합은 정확히 shared edge image다. local floating trace 차이가 작다는 검사만으로 exact C0 전제를 대체하지 않는다. union이 nonconvex이면 이 증명을 적용하지 않는다.

### 공유 vertex: 정확 cone separator

공유 reference vertex V와 동일 mapped vertex chi(V)를 기준으로 각 triangle의 두 ray a_i,b_i를 사용한다. ||n||_1=1일 때 제안한

    n·a_i + L_A ||a_i||_infinity < 0
    -n·b_i + L_B ||b_i||_infinity < 0

는 충분하다. triangle 안의 점은 두 ray의 비음수 결합이며, infinity norm의 삼각부등식으로 모든 비영 결합에서도 strict separator가 유지된다. mapped error는 vertex에서 triangle 내부 직선을 따라 적분하므로 L_A,L_B를 그대로 사용할 수 있다. 결과는 shared vertex 외 image 교집합이 없다는 것이다.

각 sign quadrant에서 n=(sx*p,sy*(1-p)), p∈[0,1]로 쓰면 제약이 p에 선형이므로 exact rational interval/LP로 유효 p와 strict margin을 저장한다. L 상계도 실제 outward binary 값을 exact rational로 해석한다. p=0,1도 허용하되 모든 실제 ray에 strict margin이 있어야 한다. 접하는 reference ray에 이 충분조건이 실패하면 pending으로 둔다.

### Nonshared pair: SAT와 displacement

reference triangle edge normal들의 SAT 축에 대해

    gap_n = min_B n·Y - max_A n·X
    gap_n > sum_k |n_k|*(dA_k+dB_k)

이면 제외할 수 있다. dA,dB는 각각 chi-X의 성분별 absolute displacement 상계다. 상대 중심 이동이나 방향성 상계를 더 정확히 쓰려면 실제 기존 enclosure로부터 유도해야 한다. global L만으로 서로 다른 component의 translation을 상쇄하지 않는다. strict gap이 없으면 pending이다.

각 쌍에는 method, exact geometry inputs, canonical shared entity, L/displacement bounds, inequality margin, 성공/미해결 이유를 저장한다. 전체73536 coverage와 중복 없는 ownership을 확인한다. 모든쌍 통과 시 global image interior-disjoint를 인정하며, 미해결쌍은 cap/정밀도/field scan으로 자동 보충하지 않는다.

## 2. 실제 atlas와 비용 manifest

192 whole 및367 clipped domain 중 이미 target68526에 속한 완성 자료를 재사용한다. 나머지는 완료된 실제 loop/event를 출발점으로 기존 fixed chart/regular-q/potential-P2 규칙의 **기하 atlas와 typed panel 계획**을 만든다. 새 root 전체 재생이나 적분값 계산은 하지 않는다.

각 target/source/component마다 panel 수, fixed Gauss16 query 수, analytic proof 수, 필요한 arithmetic backend, packet bytes, 기존 증거 재사용 여부를 기록한다. 계획하지 못한 domain을0비용으로 놓지 않는다. global 미해결 pair가 있으면 독립적인 atlas 준비는 가능하지만 integral-ready는 false다.

한 target당 기존817152 예산을 자동으로24배 실행 승인으로 바꾸지 않는다. 이번 단계는 실제 unique panel/query 총수·최대 rank memory·증거 출력량을 산출하고, 한 target에서 측정된 비용과 함께 추정 범위를 제시하는 것이다. 아직 해석적 나머지 자체가 새 target에서 검증되지 않았다면 analytic-ready도 false로 구별한다. 비용이 큰 경우 먼저 구조적 재사용/streaming 여지를 기록하며 guard/order/precision/허용오차 조정은 하지 않는다.

유한 patch Omega와 unbounded power cell C_j는 구별한다. 적분 domain은 D_j=C_j∩Omega이며 그 유한 overlap들만 계획한다. source support 외부는 filler로 채우지 않는다. 96 external fragment는 대상 target 및 방향을 유지한다.

## 3. 다른23 target cold normalizer의 출처

target68526의 원래 수치와 tau는 그대로 고정한다. 나머지 target은 같은 물리 계약의 계수(질량/기하 기준 및 cold K,U,potential 기준 등)를 사용하되, **각 실제 D_j의** 양의 normalizer를 출처와 함께 고정해야 한다. legacy full-cell 질량·열 또는 다른 target의 U를 부분 patch domain의 분모로 쓰지 않는다.

권장 구성은 이미 whole로 인증된 source 부분 또는 D_sj 내부에 포함됨을 인증한 유한 reference subregion들의 합이다. mutually-disjoint reference subset B_l에 대해

    m_lower = sum area(B_l)*inf_B_l rho0
    U_lower = sum area(B_l)*inf_B_l u_ref_frozen
    V_lower = sum area(B_l)*inf_B_l J

를 downward로 계산한다. rho0는 원래 연속 밀도, u_ref_frozen은 이번 고정 generating heat의 reference-volume 열, J는 같은 current map이다. 물리U를 새 EOS로 초기화하지 않는다. source images의 global disjointness가 확인돼야 서로 다른 source 하한을 중복 없이 더할 수 있다.

이미 존재하는 유효 interior witnesses를 우선 재사용한다. whole-source는 그 자체나 기존 고정 subdivision을 사용한다. clipped-only target은 기존 atlas의 인증된 내부 subregion을 사용하되 한 번의 고정 witness 규칙/예산을 manifest에 먼저 명시한다. 양의 하한을 얻지 못하면 normalizer_pending으로 남긴다. 실제 열 적분을 보거나 실패한 bound에 맞춰 분모를 늘리는 반복 최적화는 금지다. 하한이 너무 보수적이라는 사실은 정확도 기준을 낮추는 이유가 아니다.

moment별 추가 velocity/radius scale이 필요한 경우 원래 계약과 동일한 정의를 명시하고 certification/단위/zero-scale 처리를 기록한다. 샘플 최대속도나 유리한 대표점을 임의로 선택하지 않는다. empty domain은 empty record이며 양의 normalizer나 EOS를 요구하지 않는다. full power cell이 아니라 partial D_j를 기준으로 한 normalizer임을 표시한다.

이번 단계에서는 normalizer 공식·값·증명·SHA와 모든12tau를 immutable manifest로 만든다. 이후 quadrature 결과를 본 뒤 바꾸지 않는다. 미확정 normalizer는 integral-ready를 차단한다.

## 4. 종료 관문과 후속

이번 묶음의 산출물은 (a)73536 pair proof ledger, (b)24target actual atlas/cost manifest, (c)23new normalizer/tau manifest, (d)68526 기존 증거 재사용과 source/target incidence audit다. 새 물리 quadrature는0이다.

다음 조건을 별도 flags로 둔다: global_disjoint_complete, actual_atlas_complete, cold_normalizers_complete, cost_manifest_complete, analytic_certificates_complete. 앞 네 가지가 통과해도 아직 계산하지 않은 analytic certificate를 true로 두지 않는다. 하나라도 pending이면 해당 witness를 보존하고 이 묶음을 종료한다. 조용한 geometry fallback이나 전체 scan을 이어가지 않는다.

전제들이 완결되면 다음 큰 경계는 실제 전 target analytic proofs와 shared overlap의 보존 수치조립이다. 동일 overlap packet을 outgoing/incoming에 한 번씩 참조하고 전체 source mass,p,L,K,U와 각 target cold 오류를 함께 요구한다. 일반 evolved heat/restart와 실제 시간진화는 이 CONTROL recipe로 대체할 수 없다. 원래 CUDA/MPI/OpenMP t6/8/10 및 GIZMO 최장공전 목표는 그대로 남는다.
