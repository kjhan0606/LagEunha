# Astra: oblique slab 종료 평가와 실제 atlas 연결 계획

2026-10-21. SHARED_CONTEXT 전체, 원본 origin LagEunha/HEAD 8e25dfe, v5 adoption, pending29 projection, internal-event 및 slab 장부, slab 정책·구현과 실제 CSV를 확인했다. SHARED_OBLIQUE_SLABS_20261021.json의 완료된 독립 감사도 읽었다. 검토자는 새 계산·폴링·driver 코드 변경 없이 이 문서만 작성한다.

## 1. 이번 결과의 채택 범위

**14개 실패 chart의 실제 실수 oblique arrangement를 채택한다.** 112개 사건, 98개 slab 중40 nonempty/58 empty와294개 표본 증명은 기존 complete domain/root·전역 transverse derivative 전제와 결합하여 유효하다. 이는 단순히98개 표본에서 성공했다는 판단과 다르다. typed panel, complex implicit-function/density clearance, all12 remainder는 여전히 미완료이며 atlas pending29는 그대로다.

cut=2^-12 아래 apex 영역의 엄격한 제외 증거를 보존해야 한다. cut은 임의로 작은 영역을 삭제하는 물리 cutoff가 아니다. 실제 영역이 그 밖에 있다는 기존 증명이 이번 truncated chart 전체의 전제다. 해당 증명이 없는 다른 chart에 같은 cut을 적용하지 않는다.

## 2. Slab 내부에서 branch가 바뀌지 않는 이유와 남는 조건

X=A+u e+v d, e=(dy,−dx)이므로 det[e,d]=dx²+dy²>0이다. s∈[cut,1], t∈[0,1]인 기존 affine-ray chart의 reference image는 convex trapezoid다. 따라서 corner projection 사이에서 u 단면은 정확히 한 geometric interval이며, 그 양 끝 geometric side는 corner를 지나지 않고 바뀌지 않는다.

각 active target-plane pullback f_a는 전체 해당 영역에서 ∂f_a/∂v의 부호가 엄격하다. 따라서 그 제약은 단면에서 lower 또는 upper half-interval을 만들며, 전체 허용 단면은 하나의 interval 또는 empty다. 각 root는 존재하는 구간에서 연속 graph다. 선택된 경계 또는 feasibility가 바뀌려면 다음 중 하나가 필요하다.

1. graph가 geometric side와 만난다: canonical source-edge 또는 새 internal-side 사건이다.
2. 선택된 두 target 경계가 만난다: 해당 점이 실제 허용 영역의 closure에 있다면 feasible power vertex다.
3. geometric side가 바뀐다: trapezoid corner다.

이 사건들이 모두 포함되고 strict projection order가 확인되면 open slab의 한 내부 표본에서 인증된 순서는 slab 전체로 이어진다. 이번 감사의 사건 포함/제외, 순서, 공유 endpoint 및 실제 section/branch 회계가 그 결합을 지지한다. 실제 feasible power vertex가 해당 chart 밖이라는 기존 증명도 필수다.

**모든 infeasible plane-pair intersection을 추가할 필요는 없다.** 이미 다른 제약에 의해 가려진 branch끼리의 교차는 선택된 경계나 feasible set을 바꾸지 않는다. 반대로 실제 feasible vertex가 used-node 목록에서 누락되면 중간에서 얇은 영역이 생겼다 사라질 수 있으므로 표본 검사만으로는 충분하지 않다. 이번 채택은 retained actual-domain completeness를 전제로 한다. 앞으로 inactive plane 제외도 whole-domain 증명에 연결하며 표본 부호로 대체하지 않는다.

코드의 monotone graph root는 endpoint 반대 부호로 존재를, derivative 분리로 유일성을 확보한 뒤 interval Newton으로 감싼다. 고정256 guard나1e−27 폭 자체가 존재 증명은 아니다. 현재 구현은 그 필요한 부호 검사를 함께 수행한다.

### Event enclosure와 half-open 주의

서로 분리된 event **구간**의 순서가 증명된 것이지 중심값이 정확한 endpoint라는 뜻은 아니다. 후속 atlas는 실제 event를 parameter/implicit endpoint로 유지하거나 endpoint 불확실 띠의 모멘트 상계를 장부에 포함해야 한다. [upper(u_i),lower(u_i+1)]만 적분하고 event 띠를 버리면 안 된다. 반대로 두 이웃 slab에서 동일 띠를 물리 면적으로 이중 계산하면 안 된다. 공유 event ID와 half-open ownership은 위상/회계 규약이고, 유한 폭 enclosure의 오차 관리를 대신하지 않는다.

## 3. 다음 한 구현 묶음: oblique40 strips + origin15 실제 typed atlas

현재 완료 자료를 출발점으로 **남은29개 실제 atlas의 geometry/density-clearance 및 전체 비용 확정**을 수행한다. 새 PDE/상태/차수 실험은 추가하지 않는다. 기존 정상 panel, source133 완료 자료, canonical roots와 frozen field proof를 재사용한다.

### A. Oblique40 strips의 native 좌표와 implicit predictor

각 nonempty slab에서 저장된 branch ID로 v_lo(u),v_hi(u)를 정의하고

    v(u,eta)=(1−eta)v_lo(u)+eta v_hi(u),  0≤eta≤1,
    X(u,eta)=A+u e+v(u,eta)d,
    det(dX/d(u,eta))=det[e,d]*(v_hi−v_lo).

를 사용한다. 별도 정규화 없는 실제 d/e와 determinant를 그대로 보존한다. physical area에는 동일 current-map J를 추가하며 reference-mass moment에는 기존 measure 계약을 따른다. width derivative 항은 determinant에서 상쇄되지만 field/implicit derivative 계산에서는 버리지 않는다.

직선 geometric boundary는 정확 affine 식으로, target-plane graph는 f(u,v)=0으로 평가한다. 중심 predictor는 v'=−f_u/f_v, 곡률은

    v''=−(f_uu+2 f_uv v'+f_vv v'^2)/f_v

이다. 같은 branch의 실수 root enclosure를 anchor로 하여 기존 fixed complex-domain/invariance 방식으로 correction을 인증한다. 실수 transversality만으로 복소 분모의 비영성을 선언하지 않는다. 저장된 실수 root를 중심으로 한 predictor/metric 상계를 먼저 만들고 같은 고정 policy로 panel을 배치한다.

패널 크기는 실제 f_v, source barycentric/edge-extension 분모, radius sqrt germ, density r+epsrho/logistic 및 map J의 complex clearance에서 정한다. 실제 event 끝에서 width가0이 되는 것은 허용하지만 내부의 부호와 branch 순서는 유지되어야 한다. fold나 추가 전환이 발견되면 해당 기존 사건 장부에 원인을 결합하고 pending으로 남긴다. 새로운 angle/guard/차수 후보 반복은 하지 않는다.

### B. Origin15

원래 내부7/외부8 분류, 작은 제외 간격과 원래 profile을 보존한다. 원점 fan X=s R(t)에서 density 반경은 s sqrt(R·R)이며 실제 ray별 허용 s 구간을 target-plane boundary로 정한다. tiny gap을0으로 snap하거나 작은 영역을 filler로 채우지 않는다.

같은 함수에 대해 이미 검증된 collected ray, weighted angular 및 positive-real sqrt germ을 사용한다. 새 chart에는 새 실제 계수/root-gap/branch anchor 증명을 부착한다. source133의 수치를 복사하지 않는다. s=0에서 q/s 식을 직접 적용하지 않고 Duffy blow-up을 유지한다. epsrho=1e−6와 작지만 비영인 logistic 항을 그대로 포함한다.

density pole의 epsrho 규모와 실제 clipping event의 작은 기하 규모를 분리한다. 전자는 고정 clearance에 따른 graded radial 계획으로, 후자는 동일 실제 boundary graph와 event enclosure로 처리한다. 어느 한쪽의 인증도 다른 쪽 성공으로 생략하지 않는다. 직전 계획의 조건과 sqrt branch 계약을 그대로 따른다.

### C. 누락된 정상 arrangement metadata의 최소 회복

추가로 읽은 shared_remaining29_recovery_inventory_v1/contract.json은174 chart를36 excluded,36 whole,14 대체 oblique,88 completed clipped metadata 누락으로 구분한다. 누락 strip은169개이고 옛 root-window 수1970은 재실행 대상 목록이 아니다. 이 inventory를 회복 manifest의 상한 목록으로 채택한다. excluded/whole의 기존 증거는 재사용하고,14 oblique는 이번 새 packet을 사용한다. 88 clipped 중 아래 직접 결합으로 회복되지 않는 것만 geometry-only 재구성 대상으로 남긴다. 옛 status에 없는 reference frame/radial trim을 현재 기본값으로 추측하지 않는다. 원래 frozen 함수·입력과 retained boundary 증거로 동일 frame/trim을 먼저 복원해야 한다.

먼저 새29 atlas가 실제로 참조해야 하는 누락 chart/strip 목록을 만든다. 완료 CSV가 있다는 이유만으로 root/strip 상세까지 저장됐다고 가정하지 않는다. 다음 순서로 최소 자료만 회복한다.

1. 기존 typed panel의 domain/branch/event/range 자료와 retained original domain nodes/boundaries/loops, canonical root manifest에서 직접 결합 가능한 것을 읽어 재사용한다.
2. 상세가 없지만 새 조립에 꼭 필요한 **해당 chart의 arrangement만** 원래 frozen geometry 함수와 동일 저장 roots/events로 한 번 재구성하여 직후 typed snapshot을 남긴다. 이것은 metadata recovery라고 실행 manifest에 별도 기록한다. 원래 field/velocity/heat quadrature나 전체 canonical root 격리는 실행하지 않는다.
3. 회복한 domain/경계·순서·기존 panel 범위/개수와 저장된 완료 산출물을 대조한다. 가능한 부분은 byte 동일성을, representation 차이가 있는 부분은 정확한 event/domain 동치와 이유를 기록한다. 상세가 유일하게 결정되지 않거나 기존 자료와 모순되면 원래 완료 판정을 조용히 덮지 않고 해당 연결을 pending으로 둔다.

이 예외를 정상 전체 chart 재생으로 확대하지 않는다. 실제 필요한 목록을 먼저 고정하고, 추가 누락은 이유와 최소 의존성을 기록한다. 복구 과정에서 새 root 증명이 정말 필요하면 그것은 기존 root 재사용이라고 기술하지 않는다. 해당 작업의 범위를 명시하고 기존 고정 guard에서만 처리한다.

## 4. 산출물과 종료 경계

각 typed panel은 실제 source/target/chart/branch/event ID, 정확한 좌표변환과 measure, endpoint enclosure 처리, 요구 radius evaluator/branch anchor, complex geometry/density proof ID를 포함한다. 기존34955 새 패널 및1047 retained 패널과의 중복·누락도 검사한다. target당 승인된13074432 queries/51072 panels, Gauss16×16을 적용하며 tau/normalizer/eps/정밀도/원본 상태는 불변이다. 비용에는 각 target의 기존 완료 panel과 새 panel을 함께 더한다.

이번 묶음은 남은29 atlas를 완료하거나 각 최초 실패 witness를 보존하고 종료한다. 이미 통과한 축/패널은 재사용하며 cap 증가나 새 fallback 연쇄를 하지 않는다. flags는 actual_atlas、geometry_density_clearance、all12_analytic_remainders、full_quadrature를 분리한다. 첫 두 항이 통과해도 아직 실행하지 않은 all12 remainder를 true로 기록하지 않는다.

이 기하 연결이 닫히면 다음 substantive 경계는 native batch streaming의 전 target analytic certificate와 보존적 outgoing/incoming packet 조립이다. 일반 evolved heat와 native CUDA/MPI/OpenMP 시간 갱신을 실제로 연결해야 원래 t6/8/10 및 GIZMO 최장공전 목표로 진행할 수 있다. 이번 slab 결과는 그 연결을 위한 유효한 실제 arrangement 증거이며, 시간진화나 장시간 물리 안정성의 대체 성공 기준이 아니다.
