# Astra: source 분리·normalizer·atlas 종료 평가

2026-10-21. 계획·평가 담당 Astra. SHARED_CONTEXT 전체, 직전 ASTRA_SHARED_REGISTRY_COMPLETION_PLAN, source-pair/normalizer/atlas 독립 장부와 test_shared_normalizers_v1.cpp를 읽었다. 원본 worktree HEAD 8e25dfe를 확인했다. 새 계산·테스트·폴링은 하지 않았다. 이 문서만 작성한다.

## 1. 이번 묶음 판정

**Source image 내부 분리와 23개 새 normalizer는 완료로 채택한다. 실제 atlas는 35개 미완료로 이번 묶음을 종료한다.** 실패를 실제 교차, 특이점 또는 물리 오류로 해석하지 않는다.

- 73536 pair는 bbox 71220, convex shared-C0 536, vertex cone LP 1476, SAT 304로 모두 회계가 닫혔다. 정확 부등식 검사는 기존 canonical C0와 whole-source derivative/displacement 증명을 전제로 한다. 공유 경계 접촉은 보존하며, patch 외부를 채우는 증거는 아니다.
- 188개 whole subset을 이용한 다른 23개 target의 양의 하한과 원래 계수의 12tau를 고정한다. 68526의 원래 tau는 불변이다. 기준 영역은 power cell과 실제 source patch의 교집합이다.
- 새 37107 typed panel과 재사용 1047 panel/268032 query의 출처를 구분한 비용 장부를 채택한다. incomplete target의 비용은 UNKNOWN이며, 완료 target만 예산 이내라는 사실을 전체 예산 통과로 확대하지 않는다.
- atlas JSON의 `other23_normalizers_complete=false`는 해당 실행 시점의 기록이다. 최종 결합 manifest에서는 별도 완료된 normalizer 장부를 참조하여 상태를 합성한다. 원본 JSON을 소급 수정하지 않는다.
- 전 target analytic certificate, 수치 source closure, 전체 전달 및 일반 evolved heat는 미완료다. 이번 결과는 실제 시간진화나 GPU 전달의 통과가 아니다.

### Normalizer 수학과 한계

같은 좌표계에서 F=I+Dd, ||Dd||∞≤L<1이면 F의 고유값은 중심1, 반경L 원판에 있고 I+tDd가 계속 가역이므로 J>0이다. 2차원에서 (1−L)^2≤J≤(1+L)^2가 성립한다. 따라서 gamma>1일 때 원래 frozen reference heat

    u_ref = P0/(gamma−1) * J_gen^(1−gamma)

의 하한에 Jupper를 쓰는 부호는 맞다. 현재 구현의 generating/current 계수 동일성은 이 CONTROL에서 current J 증명을 generating J에 쓸 수 있게 한다. 이후 다른 generating geometry 또는 일반 진화 상태에는 이 추론을 재사용할 수 없다. gamma>1, 원래 양의 floor/P0, geometry·basis·heat SHA도 manifest 전제로 명시한다.

Reference area 곱셈은 downward로 수행하며, whole 포함과 source image 분리 덕분에 다른 source의 양의 하한 합산이 정당하다. 고정 reference anchor 속도 및 reference vertex 반경은 기존 normalization 정의이지 영역 전체 속도의 infimum이나 현재 물리 반경 하한이라는 주장이 아니다. 이 정의와 coarse nonexcluded 집합은 이후 비용·오차 결과와 무관하게 동결한다.

## 2. 다음 한 묶음: 35개 실패의 기하 원인과 대체 atlas 계획

범위는 저장 실패 목록의 15+14+4+2개 slot이다. 정상 domain, root, field, quadrature를 다시 계산하지 않는다. 동일 source/chart가 여러 target에 쓰이면 공통 증거는 한 번 만들고 target별 clipping 조건만 구별한다.

현재 실패 status CSV에 없는 Strip/event/prior 내부 자료를 이미 저장된 것처럼 취급하지 않는다. **실패한 원래 chart만** 원래 함수와 입력으로 한 번 재구성하여 수정 전 typed snapshot을 먼저 남긴다. 원래 이벤트 ID/구간, 활성 면, source/target/epoch, chart 좌표, partial strip, 방향 미분 구간, inverse 입력과 실패 잔차, 예상 panel count를 포함한다. 원래 실패 이유가 재현되지 않으면 새 정책 적용을 멈추고 provenance 차이를 보고한다. 정상 chart를 이 작업에 섞지 않는다.

이번에는 기하 증거와 fixed-rule atlas 계획까지만 구현한다. 새 유체 적분, 전 target analytic proof 및 시간진화는 실행하지 않는다.

### A. 15개 density-origin separation

먼저 원점이 reference source triangle에 속하는지 exact predicate로 판정한다. 속한다면 동일 canonical chi(0)와 실제 target half-plane, 이미 닫힌 domain/event를 이용하여 원점이 실제 pullback domain의 내부·경계·외부 중 어디인지 인증한다. 원점이 source 밖이면 reference triangle 자체의 양의 거리도 사용할 수 있다.

외부인데 axis-aligned bbox만 원점을 포함하는 경우, 실제 경계 fragment의 인증된 차트별 거리 또는 분리 half-plane으로 양의 r 하한을 얻는다. 전체 bbox를 반복 세분화하여 원점을 억지로 빼는 방식은 선택하지 않는다. 원점이 실제 domain 경계/내부이면 rmin=0을 그대로 기록한다. 이 경우 기존 separation 관문의 가정이 적용되지 않는 것이며 물리 밀도 발산이 아니다.

원점을 실제로 포함하는 부분에는 원점에서 출발하는 Duffy chart X=s R(t)를 사용한다. 실수 s≥0에서 r=s sqrt(R·R)로 원래 밀도 전체를 평가한다. epsrho=1e−6 및 작은 inner logistic 항을 모두 유지한다. finite rho(0)는 Cartesian 해석성을 뜻하지 않으며, chart의 sqrt 가지, r+epsrho의 pole 및 logistic 복소 pole 거리는 따로 인증해야 한다. 원점으로부터 멀리 있는 반대 꼭짓점의 기존 Duffy singularity 처리와 혼동하지 않는다.

### B. 14개 s/t 축 부적합

실제 활성 target 면만 사용해 공통 oblique 방향 d의 충분조건을 구성한다. source map의 기존 L 상계로

    |n·d| > ||n||1 * L * ||d||∞

이면 n·F d의 부호가 고정된다. 더 좁은 기존 gradient cone을 사용할 때도 같은 방향의 whole-region 부호 증명이 필요하다. 이 조건은 충분조건이며 실패가 곧 접선이라는 뜻은 아니다.

방향은 실제 normal cone의 금지각 구간 사이에서 deterministic 최대 여유 방향으로 한 번 선택한다. 촘촘한 angle scan은 하지 않는다. 공통 방향이 없으면 이미 저장된 vertex/tangency/event로 영역을 나누어 각 subchart의 활성 면에 적용한다. 한 경계가 implicit graph라는 사실만으로 전체 polygon이 단일 strip이라는 결론은 내리지 말고, 교차 순서·빈 구간·다중 연결 strip을 사건 장부에 결합한다.

원점 또는 Duffy apex에 닿는 subchart는 blow-up 좌표를 유지한다. 새로운 oblique Cartesian complex inverse가 apex를 가로지르도록 만들지 않는다. 차트 전환은 같은 canonical 함수와 적분 domain의 좌표변환이며 상태 재표현/refit이 아니다.

### C. source133의 4개 inverse 실패

동일한 실제 target node가 source edge/vertex 사건인지, 엄격한 interior인지 먼저 저장 event와 barycentric enclosure로 구별한다. 경계 사건이면 원래 canonical edge parameter와 shared event ID를 사용한다. 2D Newton을 경계 밖으로 확장하거나 좌표를 epsilon만큼 밀지 않는다.

엄격한 interior라면 그 점의 원래 잔차, contraction norm, prior box 및 각 barycentric margin을 저장하고, source 내부의 인증된 좌표 box에서 같은 inverse 정리를 적용한다. 사건 결합이 먼저 해결되어야 한다. Newton 실패는 outside 또는 empty의 증거가 아니다. 네 target이 공유하는 동일 source-node 증명은 재사용한다.

### D. source135/target67943 및 source187/target68000의 비용 초과

운전자가 이번 검토 중 전달한 exact reference geometry 결과에서 두 source triangle은 원점을 포함하지 않고 거리제곱이 9/524288이다. 실제 pullback overlap은 source triangle의 부분집합이므로 이 거리의 downward sqrt는 overlap 반경의 유효 하한이다. **이 두 건의 첫 선택은 기존 curved-boundary bbox 하한과 이 exact source-distance 하한의 max를 취하여 원래 atlas 규칙을 그대로 계획하는 것**이다. 서로 다른 좌표의 거리를 섞지 않고 둘 다 reference |X| 하한임을 확인한다. 이는 물리 분모/tau를 키우는 조정이 아니라 동일 domain의 더 강한 기하 증명이다. 이 선택으로 기존 예산에 들어오면 아래 graded 변경은 수행하지 않는다.

두 실패의 실제 최소 반경·pole 거리·방향 stretch와 기존 균일 panel 배치 비용을 먼저 분리한다. 그 뒤 **같은 차수와 기존 complex clearance를 만족하는 기하 기반 graded radial atlas** 하나를 산출한다.

원점 인접 density scale만 지배할 때 panel 폭은 local r+epsrho에 비례하도록 하고, 비례상수는 기존 complex-domain clearance 부등식에서 정한다. 첫 panel은 원점을 포함하고 epsrho 규모의 유한 폭을 가진다. 이후 경계는 해당 고정 recurrence와 실제 event에서 결정한다. 다른 pole, source 경계, implicit-graph contraction 또는 tangency가 더 제한적이면 그 거리를 사용한다. 임의의 새 비율을 여러 번 시험하지 않는다.

Angular panel 수도 같은 local chart의 실제 stretch로 산출한다. 전체 source의 최소 반경을 모든 각도/반경에 일괄 적용하지 않는다. 전체 계획 개수를 먼저 계산하여 **기존 combined panel/query 예산**과 비교한다. 초과하면 새 계획도 pending으로 끝낸다. normalization·tau·차수·cap을 바꾸지 않는다.

추가 전달된 exact reference 분류에서 origin-separation 실패 source130/182/185/257/270/314/317과 inverse 실패 source133은 원점을 closed triangle에 포함한다. 따라서 이들에게 source 전체 거리 하한은0이다. 이것만으로 실제 clipped overlap도 원점을 포함한다고 판단하지 않는다. B의 oblique 실패 source들은 양의 source 거리 하한이 있다는 증거도 재사용할 수 있으나, 그 사실은 투영방향의 transversality 증명을 대체하지 않는다. 이 추가 수치는 운전자의 저장 exact predicate 증거를 참조하며 검토자가 재계산한 결과는 아니다.

## 3. 종료 조건과 다음 실제 전달 경계

한 번의 묶음 결과는 35개 각각의 (i) 원래 실패 typed witness, (ii) 실제 원인 판정과 충분조건, (iii) 수정 전후 domain/event 동일성, (iv) 새 계획 또는 명시적 pending, (v) 비용·새 증거 범위다. 기존 정상 typed panel, 68526 출판 packet, current/generating state와 heat, 모든 tau의 hash가 유지되어야 한다.

성공은 같은 실제 domain을 누락·중복 없이 덮고, branch/order/half-open ownership 및 기존 고정 예산을 만족하는 atlas 계획이다. 아직 계산하지 않은 analytic remainder를 통과로 기록하지 않는다. 35개 중 하나라도 미해결이면 그 witness를 보존하고 이번 제한 묶음을 종료한다. 자동 세분화·cap 증가·전체 재실행을 이어가지 않는다.

모두 해결되면 다음 substantive 경계는 전 target의 실제 analytic certificate와 재사용 가능한 overlap packet 수치조립이다. source outgoing과 target incoming은 동일 packet을 참조하고, retained whole4의 group packet을 중복 합산하지 않는다. 전 source closure와 각 target의 원래 cold tau 및 최종 R80 변환을 함께 요구한다. 실제 새 총 query/메모리/증거 I/O 비용이 확정되기 전에는 전체 quadrature 실행을 열지 않는다.

이 순서는 현재 실패의 원인을 확정하고 실제 transactional 전달로 연결하기 위한 한정 작업이다. 일반 evolved thermal state, native CUDA/MPI/OpenMP 시간 갱신과 t6/8/10 및 논문 최장공전 안정성은 아직 검증되지 않았다. W0 개별 고정질량 환원은 승인된 coupling law가 해당 상태에서 만족될 때의 조건이며, 이 atlas 완료로 자동 성립하지 않는다.
