# Astra — compact 공간 기각과 원래 field의 정확 cell condensation

2026-10-22. SHARED_CONTEXT 전체, 실제 compact solve/curvature JSON, 네 rank worst witness, 71,080셀 failure analysis, compact native sparse GPU/solver 구현을 읽었다. 원본 HEAD8e25dfe/연구 ROOT 분리 계약을 유지한다. 계산·폴링·운전자 코드 변경 없이 이 문서만 작성한다.

## 1. 종료 판정

compact 후보는 종료한다. max chi=.00879341은 .002의4.397배이며, 160개 실패 중112개가 interface 비인접이다. 평균 gate 통과와 좋은 energy/torque ledger는 이 공간 오류를 상쇄하지 않는다. radius/order/residual 수를 조정해 재시도하지 않는다. CPU binary80 대조가 일치하므로 단순 GPU 산술 결함으로 돌리지 않는다. compact 공간의 진화는 실행하지 않는다.

새 assembled sparse DD backend와 작은 Schur 구현은 유효한 성과다. 1.58M DOF의 약2.07초 solve는 다른 연산자의 결과이므로 원래10.8M 공간의 속도향상으로 보고하지 않는다. 그러나 **원래 함수공간을 바꾸지 않고 실제 global linear solve를 작게 만드는 구현 경로**를 가능하게 하는 증거다.

이번 선택은 원래 검증된 step5 state의 **physical source cell 내부 가속도 exact static condensation + sparse DD skeleton solve**, 그리고 동일 selected rule의 불변 enriched reference evaluation 캐시다. force/heat/기저/quad/state를 바꾸지 않는다. 새 IC나 compact trajectory를 시작하지 않고 원래 step5에서 이어간다.

이는 이전 “10.8M 공간의 무한 최적화 종료”를 모호하게 무시하는 것이 아니다. 새 sparse backend 증거에 근거해 **global solve 차원을 정확하게 줄이는 단 한 번의 제한된 전략 개정**이다. 실패하면 matrix/cache 방식의 추가 변형을 자동 탐색하지 않는다. 여전히10.8M physical state DOF를 저장하므로 field 평가와 material 수명 문제는 남는다.

## 2. 명확한 수학: 내부 state는 없애지 않는다

기존 canonical basis를 유지하고 DOF를 다음처럼 순열로 나눈다.

- Gamma: physical source-cell 경계에 nonzero trace가 있는 DOF 및 여러 source가 공유하는 DOF.
- I_c: source c 내부에만 지지되고 그 source의 전체 외부 trace가 정확0인 DOF. 내부 triangulation의 node/edge/residual과 bubble이 포함될 수 있다.

값이 작다는 이유로 I_c로 분류하지 않는다. 실제 basis support/canonical trace identity로 판정한다. 서로 다른 cell의 I_c 사이 M coupling은0이어야 한다. exterior-boundary trace DOF도 Gamma에 둔다. 단순 triangle interior만 제거하는 것이 아니라 **physical cell 전체 내부**를 제거해야 의미 있는 규모 감소가 가능하다.

원래 매 RHS 방정식은

    [M_GG M_GI; M_IG M_II] [a_G; a_I] = [f_G; f_I].

M_II는 cell별 block diagonal이다. 한 번 factor하여

    S_c = M_GG,c - M_GI,c M_II,c^-1 M_IG,c,
    b_G = sum_c (f_G,c - M_GI,c M_II,c^-1 f_I,c),
    (sum_c S_c) a_G = b_G,
    a_I,c = M_II,c^-1 (f_I,c - M_IG,c a_G)

를 수행한다. 큰 inverse를 명시적으로 만들지 않고 Cholesky/triangular solve를 쓴다. S_c는 binary80에서 구성해 DD로 옮기며 원래 M을 authoritative하게 유지한다.

**x_I와 v_I는 매 RK stage에서 독립적인 원래 상태로 유지하고, 위에서 복구한 a_I로 갱신한다.** x_I를 x_G의 정적 함수로 치환하거나 interior mode를0으로 만들지 않는다. explicit RK에서 acceleration의 linear solve를 정확히 소거하는 것이므로 동적 field/주파수/에너지 저장을 없애는 reduced model이 아니다. 별도 spin/kinetic reservoir도 없다.

수치 Schur와 local factor의 rounding은 존재한다. 복구한 전체 a에 대해 기존 binary80 M의 `||Ma-f||/||f||<1e-16`을 직접 확인한다. Gamma residual만으로 종료하지 않는다. 필요하면 동일 block inverse를 사용한 원래 full residual defect correction을 기존 유한 cap 안에서 수행한다. diagonal jitter/negative pivot repair/tolerance 완화는 없다. 기존 CPU/reference apply와 저장 mass block을 재사용한다.

## 3. 구현과 데이터 재사용

### 한 번의 구조 구축

원래 full owned topology, source IDs, conditioned basis, fixed M와 step5 checkpoint를 읽는다. compact T나 실패한 compact acceleration을 원래 상태에 적용하지 않는다. parent midpoint의 비공선 rounding 문제도 이 방식에는 새 주입을 만들지 않으므로 도입되지 않는다.

각 source의 I/Gamma list, local dense M_II factor, M_IG, S_c를 구축한다. owner가 source를 소유하고 skeleton shared DOF만 halo한다. source-local RHS 축약과 내부 복구는 GPU batched DD, global skeleton은 이번 compact에서 검증된 assembled CSR/DD owned-dot/ghost-pair 구조를 재사용한다. 원래 local CSR를 전 rank 복제하지 않는다. rank2 초기 enrichment는 원래 field에 없는 것이므로 compact의 두 residual border를 옮기지 않는다.

구축 전에 정확한 skeleton DOF/CSR upper pattern, cell별 interior/trace dimension, factor/CSR/cache bytes를 산출한다. 예상치로 통과시키지 않는다. 목표는 skeleton<=2.5M scalar DOF이다. 이 상한을 못 맞추거나 요구 메모리가 실제 장치/host 예산을 넘으면 이 묶음의 비용 실패로 종료한다. 임의 cell merging이나 field 삭제로 맞추지 않는다.

새 mass 적분, 초기 curvature, 원래20억 점, 완료된 transfer/atlas는 재실행하지 않는다. 실제 source-block factor와 Schur 구축은 저장 M에 대한 대수 작업이다. 고정 setup wall 한도8시간을 두고 부분 산출물/원래 state를 보존한다. source factor가 양의 definite가 아니거나 원래 support partition을 증명하지 못하면 실패다.

### Stage 간 불변 reference cache

현재 stage-local cache는 매 stage마다 같은 비다항 basis의 radial tails/edge extensions와 generating field를 다시 계산한다. 다음 캐시는 **원래 selected point ID에서 불변인 enriched N, reference gradients, 원래 density/weight, generating j0**만 한 번 저장한다. current x/v/J/P/force는 저장해 재사용하지 않고 실제 stage마다 계산한다.

P3 부분은 동일 scalar arithmetic으로 저렴하게 평가할 수 있고, 비다항 enrichment 부분을 double 그대로 캐시한다. 낮은 precision/근사압축/새 query rule은 없다. 원래 point ordering과 reductions를 유지한다. cache key는 source/basis/rule/생성열 SHA와 ABI를 포함한다. 다른 epoch/rule로 잘못 사용하는 경우 즉시 거부한다.

대략 수억 selected points에 enriched value/gradient를 저장하면 수십 GB가 필요하다. 실제 n_enriched와 point count로 bytes를 계산한다. H200 용량과 MPI4 합산 memory를 기준으로 하며 rank마다 global cache를 복제하지 않는다. device 사용량의20% headroom을 남긴다. Ada에 들어간다고 가정하지 않는다. resident가 불가능하면 고정 source-tile의 pinned-host cache를 사용할 수 있으나 실제 PCIe 비용을 포함해 판정하고 cache partition sweep은 하지 않는다.

캐시를 만드는 한 번의 순회는 기존 selected reference 함수의 저장이며 새 physical trajectory나 새로운 quadrature는 아니다. immutable 함수 결과와 동일한 double bytes를 저장한다. 실제 다음 RHS의 제한된 고정 tile에서 cache/direct 동일 reduced 결과를 확인한다. 원래 accepted 구간을 다시 실행하지 않는다.

## 4. 관문은 새로운 작은 물리 문제가 아니라 다음 실제 RHS다

복구한 전체 가속도가 원래 M 방정식을 만족하면 초기 공간 정확도를 별도로 다시 검사할 이유가 없다. field/state/force가 동일하므로 원래 pressure curvature를 새 공간에 상속하는 문제가 아니다. 다음 실제 step6의 첫 RHS에서:

1. source-local 축약/전체 복구와 원래 full residual, 양의 local/skeleton solve, MPI owner 회계를 확인한다.
2. 동일 actual f에 대해 기존 resident full solve를 한 번만 shadow하여 work/U, torque/2K, acceleration difference를 기존 관문으로 비교한다. 초기 operator/force를 다시 만들지 않는다.
3. cell energy/torque/hybrid work에는 복구된 모든 interior acceleration을 사용한다. 축약된 Gamma만으로 spin/K를 계산하지 않는다. interior의 virtual reaction이 roundoff 수준인지 확인한다.
4. field cache/direct 비교와 full-state positivity/finite 검사를 수행한다. 실패하면 step5와 새 실패 witness를 보존하고 종료한다. h/quad/order/space를 바꾸지 않는다.

통과한 동일 step을 commit하고 실제 restart readback을 확인한다. 연속 경로는 원래 t=.00011263377801805954에서 이어지며 새 t0 branch가 아니다. checkpoint에는 원래 전체 x/v/entropy/M provenance와 optional factor/cache references가 들어간다. factor/cache는 재생 가능한 파생 자료이지 물리 state의 대체물이 아니다.

## 5. 성능 관문과 한 번의 연속 실행

기존 selected step의 mass 약70.6초와 prepare 약46.1초를 모두 줄여야 한다. mass만2초가 되어도 t10 경로가 자동 해결되는 것은 아니다. 최초 실측에 setup/one-time shadow를 분리하되 steady stage 준비, halo, solve, 내부 복구, 원래 true residual, endpoint 및 필요한 map guard 비용을 포함한다.

기존 sustainable h≈.000345838을 비용 사전 추정에 쓰고 개선을 가정하지 않는다. 이 h에서 t10<=7일을 만족하려면 지속 step 약20.9초 이하여야 한다. 이는 달성 예측이 아니라 채택해야 할 비용 문턱이다. 실제 h가 더 작으면 그 값으로 문턱을 갱신한다. compact의2초를 그대로 새 full-field solve 시간으로 복사하지 않는다.

- 첫 실제 step과 이어지는 정상 몇 step에서 t1<=24시간 및 t10<=7일을 모두 만족해야 장기 승격한다. 만족하지 못하면 그 accepted endpoint에서 종료한다. “조금 빨라졌으니 일단 며칠 돌림”은 하지 않는다.
- 통과하면 같은 trajectory를 .1/.2/1/6/8/10의 물리 출력과24시간 단위 durable checkpoint로 진행한다. 누적 wall 상한7일이며 실제 dt/비용 악화로 예산을 넘으면 종료한다. 작은 시간별로 승인 재질문하지 않는다.
- 이전 문서의 pressure/density/velocity/phase engineering ceiling과 더 엄격한 기존 물리 gate를 유지한다. E/L 보존이나 고정 entropy만으로 가짜 압축 가열 제거를 선언하지 않는다. Jensen bulk/variance와 실제 rho/P/궤도를 함께 기록한다.
- `global_map_injectivity_guard_pending`을 이 대수적 성공으로 지우지 않는다. 연속 장기 승격에는 기존 shear-lifetime 문서의 전체 enriched-map/trace/vertex-fan/외부 boundary 계약을 적용해야 한다. 명시적 unresolved와 실제 fold를 구별하고, guard 비용도 위 비용에 포함한다. 보호되지 않은 long run을 허용하는 예외는 아니다.

setup·정확성·비용 중 하나라도 실패하면 이번 정확 condensation 후보도 종료한다. 새로운 공간 반경/차수/enrichment 또는 cache/solver 조합을 자동 추가하지 않는다. 이번 한 번의 구조 변경이 비용을 못 닫으면 이 고정 material 방법을 t6–10의 production 경로로 더 승격할 근거가 없다.

## 6. 장기 물리와 원래 본체에 관한 판단

정확 condensation은 기존 공간의 곡률/회전 표현을 손상시키지 않지만 material shear 수명을 늘리지 않는다. 실제 t8 near-core condition≈1e5는 여전히 위험 근거다. 시간·quadrature 조정으로 감긴 map의 공간 표현 한계를 없앨 수 없다. 다음 실제 trajectory가 그 한계에서 막히면 마지막 상태를 보존하고 원인을 geometry/accuracy/cost로 분류한다. 그 시점에는 실제 conservative relabeling/remeshing과 recipient thermal functional이 필요하며, 현재 source-context fragment를 새 reference reset이라고 부르지 않는다.

원래 t8 CUDA snapshot은 역사적 실패의 공간·열·궤도 기준이다. 이미 생긴 가열을 새로운 spin으로 재분류하거나 그 snapshot의 U를 식혀서 재시작하지 않는다. 원본 CUDA/MPI/OpenMP 진입·실행·geometry 구현을 쓰되, 현재 material/Pext 경계와 원래 periodic/Voronoi 모델의 차이는 계속 공개한다. W0 dynamic Voronoi/개별 고정질량 환원과 weight-linked mass relaxation, shock entropy는 아직 해결되지 않았다.

force/acceleration 또는 알려진 정적 함수의 basis 사용을 일괄 금지한 것은 과거 방법 선택이지 사용자 직접 금지조건이 아니다. 자유 coefficient와 정확한 variational/time-basis 항을 갖는 enrichment 자체가 곧 steering인 것은 아니다. 그러나 이번 실패160셀을 지우기 위해 force basis를 하나 더 넣는 선택은 하지 않는다. 초기 가속도 일치만으로 전단 수명과 최장 공전이 해결되지 않기 때문이다. 사용자의 직접 조건인 no steering/Re∞/열·독립 spin 저장고 금지/원본 native 경로는 그대로 유지한다.

## 7. 이번 단일 실행 계약

이 문서는 실패한 compact fresh branch의 실행 계획을 대체한다. 원래 검증된 step5를 입력으로 exact cell condensation과 불변 reference cache 한 묶음만 구현·실행한다. 성공 시 실제 연속 physical trajectory로 즉시 연결하고, 실패 시 종료한다. 최종 보고는 actual skeleton/CSR/bytes, 원래 full residual과 cell work/torque, 실제 지속 비용·dt, 도달한 시간과 종료 이유를 한 번에 묶는다. t6/8/10과 GIZMO 최장 공전 목표는 미달이며 어떠한 local pass도 이를 대신하지 않는다.
