# Astra: native runtime 전달 종료 평가와 다음 실행 묶음

2026-10-22. SHARED_CONTEXT 전체와 원본 HEAD 8e25dfe를 확인했다. 지정 endpoint/input 정책, runtime reference·binary80 field/stream 경로 및 whole/real-root/input 감사 결과를 읽었다. 새 계산·폴링·driver 코드 변경 없이 이 문서만 작성한다. 원래 geometry bridge와73536 source-pair 완료는 유효한 선행 증거로 유지한다.

## 1. 이번 묶음 판정

실제1441792 binary80 queries,64 source-overlap,4 target의48 tau·양성·최종 변환 검사를 통과한 whole batch를 채택한다. retained68526을 포함한5 target는 완료 자료로 재사용한다. 이들을 다시 계산하지 않는다.

Regular6521 groups/104336 nodes와 oblique274 groups/4384 nodes의 실제 branch cache는 해당 interval-Newton/유일성·함수 전제를 유지한 저장 runtime 입력으로 채택한다. 최종 outward 확장263건은 숨기지 않고 실제 입력 radius에 포함한다. 단지 이전 구간과 겹친다는 사실로 root identity를 정하지 않는 기존 계약을 유지한다.

438120개 input-only 검사에서 legacy barycentric center 경로의4272개 edge sum 미분리가 direct reference 경로에서 해소된 것은 **동일 좌표의 안정적 입력 표현을 선택할 충분한 근거**다. 전체21629184 남은 query의 field/all12 통과나 모든 입력의 양성 증거는 아니다. 입력 표본에 포함되지 않은 no-implicit clipped 부분도 실제 실행에서 똑같이 검사해야 한다.

## 2. Endpoint uncertainty를 runtime 오차에 넣는 계약

실제 root/event parameter theta를 포함하는 원래 구간과 같은 branch identity가 있고, analytic remainder가 그 parameter family에 균일하다면,

    |I_theta−Q_theta|≤Ranalytic,
    |Q_theta−Qcomputed|≤Σ Enode+Esum+Econversion

로 합산하는 방식은 타당하다. Enode에 원래 좌표·root·Gauss·measure의 uncertainty와 field/native forward 오차를 모두 포함하면 별도의 같은 endpoint 띠 오차를 다시 더할 필요는 없다.

단, 이는 빠진 domain, 틀린 root identity 또는 discarded cap을 교정하지 않는다. 실제 topology/coverage는 기존 완료 증거에 의존한다. 새 runtime cache는 실제 theta를 포함해야 하고, analytic proof가 사용한 동일 함수/parameter family와 연결되어야 한다. 실제 root를 포함하는 원래 집합보다 최종 serialized ball이 아주 조금 커진 경우, analytic proof를 자동 확대했다고 주장하지 않는다. 실제 root의 원래 family 포함 증명은 유지하고, 더 큰 runtime ball에 대한 native denominator/forward guard는 실제로 검사한다. 그 ball의 불필요한 부분까지 analytic remainder가 성립해야 하는 것은 아니지만 native error graph가 평가에 사용한 전체 ball을 유효하게 감싸야 한다.

좌표와 measure에서 같은 theta를 각각 독립 interval로 감싸는 것은 상관을 잃는 보수적 과대평가다. 불확실성을0으로 만들거나 enclosure 중심만 넣는 것은 허용되지 않는다. Source outgoing/target incoming은 같은 overlap 값·오차 record를 사용한다.

## 3. Direct reference lambda 입력: 구현을 채택할 조건

선택은 **legacy default를 보존하는 명시적 direct-reference 입력 모드**다. Duffy/regular-q/oblique에서 이미 계산한 같은 reference point의 lambda를 binary80 center로 직접 변환한다. X,Y,lambda,weight는 같은 source/chart/branch/Gauss input record에서 함께 얻는다. 단순히 세 lambda를 정규화하거나 마지막 하나를 임의로1−나머지로 덮지 않는다.

현재 r80::evaluate는 q.l의 center뿐 아니라 내부에서 sf::barycentric(m,X,Y)를 다시 계산해 오차 enclosure를 정한다. **center만 direct로 바꾸고 이 낡은 enclosure 경로를 그대로 두는 수정은 충분하지 않다.** Optional 모드는 원래 chart lambda enclosure도 명시적 인자로 받아 cache(q.l[i], original_lambda[i])에 사용해야 한다. legacy 모드의 원래 경로는 남긴다.

원래 Cartesian gradient cache와 qnode/edge caches 및 동일 analytic basis는 그대로 사용한다. 변화는 좌표를 거쳐 되찾던 affine lambda의 계산 순서와 forward error 모델뿐이다. exact lambda와 exact X,Y는 같은 reference point에서 정의되며, 각각의 ball이 그 joint assignment를 포함해야 한다. 둥근 center들이 sum(lambda)=1 또는 X=Σlambda vertex를 마지막 비트까지 만족하지 않아도 forward graph가 동일 정확 함수까지의 차이를 포함하면 허용된다. 이 사실은 state refit 또는 서로 다른 물리 point 사용의 허가가 아니다.

실제 모든 query에서 edge sum, radical/tail 분모, J 및 물리U 등 기존 guard를 그대로 적용한다. 입력 모드만으로 나머지 실패를 통과 처리하지 않는다. 입력 감사4272 witness는 기존 old/new 차이와 전체 radius를 보존한다.

## 4. 다음 한 묶음: direct 모드 검증에서 남은19 실제 전달까지

아래 조건을 차례로 만족하면 같은 묶음 안에서 driver가 재승인 없이 진행한다. 새 atlas/차수/정밀도/tau 변경은 하지 않는다.

### A. 입력 모드의 최소 통합 증거

이미 저장된4272 실패 입력의 exact 감사와 optional mode 구현을 마무리한다. 추가 광범위 입력 scan은 하지 않는다. 같은 원래 reference point에 대해 direct lambda balls, native X/Y cache와 gradient provenance를 연결하고 sum·barycentric reconstruction 관계가 원래 analytic point를 포함함을 검사한다.

이후 아직 계산하지 않은 실제 mixed-query stream 안에서, 실패 입력 유형과 whole/regular-radial/regular-angular/oblique 각 경로를 포괄하는 고정 소수 Q128 oracle 위치를 **실행 전에** 정한다. 같은 query에서 native field17 및 all12 값을 oracle와 bound로 대조한다. Oracle는 per-query forward graph를 대신하지 않으며 실패점만 사후 선택해 검증하는 방식도 아니다. 완료된5 target의 성공 queries는 이 검사를 위해 재생하지 않는다.

### B. 저장 실제 root/cache를 소비하는 streaming

남은19 target/21629184 queries의 원래 순서와 panel 계획을 사용한다. Regular/oblique 실제 branch cache를 읽고 root·inverse·geometry를 다시 풀지 않는다. Cached node/branch/parameter ID와 source/target/epoch/hash를 묶고, Gauss weight 및 full endpoint measure radius를 native 입력까지 전달한다.

한 query가 성공하면 실제 all12 중심/오차, reduction-tree 상태와 committed 범위를 보존한다. Prefix 성공을 다음 실패 때문에 다시 실행하지 않는다. 재개 시 input/evaluator/arithmetic-policy hash와 정확한 next query index를 확인한다. 단순히 마지막 파일 길이로 완료를 추정하지 않는다. 실제 point record와 accumulator의 durable commit 경계를 맞춰 중복 합산·누락을 막는다. 새 input mode는 manifest에 기록하고 legacy 성공 target의 기존 모드를 소급 변경하지 않는다.

Source-overlap packet은 같은 값과 오차로 outgoing/incoming에 연결한다. Analytic remainder, 실제 native forward, 실제 합산 및 final R80 conversion을 원래 target tau에 합산한다. 전체 actual m,U 양성 및 source m,p,L,K,U closure를 검사하고, 모든 target가 닫히기 전 전체 caller 상태는 commit하지 않는다. retained whole4의 group packet을 source별로 임의 배분하지 않는다. 필요한 source 회계는 저장 positive records에서만 복구한다.

누적된 이미 확정된 오차 하한/subtotal만으로 tau를 넘으면 해당 target의 불필요한 나머지 query를 중지할 수 있다. 실패는 실제 point/branch/input/forward witness로 보존한다. Guard를 완화하거나 자동 precision/degree/atlas 변경을 하지 않는다. 정상5 target와 다른 성공 prefix는 유지한다.

### C. 묶음의 종료 결과

성공은 실제24 target와 전 source의 전달·양성·원래 오류 budget 및 transactional commit을 모두 닫은 결과다. Whole-only 성공 또는 direct-input 양성만으로 이 결론을 내리지 않는다. 실패 시에는 남은 target/query 수, 지배 오류 항, caller 보존 및 재사용 가능한 prefix를 한 번에 보고한다. primitive별 반복 review 대신 이 한 묶음의 실제 결과로 다음 판단을 한다.

## 5. 그다음 물리 진화 연결 순서

전달을 닫은 뒤 바로 일반 evolved heat/state와 native timestep 연결로 진행해야 한다. 현재 same_control frozen recipe를 다음 step의 열 초기화로 쓰지 않는다. 물리U의 stress-work/수송 이력과 실제 resolved velocity/shape DOF 또는 유지 fragment를 restart payload에 함께 담는다. 12 moment만으로 고차 속도장을 유일하게 결정할 수 없으므로 mean compression이나 kinetic reservoir로 빈 부분을 채우지 않는다.

Native RK 단계마다 geometry/state/heat가 같은 stage를 나타내고, 재격자/전달·MPI ownership·rollback까지 함께 복원되는 경로를 구현한다. 기존 CUDA/OpenMP/MPI 본체와 성공한 연산자들을 재사용하되 binary80 reference 성공을 production GPU field 구현 성공과 혼동하지 않는다. 실제 한 accepted step 및 restart 연속성, 보존/양성/원래 물리 비교를 통합 검증한 뒤 연속 실행에 연결한다.

그 경로에서 기존 전단 가짜가열·열고갈이 줄어드는지를 실제 density/entropy/물리U, 궤도·phase, E/L 및 비용으로 판단하고 t6/8/10 checkpoint를 확보한다. 양성 또는 전역E만으로 성공을 선언하지 않는다. 최종 GIZMO 최장공전과 shock 열역학 검증 목표도 유지한다. 현재 결과는 실제 전달의5 target 완료이며, 아직 장시간 물리 안정성 증거는 아니다.
