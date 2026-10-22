# Astra: shared Duffy 전달 거절 평가와 다음 단일 묶음

2026-10-21. 계획·평가 담당 Astra. 운전자가 구현·실행한다. 지정 소스와 v2/v3 계약을 정적으로 읽었으며 계산·작업 조회는 하지 않았다. 원본 native 저장소 origin LagEunha, HEAD 8e25dfe 및 운전자 변경을 확인했다.

## 1. 종료 판정

v3도 정당한 REJECT이다. 같은 384 source에서 344 제외, 4 전체 포함, 36 교차 후보. canonical edge/plane 실패는 11→0, 미확정 실제 영역은 30→26으로 줄었다. 그러나 40 활성 source 모두 all12 미완료, 인증 적분점 0이다. 2304 witness 점은 적분 완료 증거가 아니다. 원본 state/heat와 caller 보존은 통과했다. 이 결과는 물리 해의 부재나 remap 불가능성의 증거가 아니다.

295296 AD2 비교와 Duffy의 실수 J/분모 검사는 유지한다. 817152점 원본 증거도 재생하지 않는다. 평면을 canonical edge에 먼저 합성해 사건을 해결한 변경은 같은 함수를 안정적으로 계산한 것이다. v3의 완료 사건·영역 기록은 입력/소스 해시가 맞는 범위에서 재사용한다.

현재 코드의 `build()`는 아직 physical fan의 star anchor를 영역 완결 조건에 묶고, 적분기는 교차 영역을 모두 Cartesian 복소 역상으로 보낸다. 따라서 앞 계획의 **clipped reference-Duffy 경로가 아직 구현되지 않았다.** 다음 묶음은 이 누락을 연결한다. 새 PDE, 상태, 기저 또는 적분 정확도 기준을 선택하지 않는다.

## 2. 수학·코드 점검

### 방향과 영역

읽은 코드에서 즉시 확인되는 뒤집힘 결함은 없다. source 방향을 뒤집을 때 canonical parameter 폭은 음수이며 tangent에는 폭, second derivative에는 폭의 제곱을 곱한다. target의 외향법선 n에 대해 정렬 접선 (-ny,nx)는 허용영역을 왼쪽에 둔다. 이 계약을 유지하고 각 loop의 방향·공유 경계 상쇄를 실제 기록으로 확인한다. 임의 abs(J)로 뒤집힘을 감추지 않는다.

`membership()`의 실패는 outside가 아니다. whole-displacement 제외가 안 되면 곧바로 strict-interior Newton을 요구하는 현재 방식은 경계·꼭짓점·가까운 외부점을 처리하지 못한다. `empty arrangement`에서 같은 coarse single-plane 제외를 다시 요구하는 것도 완전한 공집합 판정이 아니다.

### 복소 실패의 의미

out.X를 독립적인 세 lambda*vertex의 합 대신 `A+s*((B-A)+t*(C-B))`로 계산하고, displacement/derivative도 기준 affine 부분과 작은 잔여항을 먼저 합성하는 변경을 채택한다. 상태 변경이 아니다. sigma<.75 및 positive-half-plane 검사의 실패는 일반적으로 상계/좌표계 실패이며 실제 특이점의 증명이 아니다. 실제 J=0, 반대 꼭짓점의 Cartesian 비해석성, 밀도·속도 함수의 복소 pole/branch와 구별해 기록한다.

특히 원점 chart에서 r=s R(t)이면 원래 밀도의 r+epsilon 분모는 복소 s=-epsilon/R(t)에 실제 pole을 갖는다(실제 numerator 소거가 증명된 경우 제외). 고정 whole radial Cauchy 근방이 이를 포함하면 affine cancellation만으로 해결할 수 없다. 실수 적분은 유한해도 해당 해석적 나머지 증명은 불가능하다. `sqrt` 인수가 계측되지 않은 행의 0은 관측값으로 해석하지 않는다.

### 열·native 오차 범위

현재 physical_query는 frozen CONTROL과 current map이 같은 경우에만 Uref=P0/(gamma-1)*J_generator^(1-gamma)를 current f.J로 평가해도 같다. 다음 API에서는 generating map과 current map을 구별하거나 두 fingerprint의 동일성을 assert한다. 일반 진화에서 이 함수를 그대로 사용해 U를 재생성하면 안 된다.

Udot=-P*(Jdot+delta_dot), P=(gamma-1)*Uref/(J+delta)는 이번 고정 geometric-EOS context에 맞다. rho0/J의 수송 밀도와 rho0/(J+delta)의 EOS 밀도를 혼동하지 않는다. v2/v3 context byte 일치를 유지하고 delta/delta_dot의 binary 값은 정의된 고정 입력으로 보존한다. 이는 원본 이산 context의 구현 정의이지 연속 적분 오차 0이라는 주장이 아니다.

새 좌표계의 실제 query에서도 cache/coordinate/native 연산 그래프 bound를 적용한다. 기존 점 증거는 임의 offnode의 실측 증명이 아니다. libm 8ulp는 명시된 조건부 engineering model이다. all12에는 동일 양의 measure와 generating heat, Udot context를 사용한다.

## 3. 선택한 다음 묶음: 경계로 소유권을 확정하고 clipped reference chart로 적분

### A. 소유권과 coverage를 적분법에서 분리

1. canonical edge의 exact-composed plane 함수/endpoint/접선 사건을 재사용한다. 인증된 ray crossing 또는 winding으로 source 이미지 안/밖/경계를 반환한다. ray를 지나는 꼭짓점은 half-open endpoint 규칙으로 한 번 계산하고, tangency는 방향 변화 여부를 구별한다. source별 L<1과 J>0에 의한 injectivity를 전제로 한다. Newton 실패로 밖이라고 분류하지 않는다.
2. source-edge fragment, target segment, shared endpoint를 같은 사건 ID로 연결한다. 교차가 없을 때도 source의 내부 증인과 target 내부 증인의 containment를 검사해 disjoint/source-contained/target-contained를 구별한다. 일치 경계와 접촉만인 경우는 별도 상태로 남긴다. 불명확하면 reject한다.
3. `Domain`의 완결 조건은 방향 있는 boundary/containment/coverage이다. 단일 star-shaped fan 존재를 필수 조건으로 두지 않는다. 물리 fan은 인증 가능한 부분에서만 선택하는 적분 경로이다.

### B. 실제 clipped reference-Duffy 영역

기존 고정 6개 chart에서 f_j(s,t)=n_j·chi(X(s,t))-b_j<=0을 직접 사용한다. 원래 chart 경계와 target pullback curve의 사건으로 양의 영역을 나눈다. 각 부분에서 transversal 조건을 인증한 좌표로 s_minus(t),s_plus(t) 또는 t_minus(s),t_plus(s)를 추적한다. fold/꼭짓점 사건에는 이미 구현한 사건 분할·제곱 좌표화를 연결하며, 단순히 분모에 epsilon을 넣지 않는다.

각 radial graph 부분은 s=s_minus(t)+eta*(s_plus(t)-s_minus(t))로 [0,1]^2에 옮긴다. Jacobian은 reference Duffy Jacobian과 양의 구간 폭을 포함한다. 두 경계가 붙는 사건은 event chart가 처리한다. s=0에서 2차원 Duffy 역함수를 요구하지 않는다. 모든 target plane 및 chart 소유권이 해당 영역 전체에서 증명되어야 한다.

새 물리 family/큰 증명 체계를 만들지 말고 기존 canonical 사건, fold, interval implicit-function 부품을 연결한다. source118/254의 membership과 source246의 apex 포함 fan이 이 경로의 실제 필수 사례다. 새 toy target은 만들지 않는다.

### C. 실제 singularity에 맞춘 유한 적분 atlas

Gauss16/기존 나머지 공식·전체 오차 기준을 유지한다. 같은 좌표계에서 먼저 cancellation-preserving affine/잔여항 평가를 적용한다. 그래도 실제 pole/branch가 Cauchy 근방에 들어오면 **입력 함수의 특이점 거리와 사건 위치로 정한 유한 subpanel atlas**를 계산 전에 만든다. 이는 차수·정밀도·cap·guard scan이 아니다. 각 panel의 고정 Cauchy 확장 전체가 유효한지를 증명한다. density r+epsilon, softened velocity/potential, logistic, J 및 extension 분모를 모두 포함한다.

원점의 density scale은 실제 함수의 scale이므로 radial chart를 그 scale에서 분리할 필요가 있을 수 있다. 오차가 작아질 때까지 임의 반복 분할하지 않는다. deterministic atlas가 기존 총 작업 예산을 초과하면 실행 전에 cost/domain reject한다. 실제 pole이 포함된 box를 과대평가 문제로 표기하고 무한 정련하지 않는다.

4개 whole-inside 영역은 이 동일 atlas의 unclipped 경우다. 현재 192 axis-box 실패 witness를 사용해 변경 범위를 고정한다. 루트/함수 preflight를 다시 전부 돌리지 않고, 바뀐 evaluator와 실제 panel에 필요한 신규 증거만 만든다.

## 4. native 호출·종료 관문

동일 all384→target68526 호출 하나가 목표다. 완료된 v3 사건·영역과 기존 원본 증거는 해시 확인 후 재사용한다. 분류·영역·적분·R80 직렬화는 별도 완료 플래그를 갖는다.

- 원본 state/mass/heat fingerprint 불변. 새 EOS 초기화, trace averaging, refit 없음.
- all384 coverage, 12개 모멘트의 양의 공통 packet, source/target ownership 및 shared-edge 취소 확인.
- 모든 실제 적분영역에서 branch/역상 또는 reference chart와 analytic remainder 완결.
- native/geometry/analytic/sum 오차를 기존 whole-target tau와 비교. **R80 value 변환 오차와 최종 outward error 변환을 먼저 더한 뒤 최종 gate**를 적용한다. 그전 gate만 통과한 packet은 publish하지 않는다. 이전 v2/v3는 전부 reject라 이 누락으로 잘못 publish된 결과는 없다.
- 기존 모든 tolerance/normalizer 값을 manifest에 고정한다. 새 분류나 sampled velocity scale로 유리하게 재정의하지 않는다.
- 누락·범위 불명확·예산 초과·오차 초과 하나라도 있으면 caller byte 불변으로 종료한다. 성공한 일부 source를 전체 target 성공으로 표시하지 않는다.

종료 후 한 번만 평가한다. 실패하면 첫 불완전 domain과 실제 장애 종류(coverage, 진짜 특이점, enclosure 과대, 적분 bound, 비용)를 보존한다. 같은 whole call의 guard/차수/precision/cap scan을 이어가지 않는다. 통과하면 full-target 보존 조립·가중치 법칙 연결로 진행하며, 한 target 성공은 시간진화나 t6/8/10 또는 최장 공전의 성공이 아니다. no steering, Re=infinity, 물리 U/회전 보존 및 조건부 W0 계약을 유지한다.
