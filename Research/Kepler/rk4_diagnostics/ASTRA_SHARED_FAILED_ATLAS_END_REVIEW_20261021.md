# Astra: 실패 atlas 제한 묶음 종료 평가와 다음 구현

2026-10-21. SHARED_CONTEXT 전체와 지정 witness/repairs/origin-distance 장부, inverse 및 direction 구현을 읽었다. 원본 repository는 LagEunha, HEAD 8e25dfe임을 확인했다. 새 계산·폴링·코드 변경 없이 이 계획 문서만 작성한다.

## 판정

이번 묶음은 계획한 범위 안에서 유효한 진전을 만들었다. 비용 실패 2개는 완료로 채택한다. 나머지 33개는 방향·inverse·원점 분류 증거를 얻었으나 실제 arrangement/atlas는 아직 미완료다. 39/39 witness 일치와 정상 chart 재생0, 원래 cost prefix 보존은 원인 귀속의 근거다. 독립 감사의 범위는 저장 부등식/ABI/회계이며 기존 전체 함수 증명을 새로 대신하지 않는다.

검토한 수식에서 현재 결과를 뒤집는 치명적 결함은 찾지 못했다. 다음 세 구별을 유지해야 한다.

1. **물리 power vertex 0과 canonical source vertex의 image는 다르다.** chi(reference0)는 약 (2.36e−21,1.05e−22)이므로 power vertex를 source vertex로 snap하거나 사건 ID를 합칠 수 없다. 새 inverse는 다른 초기 추정으로 같은 방정식을 풀었고, 잔차·L<1·반경·양의 barycentric margin으로 해를 인증했다. 반경이 약4.8e−35인 실제 내부 preimage를 그대로 사용한다. 16회 반복이 끝났다는 사실이 아니라 이 최종 인증이 통과 근거다. 기존 barycenter 초기값 실패는 비존재 증거가 아니다.
2. Normal-gap 선택은 후보 생성법이다. longdouble 각도 계산의 최적성이나 완전성이 없어도 선택된 실제 방향에 대한 Q128 derivative interval 및 충분조건의 독립 정확 검사가 통과하면 유효하다. 방향이 있다는 사실만으로 전체 domain이 단일 strip이라는 결론은 나오지 않는다.
3. 원점 거리 정리는 유효하나 작은 양의 하한 자체가 효율적인 atlas를 보장하지 않는다. 내부7개와 외부8개는 정확한 domain 분류로 유지한다. tiny offset을 반올림해서 경계 접촉으로 바꾸지 않는다.

### 원점 하한의 정확한 의미

target half-plane을 h(y)=n·y−b≤0으로 두고 h(chi(0))=a>0이면, 같은 convex reference source 안의 X∈D에 대해

    a ≤ |n·(chi(0)−chi(X))|
      ≤ ||n||1 (1+L) ||X||∞.

따라서 ||X||∞≥a/[||n||1(1+L)]이며 Euclidean r에도 같은 하한이 유효하다. a는 엄격한 outward 하한, 분모는 상한이어야 한다. 이는 source가 reference0를 포함하고 segment가 source 안에 있는 경우의 증명이다. 저장8개가 이 조건을 만족한다는 결과를 채택한다. source 전체 L의 실제 nextafter inflation 재현도 유지한다.

## 선택: 33개에 대한 실제 event arrangement와 fixed atlas를 완성하는 한 묶음

새 방향/seed를 다시 탐색하지 않는다. 이미 저장된 증거를 소비하는 native geometry 경로를 구현한다. 정상 panel과 완료된 비용2의 256 typed panel은 재사용한다. 기존 all-target 원본 registry, 실제 loop, canonical event/endpoint, 입력 state/heat/tau는 불변이다. 새 물리 quadrature와 전체 analytic remainder 실행은 이번 범위에 포함하지 않는다.

### 1. 공통 사건과 inverse binding

8개 원래 실패 호출은 동일한 source133/power-node 입력 enclosure와 함수 hash에 묶어 인증된 inverse packet을 한 번 참조한다. 각 target의 incident plane 및 방향은 별도 보존한다. 작은 양의 barycentric 좌표를0으로 만들지 않는다. source vertex118의 ID를 재사용하지 않는다.

이 preimage로 실제 chart event를 구성하고, 원래 node/boundary/loop의 연결과 half-open ownership을 대조한다. chart 좌표 변환의 enclosure가 새 prior 안에 포함되는지도 요구한다. inverse 존재 인증은 뒤따르는 implicit graph의 복소 영역 인증을 자동으로 보장하지 않는다.

### 2. 14개 oblique arrangement

인증된 d를 graph의 종속 방향으로, 독립인 e를 자유 방향으로 사용하여 X=Xc+u e+v d로 표현한다. det[e,d]의 비영 인증과 실제 Jacobian을 저장한다. d를 정확히 단위벡터라고 가정하거나 정규화 오차를 버리지 않는다.

실제 target-plane pullback의 v 미분 부호는 이미 증명된 방향 자료를 사용한다. 기존 canonical node와 source-boundary 사건을 u에 투영하여 slab을 구성한다. u-extremum이 source 경계에 존재하면 그 투영의 사건만 같은 원래 고정 예산으로 격리한다. 원래 power-plane root 전체를 다시 구하지 않는다. 동일 u 사건의 tie, 중복 endpoint 및 다중 strip을 명시적으로 다룬다.

각 slab에서 경계 branch ID, 상하 순서, 대상 source/target half-plane 부호를 인증한다. 단조 방향이 모든 target 면에서 성립해도 source 경계의 fold나 비볼록 pullback 때문에 별도 strip이 필요할 수 있다. 배열 순서나 표본 한 점으로 그 증명을 대신하지 않는다.

Duffy apex를 포함하는 부분은 기존 blow-up 좌표와 원래 analytic 대표를 유지한다. 새로운 Cartesian complex inverse를 apex를 가로질러 적용하지 않는다. 좌표 전환은 동일 함수의 평가이며 coefficient refit이 아니다.

### 3. 원점 관련 15개: density 규모와 tiny geometry 규모를 분리

원래 rho=floor+2/(r+epsrho)×두 logistic, epsrho=1e−6를 그대로 사용한다. rho(0)가 유한해도 Cartesian |X|는 해석적이지 않다. 원점을 기준으로 한 source fan X=s R(t), s≥0에서 r=s sqrt(R·R)로 평가하는 경로를 구현한다. sqrt의 가지, r+epsrho pole, logistic pole 및 기존 source enrichment의 각도 분모를 각각 인증한다.

내부7개는 s=0을 포함하는 실제 영역을 덮는다. 외부8개는 원점과의 실제 작은 간격을 보존하고 각 ray의 실제 하한 s_in(t)>0를 target-plane event로 결정한다. 원점을 채워 넣거나 작은 제외 영역을 삭제하지 않는다. 기존 정확 거리 하한은 membership/separation 증거로 유지하되, 그 1e−23 규모를 모든 smooth density panel의 폭으로 강제하지 않는다.

밀도 pole clearance가 지배하는 radial panel은 기존 고정 복소 여유에서 유도한 local s·|R|+epsrho 규모를 사용한다. Tiny clipping boundary의 복소 branch가 더 제한적이면 그것은 별도의 event-chart 제약으로 처리한다. 물리 density의 epsrho 규모와 기하 branch의 작은 간격을 혼동해 둘 중 하나를 무시하지 않는다. 실제 branch가 기존 예산에서 인증되지 않으면 pending으로 종료한다.

원점 fan 경계, power-vertex preimage, 원래 source-edge 사건을 연결하는 과정에서 새 projection 사건이 필요하면 해당 실패 domain 안에서만 기록한다. 기존 root/node ID와 구간은 변경하지 않고 새 좌표계 사건을 참조 관계로 추가한다. full source 또는 passing chart의 재생은 하지 않는다.

### 4. 고정 계획과 비용

이번 산출물은 33개 각각의 실제 typed atlas, 미완료이면 그 최초 증명 실패 witness다. 같은 fixed Gauss 차수, cap, query 예산 및 complex-clearance 정책을 사용한다. 새 방향이나 panel 비율을 여러 후보로 반복하지 않는다. 각 source-target에서 ordinary/clipped 부분을 합쳐 원래 예산을 검사하고, 전체 target 총합도 대조한다. Slot 하나의 통과를 target 예산 통과로 대신하지 않는다.

패널에는 domain/branch/event provenance, 변환 determinant, 원점/apex 유형, 실제 radius와 pole clearance, 필요한 analytic proof 유형, arithmetic backend, 예상 query/bytes를 남긴다. 지금 계산하지 않은 analytic remainder를 true로 두지 않는다. 빈 영역은 실제 arrangement 증명으로만0비용 처리한다.

## 종료 및 다음 경계

한 묶음에서 기존 33개의 동일 domain을 빠짐없이 다루고, 사건 순서·coverage·half-open ownership·고정 비용을 확인한다. 기존 passing packet/typed panel, caller와 m/U/x/v, generating heat 및 원래 모든 tau의 hash 불변을 요구한다. 정상 함수817152점의 재검사는 필요 없다.

한 개라도 남으면 해당 이유와 부분 atlas를 저장하고 종료한다. cap/정밀도/차수/허용오차 증가나 자동 fallback 연쇄로 이어가지 않는다. 모두 닫히면 다음 substantive 단계는 전 target analytic certificate와 source-outgoing/target-incoming 동일 packet의 보존적 수치조립이다. 실제 비용을 고정한 뒤 그 단계의 실행 범위를 정한다. 현재는 전체 quadrature를 열지 않는다.

이는 native transactional transfer를 완성하기 위한 기하 연결 단계다. 일반 evolved heat, weight coupling, native CUDA/OpenMP/MPI 시간 갱신과 실제 t6/8/10 및 논문 최장공전 목표는 여전히 미검증이다. 이번 기하 성공으로 조건부 W0 개별 고정질량이나 장시간 물리 안정성을 대신 증명하지 않는다.
