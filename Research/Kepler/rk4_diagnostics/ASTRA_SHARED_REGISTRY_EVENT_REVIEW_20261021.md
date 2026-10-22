# Astra: retained domain과 새 canonical event 연결

2026-10-21. 계획·평가 Astra. SHARED_CONTEXT, roots audit, domains 계약 및 `test_shared_registry_domains_v1.cpp::node_events`를 읽었다. 원본 HEAD8e25dfe 확인. 계산·root 재평가·작업 폴링은 하지 않았다.

## 판정

제안한 **저장 증거만 사용하는 5개 domain bridge를 채택한다.** 현재 실패 조건은 old boundary의 t.n/t.e와 새 root record의 t.n/t.e가 정확히 같아야 한다는 metadata 요구다. 같은 물리적 root에 대한 서로 다른 유효 enclosure는 비트 일치할 필요가 없다. 새 근 계산이나 원래 geometry 교체가 필요한 실패로 해석하지 않는다.

단, '구간이 겹친다'는 조건만 단독으로는 root 동일성 증명이 아니다. 아래 존재·완전성·동일 함수 bindings를 함께 요구한다. 이번 범위는 기존 geometry를 그대로 두는 event identity bridge이며 새 domain 적분/전체 geometry 재생이 아니다.

## 동일 event의 충분조건과 증명

같은 고정 함수 f(t), 같은 canonical parameter 범위 [0,1]에 대해:

1. old interval I가 실제 root a를 포함한다는 저장 존재증명이 있다.
2. 새 complete set J_0,...,J_n이 f의 해당 범위 모든 root를 포함하며 서로 분리된 isolated-root records다.
3. I와 교차하는 J_k가 정확히 하나다.

그러면 a는 complete set의 어떤 J에 속하며 동시에 I에 속해야 하므로 그 J는 J_k다. 새 record의 유일성까지 있으면 동일 root/event ordinal임이 결론난다. I가 J_k에 완전히 포함될 필요는 없다. I와 J_k의 교집합은 동일 root를 포함하지만 **원래 boundary 구간을 그 교집합으로 교체하지 않는다.** bridge가 보존하는 것은 event identity다.

여러 새 구간과 교차하거나 하나도 교차하지 않으면 bridge 실패다. epsilon, nearest ordinal, 가장 가까운 midpoint로 해결하지 않는다. endpoint root는 exact endpoint 증명과 별도 vertex identity가 있어야 한다. 지금 root 감사의 endpoint_roots=0을 모든 미래 사례의 가정으로 고정하지 않는다.

## 필수 saved-proof bindings

각 bridge 행에 다음을 묶는다.

- epoch1/state/heat SHA, 실제 target registry 및 weight SHA.
- canonical edge의 ID, 실제 endpoint 순서, shared analytic function/basis/계수 SHA. 방향이 반대면 정확히 t↦1-t로 interval을 변환하고 orientation을 기록한다.
- unordered target pair와 정확 power plane. own-target가 달라 평면 부호가 반대인 경우, 정확한 전체 부호 반전임을 확인한다. 일반적인 '비슷한 normal'은 동일 함수 근거가 아니다.
- old 존재증명 생성 run/source/manifest, old root index 또는 그 boundary에 연결된 사건 provenance, 원래 serialized interval. old geometry CSV에 숫자만 있다는 이유로 존재증명을 새로 가정하지 않는다. 당시 complete-root 결과와 그 구간의 연결이 있어야 한다.
- new complete-root packet/manifest, new isolated-root ordinal, 분리·완전성의 증명 상태, serialized interval.
- exact-decoded interval intersection 및 유일 교차 개수, old→new eventID mapping.

interval endpoint는 binary payload의 실제 중심/radius로부터 정확하게 복원한다. CSV가 원래 값을 충분히 roundtrip하지 못하면 그 파일만으로 연결을 승인하지 말고 원래 packet/outward interval을 사용한다. 반올림한 decimal midpoint끼리 비교하지 않는다.

고정 exact 함수가 아니라 입력 불확실성 family에 대한 root enclosure이면 old/new가 **같은 family와 parameter에 대해 존재·완전성을 갖는지**도 확인한다. 서로 다른 family의 구간 겹침은 동일 사건 증명이 아니다. 현재 exact coefficients를 중심으로 한 계산 오차 enclosure와 물리입력 family를 구분한다.

## 한정 실행 범위와 구조 검사

bridge 대상은 target68526 source246/253/255/285/287의 실제 사용된 crossing node뿐이다. 원래 v5 domain node/boundary/loop CSV는 byte-preserve하고, 별도 immutable bridge 파일에 새 event ID만 기록한다. v1 실패자료는 보존한다. 다른9211 완료 record는 그대로 참조한다. 전체 원래 target 출판 packet은 변경하지 않는다.

같은 node가 두 boundary에서 사용되면 두 참조가 같은 eventID에 결합돼야 한다. 서로 다른 old root가 하나의 new ordinal에 결합되는 경우 단순 허용하지 말고 동일한 기존 사건의 복수 참조인지 확인한다. 원래 분리된 두 사건이면 실패다. multi-plane vertex는 실제 active-plane set와 rank2 교점 증거를 유지한다.

source-edge chain은 canonical t 방향의 event 순서와 segment 소유권을 검사한다. 같은 event의 여러 enclosure는 **공통 물리 root를 포함한다는 위 증명**으로 연결하며, 단순한 pairwise-overlap 사슬을 equivalence relation으로 사용하지 않는다. 정렬은 인증된 new isolated-root 순서에 따른다. 기존 segment interval은 그대로 보존하고, event ID만으로 양의 길이·누락 없는 coverage까지 자동 승인하지 않는다.

power-interface의 두 target-side fragment는 unordered plane/event endpoints가 같고 방향이 반대여야 한다. 외부 support 경계는 짝 없는 것이 정상이며 실제 external 역할로 기록한다. zero-measure tie와 half-open owner 규칙도 기존대로 유지한다. global source-image disjointness 및 source partition cancellation은 별도 미완료 관문이다.

## 완료·중단 기준

5개 bridge가 위 조건을 모두 통과하면 원래9211개 완료 record와 함께 **9216 domain arrangements complete**로 조립할 수 있다. 이는 source partition/global image disjointness/all-target numerical transfer 통과가 아니다. 해당 flags는 실제 후속 검사 전 false/pending으로 둔다.

하나라도 provenance·동일함수·unique intersection·event order를 증명하지 못하면 해당 bridge만 미완료로 남긴다. root/geometry/540개 새 candidate를 자동 재계산하지 않는다. 필요 증거가 없는 이유를 정확히 기록하고 이 저장-record 묶음을 종료한다.

새 quadrature/열·상태 변경/시간진화는 없다. 전체 quadrature 차단은 유지한다. 목표는 원래 native CUDA/MPI/OpenMP의 실제 장시간 Kepler이며 metadata 연결 통과를 물리 진화 통과로 확대하지 않는다.
