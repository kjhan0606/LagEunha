# Astra: 동일 종료 평가의 target 예산 개정

2026-10-21. 사용자의 명시적 최신 지시 “예산을 충분히 늘려요”를 반영한다. 원본 ASTRA_SHARED_FAILED_ATLAS_END_REVIEW_20261021.md 및 그 SHA 참조는 변경하지 않는다. 이 부록은 해당 평가의 **target query/panel 비용 제한만** 개정하며 나머지 수학·물리·출처·종료 조건을 유지한다. SHARED_ATLAS_COST_FLOOR_20261021.json과 대응 native contract의 저장 회계를 읽었다. 새 계산·폴링은 하지 않았다.

## 적용 정책

운전자가 제시한 target당 16배 예산을 채택한다.

- 기존: 817152 queries = 3192 panels.
- 새 정책: **13074432 queries = 51072 panels**, 각 panel의 16×16 규칙 유지.
- 새로운 immutable budget-policy ID/SHA와 적용 run/manifest를 기록한다. 기존 실패와 그때의 cap을 소급 변경하지 않는다.
- root isolation/iteration guard, analytic clearance, 차수, 정밀도, 물리 tau/normalizer, density profile, 상태와 heat는 그대로다. target query 예산 확대가 다른 수치 guard의 자동 확대를 허용하지 않는다.

원본 계획에서 “기존 고정 query/panel 예산”이라고 한 부분은 새 작업에 한하여 이 정책을 뜻한다. 이미 완료한 원래 예산 검사는 그대로 유효한 기록이다. 별도 사용자 승인 재요청은 필요 없다.

## 두 target의 완료 계획은 재사용한다

저장된 완료 부분만 합해도 target67943은 4621 panels/1182976 queries, target68527은 5043 panels/1291008 queries다. 원래 cap 초과는 실제 회계 사실이다. 두 값은 새 cap 미만이지만 **미완료 slot을 제외한 하한**이므로 전체 target 비용 통과로 기록하지 않는다.

source131/slot3152 및 source271/slot6520의 각3078 panels와 나머지 완료 typed plan을 그대로 재사용한다. 두 source의 chart0/1/4/5 whole, chart2/3 clipped 구성은 모든6 chart가 whole인 앞선 비용수정 사례와 다르다. 예산 확대 후 이 정상 계획들을 줄이기 위한 추가 replan/geometry replay는 필요 없다.

기존 contract의 `two_complete_target_plan_families_need_replanning=true`는 원래 cap 아래의 판단으로 보존한다. 새 결합 manifest에서는 “새 예산에서는 이 이유만으로 replan 불필요, 전체 비용은 아직 UNKNOWN”이라고 명시한다. 기존 문자열을 바꾸거나 원래 실행이 새 정책으로 통과했다고 쓰지 않는다.

## 다음 실행 경계

원래33 pending의 실제 arrangement/atlas 구현 계획을 계속한다. 예산 증가는 oblique arrangement, source133 inverse binding, 원점 및 pole clearance, analytic certificate를 해결한 증거가 아니다. 완료된 typed plan을 재사용하고 새33 계획을 더한 후 target별 전체 비용 및 전 target 고유 계산량을 확정한다.

실제 적분은 actual atlas, 필요한 analytic proof 및 전체 비용 계획이 확정된 후 진행한다. 큰 point packet 전체를 메모리에 올리는 대신 native batch streaming을 사용한다. batch별 domain/panel/query ID, 출처 hash, 처리 개수, 중복·누락 및 재시작 상태를 기록하고, 기존 결정적 합산 순서/오차 장부와 최종 R80 변환 검사를 보존한다. batch 크기 변경을 물리 적분 규칙이나 오차 허용치 변경으로 연결하지 않는다. 총 저장량·최대 상주 메모리·예상 walltime도 비용 manifest에 포함한다.

새 cap에서도 초과하거나 기하/analytic 증명이 실패하면 해당 이유를 보존한다. 자동으로 추가 cap을 늘리거나 불완전 packet을 publish하지 않는다. 이번 개정은 더 충분한 작업 예산으로 기존 과학적 관문을 수행하기 위한 것이며, native CUDA/MPI/OpenMP 장시간 Kepler 또는 t6/8/10 및 논문 최장공전 성공 판정은 아니다.
