# 국소 전달 연산 검증 결과 — 2026-10-22

## 완료한 범위

같은 CONTROL 상태의 **384개 소스 삼각형·24개 대상 셀**을 검증했다.
128²/256² 원반의 시간진화 결과는 아니다.

| 항목 | 실제 결과 |
|---|---|
| 전체 적분점 | 기존 결과를 포함한 23,339,008점 완료 |
| 이번 새 적분 | 21,629,184점, 465/465개 겹침 구간 |
| 기준 함수 대조 | 같은 계산점의 2,241개 field17·물리량12 비교 통과 |
| 독립 정수 감사 | 모든 새 점의 합·오차와 84,489개 체크섬·합산 상태 통과 |
| 대상 셀 | 24/24개 기존 총 오차 한도·질량/내부에너지 양수성 통과 |
| 보내는 쪽/받는 쪽 회계 | 559개 같은 겹침 레코드 사용, 물리량12의 정확한 합계 차이 0 |
| 기존 결과 | 기존 5개 대상 셀 재계산·교체 없음 |

새 19개 셀에서 총 오차 상한/허용치의 최댓값은 **0.28694**였다.
이는 물리해의 상대오차가 28.7%라는 뜻이 아니다. 보존된 기존 68526번
셀까지 포함한 24개 전체의 최댓값은 **0.71296**으로, 역시 허용치 안이다.
오차 장부에는 적분 나머지, 경계·좌표 입력 불확실성, 네이티브 계산,
실제 합산 및 최종 저장 형식 변환을 포함했다. 네이티브 함수/libm 오차
모델의 기존 가정은 유지한다.

실제 적분은 lageunha의 C++/OpenMP16·MPI1에서 약23분 걸렸다.
저장된 점의 독립 정수 감사는 약14.2초였다. 이 binary80 참조 계산과
약23.89GB 진단 자료를 production GPU 시간 갱신 방식으로 쓰지는 않는다.

## 발견하고 수정한 입력 문제

438,120개 실제 적분점 좌표 표본 중 기존 Cartesian 역변환은 4,272개에서
모서리 분모의 양수성을 보장하지 못했다. 원래 셀 좌표의 중심과 오차
구간을 함께 직접 전달하는 경로로 수정했다. 같은 표본의 미분리 건수는
0이 되었고, 이후 위의 실제 적분·오차 감사를 완료했다. 이것을 t>6
가짜 가열의 물리적 원인이라고 판정한 것은 아니다.

## 아직 증명하지 않은 것

보내는/받는 레코드 합계의 일치는 전달 회계의 검증이다. 실제 유체의
시간진화에서 힘이나 가짜 가열 오차가 0이라는 뜻은 아니다. 소스 분할
적분은 기존 coverage 증명과 오차 장부로 원래 연속 적분을 감싸지만,
별도의 거대한 whole-source 재적분을 수행한 결과는 아니다.

현재 일반 evolved heat·full-field fragment caller state·accepted RK4
step·장시간 안정성은 아직 완료되지 않았다. 다음 물리 진화 연결은
[Astra의 실행 계획](ASTRA_SHARED_STATE_EVOLUTION_PLAN_20261022.md)에 따라
source field/열역학 context를 보존하고 기존 CUDA 힘 연산자·RK4·실제
restart로 진행한다. 12개 모멘트로 고차 속도/회전을 임의 복원하거나
현재 J로 열을 다시 초기화하지 않는다. 작은 연결 검증 이후 실제
t=6·8·10과 GIZMO 최장 공전 비교를 계속해야 한다.

## 근거

- `shared_mixed_native_v2_mpi1/contract.json`
- `shared_mixed_native_points_audit_v2/contract.json`
- `shared_mixed_native_summary_audit_v1/contract.json`
- `shared_original_whole4_source_accounting_v1/contract.json`
- `shared_transfer_incidence_v1/contract.json`
- `SHARED_MISSING_SOURCE_LEDGER_20261021.json`
- `shared_runtime_final_integrity_v1/integrity.txt`
