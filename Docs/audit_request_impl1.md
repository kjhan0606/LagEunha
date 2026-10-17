# 감사 요청 (Astra, 구현 1차): 설계 1–5단계 구현과 관문 결과

작성: 2026-10-02. 사용자 지시로 2026-10-07 08:00까지 감사는 Astra가 맡는다.
설계: `Docs/design_astra.md`. Fable 감사: `Docs/design_astra_audit_fable.md`. 운전자 기록: `Docs/design_astra_audit.md`(관문 결과).

## 구현 (`Exam/exam.c`, 모두 기본 꺼짐)
- `GFS_SHARED_GEOM`(df3c796): `shf_record`, `shf_common`. 단계 기하 패스에서 쌍 면 끝점을 기록하고, TCF 면 루프에서 공통 구간을 쓴다. 나머지는 진공 면이다.
- `GFS_PERIODIC_COMPLETE`, `GFS_TCF_FACE_Z`(0ae13b6): `cen_poly`, `pc_ncap`. 면 상태 임피던스 Z = max(μ(c_M + max(m,0)), μ[a + √(a² + c_M²)]), ρ_if = ρ_i(1 − 2δ_i·e/d), P_if = P_i(ρ_if/ρ_i)^γ.
- `GFS_HG_AFFINE`(10d4734): `hg_cell`, `hg_vjac`, `hg_ensure`, 래퍼의 가속도·에너지 추가. 기준 구현은 `Docs/subzonal_checks/hg2d_ref.py`(+`_test`, `_test2`).
- Q 법칙은 minmod를 기본으로 두었다(QLAW=2는 Sod L1 초과).

## 결과 (`Docs/design_astra_audit.md`)
- 관문 H(공유 면): 자유 표면 판 에너지 −3.8e-3(dt와 무관)에서 RK4 4차 수렴(−2.6e-13)으로 바뀜.
- 관문 C(2D Sod, m∝ρ): 면 Z + minmod는 통과. HG를 더하면 L1이 더 좋아짐(0.0151/0.0391/0.0239).
- HG 장부·dt 수렴은 위 문서 참조. Gresho는 진행 중.

## 요청
1. 구현이 설계·감사 정정(면 Z 우선, 법선 회전 항, 기하 dt, minmod 기본)과 일치하는지 코드 수준에서 점검하라. 특히 `hg_cell`의 기울기·필터·계수, `shf_common`의 공통 구간, 면 상태 Z, 래퍼의 에너지 장부(`die`가 dE/dt인 경로).
2. 관문 판정이 타당한가? 빠진 관문 항목(예: Noh 유입 경계, Sedov, Hill, Kepler)의 순서와 기준을 지금 결과에 맞게 확정하라.
3. 오류나 위험이 있으면 목록으로, 없으면 "감사 통과"와 다음 단계 지시를 달라.
읽기 전용, 한국어, 결론 먼저.
