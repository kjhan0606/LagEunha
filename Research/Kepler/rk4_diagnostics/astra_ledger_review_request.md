# Astra 묶음 판정 요청: 실제 cold 실패 장부 + 과거 가역 불안정

한국어로 답해 주십시오. 코딩은 운전자가 합니다. 이번 완성된 원인 판별 묶음에서 다음 한 개의 수정/판별 후보를 선택해 주십시오. 사용자 최종 목표는 native CUDA+MPI의 GIZMO 최장 약1000공전 안정성이며 짧은 양성 통과가 아닙니다. 조향 금지, Re∞, W0→Voronoi/개별 고정질량, 임의 열 보충·감쇠 금지. 시간 의존 weight/내부 spin은 이 계약을 닫는 경우 허용됩니다.

정확한 근거:

- `rk4_diagnostics/fine_cold_interval/ledger/results.json`, `selected_faces.tsv`
- `rk4_diagnostics/RESULTS.md`의 새 fine ledger 절
- `rk4_diagnostics/fine_cold_interval/manifest.json`
- `/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/HANDOVER_2026-10-09.md` 전체. 이력의 TCF/ring과 현재 Plummer+기존HLLC native IC는 동일 모델/정의가 아니지만 단순 adjoint/VPH 복귀가 근본 가역 미끄러짐 불안정을 다시 겪을 위험을 평가해야 합니다.
- Sedov 실제 JSON: `rk4_diagnostics/failure_replay/results.json`, `rk4_diagnostics/failure_geometry_v4/results.json` (이전 전달에서 rk4_diagnostics prefix가 누락됐습니다).

새 증거: full/half 실패, quarter4 동일 endpoint 통과. full의 U/V/m은 원본 실패와 비트 일치; 네 ID 모두 일치합니다. quarter 통과 후에도 U/Uad≈.0406으로 소진 상태가 지속돼 장시간은 재개하지 않았습니다. native 첫 RHS를 실제 면장부로 재구성했고 force오차1e-22, Udot오차7e-15입니다. 실패 셀의 Udot≈-2.369e-8 대비 -Pi Vdot≈-2.93e-11. 추가 냉각은 모두 (Pi-Pstar)s_geom이며 pair=0, work-geom 차이≈1e-21입니다.

실제 native 면에서 fixed-W volume Jacobian adjoint force를 **읽기 전용 기준으로만** 계산했습니다. 총 power잔차1.84e-20/sum|PdV|2.89e-6. 독립 Qhull x±delta v의 Vdot오차는 delta반감 때4배 감소해 최종L2relative8.73e-9였습니다. 이 기준으로 진화하지 않았습니다.

CUDA `Exam/exam_gpu.cu:625` 부근은 cellP ratio100에서 pressure뿐 아니라 work face speed도 HLL vstar로 갑자기 바꿉니다. quarter stage의 Udot은 -2.3e-8/+4.7e-9로 오갑니다. 임계 분기가 chatter를 만드는지 아직 직접 계측하지 않았으므로 단정하지 않습니다. Pi/Pj와 그 분기 flag를 목표 셀/몇 번의force call에만 기록하는 작은 native 계측이 한 후보입니다.

물어볼 사항:

1. 이 signed ledger로 현재 fixed-mass pressure/work의 어떤 수학적 결함까지 확정할 수 있습니까? 단순 adiabatic값 비교만으로 일반 Riemann 교환을 오판하지 않게 범위를 정해 주십시오.
2. 다음 실행은 branch 기록/경계 U→0 operator 판별을 먼저 해야 합니까, 아니면 한 개의 에너지/엔트로피 호환 force-work 후보를 구현해야 합니까? 한 개만 우선 선택하고 종료 기준을 주세요.
3. adjoint 가역 폐쇄가 과거의 sliding-cell 선형 불안정까지 해결하지 못함을 감안하여, 그 후보가 최종 장시간 목표로 가는 데 어떤 역할인지 명확히 해 주십시오. 막연한 limiter/해상도/CFL scan, 과거실패 DAE 재활성, spin감쇠는 피하십시오.

운전자는 최신 native 소스와 frozenv5를 유지합니다. 모든 기존 결과를 보존하고 읽기 전용 분석 외 새 물리 수정은 아직 하지 않았습니다. H200분석은QOSGrpGRES 대기로 취소했고, 사용자 허용 CPU서버lageunha에서 Python 읽기 전용 장부만 끝냈습니다.
