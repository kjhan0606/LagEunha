# 중심부 세분화 이후 Astra 묶음 평가

2026-10-20. 기존 Astra 계획/평가 에이전트의 두 번째 묶음 검토를 한국어로 정리했다. 운전자가 코드를 작성하며 매 단계 Astra를 호출하지 않는다.

## 판단

고정 eps 세분화는 초기 양의 열률을 약973배 줄여 공간 부족의 영향을 보여 줬다. 하지만 fine의 과도 냉각/단계 음의 U는 남아 있다. 새 vector limiter는 실제 cold/shock에서 악화된 후보로 보류한다. 실제 구현은 Udot=die-m v·a를 RK4로 적분하므로 보수적 총에너지 FE의 양성 이차식을 현재 stage의 실제 양성 조건으로 사용하지 않는다.

## 제한된 판별 묶음

1. 정상 fine step-start에서 같은 실패 구간의 full1/half2/quarter4를 최대7시도 단계로 재생한다. 모든 stage를 EOS floor 전에 검사한다. quarter 통과도 시간간격 민감성만 뜻한다. 실패 stage를 재시작 입력으로 쓰지 않는다.
2. 실패 셀·직접 이웃·같은 반경 대조 셀의 Udot=Edot-v·F 및 R_U=Udot+P Vdot_geom를 signed 항으로 기록한다. 실제 Pstar, 양쪽 압력, pair pressure, 실제/기하 면속도, force/work, 체적과 근접쌍을 포함한다. PdV≈0일 때 비율로 원인을 판정하지 않는다.
3. 실제 체적 Jacobian과 압력력의 adjoint 관계를 우선 검증한다. F_i=sum_j P_j dV_j/dx_i, Udot_i=-P_i Vdot_i의 가역 기준은 전체 압력 에너지 교환을 닫는다. 현재 공유면 기하에 독립 체적 미분을 대조한 후 현재 force/work와 차이를 확인한다. 옛 VPH flag를 켜는 방식은 아니다.
4. 이후에만 호환되는 유한 시간 에너지/양성 갱신 후보를 구현한다. 보수적 E 적분으로 바꾸는 것만으로 잘못된 반이산 냉각을 고칠 수는 없다. 동적 weight/material volume와 고차 moment/spin은 장부가 필요성을 보여 준 뒤에 진행한다.

고정 엔트로피의 U=K m^gamma V^(1-gamma)/(gamma-1)는 K,V>0에서 양수인 가역 기준이다. 일반 Riemann 교환에서 이 값을 강제로 보충하는 해법은 아니다. Sedov에는 일치하는 수치 충격파 엔트로피 생성이 필요하며 Re∞가 이를 금지하지 않는다. Kepler의 임의 물리 점성·열 보충·spin 감쇠는 넣지 않는다.

## 근거의 한계와 정정

검토 시 agent가 직접 읽은 core 결과/보고서와 달리 Sedov JSON 두 경로는 잘못 전달돼 독립 열람하지 못했다. 실제 위치는 이 문서와 같은 `rk4_diagnostics/` 아래 `failure_replay/results.json`, `failure_geometry_v4/results.json`이다. 운전자의 이후 직접 열람 및 native 보존 자료로 추적할 수 있다.

앞서 실패 셀 두 개에 초점을 맞춘 보고는 전체 실패 개수를 뜻하지 않는다. 원래 fine old/vector 로그 모두 네 셀67935,68001,68526,68592를 기록한다. full 재생의 U/V/m은 정확 일치하므로 ID가 달라진 재현 문제가 아니다.

짧은 구간 통과는 t6–10이나 GIZMO 약1000공전 안정성의 대체 증거가 아니다. 최종 목표는 `../GIZMO_COMPLETION_CONTRACT.md` 그대로다.
