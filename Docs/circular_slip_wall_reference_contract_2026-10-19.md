# 원형 slip wall 참조 구현 계약

계획 근거: `astra_kepler_plan_review_2026-10-18.md`의 2026-10-19 Ψ-coupled 묶음. 이 문서는 구현 내용과 현재 검증 범위를 기록한다.

## 기하와 경계

원점 중심 R=2.5의 고정 disk와 power halfplane의 정확한 교집합을 표현한다. 직선 공유면에는 이웃 ID를 저장하고, 벽에는 각도 구간과 반경을 저장한다. 원호를 chord로 치환하지 않는다. 원호의 quadrature는 x=R(cosθ,sinθ), n=x/R, ds=R dθ, wn=0이다.

벽에서는 normal velocity만 mirror하고 tangential velocity는 보존한다. 내부 analytic reset, 생성자 clamp/projection, 점성 또는 조향을 추가하지 않는다. 이번 static fixture는 채택된 native exact Euler star와 기존 kinetic transport 함수를 실제 호출한다. 질량/E/Phi flux와 원점 torque가 반올림 오차 범위에서 0인지 각 quadrature node에서 검사한다.

## 고정질량과 가중치

모든 내부 유체와 r=2..2.5 buffer가 동적이어야 한다. 닫힌 계에서 constant W는 기하를 바꾸지 않으므로 massroot의 nullspace를 명시적으로 제거해야 한다. 작은 정적 reference는 augmented system과 sum(z)=0을 사용한다. 전체 residual이 부적합하면 명시적으로 실패한다. 기존 frozen Dirichlet root를 전체계에 그대로 쓰지 않는다.

이 augmented dense solve는 향후 production root의 구현을 의미하지 않는다. 아직 nonlinear closed-system massroot 또는 동적 원형 경계 RHS는 실행하지 않았다.

## 정적 검증

single/pair/weighted 5-cell에서 partition 면적 πR², q=0, Q=(πR⁴/4)I, 각 cell의 ∮n ds=0 및 ∮x⊗n ds=VI를 검사했다. 공유면의 뒤집힌 endpoint/normal/measure와 flux의 쌍별 반대 부호 배분을 보장하는 기하를 확인했다. constant W shift invariant, empty cell 및 incompatible residual 거부도 검사했다.

결과: `canonical_circular_power_reference_2026-10-19.json`. 모든 fixture exit0 통과. 이것은 작은 all-constraint reference이며 약2만 셀 production geometry의 비용·robustness 또는 장시간 안정성 검증은 아니다.

## 후속 초기화 시 지켜야 할 사항

정상 Kepler장을 실제 disk-clipped cell에 적분해 m/raw/Psi를 초기화한다. 원래 square/frozen 경계의 mass를 새 경계셀에 재사용하지 않는다. 바뀐 명시적 경계 benchmark임을 기록하고 이후 buffer/반경 민감도를 확인한다. 빈 셀과 해상되지 않은 기하는 조용히 잘라내지 않고 실패시킨다.

공유면 유량의 추가 정적 검증은 `canonical_circular_shared_flux_2026-10-19.json`에 저장했다. weighted 5-cell의 실제 128개 공유면 node에서 native face-record, 기존 kinetic transport 및 vector moment accumulation을 호출했다. raw18의 전체 합과 면적/Phi rate 합은 각각의 반올림 bound 내에서 상쇄됐고 exterior ledger는 정확히0이다. 테스트 Phi=-1/r는 쌍별 상쇄를 확인하는 공통 시험함수이며 softened Kepler potential을 사용한 동적 경계 검증은 아니다.
