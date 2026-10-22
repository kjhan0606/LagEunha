# Astra: 공통 경계 대표로 적분 오차와 경계 불확실성 분리

2026-10-21. 기존 astra_planning_review의 한 번의 묶음 평가. 코딩/실행은
운전자가 한다. 이전 계획의 사건 격리와 양의 보완영역 조립은 유지한다.

## 판단과 실행 범위

다음 구현은 같은 매끄러운 경계 대표를 두 적분 규칙에서 공통 사용하는
branch-tracked 적분이다. triangle225와137만 먼저 수행하고, 같은 작은
닫힌 라게르 셀 fixture를 한 번 회귀 확인한다. 이 제한 묶음이 끝나기 전
전체384개, weight Newton, 압축, 시간 진화를 시작하지 않는다.

복구 자료는 유효한 실패 분석 자료이다. 첫 native 종료코드는 unknown이다.
384개 중344개 적분 추정 기준 통과/40개 실패와 별개로 실행 전체 성공을
기록하지 않는다. 원래 cap, 두 양의 규칙, 허용오차는 유지한다.

## 기존 추정량의 문제 가능성

현재 비교는 Q8,12(target-unknown)-Q16,16(target-unknown)이다. 그러나
radius마다 다시 격리한 root window의 끝점이 known 영역 경계로 사용된다.
두 규칙은 물리적 적분함수와 수치적인 window 절단 변화를 함께 비교한다.
이는 실패 가능한 기전이지 ratio42의 원인 확정은 아니다. 실제 원점 밀도
변화나 endpoint 특이성에 따른 적분 부족도 구분해야 한다.

## 구현

사건 사이 각 simple branch에 고정 ID/active plane을 부여한다.

    f(r,theta_star(r))=0
    d theta_star/dr = -f_r/f_theta

기존 compiled map의 고정 계수가 같은 해석적 대표를 정의한다. Continuation
predictor와 safeguarded root solve로 연결하고, interval 검사는 근의 존재,
유일성, 오차를 확인한다. Source 경계에도 같은 원칙을 적용한다.
Interval window의 양 끝점을 물리 적분 경계로 쓰지 않는다.
f_theta가0을 포함하면 simple-branch 공식을 강행하지 않는다. 접선/교차/
chart 변경은 기존 event band에서 처리한다.

두 규칙은 같은 theta_star_L(r),theta_star_R(r) 사이의 양의 영역을 적분한다.

    |I-Ihat| <= E_quad + E_root/representation + E_event + E_gluing

대표를 선택해도 계수/평가/root 불확실성이 사라지지 않는다. Unknown을 뺀
값을 임의로 매끄럽게 만드는 방식이나 허용오차 변경으로 대체하지 않는다.

## 두 삼각형의 판별 witness

각 실패 panel/target/moment에 대해 다음을 같이 기록한다.

1. 기존 known-only coarse/fine 차이.
2. 공통 대표 coarse/fine 차이.
3. 기존 known 영역과 대표 영역 사이의 signed 적분 차이.
4. 그 차이를 덮는 경계 띠의 절대 상계.
5. Branch ID, active plane, 근 폭, inf|f_theta|, event와의 거리.
6. Refinement 지배 moment.

대표 적분만 통과하면 window 절단의 영향이라는 증거이다. 대표도 실패하면
실제 적분 부족/평가 잡음을 조사한다. 경계 띠 상계만 실패하면 표현/root
문제이며 radial subdivision을 더 하지 않는다. 독립 확인은 지배 panel의
한 번의 이분 적분으로 제한한다.

## 필수 오차 장부

Global U로 나눈 node별 angular K bound는 target cold gate가 아니다.
Radial panel 전체의 경계 띠 상계, 가능한 모든 수신 target, source outgoing
차감과 target incoming 오차를 각각 계상하고 target U_j로 K 오차를 나눈다.
Gauss node 상계의 합만으로 연속 radial 상한을 주장하지 않는다.

기존 unknown 장부에는 THERMAL_RATE 상계가 빠져 있다. 실제 integrand에
이 항이 있으므로 E_Udot_band <= integral_band P J |div v| dX 상계를
추가해야 한다. 현재0을 무오차로 해석하지 않는다.

이번 제한 묶음에서 물리 상태를 다시 맞추는 trace projection을 추가하지
않는다. 고정 대표와 기존 whole-trace mismatch bound를 사용하고 영향받는
target별 gluing 오차를 보수적으로 계상한다. 이것이 cold 예산을 넘을 때
공통 연속 표현의 필요성을 구체적 실패로 남긴다. 경계를 평균해서 오차를
지우지 않는다. Target m,U,I 양수만으로 complement 질량행렬 양성을
증명하지 않으며, 그 검사는 weight solve/압축 전에 남아 있다.

두 삼각형의 기존 m,p,L,K,U 기준과 오차 귀속이 모두 완결돼야 통과이다.
실패하면 지배 항 하나와 수치 witness를 남기고 묶음을 닫는다. 새 PDE,
IC/궤도 reset, 열 보충, 독립 spin/kinetic reservoir를 도입하지 않는다.
W0 고정질량은 유한 고정 mu에서 존재하는 해군에 대한 조건부 관계이다.
임의 W0 overlap 또는 mu→infinity로 환원을 보장하지 않는다.
