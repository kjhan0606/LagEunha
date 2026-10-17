# 양의 평활 kinetic splitting 후보

## 목적과 명시적 모델 변경

기존 유량의 max(s,0), min(s,0)은 표본의 상대 법선 속도 s=0에서 미분이 꺾인다.
주어진 모멘트 경로의 기하 ODE에서 많은 내부 단계가 필요해, 다음 유량을 별도
후보로 시험한다. 이는 root 허용오차 변경이 아니라 **수치 수송 플럭스 변경**이다.
원 저장 m,q,Q,p,D,Pi,E, 원 질량, material generator 법칙은 바꾸지 않는다.

\[
s_\pm=\tfrac12(s\pm\sqrt{s^2+\epsilon^2}),\qquad
\epsilon_{ij}=\kappa\,\bar c_{ij}\,\frac{\bar A_{ij}}{\ell_0^2}.
\]

c는 수치 평균 열압력과 밀도로 구한 음속, A는 실제 셀 면적의 산술평균이다.
κ=0이면 종전 유량과 동일하다. 비교 시험의 길이 척도 ℓ0=1이다.
고정 κ에서 ε는 격자 길이의 제곱과 함께0으로 간다.

왼쪽 표본 질량 가중치는 rhoL*s_plus/4>0, 오른쪽은 rhoR*s_minus/4<0이다.
따라서 원래 유출 K와 유입 Gram의 양성 논거가 유지된다. 공유 면의 수송은
양쪽 셀에 같은 값/반대 부호로 들어가 질량·운동량·각운동량·에너지를 보존한다.
표본 속도와 1,x 직교 이차 복원의 모멘트 예산은 그대로 사용한다.
새 좌우 trace가 같으면 s_plus+s_minus=s여서 affine/균일밀도 patch를 재현한다.

## 고정된 수치 상태에서 weightdot의 유일성에 대한 논거

양의 rho와 ε>0에서 d(s_plus)/ds와 d(s_minus)/ds가 모두 양수다.
따라서 질량 잔차의 음의 weightdot Jacobian은 모든 실제 내부 면에 양의 계수를 갖는 graph Laplacian이다.
비어 있지 않은 convex 셀의 연결된 partition에서 전역 gauge를 제거하면 양의 정부호다.
면 유량은 weightdot 차이에 대해 엄격히 단조이며, 그 음의 원시함수로 이루어진
퍼텐셜은 gauge 공간에서 엄격 볼록하다. 큰 차이에서는 한쪽 유입이 선형으로 커져
퍼텐셜이 이차로 증가하므로 coercive이다. 이 조건 아래 각 고정 상태의 mass root는
존재하고 gauge를 제외하고 유일하다.

이는 이번 구성의 유한차원 논거다. 밀도 재구성의 실현 가능성, 셀 붕괴,
전체 결합 ODE/PDE의 전역 해 존재, shock entropy는 이 논거로 해결되지 않는다.
κ>0의 유일한 weightdot과 원래 κ=0 문제의 유일성을 구분한다.

## 구현과 검증

- `kinetic_speed_smoothing.py`는 큰 음수 속도의 작은 양의 tail도 차감 소실 없이 구한다.
- `kinetic_factor_rhs.py`, vector mass root, 질량 전용 기하 RHS와 독립 감사에서
  같은 명시적 플럭스를 사용한다. 감사는 여전히 실제 저장 다항식 미분을 사용한다.
- 저장 source RHS의 κ와 길이 척도가 runtime/config와 다르면 재사용을 거절한다.
- 연속 positive Density quadrature와 조합한다. exact-bin/독립 BoundaryDensity와
  조합해 다른 유량을 묵시적으로 섞지 않는다.
- 부호/분할 항등식/미분, graph Jacobian 유한차분, 양의 Rdot 항등식,
  전체 RHS 대 독립 질량 감사, affine patch와 잘못된 source 거절을 포함151개 검사 통과.

## 현재의 성능 증거와 한계

보존된 FT step320의 **같은 원시 모멘트/인자 경로**를 사용했다. 종전 기하 ODE는
425회 평가였다. 확장 중간 산술만 쓰면356회/153.8초이고, κ=1 평활 유량까지 쓰면
164회/71.2초였다. 평활 후보의 독립9점 감사 비는.00010512,
생성점 결함5.08e-13이었다. 원래 κ=0 확장 산술 경로는257점에서도 감사가
통과하지 않아, 계산 비용과 시간 적분의 규칙성 문제가 남았다.

이는 모멘트 경로를 고정한 조건부 기하 비교다. 새 유량의 완전 유체 진화나
정확 물리장 통과 결과가 아니다. GB에서 같은 물리 vortex 초기 배열의 고유
RHS를 구하고, 종료 물리 오차·균일 시간 반감·κ 민감도를 별도로 검사한다.

확장 산술은 이미 저장된 초기 raw moment의 잃은 자릿수를 복구하지 못한다.
밀도 내부 Newton의 중단 잔차가 시간 함수에 만드는 작은 점프도 확인 중이다.
FZ/GA는 중간 상태의 density fit을1e-15까지 풀고 κ=0/1을 비교하는 후속 실행이다.
기존 native stage/승인 state의 float64 재구성과 source RHS는 그대로 유지한다.
