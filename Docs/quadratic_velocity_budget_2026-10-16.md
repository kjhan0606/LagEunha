# 공분산 예산으로 제한한 이차 속도 복원

**현재 상태:** 아래 첫 구성은 v1의 역사 기록이다. v1은 순간 모멘트 검사는
통과하지만 순위 변화의 불연속으로 진화에 실패했다. 현 코드는 아래 v2이며 별도
시간 진화 검증 중이다. 일반 Euler 정확도는 확보되지 않았다.


## 목적과 적용 범위

현재 선형 속도 복원에서 남는 R을 전부 위치와 독립인 속도 표본으로 운반하면,
실제 셀 내부 속도 곡률도 무작위 속도분산처럼 면 질량유량에 들어간다.
새 후보는 이 분산의 일부를 결정적인 이차 속도장으로 설명한다.
이는 이번 연구에서 구성한 수치 닫힘이며 일반 Euler 정확도나 충격파 수렴을 이미
증명한 방법으로 해석하지 않는다. 현재 연속 positive moment density에만 적용한다.

저장 상태 m,q,Q,p,D,Pi,E와 양의 일차 변수 R=L Lᵀ는 바꾸지 않는다.
정확한 물리해는 초기화와 채점에만 사용한다. 생성점은 복원한 유체 평균 속도로 움직인다.

## 양의 운동 분포를 구성하는 방법

질량 중심 c=q/m, 평균 속도 v=p/m와 다음을 둔다.

\[
S=Q-qq^T/m,\quad C=D-qp^T/m,\quad B=C^T S^{-1},
\qquad v_A(x)=v+B(x-c).
\]

양의 수치 밀도 rho와 같은 양의 체적 quadrature를 사용한다.
국소 좌표 z=(x-c)/sqrt(area)에서 세 이차 단항식의 1,z 성분을 빼서 psi를 만든다.

\[
\langle\psi\rangle_\rho=0,\quad
\langle (x-c)\psi^T\rangle_\rho=0,\qquad
G=\langle\psi\psi^T\rangle_\rho,
\quad\langle f\rangle_\rho=m^{-1}\int\rho f.
\]

구현에서는 quadrature의 가중 최소제곱 투영과 Gram QR을 이용한다.
따라서 아래 보존 관계도 해당 수치 적분에서 성립하며, 연속 적분의 quadrature 오차는
별도로 검증해야 한다.

이웃 셀 평균 속도와 저장된 m,q,Q에서 이웃의 <psi>를 구해 이차 계수 H0를 추정한다.
정보가 부족하면 기존 affine 복원과 전체 R을 사용한다. 평균/혼합모멘트에 영향을 주는
상수·선형 성분은 새 계수에 허용하지 않는다.

공분산 좌표에서 W0=sqrt(m)L^+ H0, C0=W0 G^(1/2)를 만든다.
일차 인수 L의 range 밖 후보 성분을 제거하고, C0의 특이값을
상한 1-1e-12인 제약 안에서 선택한다. 이는 **후보 계수**의 에너지 제약이며
저장된 R이나 음의 공분산 고유값을 잘라내는 연산이 아니다.
제약 후 W와 대칭 T를 구성하면

\[
WGW^T+TT^T=I,\qquad H=LW/\sqrt m.
\]

동일 확률1/4의 alpha_k=(±sqrt2,0),(0,±sqrt2)를 이용해

\[
v_k(x)=v_A(x)+\frac{L}{\sqrt m}\{W\psi(x)+T\alpha_k\}
\]

를 정의한다. 모든 표본 확률과 rho가 양수이다. 평균 속도는 v_A+H psi이고,
위치와 독립인 남은 속도분산은 L T Tᵀ Lᵀ/m이다.

\[
mHGH^T+LTT^TL^T=LL^T=R.
\]

psi의 두 직교 조건 때문에 p,D가 유지되고, 위 식 때문에 Pi도 유지된다.
따라서 m,q,Q,p,D,Pi,E와 열압력은 원 상태와 동일한 모멘트를 갖는다.
순수 이차 속도장의 uniform density 검사에서 이차 평균장을 재현했다.

## 공분산 갱신에서 필요한 일관성

R은 여전히 **최적 affine 투영에 대한 잔차**이다. 이차 평균장을 사용해도
cov_A의 B를 이차 속도의 점별 gradient로 바꾸면 안 된다.

각 표본의 affine 잔차 좌표 xi_k=W psi+T alpha_k를 사용한다.
유출 손실은 양의 K=Σ(scale/m) w_out xi_k xi_kᵀ로 계산하고,
유입 Gram은 donor 표본과 receiver의 v_A 사이 차이로 계산한다.
receiver의 이차 평균장을 빼면 저장 R의 정의와 달라진다.

\[
\dot R=-BR-RB^T-LKL^T+Q_{in}Q_{in}^T.
\]

기존 양의 Sylvester/Gram 갱신을 이 K와 Q_in에 적용한다.
공유 압력의 affine 가속도 일은 p,D만 사용하며, 이차 추가항의 1,x 직교 조건과 양립한다.
원시 모멘트 미분에서 직접 유도한 Rdot와 이 식의 일치를 회귀 검사했다.

## 검증과 실행 분리

- EG 관련 회귀122개 통과: 이차장 재현, m/p/D/Pi 보존, R 예산 제한,
  양의 유출 K, 원시 Rdot 대 Gram 식, 질량 전용 경로 대 full RHS, affine R=0 patch.
- 선형 속도를 전제로 한 exact-bin transport/독립 BoundaryDensity와의 조합은 명시적으로 거절한다.
- 실험 flag는 LAGFORCE_QUADRATIC_VELOCITY=1, 설정값은 velocity_reconstruction=quadratic_budget.
  기본 affine 방식과 구별하며, 이전 affine RHS를 새 닫힘의 source로 재사용하지 않는다.
- EH는 실제 초기 모멘트를 그대로 사용해 각 닫힘의 RHS를 별도로 평가했다.
  N32의 정확 Euler 질량 변화율L1/진폭은 .0799222→.0154508,
  압력 평균 변화율은 .0887228→.0460950, 속도 평균 변화율은 .0824181→.0782016.
  공분산 분할 상대 잔차~1e-15, source 질량 잔차~2e-19.
- 순간 변화율 개선을 종료시각 정확도나 시간차수 개선으로 대신하지 않는다.
  EI/EJ/EK에서 동일 새 source/RHS의 균일10/20/40스텝 대조를 먼저 수행한다.

코드: Docs/polybnd/quadratic_velocity_reconstruction.py, kinetic_factor_rhs.py,
smooth_mass_path.py. 스냅샷과 상세 실행 결과는 forty_second_2d_batch_2026-10-16.md에 기록한다.


## v1 실패와 v2의 연속 구성 (2026-10-16 19:29)

v1의 whitened SVD 제약은 R이1랭크에서2랭크로 변할 때 연속인 선택이 아니다.
거의0인 횡분산에 비해 큰 후보 횡곡률이 종방향 곡률까지 억제할 수 있다.
시간 간격 반감에도 생성점 경로 결함 비가5405→5427→5444로 남았다.
따라서 정적 실현 가능성은 동적 닫힘의 충분조건이 아니다.

v2는 L=U diag(s1,s2) Vᵀ, s1>=s2에서 A0=sqrt(m) H0 Ghalfᵀ를 만들고,
첫 물리 방향의 a=U[:,0]ᵀ A0를 길이(1-1e-12)s1 안에 제한한다.
정규화한 행을 b, w=1-(s2/s1)^2라 두면

- C_eig=[w b;0], C_std=V C_eig, W Ghalfᵀ=C_std.
- d=sqrt(1-w²|b|²), T=V diag(d,1) Uᵀ.
- H=L W/sqrt(m), L T=U diag(s1 d,s2) Uᵀ.

따라서 m H G Hᵀ+L T Tᵀ Lᵀ=R은 유지된다. L T는 물리 좌표에서
유일한 대칭 양의 잔여 제곱근이므로 직교 factor gauge가 속도 표본을 바꾸지 않는다.
s2→0에서 작은 횡성분이 종방향 배분에 영향을 주지 않는다.
s2→s1에서는 w→0으로 모호한 주고유벡터의 영향이 사라진다.
R→0이면 평균 곡률과 잔여 표본도0으로 간다. 이는 닫힘 구성의 연속성 논거이며
원시 전역 모멘트 차감의 부동소수 안정성이나 PDE 해 존재 증명과는 구분한다.

후보 fit은 반경3sqrt(area), Wendland C2 가중치(1-r)^4(1+4r),
정규행렬에 고정1e-10 ridge를 더한다. 이웃이 반경을 넘어가도 가중치가 연속적으로0이
되며, graph flip 때문에 후보가 갑자기 바뀌는 것을 피한다.
일반2D 이차장 완전 재현은 더 이상 요구하지 않는다.1방향 이차장의 재현과
일반 이차장의 weighted RMS 개선, 전체 저장 모멘트 일치, 양의 K/Rdot 항등식을 검사한다.
새 닫힘 tag는 `quadratic_budget_continuous`이다. EN127개 회귀 통과 후 EO에서
순간 RHS, 이어서 균일 시간 간격 시험을 별도로 검증한다.

## v3: 물리 오차 목적함수의 고정횟수 projected gradient

v2는 작은 횡분산이 종속도를 억제하는 불연속을 없앴지만 일반2D의 두 번째
곡률 방향을 남은 독립 분산으로 운반한다. 이는 회전 유동의 밀도 오차에서
불리할 수 있다. v3은 이 방향도 사용하도록 별도의 연속 후보 선택을 구성한다.

A0=sqrt(m) H0 Ghalfᵀ, A=L C로 두고 다음 볼록 문제를 생각한다.

\[
\min_{\|C\|_2\le 1-10^{-12}}\frac12\|LC-A_0\|_F^2.
\]

목적함수는 **물리 속도 곡률**의 weighted L2 차이다. L⁻¹로 후보를 whitening하지 않는다.
실제 구현은 ||L||₂로 L,A0를 정규화하고 고정64회의 accelerated projected-gradient를
수행한다. 각 회는 gradient step 후 C의 특이값을 후보 제약 안에 선택한다.
이는 저장 공분산을 바꾸는 연산이 아니다. 반복 횟수를 상태별 중단 판정으로
바꾸지 않아 유한 반복 선택도 연속인 합성으로 남는다.

정확한 볼록 최적해를 이미 구했다고 주장하지 않는다.64회 뒤의 가용 후보를 사용한다.
L≠0에서 gradient와 spectral-ball의 유일한 Euclidean projection이 연속이고,
L→0에서는 ||A||₂≤||L||₂이므로 물리 곡률도0으로 간다. 직교 factor gauge에 대해
공변인 반복이어서 LC는 gauge와 무관하다.

최종 C=U diag(s) Vᵀ에서 T0=U diag(sqrt(1-s²))Uᵀ를 만든다.
L T0의 polar 회전으로 좌표 인자를 바꿔 L T가 물리 대칭 양의 제곱근이 되게 한다.
거의 같은 두 공분산을 빼서 잔여 고유값을 계산하지 않는다.

W Ghalfᵀ=C, H=L W/sqrt(m)이면 m H G Hᵀ+L T Tᵀ Lᵀ=R이다.
이전 p,D,Pi/양의 K/원시 Rdot 항등식은 그대로 적용한다.
일반2D 이차장 두 방향의 재현(검사 절대오차2e-8), 순위 경계, factor gauge,
기하 질량 root와의 일치 검사를 포함해 FN의139개 관련 회귀를 통과했다.

별도 tag는 `quadratic_budget_projected`,
flag는 `LAGFORCE_QUADRATIC_BUDGET=projected`이다. 기본은 검증된 v2의 dominant 선택이다.
velocity flag 자체가0이면 affine 기준 계산을 유지한다. 각 모델의 고유 RHS만 재사용한다.

FN은 ES의 **같은 실제 vortex 초기 배열**에서 v3 RHS를 새로 구했다.
N32 질량 변화율L1은 .01372939→.00182687, 압력 평균 변화율/drop은
.05446088→.01402300, 속도 평균 변화율/peak은 .02303803→.01683917.
초기 속도장 오차/peak은 .00141289→.000161158로 감소했다.
공분산 분할 상대 잔차~2e-15, source 질량 잔차~3e-18이다.
이는 순간 일관성 결과이며 종료 물리 정확도·시간차수·shock entropy 검증은 아니다.
FU는 첫 v3 N16 vortex 장기 진화다. 기존 affine/v2의 승인 계산은 별도로 보존한다.
