# Astra: 실제 query와 전체 오차 예산 종료 평가

2026-10-21. 계획·평가 Astra, 구현·실행 driver. SHARED_CONTEXT, v8 독립 장부, source301 실제 첫 점, 현재 R80 radial 코드와 실제 지배 panel CSV를 읽었다. HEAD 8e25dfe 확인. 계산·작업 폴링은 하지 않았다.

## 판정: 두 독립 관문이 남았다

전체871 analytic certificate 통과 및 clipped25 source 실제 적분 완료는 유효하다. source301의 실제 첫 점에서 멈췄으므로 2816점이 빠져 있고 전체 target은 미완료다. caller와 원본 열/상태 보존도 확인된 범위다.

더 중요한 것은 source301과 무관한 **현재 오차 장부의 확정적인 실패**다. potential 알려진 항만 원래 tau의 1.2992418095294642배다. 이는 실제 물리오차가 1.299배라는 뜻이 아니라 현재 certificate로 publish가 불가능하다는 뜻이다. 구성은 retained4 forward .6517281869, 새25 forward .0016426536, 전체 analytic .6458709691이다. source301을 고쳐도 그대로는 통과하지 못한다.

whole4는 이미 Q128 합산이고 .65tau 대부분은 point-forward 항이다. 합산 순서를 바꾸거나 기존 group을 다시 더해서 해결할 수 없다. old4 R80 재적분과 old25 재적분은 하지 않는다.

## 다음 한 묶음: 301의 제거가능 특이점 처리 + 저장 적분의 potential 나머지 강화

### A. 실제 source301 첫 점

준비 중인 단일 query trace에서 실패한 정확한 나눗셈의 numerator, denominator nominal/radius, caller 경로를 먼저 보존한다. 실제 Q128 mapping/J/23planes가 통과했다는 사실은 R80 분모 실패의 원인을 확정하지 않는다. 원래 값 계산을 바꾸지 않은 trace까지만 원인 확인으로 인정한다.

만약 실제 실패가 radial의 term/s, dtail=4term/s라면, |s.n|<=s.e인 구간에서 나타나는 제거가능 특이점이다. 함수를 바꾸지 않는 다음 고정 recurrence를 채택한다:

    T4=c4*s^4, D4=4*c4*s^3, c4=.56396484375
    a_n=(-.75-n)/(n+1)
    T_(n+1)=a_n*s*T_n
    D_(n+1)=a_n*(T_n+s*D_n)

tail과 derivative tail을 각각 합산한다. s로 나누지 않는다. 현재 고정256 cap, R80, 원래 함수/series-domain 조건을 유지한다. 남은 항도 r=|s|의 상계로 나눗셈 없는 형태로 제한한다. 예를 들어 N까지 합산했을 때 |c_(N+1)|를 사용한

    value_tail <= |c_(N+1)| r^(N+1)/(1-r)
    derivative_tail <= |c_(N+1)| r^N*((N+1)/(1-r)+r/(1-r)^2)

가 r<1에서 유효하다. 상계/계수 계산은 outward로 한다. s=0도 같은 연속 한계로 처리하며 epsilon 대입이나 numerator/denominator radius 삭제는 금지다.

이것은 같은 해석적 함수의 새로운 연산 그래프다. 기존 R80 graph와 모든 bit가 같다고 주장하지 않는다. 분모가 분리된 기존 경로는 유지하고, 정확한 분모 미분리 조건에서만 새 경로를 쓰면 기존25 증거를 그대로 유지할 수 있다. 실제 원인이 물리 J/edge denominator라면 이 수정은 적용하지 말고 해당 실패를 보존한다.

원래 실패한 한 점에서 same-function oracle 및 field/all12 bound를 확인한다. 그 후 source301의 기존11 typed panel/2816점만 누락분으로 완결한다. 기존871 certificate, old25/whole4 값은 재계산하지 않는다. 새로운 오류가 나타나면 자동 구간 분할이나 cap 상향 없이 종료한다. 준비된 real-slice helper를 같이 끼워 넣지 않는다. 그것을 쓰려면 별도로 동일 전체 S와 점·branch의 포함 관계를 증명해야 한다.

### B. 저장 적분값을 바꾸지 않는 고정 이차 polynomial 제거

가장 직접적인 선택은 **정규화된 적분 변수에서 전체 potential integrand의 고정 이차 Taylor polynomial을 제거해 Cauchy majorant를 계산**하는 것이다. 적분 rule 차수는 그대로 Gauss16이다. polynomial 차수를 탐색하지 않는다.

실제 CSV에서 지배적인 251/17,251/20,249/23,249/27은 모두 whole=0, policy0, implicit_observed=0이다. 따라서 q-repaired27 또는 whole4를 먼저 손보는 것이 아니다. 최종 typed packet의 chart, 구간/radius, lower/upper face와 metric을 그대로 사용한다. implicit_observed=0이라는 이유만으로 Jacobian이나 원래 clipped geometry를 생략하지 않는다.

각 기존 panel의 [0,1]^2에서, 물리적 integrand는

    g(u,v) = -rho0(X(u,v))*Phi_abs(chi(X(u,v)))*J_ref_panel(u,v)

이다. 여기 Phi_abs=1/sqrt(|x|^2+eps_gravity^2). 기존 potential 부호와 같은 식이다. reference pullback에서는 rho0/J와 physical J가 정확히 소거되므로 그 상쇄를 먼저 유지한다. 압력·열·운동에너지 재정의가 아니다.

중심 c=(1/2,1/2)에서 AD2로 하나의 고정 polynomial

    p2(z)=a0+a·(z-c)+.5*(z-c)^T A*(z-c)

를 만든다. 계수는 저장된 정확 binary nominal 값으로 정의하고, true derivative와 계수의 차이를 bound에 포함한다. potential 값 또는 potential만의 공간 polynomial을 빼는 것이 아니라 **Jacobian/원래 density/implicit graph까지 포함한 g의 polynomial**이어야 한다. 물리좌표에서 이차라는 이유만으로 Gauss-exact라고 주장하지 않는다.

Gauss16 tensor rule은 각 변수31차까지 정확하므로 E(g)=E(g-p2)이다. 실제 저장 quadrature sum에는 p2를 빼거나 다시 더하지 않는다. 이 항등식은 나머지 증명에서만 사용한다. 기존 finite nodes/weights의 오차는 원래 point-forward 장부에 남는다.

기존 두 축 complex 근방과 같은 branch에서 g-p2의 absolute majorant를 계산한다. 큰 g와 p2를 각각 interval 계산해 마지막에 빼지 말고, affine dependence/이차 계수를 공유하는 Taylor-model arithmetic으로 cancellation을 유지한다. 구현 대안으로 같은 convex analytic 영역의 Hessian에 대해

    g(z)-p2(z) = coefficient_errors
        + integral_0^1 (1-t) dz^T*(H_g(c+t*dz)-A)*dz dt

를 outward로 감싸도 된다. 이 식은 complex line segment 전체의 해석성과 Hessian 포함을 요구한다. center 표본과 관측 Hessian 변화만으로 bound를 만들지 않는다. implicit graph가 있으면 그 미분도 같은 인증된 branch에서 계산한다.

적용 범위는 기존871 panel의 **potential certificate만 한 번**이다. 나머지11 모멘트 certificate와 모든 저장 적분값은 불변이다. 각 panel에 oldR, newR, 계수, domain/branch hash를 남긴다. 둘 다 유효한 증명이면 min(oldR,newR)을 사용할 수 있으나 작은 관측차이를 bound로 대체하지 않는다. 실제 특이점이 나타나거나 새 증명이 더 넓으면 기존R을 유지한다. 다른 차수·정밀도·guard·panel topology로 자동 재시도하지 않는다.

이는 analytic certificate만 재평가하는 범위다. 220160개 완료 quadrature 점의 profile/field/all12를 다시 계산하지 않는다. 필요한 중심/complex AD2 평가는 certificate 비용으로 별도 기록한다.

## 조립 순서와 종료 기준

1. 원래 snapshot/SHA, typed871 packet, 기존25 source 합 및 whole4 group을 고정한다.
2. source301 실제 나눗셈 원인에 맞는 최소 graph 수정과 첫 점 검사를 완료한다. 별개로 위 고정P2 potential certificate pass를 한 번 수행한다.
3. 새 analytic subtotal이 기존 known forward와 합쳐 이미 tau 이상이면, source301 적분을 굳이 실행해도 publish 불가능하므로 그 상태에서 reject하고 종료할 수 있다. 이는 missing301을 성공 처리하는 것이 아니다.
4. 예산 가능성이 남으면 missing301의2816점만 수행해 source별 all12를 보존한다. 새26+retained4를 각각 한 번 합산한다. source301의 실제 forward오차, 전체 addition, 최종 R80 serialization까지 넣은 **모든 원래12tau**가 유일한 publish 관문이다.
5. 여전히 실패하면 caller를 보존하고 종료한다. no old4 R80 재실행, no old25 재적분, no normalization/heat/state 변경. 두 고정 구현의 결과를 모아 한 번에 다음 판단을 한다.

성공하면 기존 계획대로 all-target registry/incidence/비용 및 source별 보존 조립으로 넘어간다. 일반 evolved heat transport, 실제 native 시간적분 및 GPU 전달 인증은 아직 남는다. 이번 4880초 CONTROL 호출을 장시간 효율의 근거로 삼지 않는다. 원래 CUDA/MPI/OpenMP Kepler t6/8/10과 GIZMO 최장공전 목표는 그대로 미달이며, 단일 target 인증을 최종 목표로 축소하지 않는다.
