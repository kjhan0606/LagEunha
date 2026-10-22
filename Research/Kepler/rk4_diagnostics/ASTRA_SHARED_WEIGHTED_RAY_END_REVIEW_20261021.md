# Astra: weighted angular / sqrt germ 종료 평가

2026-10-21. SHARED_CONTEXT 전체를 읽고 원본 repository origin LagEunha / HEAD 8e25dfe를 확인했다. 지정 weighted-ray/coverage 장부와 collected-radius, weighted-angular, origin-quadratic 구현을 읽었다. 새 계산·폴링·driver 코드 변경 없이 이 문서만 작성한다.

## 결정

**slot3203의 304-panel 실제 domain 및 geometry/density-clearance 완료를 채택할 수 있다.** 기존381+4축과 새223축, 17개 root-gap certificate를 합친 회계 및 동일 domain 분할을 근거로 한다. 이 채택을 새 manifest에 기록하고 원래 `mathematical_branch_review_pending` JSON은 역사 기록으로 보존한다. v4의30 pending에서29 pending으로 이동한다. 다른 context의 remaining31은 현황으로 사용하지 않는다.

현재 수식과 구현에서 채택을 뒤집는 치명적인 결함은 발견하지 못했다. 이것은 all12 analytic remainder나 실제 적분의 통과가 아니다. 기존 canonical/native field 및 libm engineering bound 전제도 유지된다.

## 1. weighted angular 항등식

q_lo=beta_lo*s+alpha_lo+delta_lo, q_hi도 같은 형태라고 하면 eta=eta0+z에 대해

    (q_lo+eta*(q_hi-q_lo))/s
    = (1-eta0)*beta_lo+eta0*beta_hi
      +[(1-eta0)*(alpha_lo+delta_lo)+eta0*(alpha_hi+delta_hi)]/s
      +z*[beta_hi-beta_lo+(alpha_hi-alpha_lo+delta_hi-delta_lo)/s].

구현은 이 동일 항등식을 eta 중심에서 먼저 모은 것이다. delta_lo/hi는 상수가 아니라 원래 implicit correction 함수의 enclosure여도 된다. 같은 s-domain에서 두 함수를 동시에 포함해야 하며, 공통 branch/parameter 및 correction의 유효 영역을 그대로 묶어야 한다. interval 합산이 남는 상관을 잃으면 보수적으로 넓어질 뿐 항등식은 바뀌지 않는다.

`S(0,eta.e)`는 z의 complex disk이다. eta 중심이 실수이고 입력 disk가 원래 panel을 포함한다는 전제하에 적절하다. s가0을 포함하는 경우 이 나눗셈을 그대로 정당화할 수 없다. 현재 regular-q 인증의 s 분리 전제를 유지하고 apex는 별도 Duffy 표현을 쓴다. 더 좁은 t enclosure를 사용한다고 actual field operation graph나 그 forward error가 저절로 바뀌지는 않는다.

## 2. origin quadratic의 analytic continuation

실수 U,V에서 P(t)=(U+tV)·(U+tV)=a+bt+ct², c>0, Delta=4ac-b²>0이다. 이는 complex t에 대한 holomorphic polynomial이며 conjugate t를 사용하는 norm이 아니다. 근은 alpha±i beta, alpha=−b/(2c), beta²=Delta/(4c²)다.

실수 중심 n, 반경 e인 disk에서 두 근까지의 동일 거리 d=sqrt((n-alpha)²+beta²)에 대해 인증된 g=d_lower−e>0이면 disk 안에 P의 영점이 없다. disk가 단순 연결이므로 P(n)의 **양의 실수 제곱근에서 시작하는 유일한 analytic sqrt germ**이 존재한다. P의 자연 interval image가 음의 실수축을 가로지르는 것과 이 영점 회피 조건은 다르다. principal sqrt의 임의 branch-cut 검사를 통과시킬 필요는 없다.

disk에서 |P|=c|t-z+||t-z−|≥c_lower*g²이므로

    |sqrt(P)| ≥ sqrt(c)_lower*g,
    |(sqrt(P))'| ≤ (|b+2cn|+2|c|e)/(2 sqrt(c)_lower*g).

중심에서 disk 안의 직선을 따라 적분하면 중심 sqrt enclosure에 e배 derivative 상계를 더하는 코드의 variation bound를 얻는다. 계수·거리·분모·중심값의 outward 오차를 포함한 구현 순서는 이 논리에 맞는다. Delta≤0, g≤0 또는 분모 인증 실패에는 이 방법을 적용하지 않는다.

실수축에서 P(t)>0이므로 이 germ은 실제 ray 길이의 양의 sqrt와 연결된다. 따라서 실수 s≥0에서 r=s sqrt(P(t))는 원래 reference 반경이다. 복소 extension에서는 |s|나 conjugation을 삽입하지 않는다. 겹치는 chart/증명 disk의 branch는 같은 실수 실제 영역과 양의 germ으로 결합해야 한다. complex overlap만 있고 실수 연결의 출처가 없는 경우 별도의 branch continuation 기록이 필요하다.

이 변경은 density 입력 radius의 analytic enclosure를 강화한다. 이후 r+epsrho, logistic, J 및 기타 분모의 clearance는 별도로 계속 검사한다. sqrt germ 성공이 이 나머지 관문을 대체하지 않는다. 향후 remainder 코드에서도 이 evaluator/branch ID를 실제로 소비해야 하며, 낡은 principal-sqrt 검사에 새 인증만 붙여서는 안 된다.

## 3. 다음 한 묶음: 남은29 실제 atlas

범위는 **14 oblique +15 origin**이다. source133의 네 영역, slot3200의3128 panels, slot3203의304 panels와 기존 완료 자료는 재사용한다. 새 target 예산은 승인된13074432 queries/51072 panels이며 차수16×16, tau/normalizer, profile, eps, 정밀도와 root guard는 불변이다.

### Oblique14

이미 인증된 방향을 그대로 사용한다. 실제 node/edge/loop의 투영 사건으로 slab과 다중 strip을 만들고 각 branch의 순서·활성 면·half-open ownership을 기록한다. 방향의 전역 최적성을 증명하거나 새 각도 후보를 탐색하지 않는다. projection-extremum이 추가로 필요한 경우 해당 실패 chart의 실제 경계에서만 같은 guard로 처리하고 canonical root 전체를 재생하지 않는다.

변환 determinant와 actual Jacobian을 포함한다. apex에 걸치는 부분은 기존 blow-up을 유지하며 Cartesian complex inverse로 대체하지 않는다. 같은 s,t 또는 regular-q가 생기면 이번 collected/weighted/germ 증명 경로를 사용할 수 있으나 그 chart의 실수 계수·root gap·branch 연결을 새로 인증해야 한다. source133의 수치 gap을 다른 chart에 복사하지 않는다.

### Origin15

기존 내부7/외부8 분류와 작은 양의 제외 간격을 보존한다. 원점 fan X=sR(t)에서 기존 전체 density를 사용하고 실제 target-plane boundary로 각 ray의 허용 s 구간을 정한다. 외부8의 tiny positive 하한을0으로 바꾸거나 작은 제외 영역을 채우지 않는다.

새 atlas는 밀도의 epsrho 규모와 실제 clipping branch 간격을 분리해서 다룬다. 기존 고정 pole clearance로 정한 radial 배치와 실제 사건을 사용한다. tiny 기하 간격을 smooth 전체 영역의 균일 panel 폭으로 강제하지 않되, 그 경계의 branch 증명을 생략해서도 안 된다. origin quadratic을 쓸 수 있는 실수 U,V chart에는 동일 root-gap 규칙을 적용하고, 일반 apex/비원점 chart에는 자동 확장하지 않는다.

### Native 저장·재사용 계약

각 typed panel에 요구 radius evaluator, weighted coordinate 여부, sqrt germ의 실수 anchor와 계수/증명 ID, 실제 branch/domain ID를 저장한다. native 값 계산은 원래 graph를 유지하고 proof evaluator 차이를 metadata로 명시한다. 모든12 moment에 대한 새 analytic remainder는 이번 geometry/density-clearance 단계와 별도 flag다.

이미 통과한 prefix/축/packet은 hash와 영역 동일성을 확인해 재사용한다. 미검증 suffix만 평가한다. 수정 전후 panel geometry가 동일하면 byte equality를 기록하고, 새29 atlas에는 원래 domain의 정확 coverage/무중복 증명을 요구한다. 정상817152점 field proof나 완료 quadrature를 다시 실행하지 않는다.

## 종료 조건

29개 각각에서 실제 domain coverage, 이벤트/branch 연결, geometric Jacobian, density pole clearance 및 고정 target 합산 비용까지 완료하거나 최초 미해결 witness를 보존한다. 실행 전에 policy와 pending 목록을 동결하고, 실패 후 새 guard/order/precision/cap 후보의 자동 반복을 하지 않는다. 큰 packet은 승인된 native streaming 정책으로 계획한다.

모두 닫히면 전체 target별 실제 비용과 reuse manifest를 확정하여 다음 all12 analytic certificate/보존적 수치조립 단계로 연결한다. 하나라도 남으면 해당 이유와 partial atlas를 보존하고 이 묶음을 종료한다. 현재 단계에서 full quadrature, RK4 또는 장시간 실행을 통과로 선언하지 않는다.

일반 evolved heat 및 native CUDA/MPI/OpenMP 시간 갱신, 실제 t6/8/10와 GIZMO 최장공전 안정성은 계속 미달이다. 이번 채택은 동일 상태의 실제 전달 영역을 인증하는 구체적인 진전이며 물리 상태·모델·목표 변경은 아니다.
