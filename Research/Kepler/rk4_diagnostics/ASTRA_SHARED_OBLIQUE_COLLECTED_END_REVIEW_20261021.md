# Astra: oblique 공통좌표 종료 평가

2026-10-21. SHARED_CONTEXT 전체와 원본 origin LagEunha/HEAD 8e25dfe를 확인했다. 지정 정책, oblique point/v-domain 및 Duffy 코드, 두 suffix native contract, 실제 clearance CSV와 저장 source/input/output hash 기록, 독립 감사 코드를 읽었다. 새 계산·폴링·driver 코드 변경은 하지 않았다. 이 문서는 정적 수학 검토이며 저장 hash를 새로 계산한 검증은 아니다.

## 판정과 채택 경계

**동일 함수의 enclosure를 강화한 변경으로 판단하며,14개 oblique chart의 geometry/density clearance 완료를 채택하는 수학적 근거는 충분하다.** 기존2805 passing axes와 새5963을 합친4384 panels/8768 axes가 native에서 모두 통과했다. 같은 상태·domain·panel을 유지한 suffix 처리라는 계약도 적절하다.

초기 검토 후 같은 요청의 독립 감사가 완료되어 SHARED_OBLIQUE_COLLECTED_SUFFIX_20261021.json을 추가로 읽었다.348 domain의 correlated quantities 및 caps,8768 axes, 기존 packet 보존 검사가 PASS여서 이번 subchart 채택의 회계 관문도 닫혔다. 별도 새 수치 실행이나 재평가는 하지 않았다.

과거 phase66 실패13개 중 실제 denominator guard 실패는7개, binomial guard 실패는6개이며, 별도7169 final-point binomial 실패도 있었다는 운전자의 구분을 유지한다. 공통좌표 수집이 모두 해소했더라도 이들을 모두 실제 denominator zero나 동일 guard 실패로 설명하지 않는다. 작은 상계에서 성공한 것은 이전 실패가 실제 함수 특이점을 증명하지 않았다는 근거이며, 이번 유효성은 각각의 새 denominator/binomial/domain 인증에서 나온다.

이 판정으로 전체29 slot을 완료로 바꾸지 않는다. slot의 다른 chart/strip 연결, origin15, 전체 atlas 회계와 all12 remainder는 미완료다. 이번 복소 enclosure 계산은 physical quadrature나 시간진화가 아니다.

## 1. 공통 u/gamma 수집은 같은 함수인가

원래 X=A+(B−A)s+(C−B)q에

    v=beta*u+gamma,
    s=su*u+sv*v,  q=qu*u+qv*v

를 대입하면

    X=A+[du+beta*dv]u+dv*gamma,
    du=(B−A)su+(C−B)qu,
    dv=(B−A)sv+(C−B)qv

이다. 새 affine_coordinate는 이 항등식을 interval magnitude 이전에 모은다. barycentric combination, edge denominator/parameter 및 좌표의 동일 선형 결합에도 적용할 수 있다. beta/gamma가 disk이면 이들 사이의 남은 상관을 잃어 상계가 넓어질 수 있지만 같은 실제 변수값을 포함하는 한 enclosure의 유효성은 유지된다. 이는 새로운 DOF나 상태 refit이 아니다.

point의 eta 중심과 편차를 beta/gamma에 분배하는 식도 원래 v_lo+eta(v_hi−v_lo)와 같다. width는 독립적으로 hi−lo의 같은 식으로 남겨야 하며, 좌표 enclosure가 좋아졌다는 이유로 measure의 width 오차를 삭제하지 않는다. 실제 면적은 기존 det[e,d]·width와 current-map J 및 panel scale을 각각 필요한 moment 계약에 맞춰 사용한다.

## 2. Cartesian derivative 및 implicit neighborhood

Duffy basis에서 barycentric Jet의 gradient는 원래 Cartesian grad(lambda)이며, 새 좌표식은 그 Jet를 평가할 점 집합을 감싼다. beta/gamma를 Cartesian 미분 변수처럼 추가 미분하는 방법이 아니다. oblique face의 (s,q)→(u,v) 변환 행렬은 원래 고정 선형변환이므로 H_uv=D H_sq D^T도 그대로다. 이 구별이 유지된 현재 구현은 같은 field 값·Cartesian gradient·Hessian의 enclosure다.

implicit predictor v=alpha+beta*u+delta의 Hessian은 delta의 **전체 prior disk**에서 먼저 감싸야 한다. geometric_v_prior는 radial/edge denominator를 원래 predictor에서 양수로 분리한 뒤 delta가 이를 절반 이상 훼손하지 않도록 제한하고, 실제 Möbius parameter 변화를 binomial clearance에 연결한다. 이어 실제 코드가 R≤prior, q<1 및 residual/prior+q<1 형태의 invariance를 요구하므로, 성공한 좁은 root에서만 Hessian을 평가해 순환적으로 정당화하는 구조는 아니다.

특히 parameter 변화의 분자는 gb*D−N*(ga+gb)이며 delta 방향에서 상수다. |D_new|≥D_floor/2를 사용한 2|분자|/D_floor²의 Lipschitz 상계는 타당하다. 이 증명은 전체 complex u disk의 predictor와 그 delta neighborhood를 포함해야 한다.348개 새 domain 감사는 이 연결과 사용한 branch/panel/axis를 확인해야 하며,139개 이전 감사만으로 새348개를 대신하지 않는다.

원래 physical field의 libm/engineering 전제와 finite tail은 그대로다. 공통좌표 변경은 그 전제를 강화하거나 새로운 물리 정확도 증명을 주지 않는다. 실수 branch와 complex implicit predictor의 연결도 기존 event ID/실수 root로 유지한다.

## 3. ABI 및 다음 all12에서 필요한 자료

새 AffineAngular는 private evaluator 자료다. historical RayTrace를 이 새 layout으로 읽지 않는 계약을 유지한다. 기존 Panel192/Axis656/Frame2064/Linear144/VDomain880/SavedPoint1488의 명시 width와 epoch/hash 검사는 필수다. 미래 typed panel에는 요구 evaluator/helper SHA를 명시하여 낡은 common-u만 수집하는 evaluator로 증명을 해석하지 않도록 한다.

감사 코드가 구별하듯 새 full-field snapshots는 이전2005+새5963개다. **원래800 passing axes의 complete field 값은 all12 moment bound에 아직 필요한 별도 의존성**이다. geometry-clearance pass bit만으로 all12 sup bound를 만들 수 없다. 우선 해당 원래 snapshot이 보존됐는지 읽어 연결한다. 없다면 필요한800개 complex proof query만 같은 frozen 입력/함수로 복구하는 범위를 manifest에 따로 고정한다. 이것은 물리 quadrature나 전체817152점 재생이 아니며, 새 데이터 회복임을 숨기지 않는다. 재사용 가능한 기존 field/나머지 증명은 반복하지 않는다.

## 4. 다음 한 묶음

**origin15 actual atlas + 누락88 clipped/169 strips의 최소 metadata 회복 + 전 target analytic 준비 manifest**를 하나의 연결 묶음으로 진행한다. oblique14는 이번4384 panel 증거를 재사용한다.

1. Origin15는 기존7 inside/8 outside와 tiny gap을 보존하고 실제 ray별 target 경계를 사용한다. 원래 epsrho=1e−6/logistic 전체를 유지한다. apex Duffy, collected/weighted radius 및 양의 sqrt germ의 기존 계약을 적용하되 새 chart의 계수·root gap·branch 연결을 인증한다. radial pole의 분리 실패를 상태/밀도 수정으로 우회하지 않는다.
2. 누락88 chart는 먼저 기존 typed panel/domain/root manifest에서 직접 결합한다. 없을 때 실제 필요한 chart의 arrangement/frame/radial trim만 원래 frozen geometry 경로로 회복하고 즉시 typed snapshot을 남긴다. 옛1970 root windows나 canonical global root 작업을 재생하지 않는다. 실제 기존 자료로 복원 불가능한 사건은 명시적 pending이며 추측한 기본값으로 채우지 않는다.
3. 정확 domain coverage, half-open ownership, event enclosure의 오차 귀속과 measure를 조립하고, target당13074432 queries/51072 panels 예산 안에서 전체 비용을 합산한다. completed panel은 byte/hash 또는 증명된 domain 동치로 재사용한다. 저장 없는 metadata 복구를 field/적분 재실행으로 확대하지 않는다.
4. 같은 묶음에서 all12의 입력 의존성 표를 완성한다. 각 panel의 full complex field, frozen heat context, endpoint error, analytic rule/P2-potential proof, native forward 및 final R80 합산·변환 항을 빠짐없이 연결한다. 새 analytic proof는 geometry/field 의존성이 완결된 해당 panel만 대상으로 할 수 있다. 전체 field query나 물리 적분을 선행하지 않는다. target 총 예산은 부분 source 예산보다 우선하며 원래 tau는 유지한다.

한 묶음 종료 시 actual atlas/geometry clearance/field completeness/all12 remainder/cost를 별도 flags로 보고한다. 실패가 있으면 최초 실제 witness를 남기고 cap·정밀도·차수·tau 변경이나 자동 fallback 반복을 하지 않는다. all12 및 실제 비용이 닫히기 전 full quadrature는 차단한다.

## 5. 장기 목표로의 연결

이 경로의 의미는 회전·strain을 포함한 현재 물질장의 운동에너지와 물리U를 섞어 열보정하지 않고, 기하 이동 뒤 동일 장을 native target으로 보존적으로 전달하는 것이다. 실제 전 target packet/source closure가 통과하면 다음 경계는 batch streaming의 transactional native 전달이며, 일반 evolved heat를 보존·수송하는 상태/restart 계약과 CUDA/OpenMP/MPI 시간 갱신을 연결해야 한다. 현재 frozen CONTROL heat recipe를 다음 시각의 EOS 재초기화로 재사용할 수 없다.

그 이후 실제 perturbation/shock/positivity·보존 관문과 Kepler t6/8/10 및 논문 최장공전 검증이 필요하다. geometry clearance의 성공은 그 물리 관문을 대체하지 않는다. 연구 목표를 이번 국소 인증으로 축소하지 않는다.
