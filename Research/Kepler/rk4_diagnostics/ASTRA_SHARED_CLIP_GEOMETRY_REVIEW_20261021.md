# Astra: clipped geometry 종료 평가와 다음 단일 native 묶음

2026-10-21. 계획·평가 Astra, 구현·실행 driver. SHARED_CONTEXT와 지정 clip/pullback/driver 및 원래 whole-target 산술 계약을 읽었다. 저장소 origin LagEunha, HEAD 8e25dfe 및 운전자 변경을 확인했다. 계산·폴링·추가 에이전트 호출은 하지 않았다.

## 판정

ALL384 coverage 완결은 중요한 전진이다. 344 제외, 4 whole, 10 empty, 26 실제 clipped를 그대로 유지한다. 4 whole source의 176 panel 인증/45056 양의 점도 유효한 부분 증거다. 그러나 전체 target은 아직 실패다. 최신 호출은 graph plan 단계에서 정상적으로 거절했으므로 clipped analytic 실패 26이라는 숫자는 실제 analytic 계산 26회 실패를 뜻하지 않는다.

현재 `axis()`는 각 고정 Duffy chart에서 모든 보수적 active face에 대해 s 또는 t 한 축만 공통으로 횡단해야 한다고 요구한다. 이것은 충분조건이지 실제 영역의 표현 가능성에 대한 필요조건이 아니다. 19개 실패만으로 fold나 물리적 특이점을 주장할 수 없다. 다음은 더 작은 동일 축 상자를 무작정 늘리는 작업이 아니라 **실제 경계의 법선으로 정한 oblique projection과 꼭짓점 event chart**를 연결하는 묶음이다.

## 1. 실제 경계에서 projection을 정한다

완료된 oriented domain/target vertex/shared edge 사건을 재사용한다. coarse displacement 상계로 active인 모든 평면을 같은 수준으로 요구하지 말고, 경계에 실제 참여하는 면과 전체영역에서 비활성으로 증명한 면을 구별한다. 비활성 제거는 샘플 부호가 아니라 기존 coverage/interval 증거에 근거한다. 최종 packet은 모든 원래 평면을 계속 만족해야 한다.

reference Cartesian 좌표에 고정 단위 방향 d와 횡방향 e를 두고 X=Xc+u e+v d로 쓴다. f_j(X)=n_j·chi(X)-b_j의 횡단 조건은

    |n_j·d| > ||n_j||_1 L ||d||_infinity

이다(현재 L은 ||F-I||_infinity 상계). 더 좁은 이미 인증된 gradient cone이 있으면 그것을 사용한다. 실제 방향의 반올림과 projection 연산 오차도 포함한다.

방향 선택은 법선이 정하는 금지 angular interval들의 여집합으로 한다. 가장 큰 인증 여유를 갖는 구간의 중앙을 택하고 동률은 고정 ID 규칙으로 결정한다. 각도를 여러 번 돌려 통과하는 방향을 찾는 scan이 아니다. 필요하면 해당 실제 domain의 인접 경계 사건으로 나눈 부분별 방향을 사용한다. 유효 방향이 없으면 실제 tangent/fold event를 분리하며, 예산 내에서 분리할 수 없으면 거절한다.

이 projection은 기저·저장 좌표를 변경하지 않는다. 적분 좌표만 바뀐다. Xc는 실제 영역의 고정 기준점이며 affine 부분과 작은 bend를 먼저 합성한다.

## 2. graph arrangement와 꼭짓점 정칙성

projection된 실제 source/target vertex, 경계와 chart 경계의 교차, projection critical event를 정렬하여 slab을 만든다. 각 slab에서 실제 경계 가지 v_-(u), v_+(u)를 인증한다. 필요한 것은 선택된 lower/upper 가지의 횡단성과 나머지 평면의 소유권 증명이다. 모든 후보 면이 같은 Duffy 축을 횡단할 필요는 없다. 양의 폭과 중복·누락 없는 coverage를 기존 경계 장부와 대조한다.

보통 slab에서는 v=v_-(u)+eta*(v_+(u)-v_-(u))를 사용한다. Jacobian은 선택한 oriented (e,d)와 양의 폭으로부터 유도한다. 음수 determinant를 abs로 감추지 않는다.

원래 source 꼭짓점을 포함하는 slab은 그대로 Cartesian analytic box로 보내지 않는다. 꼭짓점 (u0,v0)에서 u=u0+r, v=v0+r*zeta로 blow-up하고, 그 좌표를 기존 Duffy evaluator에 대입한다. 정확히 apex를 지나는 평면은

    f_j(X0+r(e+zeta*d)) = r * fbar_j(r,zeta)

로 **정확한 apex 상수 0을 증명한 뒤** 나눈다. fbar(0,zeta)=n_j·F(X0)(e+zeta*d), 따라서 n_j·F(X0)d가 분리되면 zeta 경계를 implicit graph로 추적할 수 있다. r=0에서 2차원 Duffy 역함수를 쓰지 않는다. 각 apex 방향 구간은 원래 6개 chart 안으로 제한하여 다른 extension 분모를 보존한다. 경계가 apex를 지나지 않으면 이 인수를 임의로 제거하지 않는다.

`scalar_roots()`에는 현재 exact endpoint와 identically-zero edge를 별도 처리하는 코드가 없다. 새 chart 경계가 기존 평면과 일치하면 전체 구간을 계속 이분하지 말고 canonical endpoint/coincident-branch 사건으로 소유권을 기록한다. 같은 target vertex에서 생긴 중복 사건도 근접점 merge가 아니라 기존 사건 ID로 결합한다. 이 보완은 앞으로 숨어 있는 cap 실패를 피하기 위한 계약 수정이다.

## 3. 나머지·비용의 유한 연결

새 slab/꼭짓점 chart에 기존 analytic singularity-distance atlas와 고정 Gauss16을 연결한다. implicit graph의 복소 continuation은 해당 slab의 실제 경계와 연결된 branch여야 한다. 실수 횡단성만으로 복소 나머지를 통과시키지 않는다. 원점 밀도/softening/logistic/J/extension 분모의 기존 검사를 유지한다.

현재 261 planned panel은 일부 source의 비용일 뿐이다. 모든 source의 atlas를 만든 뒤 retained 45056점을 포함한 총 비용을 기존 817152 positive-query budget과 비교한다. 계획·전체 analytic gate가 먼저 모두 통과해야 새 clipped quadrature를 한 번 실행한다. source별로 cap을 리셋해 총예산을 우회하지 않는다. 기존 완료 coverage, roots, whole integrals를 SHA 확인 후 재사용하고 원본817152 증거를 재실행하지 않는다.

## 4. 산술 계약: B<double>를 원래 CPU R80이라고 부르면 안 된다

`material_whole_target.h`는 명시적으로 CPU binary80 integration/Q128 geometry를 선언하며, `tg::evaluate<R>`, `mf::B<R>` all12 및 합산을 사용했다. 새 double field/double moment 경로는 같은 물리 함수를 근사하더라도 그 산술 계약과 동일하지 않다. 최종 출력만 R80으로 바꾸어도 복구되지 않는다.

다음 CPU reference의 새 clipped packet은 **명시적인 R80 field/coordinate/profile/all12/accumulation 경로**로 연결한다. binary64 저장 계수·고정 native cache는 그대로 입력이며, 그 오차를 R80이라고 지우지 않는다. 현재 shared basis를 R80로 평가하는 실제 operation graph와 Q128 canonical oracle의 bound를 새 실제 query에서 검사한다. 이는 원래 CPU reference 계약을 이어받는 구현이며 precision scan이 아니다. CUDA double 증거는 별도 backend 범위로 유지한다. 새 R80 경로가 기존 double graph와 bitwise 동일하다고 주장하지 않는다.

이미 인증된 4 whole-source double 기여는 그대로 사용할 수 있다. 이 경우 최종 packet은 명시적으로 mixed-arithmetic certified packet이며 homogeneous R80 replay라고 부르지 않는다. 해당 1.4441e-21 native/geometry bound와 4.3196e-23 remainder를 한 번만 전액 계상한다. 남은 26 source가 예산을 넘을지 아직 계산하지 않았으므로 지금 실패를 예단하지 않는다. 전체 gate가 실패하면 그 이유를 보존하고 종료한다. 이를 피하려고 자동으로 4 source를 더 높은 정밀도로 재실행하거나 normalizer를 바꾸지 않는다.

새 입력·cache·field·moment·sum 항을 구분해서 저장해야 bound의 지배 항을 알 수 있다. R80 연산은 이미 존재하는 binary64 input/cache 불확실성을 줄이지 못한다. 최종 R80 value/error 변환까지 더한 기존 publication gate를 그대로 유지한다.

## 5. 한 묶음의 끝과 다음 경계

구현 순서는 actual-boundary oblique/event plan → 전체 비용 → 전체 analytic certificate → 새 clipped positive quadrature 한 번 → retained whole4와 ALL384 조립 → 원래 12 tau 및 최종 직렬화 gate이다. 원본 열 생성자/계수/삼각형 fingerprint와 current map 동일성 확인, physical U 보존, caller rollback을 유지한다.

모든 면/branch/영역/오차가 완결되고 실제 질량·U가 양수일 때만 target68526을 publish한다. 어느 단계든 불완전하면 caller 9/1234를 바꾸지 않고 실패 region·방향 cone·branch/majorant·비용을 남긴다. 계획 미완료를 analytic 실패로, bound 초과를 측정된 물리 오차로 기술하지 않는다. 반복 cap/guard/order/precision/IC scan은 하지 않는다.

통과 다음은 full-target 보존 조립과 실제 가중치 법칙 연결이다. 이 묶음은 accepted RK4 step도, t6/8/10 혹은 최장 GIZMO 공전 검증도 아니다. no steering, Re=infinity, 물리 열/분해된 회전 에너지 보존, 조건부 W0 환원과 원래 native CUDA/MPI 목표를 유지한다.
