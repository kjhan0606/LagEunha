# Astra — compact 재격자의 곡선 외부 영역 계약 보완

2026-10-22. SHARED_CONTEXT 전체, RESULT_COMPACT_SOURCE_GEOMETRY_BUILD 및 실제 source-context import 코드를 읽었다. 현재 source identity는 보존된 원본 HEAD8e25dfe, accepted step5다. 운전자가 구현하며 이 문서는 기존 compact 계획의 **영역·공간 계약에 대한 한 번의 명시적 수정**이다. 계산·폴링·코드 변경은 하지 않았다.

## 1. 원래 명세의 결함 인정

지적이 맞다. 비영 시각의 실제 exterior가 곡선이면, 유한 straight polygon union과 identity target map으로 같은 물리 영역을 정확히 덮을 수 없다. 이전 계획의 세 조건을 동시에 요구한 것은 불완전했다. 원래 square로 돌아가거나 얇은 source strip을 버리는 것으로 해결하면 안 된다.

analytic canonical exterior trace에서 세 실제 점이 같은 직선에 없다는 인증은 그 trace가 직선이 아니라는 충분한 증거다. analytic 비직선 arc는 유한 polygon의 경계와 일치할 수 없다. 다만 그 witness가 없어도 직선이라고 결론낼 수 없고, 유한 산술의 cross가 0을 제외하는 enclosure가 필요하다. 이 확인은 이미 선택한 transaction의 domain 판정에 포함할 수 있으며 별도 작은 실행을 요청하지 않는다.

**선택 수정: 실제 curved exterior를 정확한 reference geometry로 보존하고, 내부는 compact P3로 재격자한다.** 영역 근사를 허용해 tau를 소비하는 선택은 하지 않는다. source embedding이 인증되는 조건에서 같은 물리 영역 Ω*=chi_old(Ω_ref)를 그대로 사용한다.

## 2. 실제 target partition과 reference chart

target은 전체 plane의 power cell 자체가 아니라

    T_i = Ω* ∩ PowerCell_i(q,W)

이다. 내부 power interface는 직선이며 외부는 원래 canonical curved trace의 실제 부분구간이다. convex power cell과의 교집합이라고 T_i의 convexity/star-shapedness를 자동 가정하지 않는다. disconnected component가 있으면 동일 target ID의 별도 양의 component로 남긴다. 외부 미정 gas/filler는 없다.

source/target overlap은

    D_si = chi_old(source_s) ∩ T_i

로 정의한다. 동일 overlap packet을 outgoing/incoming에 한 번씩 사용한다. 반열린 tie 규칙과 외부 canonical trace/event ID를 유지한다. source 전체와 target 전체의 union은 같은 Ω*이며, numerical enclosure/적분오차만 기존 장부에 남는다. geometry 변경에 따른 모델 오차를 새로 숨길 여지가 없다.

내부 target fan triangle은 기존 straight identity chart를 그대로 쓴다. 외부에 닿은 component만 곡선 경계 collar를 둔다. 한 curved-edge patch의 chart를 예를 들어 다음과 같이 구성한다:

    G(s,t) = A + s[(1-t)B+tC]
             + s^4 [Gamma(t)-((1-t)B+tC)],   0≤s,t≤1.

Gamma는 현재 source의 canonical exterior trace를 해당 실제 부분구간에 재매개화한 **동일 함수**다. s=1에서 원래 곡선, 나머지 두 변에서 원래 straight segment, apex에서는 안정된 zero-trace extension이다. 새 orbit 함수나 새 초기해를 넣지 않는다. s^4는 기존 안정된 edge-extension 구조를 사용할 수 있다는 뜻이며 source curve를 polynomial fitting하는 명령이 아니다.

단일 chart의 J>0 및 interior disjointness/coverage를 인증할 수 있을 때만 쓴다. chart 하나가 곡선 component를 덮지 못하면 실제 저장 boundary event/분할에 따른 유한 patch로 구성한다. 기존 예산을 넘어선 자동 refinement나 parameter scan은 없다. source proof가 unresolved이거나 필요한 chart가 구성되지 않으면 transaction을 거부한다.

## 3. 숨기면 안 되는 두 번째 문제: 회전과 geometry-completed P3

곡선 G만 두고 velocity를 reference P3로 제한하면 좌표 x가 velocity/test space에 들어가지 않을 수 있다. 그러면 정확한 rigid-rotation virtual work와 각운동 보존을 기존 증명에서 그대로 가져올 수 없다. 단순히 'geometry만 curved, 물리는 pure P3'라고 명명하고 넘어가지 않는다.

따라서 이번 수정 후보를 **geometry-completed P3**로 명시한다. 기준 straight P3 geometry를 G_P3라 하고, 전체 conforming reference patch에 조립한 고정 geometry remainder를

    h1 = G_x - (G_P3)_x,
    h2 = G_y - (G_P3)_y

로 정의한다. 내부 straight 영역에서는 h=0이다. 동일 connected component의 scalar space를

    V_h = P3_conforming + span{h1,h2}

로 둔다. 정확한 선형종속이 있으면 algebraic rank만 처리한다. 작은 singular value를 임의로 삭제하지 않는다. 추가 자유도는 **component당 최대 두 scalar geometry mode**이며, 기존 전체 IC residual DOF들을 되살리는 것이 아니다. 이것은 pure P3 명세의 공개된 최소 수정이고 polynomial order 상승은 아니다.

위 공간으로 현재 위치와 속도 양쪽을 표현한다. 시작 x는 정확 G이며, h의 위치 계수는 G를 재현하도록 결합된다. 이후 h 계수도 동일 mass/force 방정식의 동적 계수다. curved 경계 형상을 시간에 고정하거나 generator를 조향하지 않는다. constant와 현재 coordinate x,y가 공간 안에 있으므로 translation 및 x의 회전 variation이 허용된다. symmetric stress의 torque/work와 동일 mass metric으로 각운동 계약을 검사할 수 있다. h를 geometry에는 쓰고 mass/force/velocity에는 누락하면 실패다.

새 두 mode는 global support처럼 보이지만 실제 비영 support는 exterior collar다. sparse bulk와 최대2-mode/component의 작은 coupling으로 구현한다. global dense solve나 전체10.8M 공간 복원은 필요하지 않다. MPI mode owner/halo와 checkpoint에 이 geometry 함수·계수를 포함한다.

## 4. 밀도·열·mass·boundary work

Bernstein P3 positive rho/e 후보는 각 새 reference patch에서 정의하되 물리 적분에는 det(DG)를 포함한다. 따라서 degree9 mass exactness는 **straight interior triangle에만** 적용된다. curved collar에는 실제 geometry 함수와 geometry-completed basis를 포함한 양의 적분/오차 장부가 필요하다. 그 비용을 원래 compact 비용표에서 빼면 안 된다.

현재 source rho/v/e를 overlap에서 읽고 m/q/Q/Phi 및 p/L/K/U 보존 제약을 원래대로 적용한다. 새 P3 표현이 현재 source를 근사하는 오차는 이전 계획의 representation budget에 속한다. curved exterior 보존은 그 field 오차를 없애주지 않는다. entropy 생성자료는 실제 투영된 physical e와 generating jhat를 사용한다. 원래 P0나 current J만으로 열을 재초기화하지 않는다.

reference 좌표계의 density는 rho_ref=rho_physical(G)*det(DG)다. 코드에서 저장한 Bernstein 계수가 physical density인지 reference density인지 tag를 고정하고, Jacobian을 두 번 곱하거나 누락하지 않는다. geometric-pressure의 nu/delta/alpha도 이 reference measure에 맞게 일관되게 정의한다. geometric basis h가 포함된 실제 V_boundary와 그 configuration derivative를 사용한다. 기존 straight identity 전용 G/H를 그대로 넣지 않는다.

외부 조건은 변경 없이 Pext=P0다. 압력 traction과 실제 외부 work는 같은 canonical Gamma/current trace에서 평가한다:

    F_ext = -∫ Pext n N dS,
    W_ext = -∫ Pext v·n dS.

순간 재격자는 물리 시간이 흐르지 않으므로 외부 일을 더하거나 빼지 않는다. 다만 velocity projection으로 생기는 boundary normal velocity와 W_ext의 차이는 **다음 RHS의 representation error**로 평가한다. 면적/질량/E가 보존됐다고 이 차이를 면제하지 않는다. 기존 force/work/torque gate와 cold field 예산을 그대로 적용한다.

## 5. 비용·실행 범위와 끝내는 조건

추가로 허용하는 일은 같은 compact transaction 안의 exterior geometry adapter, 최대2개 geometry mode/component와 그 mass/force/열·경계 연결이다. 기존 source import와 누락 endpoint context의 CUDA pass0/1 재구성은 유효하다. 실제 gamma=1.6666669845581055와 저장 rule/생성자료를 유지한다. stage4를 endpoint로 대신 쓰지 않는다.

실행 전에 실제 curved collar patch 수, 추가 geometry 함수 저장량, 적분점·mode coupling 비용을 전체 compact 비용에 합산한다. interior P3 생성/degree9 mass/checkpoint 코드는 재사용하고, 원래 전체 원장 basis·817152/전체 초기 점·완료 trajectory를 반복하지 않는다. bulk와 collar의 비용을 구분해 곡선 경계 때문에 또 전체 원장 공간에 준하는 비용을 지불하는지 드러낸다.

동일24시간 bounded transaction 안에서 다음 중 하나로 종료한다:

1. source embedding 및 같은 Ω* coverage, geometry-completed 공간, positive field와 보존 투영, 기존 cold error/원래tau가 모두 통과하면 compact 후보를 commit하고 이전 계획의 조건부 다음 실제 native step으로 진행한다.
2. exterior chart/coordinate-mode rank·mass positivity·projection·field budget·비용 중 하나라도 실패하면 accepted step5를 그대로 보존하고 해당 구조적 실패를 보고한다. straight exterior 근사, source-strip 삭제, outside gas 생성, h mode 삭제, heat/torque correction으로 우회하지 않는다.

이번 수정은 material shear 수명을 줄이기 위한 reference 갱신의 물리 영역을 닫는 것이다. 반복 remap 수명, shock entropy, W0의 동적 Voronoi 환원, 경계가 다른 기존 periodic/논문 비교는 여전히 미완료다. t6/8/10와 올바르게 정규화한 최장 공전 목표는 유지한다. 이 부록에 따른 구현은 기존 주요 묶음 안에서 진행하며 primitive별 별도 평가 호출이나 승인 질문을 요구하지 않는다.
