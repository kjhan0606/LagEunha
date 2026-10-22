# Astra regular-q 전달 묶음 평가

2026-10-21. 계획·평가 Astra, 코딩·실행 driver. SHARED_CONTEXT와 현재 duffy/pullback/clip, reference80 integral, whole-call 정책 및 witness v5 계약을 읽었다. 원본 저장소 HEAD 8e25dfe를 확인했다. 작업 실행·폴링은 하지 않았다. 읽은 시점에 v8 contract.json이 없어 v8 결과를 판정하지 않는다.

## 1. 이번 증거의 판정

실패했던 29개만 재검사해 29개 통과, 실제 same-function 값/AD2 비교 5916개 통과는 유효한 연결 증거다. 기존 842개를 되풀이하지 않은 범위도 적절하다. 이는 전체871 나머지/새26 적분/최종 target budget 통과와 별개다. v8의 기존842 나머지 재산출은 이전 CSV에 R12가 없어서 필요한 한 번의 보완이며, 이번에는 typed panel과 나머지를 보존한다.

읽은 핵심 algebra에서 결과를 뒤집는 결함은 발견하지 못했다. 아래 조건을 전제로 전체 호출을 진행할 수 있다. 상태·기저·물리 오차 기준은 변경하지 않는다.

## 2. 수학적 확인

- q=s*t이면 X=A+s(B-A)+q(C-B)이고 det(dX/d(s,q))=detChart다. 따라서 strip의 물리 reference measure는 detChart*(q_hi-q_lo)*ds*deta이다. 현재 regular_angular_strip_point의 Jacobian은 이 식과 맞는다. s를 또 곱하거나 t 폭을 중복 적용하면 안 된다.
- 이 경로는 실제 s=0을 제외하는 regular 경로다. 꼭짓점에는 기존 별도 Duffy event 계약을 유지한다. q 좌표 도입이 apex의 Cartesian holomorphic 역함수를 만들어 주는 것은 아니다.
- eta의 nominal 계수를 먼저 합성하고 eta radius 항을 따로 더하는 AffineAngular는 동일 affine 조합의 상계다. 계수 간 cancellation을 유지하면서 eta 불확실성을 버리지 않는다.
- Mobius 보정은 N/D에 대해 (N+Eb*delta)/(D+Ed*delta)=N/D+delta*(Eb*D-Ed*N)/(D*(D+Ed*delta))이다. 코드의 선형 barycentric determinant와 sequential division이 이 항등식에 맞는다. 두 분모의 분리는 모두 필요하다.
- geometric_q_prior는 Hessian 계산 전에 실제 분모·고정 binomial domain으로 제한한다. 더 작은 prior를 선택하는 것은 정확도 기준을 낮추는 일이 아니다. 다만 최종 correction ball이 그 prior 안에 들어가야 하며 prior 밖 Hessian 사용은 금지다.
- radial fallback에서 기존 prior 자체가 q<1 및 residual/prior+q<1을 만족하면 그 prior를 불변 수축영역으로 쓰는 것은 정당하다. 2*residual/gap은 충분한 구성법 하나일 뿐 필요조건이 아니다.

복소 branch는 real center의 원래 graph root와 이어져야 한다. real-root interval이 최종 correction ball에 포함되고, predictor 전체에서 같은 implicit 함수의 수축/불변성과 모든 해석적 분모가 인증되어야 한다. 단순 값/AD2 표본 일치가 이 branch 조건을 대신하지 않는다. 현재 코드의 real_radius, R<=prior, contraction/invariance 기록을 전체871에서 유지한다.

## 3. 호출 조립·산술 계약

v7에서 미리 고정한 842 Duffy/27 q/2 radial 정책은 허용한다. 새 결과를 보고 더 작은 나머지를 골라 쓰는 방식은 금지다. 정책의 source/panel/방향/branch ID와 모든 경계 중심·radius를 typed packet/hash로 묶는다. 현재 CSV 대조가 nominal endpoint만 비교하는 부분은 전체 구간의 bitwise 동일성 주장 근거가 아니다. 실제 나머지를 새 구간에서 다시 인증한다면 안전성은 별개로 유지되지만, provenance에는 그 차이를 정확히 구분해야 한다.

새26 source는 R80 field/profile/all12 및 pairwise tree, retained4 source240/244/248/252는 기존 double 증거다. 최종 결과를 mixed-arithmetic으로 명시한다. pairwise tree는 실제 R80 덧셈과 각 연산 오차를 포함하므로 허용한다. finish_reference80의 마지막 occupied-level 합도 forward radius에 포함된다. sum_local_rounding 진단만으로 전체 tree 오차를 대신하지 않는다.

retained4 group은 한 번만 초기 합에 넣고 group source rows는 다시 더하지 않는다. 새 source는 한 MPI owner에서 완결한 뒤 고정 source 순서로 합산한다. 원본 열 generating fingerprint와 current CONTROL 동일성, delta/delta_dot, all12의 공통 양의 measure를 유지한다. 이 열 recipe는 일반 evolved thermal transport가 아니다.

source별 point bound/analytic remainder/sum 오류 및 최종 R80 conversion을 합친 뒤 원래 12 tau를 검사한다. native R80 oracle 비교는 실제 새 query마다 수행하고 기존817152 증거를 재생하지 않는다. point oracle의 관측 차이는 absolute bound를 대신하지 않는다. 최종 caller는 모든 source/원본SHA/양성/전체오차 관문이 통과한 경우만 commit한다.

## 4. v8 이후의 단일 결정 경계

### v8가 아직 없거나 실행 중인 경우

현재 승인된 한 호출을 완료하고 frozen871 panel/R12, source별12기여, 양의packet, typed manifest, 최종오류장부를 저장한다. 별도 새 계산을 추가하지 않는다. reference80 diagnostic 통과는 선행 구현검사이며 whole-call 통과로 기록하지 않는다.

### 전체871와 적분은 완료했지만 whole tau가 실패하는 경우

caller를 보존하고 종료한다. 이미 저장한 source기여와 panel별R12만으로 지배항을 분리한다: retained double group, 새field/cache/profile, 실제합산/직렬화, analytic remainder. 열 또는 normalization을 바꾸거나 자동 higher-precision/차수/guard/cap scan을 실행하지 않는다. 새 전체 적분을 반복하지 않는다. 표본오차가 작다는 이유로 bound를 삭제하지 않는다.

### v8가 실제 publish까지 통과하는 경우

다음 substantive 묶음은 **동일 frozen state의 전체 target 보존 전달 조립**이다. target68526의 완성 packet/증거를 그대로 재사용한다. 기존24-cell target registry와 실제 외부 유출 domain을 명시하여 source를 빠짐없이 분할한다. 각 canonical overlap은 한 번 소유하고 source outgoing/target incoming에 같은 packet을 반대 장부 역할로 연결한다. 전체 영역이 닫히지 않는 부분은 명시적 external flux/fragment로 보존하며 filler를 만들지 않는다.

그 묶음은 먼저 all-target incidence/예상node·저장비용/target-U 기존계약을 고정하고 예산을 확인한 뒤 native 조립을 한 번 수행한다. 각 target의 cold오류와 source별 m,p,L,K,U,potential 등12합 보존을 모두 요구한다. global합만으로 개별 target을 승인하지 않는다. 압축·weight Newton·새 물리시간 적분은 아직 섞지 않는다. 목적은 하나의 임의 target을 계속 개선하는 단계에서 실제 remap operator로 전진하는 것이다.

일반 evolved heat state와 restart/수송 정의는 그 후 실제 time driver 연결의 필수 전제다. 현재 CONTROL recipe를 시간마다 재평가해서 대체하지 않는다. resolved rotation은 같은 x,v,m의 모멘트로 보존하고 별도 에너지 reservoir를 도입하지 않는다.

## 5. 주장 범위

현재 확보한 것은 새 좌표·상계의 실패29개 해결이다. 이 문서는 v8 성공을 선취하지 않는다. 최종 target pass도 accepted RK4 step이나 production GPU remap, t6/8/10 또는 GIZMO 최장공전 성공이 아니다. native CUDA/OpenMP/MPI, no steering, Re=infinity, 원래 물리U 및 조건부 W0 환원 목표를 유지한다.
