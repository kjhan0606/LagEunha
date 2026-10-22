# Astra: 실제 shared state와 native 진화 연결

2026-10-22. SHARED_CONTEXT 전체와 원본 HEAD 8e25dfe를 확인했다. 완료된 mixed/partition 장부, Factory·FrozenThermal·shared RHS·geometric pressure의 실제 인터페이스 및 KEPLER_GIZMO_COMPLETION_CRITERIA를 읽었다. 계산·폴링·코드 변경 없이 이 문서만 작성한다.

## 1. 완료 판정과 이번 선택

24 target,559 overlap,384 source partition의 동일 packet incidence와 실제23339008 queries/원래 target tau 완료를 채택한다. 기존 점·기하·atlas·성공 target를 반복하지 않는다. 이번 결과는 물질장을 적분하여 target에 나누는 실제 보존 전달이다. 고차 속도장을 갖는 새 caller state나 accepted step까지 완성한 것은 아니다.

다음 한 묶음은 **기존 source material context를 보존하는 fragment state를 commit하고, smooth Euler entropy 상태를 명시하여 native shared-DOF RK4 한 step과 실제 restart continuation을 연결**한다. Target partition을 곧바로 새 thermodynamic cell로 바꾸지 않는다. 새 reference mesh, compression 또는 recipient pressure 재투영은 이번 선택이 아니다.

이는 아직 미해결인 recipient delta 변경을 숨기지 않으면서 실제 동적 연산자를 실행하는 최소 경로다. 향후 production 재격자에서 material thermodynamic context를 바꾸려면 별도 보존 재표현이 필요하다. 이번 성공을 그 문제 해결로 선언하지 않는다.

## 2. State 계약: full field가 원본이고12 moment는 집계

별도 native state 구조를 두며 legacy particle/RK scratch를 덮지 않는다.

- `MaterialState`: 기존 Factory의 reference topology, canonical edge/basis/cache와 전역 DOF ID, `STCoefficient{x[2],v[2]}` 배열, 고정 reference density/mass, material-cell grouping, reference quadrature/alpha, gamma·중력·epoch.
- `ThermalState`: 아래의 generating geometry/nu0와 원래 physical U0 함수·값, entropy 함수 provenance. 원래 FrozenThermal은 **시작 상태 U0를 정의하는 역사 자료**로만 남긴다.
- `FragmentView`:559개의 source/domain/target 연결, 기존 boundary/event와 packet SHA, 각 target의 fragment 목록. Fragment는 source full field와 thermal-context ID를 참조한다. 별도의 독립 velocity coefficient를 만들지 않는다.
- `TimeState`: t,h,stage, accepted epoch, owner/halo 정보, stage base와 전체 rollback payload, ledger 및 restart version.

순수 partition commit에서 x/v/basis/A/material context는 bitwise 보존하고, source와 target 소유 관계만 바뀐다. Shared DOF를 target별로 복제해 독립적으로 진화시키지 않는다. Target12 moment를 사용해 B/spin/velocity를 재구성하지 않는다. 운동에너지와 intrinsic rotation은 실제 full field 및 기존 질량 metric으로 계산한다.

기존 target geometry는 **현재 완료된 remap epoch의 view**다. 이후 material step에서 그 domain을 물질적으로 운반하는 것과 새 power tessellation의 현재 overlap은 다르다. 이번 step은 전자를 보존한다. 새 power-cell overlap이라고 거짓 표시하거나 오래된 target packet을 새 시각의 값으로 사용하지 않는다. Generator/weight는 이 material step의 힘을 조향하는 장치로 쓰지 않는다.

## 3. 선택한 evolved heat: source context를 고정한 smooth entropy 함수

기존 quadrature weight를 w_q, material-cell correction을 delta_c(x)라 하고

    jhat(X,x)=J(X,x)+delta_c(x),
    delta_c=(V_boundary−sum_q w_q J_q)/sum_q w_q,
    nu_q=w_q*jhat_q.

생성 시각의 실제 reference-volume heat u0(X)와 jhat0(X)를 저장해

    a(X)=u0(X)*jhat0(X)^(gamma−1),
    u(X,t)=a(X)*jhat(X,t)^(1−gamma),
    U_q=w_q*u(X_q,t), P_q=(gamma−1)U_q/nu_q

로 정의한다. 동등한 nodal 표현은 A_q=U_q0*nu_q0^(gamma−1)이다. 이 두 A의 단위와 weight 거듭제곱을 혼동하지 않는다.

이것은 smooth inviscid 단열식 Udot=−P nudot의 정확한 상태 표현이다. **현재 J로 generating state를 매 stage 갱신하는 열 초기화가 아니다.** 저장 CONTROL의 기존 u0 함수를 생성자료로 단 한 번 결합하고, jhat0에는 실제 geometric-pressure correction을 포함한다. nu0=J라고 임의 대체하지 않는다. 모든 generating coeff와 cell grouping/quad/alpha를 deep-copy 또는 immutable hash로 보존한다. 첫 RHS 전에 원래 저장 U와 이 생성식의 일치를 확인한다.

Pure target repartition에서는 source material-cell grouping과 delta 정의를 그대로 참조한다. Recipient별 delta를 새로 계산해 같은 a에 넣으면 U/P가 변하므로 이번 단계에서 금지한다. 이렇게 context를 보존하면 fragment 내부와 공유 trace 양쪽에서 동일한 기존 thermodynamics를 평가할 수 있다.

이 표현은 **smooth Kepler 단계에 한정**한다. Shock에서 entropy production을 금지하는 새로운 PDE로 채택하지 않는다. 나중에 실제 충격을 진화할 때는 conservative total-energy/compatible numerical stress로 생긴 열에 따라 entropy 상태가 변해야 한다. 임의 damping이나 A 고정으로 충격 검증을 통과했다고 주장할 수 없다. Shock operator가 없는 현재 묶음에서는 discontinuous shock 진화를 실행하지 않는다.

## 4. Native RHS와 에너지·각운동 계약

고정 material mass matrix M에 대해

    xdot=v,
    M vdot = sum_q P_q * d(nu_q)/dx + f_gravity + f_external,
    Udot_q=−P_q*d(nu_q)/dt.

를 기존 geometric-pressure 조립과 동일하게 사용한다. 부호는 현재 검증된 force convention을 따른다. 이때 d(sum U_q)/dx=−sum P_q dnu_q/dx이며 kinetic-pressure work와 상쇄된다. 기존 mean_deltaP 및 canonical boundary G/H 항을 함께 사용하고 volume residual만 떼어 쓰지 않는다. Entropy식으로 U를 평가하더라도 이 Udot는 매 stage 독립 work ledger로 계산한다.

`GPInput.U`에는 해당 stage의 evolved U_q를 넣는다. 헤더의 'immutable OLD U' 설명을 새 mode에서는 stage input으로 명확히 구분해야 한다. `deltaJ,deltaJdot,mean_deltaP`는 같은 stage x/v/heat에서 얻는다. Canonical boundary volume/gradient 및 fixed mass CSR/preconditioner를 재사용한다. 상태를 섞은 cache는 epoch 검사로 거부한다.

주의할 실제 연결 지점: `shared_trace_rhs::external`의 P0*J^(−gamma)는 현재 CONTROL의 특수 식이다. 새 RHS에서는 같은 immutable a와 해당 material-cell jhat로 외부 경계 압력을 평가하고 기존 Piola traction/work/torque 경로에 넣는다. 외부가 주는 물리 BC는 현재 제한 patch의 기존 계약으로 명시한다.24-cell patch의 boundary 조건으로 전체 disk Kepler 안정성을 주장하지 않는다.

CUDA는 기존 gp_cuda_evaluate의 shape/field/stress/thermal-rate를 사용하고 CPU/OpenMP는 canonical boundary·sparse M solve·전역 조립을 담당한다. 초기에는 실제 작동하는 이 혼합 경로를 우선 연결한다. 전역 dense Schur를 만들지 않는다. Zero force이면 solve가0 acceleration을 반환하도록 하며 현재 solve의 norm>0 전제를 물리적 실패로 혼동하지 않는다.

Intrinsic L/I/Krot는 기존 material-cell centroid/shape에서 진단한다. 독립 spin energy를 추가하지 않는다. 외부 torque와 중력 항을 포함해 global L을 닫고, cell-local ledger는 이미 검증한 hybrid reaction traction을 포함한다. CG projection residual을 물리 pressure traction으로 이름만 바꾸지 않는다.

## 5. 이번 실행 범위와 commit/실패 조건

### 상태 import 및 relabel

완료된 packet을 읽어24 target view와 source field를 같은 epoch로 연결한다. 새로운 물리 적분은 없다. m/U positivity, 같은 outgoing/incoming 참조, 원래 state/열 provenance 및 source-field 동일성을 확인한 후 이 **full-field fragment state**를 commit한다.

### 실제 RK4와 restart

저장된 같은 finite CONTROL state를 시작점으로 기존 검증 h=2e−4에서 native classical RK4를 한 번 실제로 수행한다. 이것은 새 형상에 대한 field/J/boundary RHS 평가이지 새 remap atlas 검색이 아니다. 기존 reference mass matrix·basis·quadrature를 재사용한다. 각 RK stage의 x/v에서 jhat/U/P를 계산하고 실제 CUDA force를 조립해 sparse solve한다. P0 원운동 forcing이나 잔차 차감은 없다.

Stage 및 final에서 J,nu,U 양성/finite, shared DOF/trace, mass 고정, stress-work·중력·boundary ledger, solve 잔차를 검사한다. 기존 물리/수치 기준을 그대로 적용하고 유한시간 energy defect는 RK 오차와 반이산 work residual을 분리한다. 고정 A를 사용한 U 변화와 실제 Udot ledger의 유한 increment도 기록한다. Energy를 강제로 맞추는 보정은 하지 않는다.

첫 accepted state를 durable restart로 저장하고 같은 h의 다음 한 step을 uninterrupted 경로와 실제 restart 경로에서 각각 수행한다. 총3 actual steps/12 RHS가 이번 진화 비교의 최대 범위다. Restart는 full coefficients뿐 아니라 generating a/nu0, material grouping, basis·mass·quad·alpha, fragments, ownership 및 ledger를 복원한다. 동일 backend 결정적 경로의 상태/entropy/provenance 일치를 요구하고 CPU/CUDA shadow는 기존 연산별 bound로 비교한다. 실패 시 마지막 accepted state와 stage 입력을 보존하고 자동 dt scan 또는 열수정을 하지 않는다.

이 구간의 시작은 비정상 finite CONTROL이므로 초기 정상 Kepler의 곡률/오차 조건을 잘못 적용하지 않는다. 결과는 실제 연산자 연결·restart 성공이며 fresh Kepler 궤도 성공이 아니다.

## 6. Source-before closure의 정확한 해석

완료된 topology/disjointness와 동일 integrand의559 overlap이 원래 source를 정확히 partition하므로, 각 source partition sum은 같은 whole-source 연속 적분을 오차 상계 안에서 감싼다. 따라서 **새로운 거대한 독립 whole-source 재적분은 필요 없다.** Incidence0 자체만으로 before/after 검증을 주장하지 않고, 이 coverage 및 적분 오차 전제를 함께 붙인다.

기존 저장 source mass·kinetic M·physical U·L 장부가 있으면 그 enclosure와 partition enclosure를 read-only로 대조한다. 이산 material mass/열 quadrature와 연속 적분은 다른 근사이므로 둘을 exact 동일로 만들거나 source 값을 target 합계로 reset하지 않는다. 향후 재표현이 없다는 이번 relabel에서 authoritative physical state가 그대로이므로 그 state-before/after 보존은 byte identity로 따로 확인된다. Source-by-source 독립 quadrature 비교가 미실행이라는 사실도 계속 명시한다.

## 7. 비용과 실제 장기 경로

23분/23.89GB의 binary80 audit 적분은 production step 전략이 아니다. 다음 step은 기존 고정 material quadrature, persistent GPU buffers, 재사용 mass matrix와 native boundary 조립만 수행한다. Point별 audit dump는 기본 off이며 stage/source reduced ledger와 실패 witness를 저장한다. 같은 frozen geometry나 완료23339008점을 다시 계산하지 않는다. 실제 RHS/solve/communication/IO 시간과 GPU kernel 비중을 보고한다.

이번 성공 후 남는 주요 경계는 full-disk material state/외부 BC, shock entropy와 compatible dissipation, 필요 시 fragment 증가·recipient thermodynamic-context 재표현, 실제 production CUDA/MPI 분산 갱신이다. 이번 source-context 유지 step을 새 power mesh의 완전 재격자 time solver라고 부르지 않는다. W0 fixed mass 환원도 기존 승인 coupling law의 조건으로 유지한다.

다음 실제 fresh full-disk 연결에서 같은 continuum IC와 일반 evolved heat를 사용해 원래 t6/8/10의 density·가짜가열·궤도·E/L를 확인한다. RK stage 또는 U positivity만으로 장기 성공을 선언하지 않는다. 이후 KEPLER_GIZMO_COMPLETION_CRITERIA의 시간/궤도 normalization과 그림120/600 및 약1000공전 목표를 따라 확장한다. 이번 문서의 조건을 만족하는 routine import/빌드/3-step 및 restart 실행에는 추가 승인 질문이 필요 없다. 물리 계약 실패는 한 묶음 결과로 보존해 보고한다.
