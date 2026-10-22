# Astra — full-disk evolved RHS와 실행 비용 연결

2026-10-22. SHARED_CONTEXT 전체와 원본 HEAD `8e25dfe`, 실제 full initial/curvature JSON, owned state, cached edge basis, distributed mass, 초기화 driver/ABI를 확인했다. 계산·폴링·코드 변경 없이 작성한다. 실행 중 독립 initial audit의 완료 여부는 여기서 추정하지 않는다.

## 1. 종료 평가와 선택

실제 71,080셀, 1,138,432 triangle, 10,816,441 scalar DOF의 분산 초기 상태를 채택한다. 2,040,256,512개의 양의 점을 스트리밍한 초기 force/mass와 원래 허용치의 binary80 true residual 9.79e-17, 초기 max chi=1.62962e-4 및 셀 평균 최대4.67988e-6은 기존 관문을 통과했다. 초기 전체장을 다시 만들거나 곡률·질량 적분을 반복할 이유가 없다. cached64 Gauss는 같은 basis의 재사용이며 새 근사공간을 도입하지 않는다.

아직 local hybrid rate, 일반 evolved force, accepted full step과 full restart는 없다. 초기 일정압력에서는 volume pressure residual이 사라지므로 초기 성공만으로 변화한 압력의 force/work 구현이 검증되지는 않는다. 현재 경로는 Pext=P0의 움직이는 외부 material 경계 연구 모드이며 원래 periodic baseline과 동일 설정이 아니다.

**다음 한 묶음은 저장된 full initial state를 직접 진화 상태로 가져와, 일반 GPU evolved RHS와 분산 mass 가속을 연결하고 첫 실제 full RK4 step부터 t1/24시간의 동일 연속 경로를 진행하는 것이다.** 질량·초기화 시험을 더 늘리지 않는다. 원본 실행파일의 opt-in main/MPI 수명을 계속 사용한다.

## 2. 열 상태를 명시적으로 버전 고정

현재 thermal v1의 `law=0`은 의미를 검증하는 import adapter를 둔다. 이를 무조건 current-J recipe로 해석하지 않는다. 새 restart에는 다음 law와 생성자료 SHA를 명시한다.

    law = FRESH_CONSTANT_PRESSURE_GENERATING_CONTEXT
    j0(X) = J(X; generating_coeff) + delta0(source)
    Uq0 = P0*wq*j0/(gamma-1)
    Aq = Uq0*(wq*j0)^(gamma-1)
    j(X,t) = J(X; current_coeff) + delta(source,t)
    Uq(t) = Uq0*(j/j0)^(1-gamma)
    Pq(t) = P0*(j0/j)^gamma.

동일 source grouping, 고정 reference quadrature/alpha, generating coefficient와 delta0를 보존한다. 생성 deltaDot은 provenance이며 entropy 생성식에 들어가는 새 시간상태가 아니다. j0를 1로 반올림하거나 매 stage 새 J로 바꾸지 않는다. Uq와 reference-volume heat density의 weight 단위를 구별한다. 생성 시점의 U 합은 저장 initial physical U와 기존 오차 관문으로 일치해야 한다. 생성식은 boundary에서도 같은 함수다.

같은 법칙의 안정적인 log1p/expm1 평가는 허용하지만 equation을 바꾸지 않는다. 압력·열 increment를 작은 차이로 계산해 큰 kinetic E의 차감으로 열을 얻지 않는다. fixed A는 smooth Euler 단열식의 상태 표현이며 shock entropy 생성 구현은 여전히 별도 미완료다.

## 3. 일반 evolved GPU RHS: 실제 필요한 세 pass

`material_full_initial_cuda.cu`의 현재 공유 basis/point 생성과 DD reduction을 재사용한다. 원래 rule, source density, coefficients와 point ownership을 유지한다. mass 적분과 conditioning은 매 stage 하지 않는다.

1. **기하 pass:** owned current coefficient를 pull한 뒤 GPU에서 J, J div(v)를 계산하고 source별 QJ/QJdot를 DD로 축약한다. canonical current boundary에서 Vb/Vbdot와 G/H를 구해 delta와 deltaDot을 만든다. 외부 Pext traction은 외부면에만, 공유면은 canonical ID와 한 번의 평가·양쪽 scatter를 따른다.
2. **EOS pass:** immutable j0와 current j로 U/P를 평가하고 기존 geometric-pressure projector에 필요한 mean deltaP를 source별 축약한다. positivity/finite 실패는 source/triangle/point/epoch를 보존하고 collective fail한다.
3. **force/work pass:** 동일 nu의 미분에 대한 residual stress와 mean-pressure boundary 항, gravity를 계산한다. 기존 geometric-pressure 식과 adjoint Udot를 함께 사용한다. volume residual만 바꾸고 boundary work를 옛 초기식으로 두는 조합은 금지한다. element force, work, gravity torque와 physical U/K/L를 축약해 회수한다.

GPU 메모리에 20억 점 전체를 보관하지 않는다. 고정4096-triangle batch 또는 메모리상한 내 동일 분할을 사용한다. 필요한 multi-pass 점 계산은 **새 시간상태의 RHS**이지 옛 초기 증거의 반복이 아니다. 원래 초기 raw mass와 초기 curvature를 재평가하는 것은 불필요하다. 전 점 CPU shadow도 금지한다. 고정된 기존20triangle의 첫 변경 stage shadow, 모든 셀 algebraic ledger와 실패 witness가 충분하다.

owner/ghost 위치·속도는 각 RK stage의 동일 epoch에서 먼저 일치시킨다. source cell 전체가 한 owner에 있는 현재 grouping을 유지하고, 공유 DOF만 halo한다. 공유 edge가 두 rank에서 평가되더라도 canonical trace/방향/epoch를 같게 하며 합산은 한 소유 규칙으로 한다. GPU UUID와 rank mapping을 남긴다. 한 GPU에 네 rank인 이번 실행을 네 GPU 성능으로 기술하지 않는다.

## 4. 질량 병목: DD GPU inner solve + binary80 실제 방정식

측정72.268초는 **local factor 생성과 두 CG solve의 합**이다. 이것만으로 매 RHS가72초라고 단정할 수 없다. 다음 실행에서 factor-once, apply, precondition, halo, dot/reduction, true residual을 분리한다. 현 CPU `Mass::apply/precondition`은 매 호출 vector와 element solved scratch를 할당하므로 재사용 버퍼로 바꾸는 것도 우선적인 의미 보존 최적화다.

선택한 GPU 경로:

- 저장 conditioned binary80 mass block과 local Cholesky/partition을 한 번 읽어, 각각 hi+lo double 쌍으로 변환한다. 변환의 복원 잔차를 기록한다. 원래 binary80 block이 authoritative M이다. 순수 double M으로 바꾸지 않는다.
- GPU에서 작은 element block matvec 및 두 번의 triangular solve를 DD 연산으로 수행한다. FMA 기반 error-free product/sum, DD 나눗셈, 정해진 row/support 순서를 사용한다. `partition * sum(local inverse * partition)`의 양쪽 partition을 모두 보존한다. RHS/solution/halo 값도 hi/lo를 함께 이동한다. unordered atomic double scatter로 대체하지 않는다.
- local force와 mass의 ghost contribution을 한 번 owner로 합산한다. global norm/dot는 owned 성분만 포함한다. ghost를 두 번 센 작은 residual은 허용하지 않는다.
- **outer defect correction은 binary80**로 `r=b-Mx`를 직접 계산한다. GPU DD PCG는 이 r의 correction을 푸는 inner solver로 사용한다. DD의 recurrence residual만으로 완료하지 않는다. GPU 경로를 기존 exact-operator CG의 수학적으로 동일한 반복이라고 표시하지 않는다.
- 종료는 기존 `||b-Mx||/||b|| < 1e-16` 그대로다. zero RHS는 정확 zero solution으로 처리한다. correction이 줄지 않거나 positive metric/denominator가 깨지면 종료 실패다. 허용치 완화나 무제한 residual restart는 없다. 최대로 허용할 outer correction 수는 구현 manifest에 작은 유한 수(3회)로 고정한다.

첫 실제 변경-stage RHS 하나에서 기존 binary80 CPU solve와 GPU correction의 **같은 RHS** 결과를 대조한다. 기존 방정식 residual 및 work/torque 관문으로 판정하고 acceleration bitmatch는 요구하지 않는다. 이것은 초기화 재실행도 trajectory 두 개도 아니다. 통과 후 stage마다 binary80 true residual은 유지하되 CPU 독립 solve는 반복하지 않는다. GPU 경로가 실패하면 원인과 비용을 보존한다. 기존 검증 CPU 경로를 명시적으로 선택해 첫 full step을 완성할 수는 있으나 이를 GPU mass 가속 성공으로 보고하지 않는다.

## 5. Local hybrid ledger와 RK4 연결

저장 초기 acceleration, fixed M와 초기 element force로 첫 local hybrid ledger를 만들 수 있다. 초기20억 점을 다시 평가할 필요 없다. 변경 stage에서는 GPU가 출력하는 각 element/source force를 사용한다.

    r_c = M_c a - f_c,
    work_c = v^T r_c,
    torque_c = x cross r_c.

실제 canonical trace dual/reaction과 위 집계를 연결해 셀 E/L ledger를 검사한다. 잔차를 열이나 torque에 더해 맞추는 correction은 없다. nonzero reaction 자체를 실패나 물리 가열이라고 부르지 않는다. 초기 결과의 global rate PASS를 local closure PASS로 대신하지 않는다.

fresh t=0, 저장 full coefficients/mass/entropy로 RK4를 시작한다. 저장 초기 acceleration은 동일 physics/BC/owner/basis SHA와 방정식 검사가 일치하면 첫 k1로 재사용한다. 이때 새 local ledger와 일반 RHS 구현의 초기 대수적 reduction이 저장 F와 같음을 확인한다. 새로운 pressure state의 k2/k3/k4는 실제 계산한다. 이미 계산한 k1을 단지 형식 때문에 재적분하지 않는다.

accepted state는 네 stage 모두 positivity/geometry/true residual/work/torque 관문을 통과하고, collective commit 후에만 바꾼다. 실패하면 모든 rank의 coefficient, time, dt controller, stage/accepted ledger 및 cache epoch를 복원하고 정상 checkpoint+실패 witness에서 종료한다. 이번 묶음은 dt/threshold scan이 아니다.

restart는 x/v뿐 아니라 생성 entropy 함수·j0, fixed M/conditioned basis, ownership/global IDs, source grouping/BC, quadrature, accepted ledger/time/dt를 포함한다. immutable payload는 SHA 참조하여 중복 저장을 피하되 누락 시 restart를 거부한다. 원래 source masses와 U0를 legacy particle reader가 재생성하지 못하게 한다. 첫 accepted step을 저장한 뒤 다음 한 실제 step에서 uninterrupted와 restart branch를 대조한다. 이는 **새 full-disk restart 구현의 한 번 검증**이며 과거 작은3step 반복이 아니다. 이후 성공 구간을 재실행하지 않고 한 경로를 계속한다.

## 6. 실행·비용·물리 종료 조건

독립 full initial 감사가 PASS한 뒤 위 native mode를 실행한다. 새 IC/Factory/mass/곡률/atlas 구축은 없다. 첫 full RK4 step의 실제 pass별 GPU/CPU/통신 시간, peak RAM/VRAM, CG/outer 수, 출력 크기를 기록한다. 초기화108초나 curvature19초를 그대로 RHS 시간으로 대입하지 않는다.

첫 batch는 이전 결정대로 t1 또는 wall24시간까지다. 초기 timestep은 실제 material 음향/변형 제한과 기존 native controller로 결정하며 큰 h로 예산을 맞추지 않는다. .1/.2/1은 중간 물리 출력, 이후6/8/10이 목표다. checkpoint는 초기·지정 시각·wall 종료/실패로 제한하고 perpoint dump는 없다. 첫 accepted step의 실측으로 t1/t10 walltime을 제시한다. 초기 비용 예측이 틀렸더라도24시간 끝에는 정상 accepted state에서 종료하고, 수십일 작업을 자동 제출하지 않는다.

매 stage 양의 J/jhat/U, source mass 고정, sparse solve, canonical 공유 trace와 E/L ledger를 검사한다. accepted 출력에는 Eulerian rho/P, annular 압력·밀도·궤도/phase, actual U/K/L, 외부 pressure work 및 gravity를 포함한 E, mesh condition과 비용을 남긴다. A 고정에 따른 entropy 일정만을 정확도 증거로 쓰지 않는다. 원래 .1 부근 열고갈 경로와 t6–10 실패를 압력·밀도·궤도 수준에서 비교해야 한다.

점 J 양성만으로 전역 비접힘을 선언하지 않는다. material 경계의 접근/교차와 local map 품질을 감시하고 미보증·folding·비정상 condition 시 마지막 accepted state에서 멈춘다. fragment labels가 새 power tessellation을 구현한 것처럼 표시하지 않는다. 재분할 시 source/recipient thermal context, W0 동적 Voronoi 환원, 충격 entropy 생성은 계속 미완료다. 이번 trajectory가 잘 가더라도 이 항목들을 자동 PASS로 바꾸지 않는다.

## 7. 한 묶음의 종료 보고

운전자는 위 구현·조건부 실행을 재승인 질문 없이 진행할 수 있다. 보고는 첫 accepted full RK4+full restart, 연속 t1/24h 결과, 명시적 full-state 실패 또는 실제 지배 비용을 한 번에 묶는다. primitive별 추가 review는 요청하지 않는다. 성공한 full initial과 curvature를 재실행하거나 더 좁은 인증문제를 목표로 삼지 않는다. 목표는 원래 native CUDA/MPI/OpenMP에서 t6/8/10와 최장 공전의 실제 물리 정확도이며, 이번 문서는 그 진화 경로를 연결하는 계획이다.
