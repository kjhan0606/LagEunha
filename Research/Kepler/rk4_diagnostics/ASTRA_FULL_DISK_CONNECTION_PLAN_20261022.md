# Astra — 실제 full-disk native 연결 결정

2026-10-22. SHARED_CONTEXT 전체, 원본 LagEunha HEAD `8e25dfe`, 실제 evolution/audit, Factory, persistent CUDA, assembly와 원본 `Exam/exam.c` 진입부를 확인했다. 계산·폴링·운전자 코드 수정 없이 작성한 종료 평가다.

## 1. 이번 결과의 채택 범위

`shared_evolution_v2_h200_mpi4`와 독립 audit의 실제 GPU force, 변경된 stage 계수, classical RK4 및 restart continuation을 채택한다. 12 RHS에서 전역 에너지율/U 최대 2.48e-15, 셀 3.06e-13, 각운동/2K 1.64e-17, CPU/CUDA force 5.24e-13이며 restart stage가 일치한다. 이는 열 생성자료를 보존한 **실제 동적 연산자**의 증거다. 성공한 세 step, 이전 23,339,008 전달점, atlas와 source partition을 반복할 필요 없다.

그러나 시작 상태는 24셀 finite CONTROL이다. 원래 t8 snapshot도, fresh full disk도 아니다. 고정 material reference를 변형한 것이며 새 power tessellation으로 재분할한 진화가 아니다. 매끄러운 단열 entropy 법칙은 충격 entropy 생성 구현을 대체하지 않는다. 이 결과로 t6/8/10, W0의 동적 Voronoi 환원, GIZMO 최장 공전 성공을 선언할 수 없다.

**다음 한 묶음은 원본 실행파일의 opt-in full-disk material evolution 연결과 동일 경로의 연속 진화다.** 또 다른 작은 frozen 시험이나 전달 인증 확장을 선행하지 않는다. 현재 test executable을 입력 개수만 늘려 실행하는 것은 승인하지 않는다. 아래 필수 연결을 구현한 뒤 조건부 자율 실행한다.

## 2. 실제 수정 위치와 state 계약

원본 `eunha2.c`의 MPI/설정/실행 수명과 `Exam/exam.c::exam2d_vph_rk4_int_blend`의 opt-in 분기를 이용한다. 현재 DG/moment 초기화 분기 옆에 별도 material mode를 연결하고, legacy particle scratch나 기존 RK4 분기로 떨어지지 않게 한다. C++ driver를 Python 진화기로 대체하지 않는다. 기존 옵션과 binary는 보존한다.

필수 일반화:

- `Factory.cell/triangle/dof` 크기로 모든 scratch를 할당한다. `sums[24][4]`, `projection[24]`, audit의 셀 배열, MPI count 96/24, restart의 3797/384/24/559/817152 검사는 runtime manifest의 실제 개수와 전역 ID 검사로 바꾼다. 64-bit 길이와 MPI count 분할을 사용한다.
- full disk는 새 reference source를 생성한다. 559개 과거 fragment를 모든 디스크에 붙이지 않는다. fresh 초기 view는 각 source의 identity view이며, 기존 24셀 전달 증거는 역사적 validation registry로 참조한다.
- 고정 reference mass, canonical shared edge/basis, x/v 계수, 생성 시각 U0와 jhat0, material-cell grouping, quadrature/alpha, boundary 계약, MPI owner와 accepted time/epoch가 restart의 물리 상태다. 12개 집계 moment는 고차 velocity나 열 상태의 원본이 아니다.
- 현재 CONTROL용 `boundary_pressure`의 `P0*Jgen^(1-gamma)/(gamma-1)`를 fresh 상태에 무조건 복사하지 않는다. fresh 상태의 실제 U0 함수와 생성 jhat0로 정의한 entropy를 volume과 boundary에서 동일하게 평가한다. 매 stage 생성자료를 갱신하지 않는다.
- CUDA `cudaSetDevice(0)`은 rank/GPU 배치 계약으로 바꾼다. GPU UUID, local rank 및 실제 device 소유자를 기록한다. 네 rank가 한 H200을 공유하는 경우 이를 네 GPU라고 표시하지 않는다.

수식은 이번에 통과한 그대로 유지한다:

    xdot = v,       M vdot = F_pressure(x,A) + F_gravity(x),
    jhat = J + delta_c,       nu_q = w_q jhat_q,
    A_q = U_q0 nu_q0^(gamma-1),       U_q = A_q nu_q^(1-gamma),
    F_pressure = sum_q P_q d(nu_q)/dx,       Udot_q = -P_q nudot_q.

Canonical boundary G/H, 실제 hybrid reaction, physical U, 동일 mass metric의 resolved K와 intrinsic L/I를 유지한다. spin은 별도 에너지 탱크가 아니다.

## 3. Full IC와 경계: 숨길 수 없는 선택

첫 연결 IC는 이미 승인·보존된 `core_refinement/fine_integrated`의 71,080 생성점, box24, eps=.047, 같은 연속 rho/P/중력을 사용한다. 새 graded IC나 새 refinement sweep을 추가하지 않는다. `prepare_core_refinement.py`의 실제 다각형 적분 질량과 `geometry.npz`/checkpoint SHA를 연결한다. 원본 native raw complete geometry에서 한 번 canonical 공유 꼭짓점/면 및 conforming triangle reference를 생성하고 Factory에 전달한다. 좌표 근접 병합으로 공유 DOF를 만들지 않는다.

Factory에는 patch와 동일한 P3/shared residual 구조와 선택된 reference subdivision 규칙을 적용한다. patch에서 통과한 수렴 결과를 모든 full-disk 셀의 정확도 증명으로 복사하지 않는다. 전체 fresh IC의 mass, 실제 연속장 p/K/U/E, 초기 curvature와 상태 양성을 초기화 중 집계한다. 기준을 못 맞추면 그 full-state 실패를 보존하고 이 묶음을 끝낸다. 자동 h/order 탐색은 하지 않는다.

기존 pointwise generator velocity를 셀 평균으로 취급한 legacy kinetic E와 새 연속장 적분 E는 다를 수 있다. 두 초기 장부를 나란히 기록하고, 차이를 열로 넣거나 빼지 않는다. 새 material 계수는 원래 연속 원운동장을 표현하며 이후 해석 궤도를 다시 대입하지 않는다.

**경계의 결정:** 원래 square periodic box의 방위 속도는 반대쪽 경계에서 주기적이지 않다. 정확 초기 v를 보존하면서 그 양쪽을 H1 동일 DOF로 묶을 수 없다. 따라서 첫 full-disk material 실행은 같은 square 초기 영역의 움직이는 외부 material 경계에 **주어진 일정 외압 Pext=P0**를 적용한다. 외부 traction은 -Pext n이며 에너지 장부에 그 실제 경계 일을 포함한다. 외부 위치/속도를 해석 원궤도로 고정하지 않는다. 내부 evolved P를 외부 prescribed P로 슬쩍 대체하는 자동 경계도 쓰지 않는다.

이것은 명시적인 경계 변경을 가진 full-disk 연구 경로다. 기존 periodic baseline과 동일 경계 대조라고 주장하지 않는다. 기존 t8/t10 결과는 실패의 공간·시간 이력을 비교하는 자료로 보존하며, 최종 원래 설정 및 논문 직접 비교에는 경계 일치가 별도 필요하다. 원래 periodic 경계 그대로라는 표시를 요구하는 모드에서는 진입을 거부한다. 경계 차이를 감추는 것보다 이 제한을 공개하고 전체 물리 진화를 진행하는 것이 타당하다.

## 4. 전체 디스크 비용을 감당하는 필수 구현

24셀 817152점과 longdouble 전역 mass solve/전 셀 CPU shadow를 그대로 복제하면 실행 가능한 production 경로가 아니다. 12.55초를 셀 수에 선형 외삽해 성능을 약속하지 않는다. 다음은 수학 변경 없이 이번 묶음에서 수행한다.

1. GPU가 shape/J/gradv, nu/U/P, stress와 gravity 및 셀별 reduction을 실제 수행한다. 지금처럼 CPU가 전 점 field를 먼저 계산하고 GPU 결과 전체를 매번 회수하는 경로를 production에서 제거한다. 셀 delta와 mean pressure에 필요한 두 pass는 같은 점 순서/함수로 유지한다. CPU shadow는 첫 실제 full-state RHS의 제한된 고정 셀 집합과 전체 scalar ledger에만 둔다. 이미 통과한 작은 문제를 재실행하지 않는다.
2. 고정 reference mass operator와 local preconditioner를 한 번 구축한다. DOF owner/ghost를 기준으로 sparse 또는 element matrix-free matvec을 수행하고 공유 DOF만 통신한다. 매 rank 전체 element force/전역 CSR를 복제·Allreduce하지 않는다. 기존 true residual/양의 metric 기준을 유지한다. 순수 double 전환이나 느슨한 CG tolerance는 이번 최적화가 아니다. longdouble 필요한 누적/잔차는 유지하되 고정 local blocks와 native 병렬화를 사용한다.
3. 고정 quadrature는 셀/triangle batch로 스트리밍하고 reference rule/geometry/entropy 생성자료를 압축 가능한 불변 구조로 저장한다. 매 stage 수십억 STPoint/GPOutput를 디스크나 host에 전개하지 않는다. rule을 임의로 낮추지 않는다. 기존 64x28이 전체에서 과도하더라도 이 묶음에서 낮은 차수 탐색을 시작하지 않는다. 최초 full-state count/bytes/실측이 production 비용의 판정 근거다.
4. 원본 MPI/OpenMP 및 persistent CUDA 경로를 사용한다. local triangle ownership, halo, CUDA batch, sparse solve와 audit의 시간/메모리를 각각 기록한다. 원본 실행파일/실제 mapped library/GPU identity를 남긴다.

초기화 전에 실제 triangle/DOF/CSR 또는 element-block/점 수와 checkpoint bytes를 산출한다. 장치와 host에 모두 들어간다는 가정으로 할당하지 않는다. 스트리밍 한도를 명시하고 동적 부족 시 caller/accepted checkpoint를 보존하며 실패한다. 첫 full-state accepted step에서 실제 초/step 및 필요한 step 수로 t1/t10 비용을 제시한다. **24셀 CPU reference 전달 23분/23.89GB 경로는 매 step에 호출하지 않는다.**

## 5. 실제 실행 순서와 종료 조건

이번 묶음은 full IC 한 개, 연속 trajectory 한 개다. 조건을 통과하면 동일 accepted state에서 진행하며 성공한 구간을 반복하지 않는다.

### 초기화와 첫 실제 full step

- 새 native main/mode에서 identity, 전체 공유 trace/양의 mass metric, 물리 U, J/jhat, 생성 entropy, 초기 pressure와 기존 초기 chi 기준을 확인한다. full source 적분은 fresh 상태 구성에 필요한 부분만 수행한다. 옛 24셀 전달/원점 proof를 다시 계산하지 않는다.
- 초기 t=0에서 시작한다. 저장 finite CONTROL의 offset를 full disk 시간으로 붙이지 않는다.
- RK4의 실제 stage에서 geometry→EOS→공유 force→mass solve를 사용한다. 단계 실패는 전 rank collective reject와 완전한 stage-base/ledger 복원을 요구한다. 임의 U clip, heat repair 또는 조용한 stage floor는 없다. 첫 full-state 실패를 보존하고 종료하며 이번 묶음에서 dt/threshold sweep은 하지 않는다.
- 첫 accepted full step은 다음 연속 실행의 첫 step이다. 전역/셀 energy-work와 torque ledger, fixed mass, entropy 생성자료 불변, GPU 실제 구동, restart payload와 runtime identity를 확인한다. 새 분산 구현이 미검증이므로 이 full-state 첫 RHS에서 owned/halo 합산과 global norm residual을 확인하는 것은 필요하다. 성공한 과거 3step 재시험은 필요 없다.

### 연속 경로

원래 native의 timestep 산정에 현재 material map의 음향/변형 제한을 명시적으로 추가하고 선택 원인을 기록한다. 24셀 control h=2e-4를 full mesh 전체 CFL이라고 가정하지 않는다. 시간을 줄여 실패를 숨기는 무제한 retry는 허용하지 않는다.

초기 체크포인트, t=.1/.2/1/6/8/10의 지정 출력과 가벼운 매-step 장부를 사용한다. .1/.2는 옛 열고갈 구간을 확인하는 중간 출력이며 목표 종료시간을 대신하지 않는다. 첫 배치는 **최대 wall 24시간 또는 t1 중 먼저 도달하는 시점**까지의 bounded full-disk 실행으로 한다. 예산 끝은 정상 restart로 종료한다. 첫 step의 비용 외삽이 t1조차 이 범위에 전혀 맞지 않으면 장시간 작업을 방치하지 말고 성능 자료와 accepted state에서 종료한다. 다음 성능 작업은 실제 병목 한 곳에 한정한다.

초기·수치 관문을 통과하고 기존 문서의 물리 오차 기준 및 실제 비용을 만족한 경우 동일 checkpoint에서 t6/8/10으로 연장할 수 있다. 재코딩 승인 질문은 필요 없다. 다만 단지 양수/보존이라는 이유로 장시간 물리 통과를 선언하지 않는다. 기존 기준이 없는 새 경계 항목은 미평가로 남기며 임의로 PASS를 생성하지 않는다.

필수 출력은 같은 Eulerian 영역의 rho/P, 밀도·압력 오차, 궤도 위치/속도/phase, 실제 U/K 및 경계·중력 일을 포함한 E, 전역/annular L, map/nu minima와 solve 비용이다. **A 고정으로 material entropy가 일정한 것은 모델식의 결과다.** 그 사실만으로 가짜 가열·밀도·압력·궤도 문제가 해결됐다고 주장하지 않는다. 과거 uniform/graded 비교의 같은 Eulerian bin을 재사용해 운동학적 압축 오류를 드러내야 한다.

## 6. 변형·재분할의 남은 실제 장애

현재 Lagrangian map에는 differential rotation이 계속 쌓인다. 점 J>0만으로 전역 단사성이나 공유 곡선의 비교차가 보장되지 않는다. 각 accepted stage의 local J/nu와 map condition에 더해 인접 trace 및 공간 broad-phase로 접근한 비인접 material 경계의 교차를 감시한다. 양의 영역/단사성 증명을 잃거나 solve가 정의되지 않으면 마지막 accepted state에서 중단한다. 임의 generator steering이나 weight로 모양을 되돌리지 않는다.

이 실패가 t6 이전에 발생하면 결과는 **material topology 수명의 한계**다. dt를 끝없이 줄이거나 frozen fragment labels를 새 Laguerre tessellation이라고 바꾸면 안 된다. 다음에는 실제 현재장에 대한 conservative repartition과 recipient thermal-context 호환성이 필요하다. 지금 완료된 transfer와 source-context fragment는 그 기반이지만 새 reference의 delta/entropy 재표현과 장기 비용을 아직 해결하지 않았다. 고정 mass는 material source 기준이며 현재 cell shape가 Voronoi라는 뜻은 아니다. W=0의 개별 고정질량+동적 Voronoi 환원은 계속 미완료로 표시한다.

Shock entropy 생성 역시 별도 미완료다. 이번 smooth Kepler 경로에 독립 kinetic reservoir나 점성 감쇠를 넣지 않는다. 실제 shock이 형성되어 매끄러운 entropy law의 적용 범위를 벗어나면 충격 해법을 통과한 것으로 해석하지 않는다.

## 7. 이 묶음의 완료 산출물

새 원본 native mode의 full-state manifest/restart, 실제 GPU/MPI/OMP 경로, 자원·비용표, 첫 accepted step을 포함한 연속 trajectory와 지정 물리 출력, 또는 구체적 full-state 실패 witness다. 전체 디스크 초기화/진화가 막히면 해당 구현·비용·물리 원인을 한 번에 보고한다. 작은 proof/primitive PASS를 더 모으는 것으로 대체하지 않는다.

최종 목표는 그대로 t6/8/10 및 `KEPLER_GIZMO_COMPLETION_CRITERIA_20261022.md`의 최장 공전이다. 이번 경로는 그 목표를 향한 실제 full-disk 동적 검증이며, 경계·topology·shock·재분할 미해결 사항을 숨긴 완성 주장도 아니다.
