# 원래 셀 내부 가속도 소거 — 구조·행렬 구축 완료

2026-10-22. 연구 ROOT: `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics`.

## 목적과 인정 범위

목적은 역사적인 t≈6–10의 가짜 압축 가열·궤도 붕괴를 해결하고, 궁극적으로 GIZMO 최장 공전 검증까지 진행하는 것이다. 이번 결과는 **그 검증을 감당할 계산 경로의 준비**이며 전단 수명 개선이나 물리해 통과를 의미하지 않는다.

`ASTRA_EXACT_CELL_CONDENSATION_DECISION_20261022.md`의 단일 exact cell condensation 경로를 구현했다. 원래 10,816,441 scalar DOF의 x/v, 원래 기저, 고정 M, 기존 selected quadrature와 열 생성 상태를 유지한다. 실패한 compact 공간의 T, 두 residual 모드, 초기 가속도는 사용하지 않는다. 새 초기조건·물리 적분·초기 curvature 재실행은 없다.

인정된 물리 상태는 여전히 원래 **step5, t=0.00011263377801805954**이다. 이번 묶음의 accepted step은 0이다. 원래 경로:

`material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r[0-3].bin`

## 1. 실제 원래 support 분류와 CSR 패턴

`material_cell_condensation_structure_v1_mpi4/run/structure_result.json` 및 네 `cell_partition_r*.bin`.

- 전체 71,080 source cell, 1,138,432 triangle을 읽었다. MPI4, 원래 conditioned basis metadata를 사용했다.
- 물리 source cell 안의 edge incidence 1을 외부 trace로 분류하고, 해당 네 P3 edge DOF와 edge residual DOF를 Gamma에 둔다. 내부 삼각형/edge/bubble의 zero trace는 원래 기저의 해석적 항등식으로 판정한다. 작은 수치값으로 분류하지 않는다.
- 모든 I DOF가 한 source에만 지지되며 그 source owner에 속함을 확인했다. 실제 원래 mass의 triangle 지원으로 다른 cell의 내부 DOF 사이 coupling이 없음을 확인했다.
- 전역 Gamma **1,355,777**, 내부 I **9,460,664**, 합계 **10,816,441**. Gamma≤2.5M 관문 통과.
- 최대 source 내부 차원 **463**, 최대 trace 차원 **140**.
- 실제 조립용 partial CSR 패턴 **96,946,492 entries**. MPI ghost row가 포함된 부분 행렬 합계이며 전역 고유 nnz라고 부르지 않는다.
- native exit 0; 구조 구축 native 8.3666초. 분류/CSR/selected point offset 파일을 실제 저장했다.

## 2. 저장된 원래 M의 실제 대수 구축

`material_cell_condensation_factor_v1_mpi4/run/factor_result.json` 및 네 factor/네 skeleton matrix 파일.

원래 binary80 triangle mass를 source별로 조립하여 `M_II=L L^T`, `Y=L^-1 M_IG`, `S=M_GG-Y^T Y`를 구성했다. 큰 역행렬은 만들지 않는다. L과 원래 M_IG를 저장하고 S를 source 순서대로 partial CSR에 합산했다. x_I/v_I를 제거하지 않으며, 매 RHS에서 복구할 a_I만 소거한다.

- 전체 71,080 source의 M_II 및 S Cholesky가 양의 pivot으로 통과했다. jitter, clipping, pivot repair 없음.
- 원래 M_II에 대한 최대 `||LL^T-M_II||_F / ||M_II||_F` = **3.722651914420799e-19**.
- 최소 II pivot / 해당 source 최대 II diagonal = **0.0099589696036080**.
- packed L **635,433,712 entries**, M_IG **379,124,224 entries**.
- 조립 partial CSR는 정확 대칭이며 MPI 합산 diagonal도 모두 양수이다.
- native **26.0184초**, MPI4×OpenMP8, native/runner exit 0. 출력은 순수 파생 행렬이며 물리 dump는 없다.

이것은 실제 f를 푼 결과가 아니다. 원래 전체 방정식 `||Ma-f||/||f||<1e-16`, cell work/torque, 같은 f의 원래 full-solve shadow 비교는 **미실행**이다.

## 3. 전체 파일 native 재읽기

`material_cell_condensation_readback_v1_mpi4/result.json`.

전체 71,080 source, 모든 L/M_IG 및 96,946,492 CSR entries를 다시 읽었다. 샘플 검사가 아니다.

- 각 partition/factor/CSR 파일의 전체 FNV payload checksum와 EOF 통과.
- 모든 저장 L diagonal 양수, 모든 partial CSR 대칭.
- 모든 binary80 L, M_IG, CSR 값이 double-double 두 성분으로 정확하게 표현됨을 확인.
- native exit 0, native **11.5009초**. actual full-equation solve와 GPU 복구의 정확도 관문은 계속 pending.

## 4. GPU 코드와 메모리

`material_cell_condensation_gpu_build_v1`의 frozen source와 `cell_condensation_gpu.o`: Ada sm89 및 H200 sm90를 포함한 CUDA 컴파일 통과. 원래 CUDA 물리 코드와 기존 DD CSR backend를 재사용하고, 원래 source별 L/M_IG로 두 가속도 성분을 축약·복구하는 커널을 추가했다. GPU 수치 실행은 아직 하지 않았다.

full f_G는 이미 전역 조립된 값이므로 source마다 중복해서 더하지 않는다. source-local `-M_GI M_II^-1 f_I`만 고정 source 순서로 모으고, f_G는 그 DOF의 owner에서 한 번 더한 뒤 MPI 합산하도록 구현했다. 내부 복구에는 모든 원래 I 모드가 포함된다.

현재 보수적인 제안 배치의 합산 GPU 예산은 **106,707,056,260 bytes**. 불변 cache만 **59,891,910,888 bytes**이다. 이는 실제 모든 backend allocation을 측정한 결과가 아니며 최종 runtime 회계가 필요하다. Ada 단일 장치에는 이 배치가 들어가지 않는다.

H200 capacity 작업 **434126**, `--partition=h200 --gres=gpu:H200:1`, 마지막 실제 squeue 확인은 **PENDING / QOSGrpGRES**. 작업을 재제출하지 않았다.

CPU 구축은 GPU 예약과 독립적으로 실행했다. 사전 예산 근거는 같은 H200 pool에서 실제 실행됐던 `material_full_initial_v2_h200_mpi4/device.txt`의 **143,771MiB = 150,754,820,096 bytes** 장치 기록이다. 제안량은 기록된 용량의 80% 이내이다. 이 사실은 **현재 free-memory 측정이나 현재 allocation 성공을 뜻하지 않는다**. `setup_capacity.json`과 factor 결과에 fresh GPU allocation pending을 명시했고, GPU 실행 전에 별도로 확인해야 한다. 최초 준비 중 읽기 전용 복사본 덮어쓰기 오류는 실행 전에 수정했다. `src_revision2`가 실제 CPU 구축 source이며 앞선 준비본으로 native 계산을 실행하지 않았다.

## 다음 연속 작업

1. 기존 step5에서 다음 실제 RHS를 계산하는 경로에 factor reader, DD local condensation/recovery, assembled skeleton solver를 연결한다. step4의 high k1을 새 step5 RHS로 재사용하지 않는다.
2. 실제 같은 f의 원래 full solve를 한 번만 shadow하고 전체 원래 M residual, cell work/torque 및 모든 내부 가속도를 확인한다. GPU 커널 컴파일 통과로 대체하지 않는다.
3. 같은 selected point의 불변 enriched basis/gradients와 generating j0 cache를 구현한다. cache/direct 같은 고정 tile 비교와 source/basis/rule/thermal SHA binding이 필요하다. 현재 x/v/J/P는 매 stage 다시 계산한다.
4. 같은 실제 RK4 step commit과 restart readback, guard를 포함한 실제 step 비용을 확인한다. h≈0.000345838 기준 t10≤7일에는 약20.9초/step이 필요하다. 이번 구축 시간이나 실패 compact의 2초 solve를 이 비용으로 대신하지 않는다.
5. 비용·물리·map 관문 통과 시 원래 경로로 t=.1/.2/1/6/8/10과 최장 공전까지 이어간다. material map 수명, 전체 injectivity guard, conservative relabeling/열 recipient, W=0 환원 및 GIZMO 조건 차이는 해결되지 않았다.

단일 exact-condensation 후보의 setup/정확성/비용 관문이 실패하면 이를 보존하고 종료한다. 반경·차수·cache partition·tolerance sweep으로 실패를 가리지 않는다.

증거 묶음: `material_cell_condensation_setup_completion_v1/`. 구조·행렬 준비 통과는 장기 Kepler 통과가 아니다.
