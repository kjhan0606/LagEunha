# 원래 step5의 실제 다음 RHS — exact cell condensation 통과

2026-10-22. 연구 ROOT `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics`.

## 이번에 실제로 계산한 것

`material_cell_condensation_actual_rhs_v1_ada_mpi4/run/actual_rhs_result.json`.

원래 accepted **step5, t=0.00011263377801805954**를 v3 체크포인트로 읽고, 같은 committed selected quadrature generation **11882225342103221522**에서 **step6의 첫 실제 RHS 하나**를 계산했다. 전 영역 71,080 cell, 1,138,432 triangle, 249,691,129 selected point/pass를 사용했다. 새 초기조건·원래 20억 점 초기 적분·초기 curvature·완료 구간 재실행은 없다.

물리적 힘은 한 번만 계산했다. 동일한 조립 F에 대해 기존 원래 전체 resident CUDA 풀이를 한 번 수행하고, 그 GPU mass/vector buffers를 해제한 뒤 새 cell-condensed 풀이를 수행했다. 같은 F와 두 가속도 파일을 저장했다. 모든 원래 내부 가속도는 복구했으며 x_I/v_I를 제거하지 않는다.

native/runner exit **0**. 이 묶음은 실제 RK4 첫 단계의 미분값을 준비한 것이며 **accepted step 추가는 0**이다. 아직 t6/8/10이나 GIZMO 최장 공전 검증은 아니다.

## 정확성 결과

`run/stage1/cell_condensation_actual_solve.json`, 원래 `rhs.json`, 네 `cells_r*.csv`.

| 항목 | 실제 값 |
|---|---:|
| 복구된 전체 원래 `||Ma-f||/||f||`, x | 1.7033805179717221e-19 |
| 복구된 전체 원래 `||Ma-f||/||f||`, y | 1.7043260912422612e-19 |
| 전체 defect correction 횟수 | 1 |
| 두 성분 skeleton DD PCG 반복 합계 | 104 |
| 기존 전체 풀이와 가속도 차이의 원래 M 상대 노름 | 3.6143852408318596e-18 |
| 가속도 계수 최대 차이 | 3.9378222904673521e-15 |
| 비교 가속도 계수 최대 규모 | 174.23592925946777 |
| 같은 F의 work 차이 / U | 2.3275571388216054e-19 |
| 같은 F의 torque 차이 / 2K | 5.5809900361340136e-23 |
| 모든 source 내부 I 모드의 최대 조립 반력 / 최대 load | 4.2812901809322347e-19 |
| 모든 셀 energy rate 오차 / U의 최대 | 7.575388533598615e-13 |
| 모든 셀 angular rate 오차 / 2K의 최대 | 3.6623927223499035e-16 |
| 질량행렬 K와 point K의 차이 / U 최대 | 1.8664348521008514e-11 |
| 질량행렬 L과 point L의 차이 / 2K 최대 | 9.597681052579663e-15 |

원래 point positivity, timestep, 각 source의 work/torque, resolved spin/inertia/omega/Erot 및 bubble/shared reaction 관문을 그대로 실행했다. 저장된 CSV **71,080행 전체**를 별도로 읽어 energy/angular 최대 관문을 재확인했다. 독립 spin 저장고, 열 리셋, 조향, 점성 추가는 없다.

## 속도 — 실제 같은 RHS 비교

Ada 한 장, MPI4×OpenMP8.

- 새 primary 풀이 **3.58819229초**. 이 값에는 축약, DD skeleton 풀이, 전체 내부 복구, 원래 binary80 M true residual이 포함된다.
- 축약 0.58049275초, 복구 0.557381686초, 원래 M true residual 적용 0.476916803초.
- skeleton inner 1.85632246초, 그 중 halo 회계 1.62402977초. 이 시간들은 포함 관계이며 더해서 총 시간이라고 보고하지 않는다.
- 기존 full shadow 17.0945575초, 그 중 한 번의 GPU 업로드 1.03352339초. 업로드를 제외한 기존 풀이 **16.06103411초**와 비교하면 같은 힘의 풀이가 **4.4761배** 빨라졌다.
- 새 solver constructor 본체 setup 15.3390278초는 primary 풀이와 분리한다. 그 timer보다 앞선 partition/CSR reader 시간은 별도이며, 기존 full shadow의 첫 resident vector 초기화 시간은 분리 측정하지 못했으므로 4.4761배를 순수 steady-state 비율이라고 단정하지 않는다.
- 캐시 없는 실제 force/geometry 준비 21.7019741초. 이 force는 이미 저장했으며 재계산해서 속도 검사를 반복하지 않는다.
- native 전체 실행 87.5580773초. 입력/출력 SHA 검사와 wrapper 시간은 별도다.

이 결과는 **전체 RK4 스텝의 속도 향상이나 t10 비용 관문 통과를 증명하지 않는다**. 불변 reference cache, 실제 후속 RK stages, endpoint, map guard 비용을 포함한 측정이 남아 있다.

## 실제 GPU 메모리

`run/stage1/cell_solver_allocation.json`.

- 새 cell/skeleton backend 전체 MPI 합산 explicit CUDA allocation **19,244,169,964 bytes**.
- 기존 원래 physical context는 이미 올라간 상태에서 측정했다.
- 실제 장치 전체 **33,794,359,296 bytes**, 새 풀이 업로드 직전 free **27,056,996,352 bytes**.
- 전 rank 업로드 후 free **7,721,254,912 bytes**. 장치 전체의 20% 이상 여유 관문 통과.
- 이 단계에는 불변 reference cache와 기존 거대한 transient cache를 할당하지 않았다. H200용 전체 resident cache 60GB가 Ada에 들어간다고 주장하지 않는다.

## 실제로 재사용할 다음 RK 첫 단계

`run/prepared_selected_k1_r[0-3].bin`.

- 실제 원래 coefficient hash, immutable manifest hash, selected rule generation, owner map, epoch, accepted time/h를 결합한 새 형식이다.
- magic `0x4d434e4558544b31`, version1, 11개 uint64 header, time/h, binary80 work/torque/heat/dt_bound 네 값, 원래 STCoefficient 형식의 k1 및 전체 checksum.
- 과거 step4 high k1을 쓰지 않았다. 이 파일은 **현재 step5에서 직접 계산한 selected k1**이다.
- 이후 loader는 모든 binding과 checksum, xdot=current velocity의 double bits를 확인해야 한다. 아직 실제 RK4 연속 단계에서 소비하지 않았다.

`stage1/actual_force_r*.bin`, `same_force_old_acceleration_r*.bin`도 보존했다. 실패 witness와 기존 체크포인트를 지우거나 덮어쓰지 않았다.

## 불변 reference cache의 다음 구현 제약

`material_cell_condensation_actual_rhs_completion_v1/fixed_tile_cache_memory_plan.json`.

기존 `point_tiles`의 complete-source partition, MF_BATCH_TRIANGLES=1024, point-cap=1,048,576을 그대로 사용해 saved metadata로 최대 실제 tile을 산출했다. rank별 tile 수 279/278/278/278, 최대 point 수 229,200 / 229,782 / 232,180 / 231,425이다. 타일 분할을 탐색하거나 물리 source를 바꾸지 않았다.

새 캐시는 fixed reference N의 비다항 value/gradients, lambda, weight/rho/j0와 current-stage point 상태를 분리해야 한다. 기존 MFCachedPoint의 18개 STJet를 매 point 중복 저장하면 Ada의 20% 여유를 넘는다. fixed point 최대240bytes + current point112bytes의 분리 형식이면 **같은 기존 타일**의 합산 추가 저장 상한은 **324,799,824 bytes**로 추산된다. 위 실제 allocation 상태 기준 free 예측 7,396,455,088 bytes로 20% 여유 안이다. 이는 **아직 구현/할당/정확성 검증이 끝난 결과가 아니다**.

H200이면 불변 전체 cache의 device resident 경로를 사용할 수 있고, resident가 불가능하면 같은 고정 source tile의 pinned-host 경로가 계획에 허용돼 있다. GPU에서 실제 사용하는 타입의 sizeof, 모든 실제 allocation, source/basis/rule/thermal/ABI binding 및 같은 새 RK state에서 cache/direct tile 비교를 확인해야 한다. current x/v/J/P는 매 stage 새로 계산한다. cache/direct 비교를 위해 전체 완료 구간을 반복하지 않는다.

## 코드와 후속 단계

- 새 CPU/driver: `material_cell_condensation_solver.h`, `material_cell_skeleton_resident.h`, `run_cell_condensation_actual_rhs_20261022.cpp`.
- 실제 실행 frozen closure: `material_cell_condensation_actual_rhs_v1_ada_mpi4/src/`.
- 원래 검증된 frozen evolution header를 별도 실행 디렉터리에 복사해 **전체 mass 풀이 지점의 hook 및 구현 보고 flag 두 부분만** 연결했다. 원래 물리 force 및 모든 cell ledger/gate 수식은 유지했다. 원본 저장소나 이전 frozen 실행 파일은 수정하지 않았다.
- `material_full_initial_cuda.cu`, evolved/cached CUDA, shared trace, geometric pressure, positive-rule reader가 step5 실행본과 동일 bytes임을 확인했다.

다음은 불변 reference cache 및 같은 새 상태의 cache/direct 관문을 구현하고, 위 실제 k1을 소비해 stage2–4와 endpoint/restart를 완료하는 것이다. 그때 실제 full-step 비용을 측정한다. 장기 승격 전 full-map injectivity guard는 별도로 필요하며 이번 대수적 성공으로 pending을 해제하지 않는다. 최종 목표는 계속 t>6 안정성과 GIZMO 최장 공전이다.

증거 묶음: `material_cell_condensation_actual_rhs_completion_v1/`.
