# Reference cache 수명 분리 구현 및 H200 재개 준비 — 2026-10-23

## 판정

승인된 단일 버퍼 수명 분리를 구현했고, 분리된 동결 빌드 `material_cell_reference_lifetime_v2_ada_mpi4`의 CUDA sm89/sm90, host, link가 모두 exit0이다. 물리 실행은 아직 하지 않았다. 원래 accepted 상태는 step5, t=0.00011263377801805954이며 저장된 실제 step6 k1을 유지한다. 이 문서는 t6/8/10 안정성이나 가짜 가열 해결의 증거가 아니다.

## 구현한 수명 계약

1. 원래 분할 metadata와 skeleton ownership으로 solver CUDA 할당량을 계산한다. 기존 원래 solver와 동일한 19,244,169,964 bytes인지 확인한다. 인수·질량행렬·trial space는 바꾸지 않는다.
2. 실제 GPU UUID를 모든 rank에서 비교한다. 원래 MPI4가 사용하는 한 카드에 네 rank의 모든 byte를 합산한다. 실제 free/total, 원래4096 tile, 1,048,576 point cap, reference buffer와 현재 point112 byte를 사용해 생성 단계와 정상 stage 단계의 free byte를 각각 기록한다. 원래20% 기준을 유지한다. Pinned-host cache 약60GB와 추가 host workspace도 검사한다.
3. Solver GPU 할당 전 원래 전체 selected ID/weight로 reference cache를 생성한다. 생성 offset을 순서대로 누적하여 전체 immutable payload가 완성된 경우에만 reference-only 전환을 허용한다.
4. 저장 k1은 그대로 사용한다. 다음 실제 stage의 원래 direct/cache witness가 완료되고 해당 stage의 F가 이미 조립된 뒤, 같은 RHS hook 안에서 생성 전용 GPU ID/weight만 해제한다. Host 원본과 GPU counts/offsets 및 generating context는 유지한다.
5. 해제 후 original direct selected-point API 및 재생성 호출을 명시적으로 거부한다. Device point 함수에도 누락 포인터 검사를 두어 해제된 ID/weight를 읽지 않도록 한다. 원래 direct 모드는 해제 전까지 보존된다.
6. 같은 이미 조립된 F에서 GPU solver를 한 번 생성한다. 실제 solver byte가 사전 계획과 일치하고 실제 할당 후20% 여유가 남는지 확인한다. 현재 stage 또는 저장 k1의 force/CG를 반복하지 않는다.

이 byte 계산은 실행 시 실제 용량을 검사하기 위한 구현이다. 현재 점유 중인 Ada에서 실제 두 단계 할당이 통과했다고 주장하지 않는다. H200에 제출된 기존 동결 후보와도 다른 빌드이다.

## 빌드 증거와 보존

- `material_cell_reference_lifetime_v2_ada_mpi4/build_result.json`
- `material_cell_reference_lifetime_v2_ada_mpi4/cuda_architectures.txt`: sm89, sm90 ELF 확인.
- `source.sha256`, `binary.sha256`, `postbuild_source_integrity.txt`.
- `material_full_evolution_host.h`, `material_full_cached_host.h`, 실제 cached physical kernel, cell solver, DD GPU kernel, skeleton backend 및 distributed M 소스는 이전 동결 후보와 byte가 같다.
- 수정된 source10개는 메모리 사전 계산, mode/ownership guard, cache 생성 완료 기록 및 할당 순서 변경이다. Rule, tile, physics 식, precision, tolerance, RK4는 바꾸지 않았다.
- 첫 빌드 `material_cell_reference_lifetime_v1_ada_mpi4`는 CUDA 빌드는 통과했지만 metadata 합산 `MPI_Allreduce`에 `MPI_SUM`을 누락해 host compile에서 거부됐다. 소스·오류 로그를 보존했다. v2는 이 호출 인자만 고친 같은 구현이다. 첫 빌드에서도 물리 실행은 없었다.

## H200 경로와 공통 실행 소유권

이전 capacity job434126은 실제로 syn104에서 COMPLETED/0:0였다. 측정 total150,110,011,392 bytes에서106,707,056,260 bytes 단일 계획 할당을 통과했다. 이는 다음 actual MPI cache/solver 실행의 성공을 대신하지 않는다.

기존 동결 continuation 바이너리를 그대로 쓰는 job435137을 `--partition=h200 --gres=gpu:H200:1`, MPI4×OMP8, host192GiB로 제출했다. 최근 관측은 PENDING/Resources이며 Slurm 예상 시작은2026-10-24T04:44:42였다. 예상은 예약 확약이 아니며 상태 확인 시 같은 job ID를 조회한다. syn104는 관측 시 CPU62/64, GPU6/8 할당 상태여서 GPU 여유만으로는 이32CPU 작업이 시작되지 않았다.

두 경로의 실제 실행은 공통 atomic 디렉터리

```text
ROOT/material_step6_reference_continuation.claim
```

를 취득해야 한다. Claim은 끝나도 남겨 accepted state를 감사하며 자동으로 삭제하거나 대체 경로를 실행하지 않는다. Claim 파일만으로 실행 중이라고 판단하지 않고 Slurm/PID의 실제 상태를 확인한다. 다른 geometry probe 및 타 프로젝트 작업은 종료하거나 변경하지 않았다.

## 다음 관문

자원을 실제로 할당받은 한 경로에서 저장 k1 이후의 원래 단계, 모든 셀 장부·양성·원래 M 잔차·restart readback을 확인한다. 그 다음 완전히 새로운 RK4 한 step의 전체 시간을 측정한다. Solver의 기존 약4.48배 비교를 전체 step 가속으로 외삽하지 않는다. 전역 map guard와 실제 t6/8/10 및 GIZMO 최장 공전 검증은 계속 미완료이다.

## 자원 배치 확인 후 최종 방침

`ASTRA_EXECUTION_RESOURCE_DECISION_20261023.md`의 정정 결정을 채택했다. H100 노드의 물리 여유가 관측됐으나, 실제 제출 전 `sbatch --test-only`가 H100/MPI4/OMP4/192G 구성의 예상 시작을2026-11-04T09:14:49로 제시했다. 이435140은 test-only 표기이며 실제 대체 작업을 제출하거나 자원을 예약한 것이 아니다. H200의2026-10-24 예상보다 늦으므로 H100 방침은 실행 전에 폐기했다.

**기존 H200435137을 CPU/OMP/메모리/바이너리까지 그대로 보존한다.** 취소·재제출·Ada 실행·추가 CPU/GPU 조합 탐색은 하지 않았다. QoS normal의 읽기 결과는 사용자 GPU 상한8이지만, 이를 이번 대기의 확정 원인으로 단정하지 않는다. 실제 scheduler reason은 Resources이다. 물리적 GPU 여유와 실제 scheduler 배정 가능성은 구별한다.

새 lifetime 코드의 host `MemAvailable`은 Slurm cgroup 한도를 증명하지 않는다. 이번에는 해당 코드를 실행하지 않았으며, 추후 제한된 job에서 pinned-host 경로를 실행할 때는 job/cgroup 남은 메모리와 host 가용량 중 작은 값에 실제 workspace를 포함한 관문이 필요하다. H100과 H200의 resident/pinned-host 비용도 동일시하지 않는다.
