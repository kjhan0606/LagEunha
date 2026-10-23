# 실행 자원 결정 — 2026-10-23

## 최종 정정: scheduler test-only 이후 기존 H200 보존

**기존 H200 job435137을 CPU/OMP/메모리/바이너리까지 변경 없이 보존한다. H100 대체 제출·기존 작업 취소를 집행하지 않는다.** 아래 H100/OMP4 방침과 그 취소·대체 절차는 실행 전에 superseded되었으며, 역사적 판단 근거로만 남긴다. 초안 H200/OMP1도 실행 지시가 아니다.

운전자가 전달한 실제 scheduler `sbatch --test-only` 결과는 H100 test-only435140의 예상 시작을 **2026-11-04T09:14:49**, 기존 H200435137은 **2026-10-24T04:44**로 제시했다. 시간은 전달된 scheduler 표기 그대로이며 이 검토에서 독립 조회하지 않았다. test-only는 물리 실행·예약·할당이 아니다. syn08의 물리 여유(CPU6/64 사용, GPU3/5 사용, AllocMem96/500G)만으로 정책/예약하의 조기 실행을 기대한 판단보다 이 배치 예측이 직접적인 근거다. 어느 시작 시각도 보장으로 부르지는 않는다.

현재 단일 실행 후보는 원래 H200 1장/MPI4/rank당8 CPU/192G 요청이다. 같은 step5·저장 step6 k1·공통 atomic claim·동결 source/binary를 유지한다. Ada lifetime v2는 build PASS/native 미실행으로 남긴다. 다른 CPU/GPU 조합 test sweep이나 다른 작업 중단은 하지 않는다. QoS typed GRES 읽기 결과는 대기 원인 기록에 추가할 수 있지만 그 자체로 새 경로 자동 실행을 허용하지 않는다.

실제 할당 후 기존 용량·cache/direct·원래 M 잔차·전셀 물리·restart·전체 비용 관문을 적용한다. 기다리는 동안 상태나 열 문맥을 바꾸지 않는다. 이번 정정은 자원 선택의 정정이며 수치 후보 실패/성공 판정이 아니다. 아래 H100 pinned-host와 H200 resident 비용을 혼동하지 말라는 한정은 계속 유효하다.

## 실행 전 폐기된 H100 방침 — 기록 보존

**추가 H100 자원 증거를 반영한 최종 방침은 H100 1장·원래 MPI4·OMP4·원래 동결 sm90 바이너리로 단일 대체 요청이다.** `--partition=h100 --gres=gpu:H100:1 --nodes=1 --ntasks=4 --cpus-per-task=4`, `OMP_NUM_THREADS=4`, 중첩 라이브러리 thread는1, 기존192G 메모리와 원래4096 source tile·20% VRAM 여유·물리 입력을 유지한다. 기존 Cache의 resident/pinned-host 선택을 그대로 사용한다. 새 Ada lifetime 코드는 실행하지 않는다. 초안의 H200/OMP1 방침은 추가 증거에 따라 이 결정으로 대체하며 두 방침을 동시에 실행하지 않는다.

운전자가 전달한 H100 두 노드 합계 CPU102개 여유는 H200 CPU62/64 사용 및 약하루 뒤 예상 시작과 비교하여 다른 기존 장치를 선택할 구체적 근거다. 합계 여유가 단일 노드의 GPU·CPU16·192G 동시 배정을 보증하지는 않는다. 제출 전 운전자의 읽기 확인에서 노드별 GPU/메모리/QoS 제약을 확인하고, 배정은 scheduler에 맡긴다. 새 요청이 대기해도 다시 CPU/GPU 조합을 순회하지 않는다.

운전자는 자기 소유 PENDING job435137이 실제 시작하지 않았고 공통 claim이 없는지 확인한 뒤 이 한 번의 자원 변경을 적용한다. scheduler가 기존 PENDING 요청 수정을 허용하지 않으면 자기 작업만 terminal 취소 확인 후 대체 하나를 제출한다. 이미 RUNNING이거나 claim을 취득했다면 변경하지 않고 기존 실행을 유지한다. 다른 프로젝트/스레드 작업은 건드리지 않는다. 이 검토자는 잡 조회·제출·취소·실행을 하지 않았다.

## 산술·물리적 정당성

읽은 `material_distributed_mass.h`의 `Ownership::dot`은 owned ID 고정 순서의 serial binary80 합 뒤 MPI 합이다. OpenMP reduction으로 dot을 바꾸지 않는다. mass apply/precondition은 출력 row 또는 triangle별 독립 작업을 병렬화하고 각 출력의 내부 합 순서는 고정이다. MPI4와 원래 ownership/halo를 유지하므로 CPU thread 감소는 그 분할을 바꾸지 않는다. 실패 수의 정수 reduction과 진단 max는 물리 force 합산 순서 변경이 아니다.

`material_full_evolution_host.h`의 stage coefficient 갱신은 DOF별 독립이고, 주 force/현재 geometry/EOS는 원래 CUDA 경로를 사용한다. 따라서 OMP1을 위해 force, geometry, MPI 순서, binary80 true residual, DD solver, generating thermal context, RK4 식을 바꿀 필요가 없다. 다만 전체 실행의 bitwise 동일성을 정적 읽기만으로 인증하지 않는다. 기존 native cache/direct witness·전셀 work/torque/열·양성·원래 M 잔차·restart 관문이 그대로 판정한다. CPU 수 비교용 별도 step이나 저장 k1 재계산은 하지 않는다.

OMP8→OMP4는 OpenMP 지원을 유지하면서 이 실행의 thread 수만 제한한다. 사용자 native CUDA+MPI/OpenMP 코드 사용 요구와 양립한다. 앞 문단의 산술 근거는 thread 수 변경에 적용되며 H100/H200 사이의 bitwise 동일성 증명은 아니다.

## H100과 lifetime 경로의 위치

동일 sm90은 코드 실행 가능성의 일부 조건이지 VRAM·host bandwidth·성능 동일성의 증거가 아니다. 기존 Ada 대안과 포괄 native GPU 실행 승인이 유효하다는 요청 문맥에서, 앞선 H200 gres는 해당 제출의 자원 지정으로 해석하며 H100 금지 조건으로 확대하지 않는다. H100 pinned-host mode는 기존 승인된 동일 함수 경로다. H200 probe의106,707,056,260 bytes 성공은 H100 용량 증거가 아니므로 H100의 실제 UUID/free memory/20% 검사를 반드시 통과해야 한다. 부족하면 cache 생성 전에 종료하며 tile/여유 기준을 바꾸지 않는다.

Ada lifetime v2는 이번에 실행하지 않는다. 읽은 host 계약의 UUID 일치, 두 단계 device 합산, 생성·direct witness 뒤 IDs/weights만 해제, counts/offsets와 host rule 보존 순서는 앞선 승인과 맞는다. 명백한 추가 주의점은 host capacity다: `/proc/meminfo`의 `MemAvailable`은 Slurm cgroup 한도가 아니다. **선택한 H100 원래 pinned-host 경로에도** job/cgroup 남은 메모리와 host 가용량 중 작은 값으로 약60GB pinned cache+현재 RSS+잔여 workspace peak를 제한한다. 이 검사/기록은 wrapper의 실행 전 자원 관문으로 두어 native 동결 바이너리를 바꾸지 않는다. 이 조건 없이 `host_capacity_pass`를 job 메모리 보증으로 해석하지 않는다. 전체 lifetime 실행·소유권 검증을 완료했다고 주장하지 않는다.

## 비용 판정과 종료

이 실행은 H100/MPI4/OMP4 및 실제 cache mode라는 구성으로 cache 생성, host-device 전송, 준비, CPU 잔차/halo, GPU solve, 전체 step 시간을 분리한다. queue 대기는 물리 step 시간에서 제외하되 실제 도달 소요시간에는 별도로 기록한다. 기존 Ada/OMP8 또는 H200/OMP8과의 차이를 cache만의 가속으로 부르지 않는다. H100의 전송/CPU 비용을 임의로 빼서 전체 비용 관문을 통과시키지 않는다. 반대로 비용 실패를 H200 resident 경로의 불가능성 정리로 확대하지도 않는다. 결과는 해당 실제 구성의 비용 실패로 기록하고 자동 CPU/GPU sweep 없이 묶음을 종료한다.

공통 `material_step6_reference_continuation.claim`은 그대로 사용한다. 원래 step5 및 실제 저장 step6 k1을 유지하며 남은 stage/endpoint와 원래 restart 검증만 수행한다. 현재 상태·source·binary·rule·thermal SHA 검사를 유지한다. 물리/용량 관문 실패 시 게시하지 않는다. 실제 전체 step 결과 전 장기 연장 성공이나 t1/t6/t10 비용 달성을 선고하지 않는다.

## 읽기 범위와 identity 한정

공유 지침을 전체 읽었고 현재 CWD는 backup이며 원본 repo로 취급하지 않았다. 지정 ROOT의 mass/evolution/cache 및 동결 lifetime v2 host/provenance를 읽었다. SSH의 원본 worktree git metadata 조회는 `.git`이 가리키는 `/tmp/lageunha-user-master/.git/worktrees/src3` 부재로 실패했다. 따라서 이번에 HEAD를 새로 검증했다고 쓰지 않으며 운전자의 최근 HEAD8e25dfe/원본 identity 증거와 frozen hashes를 계승한다. 이 문서 외의 source/결과는 수정하지 않았다.

최종 목표는 원래 native 전체 Kepler t6/8/10의 가짜가열·궤도·밀도·E/L 검증과 논문 최장 공전이다. 자원 배정 또는 한 step 성공은 그 목표를 대신하지 않는다.
