# Reference cache 용량 종료평가 — 2026-10-22

## 결정

**이번 Ada exit2는 현재 할당 순서의 장치 용량 실패이며 exact cell condensation 후보의 물리·수학적 실패가 아니다. 기존 고정 source tile과 동결 H200 빌드의 실행 자격을 유지한다.** H200 할당이 실제로 주어지면 같은 입력·원래 k1·원래 4096-triangle tile로 실제 용량 관문부터 진행한다. PENDING은 할당이나 실행 성공이 아니다. 이 검토에서 작업 조회·계산·재실행은 하지 않았다.

H200 대기 중 허용하는 구현은 **한 가지: reference 생성용 GPU 입력과 정상 stage용 GPU 데이터의 생존 기간 분리**다. tile 크기, rule, 산술, 솔버, 20% 여유 기준을 바꾸는 탐색은 하지 않는다. 이 변경은 별도 immutable 빌드로 기록하며 기존 동결 바이너리와 동일하다고 부르지 않는다. 두 경로가 동시에 실제 step6을 실행해서는 안 된다.

## 직접 확인한 근거와 정정

- `material_cell_condensation_actual_rhs_completion_v1/contract.json`: 실제 step6 k1에서 원래 전체 M 잔차 1.70e-19, 71,080셀 장부 통과. 동일 F 솔버 16.0610→3.58819초는 솔버 비교이며 전체 RK4 가속은 아직 미측정이다.
- `material_cell_reference_continuation_v1_ada_mpi4/native.log`: cache 생성 전 collective 용량 검사에서 거부했다. 새로운 accepted step이나 reference 생성 성공이 아니다. 원래 step5와 저장 k1을 유지한다.
- 실제 `material_full_initial_batch.h`는 `MF_BATCH_TRIANGLES=4096`이다. 과거 metadata의 1024와 그로부터 나온 약325MB 추정은 잘못됐으며 채택하지 않는다. 기존 파일은 보존하고 이 문서로 정정한다. 실제 `point_tiles`와 용량 검사에는 올바른 4096이 사용됐다.
- 실제 solver 할당 19,244,169,964 bytes, 장치 총량 33,794,359,296 bytes, 할당 후 여유 7,721,254,912 bytes이다. 20%를 보존하면 추가 허용량은 약962MB뿐이다. `Cache`는 모든 rank의 실제 최대 tile 크기와 payload stride를 합산하므로 이 제한의 거부가 타당하다.
- `material_reference_cache_host.h/.cu`는 전체 immutable 59,891,910,888 bytes와 별도 current-point112 bytes를 분리한다. 생성은 원래 `evolved_point`와 generating coefficient/delta를 사용한다. 정상 cached kernel은 원래 point ID/weight 대신 저장값을 읽고 `MFRuleView`에서는 count를 사용한다. 이 차이가 수명 분리의 근거다.

## 단일 구현 계약: 생성 → 입력 해제 → 정상 stage

1. **할당 전에 실제 구성으로 두 단계의 최대 동시 메모리를 계상한다.** compiled 4096, 실제 tile prefix/max-points/max-doubles, GPU UUID별 모든 rank, solver/physical/context/scratch/통신 버퍼를 포함한다. host pinned 약60GB의 실제 가용량도 확인한다. 어느 단계든 원래20% 여유를 못 지키면 cache 생성 전에 종료한다. 이후 `cudaMemGetInfo`로 실제 할당을 다시 확인한다. device 수를 rank 수로 간주하지 않는다.
2. Ada 대안에서는 **condensed solver의 GPU 할당을 늦춘다.** host의 원래 M·분할·상태·k1·rule은 유지한다. 원래 전체 selected IDs/weights로 고정 tile reference cache를 한 번 생성하고, 생성 완료·입력 binding·원래 direct first-tile witness를 확인한 다음 생성 전용 전체 GPU IDs/weights를 해제한다. host 원본 및 hashes는 유지한다. 정상 kernel에 필요한 counts/offsets와 active/generation 메타데이터는 제거하지 않는다.
3. Direct witness가 다음 실제 stage의 current field를 요구하면 **그 stage의 prepare/witness를 solver GPU 생성 전에 수행**할 수 있다. 원래 저장 k1의 force/CG는 재계산하지 않는다. witness가 끝난 뒤 전체 point IDs/weights를 해제하고 condensed solver를 한 번 생성한다. 이후 stage는 cache만 소비한다. 이 순서의 실제 peak가 맞지 않으면 중단한다. 임의로 witness를 생략하거나 원래 point 순서를 대체하지 않는다.
4. 캐시 완료는 reference point별 원래 좌표/weight/rho/j0/enriched jets의 생성 완료를 뜻한다. 현재 x/v/J/pressure/thermal은 매 stage 재평가한다. cache 생성 중 실패하면 불완전 cache를 active로 게시하지 않는다. immutable generating thermal context를 현재 상태로 재기준화하지 않는다.
5. 정상 cached mode에서 해제된 원본 GPU 포인터를 읽는 경로가 없도록 소유권·mode 계약을 명시한다. 기존 direct mode는 보존한다. 함수 시그니처에 `MFRuleView`가 남아 있다는 사실만으로 모든 멤버가 불필요하다고 가정하지 않는다. count 검사는 계속한다. 향후 direct witness가 필요하면 명시적인 별도 수명 구간으로만 복원한다.

이 경로는 cache/solver를 작게 바꾸는 것이 아니라 동시에 살아 있을 필요가 없는 버퍼를 분리한다. 물리 값을 변경하지 않지만 아직 용량·성능 성공은 아니다. 약3GB 해제가 충분한지도 실제 두 단계 회계로 결정한다.

## 실행 및 중단 조건

- H200 기존 동결 경로 또는 위 수명 분리 경로 중 **한 경로만 공통 claim을 취득**하여 step6을 계속한다. 이미 성공한 actual k1을 그대로 사용하고 남은 stage와 endpoint/restart를 원래 관문으로 처리한다. 기존 accepted 구간과 동일 F 비교를 반복하지 않는다.
- 실제 용량 거부, cache/direct 불일치, 원래 전체 M 잔차·전셀 energy/angular/thermal/양성 실패, checkpoint binding 실패는 게시 전 중단 사유다. headroom·tolerance·precision·rule·tile 완화나 자동 대체 경로 탐색으로 넘기지 않는다.
- H200의 추가 VRAM은 가능성을 주지만 같은 속도를 보장하지 않는다. Ada와 H200의 실제 prepare/transfer/solve/전체 step을 분리해 보고한다. solver 4.48배를 전체 가속으로 외삽하지 않는다.
- 기존 `ASTRA_EXACT_CELL_CONDENSATION_DECISION_20261022.md`의 전체 step 비용 관문을 그대로 적용한다. cache 생성·solver setup의 일회성 비용과 steady 비용 모두 보고한다. 실제 측정 전 t1/t6/t10 장기 연장은 승인된 것으로 해석하지 않는다.
- 이 한 번의 수명 분리로도 실제 용량 또는 전체 비용 관문을 만족하지 못하면 해당 Ada 구현 경로를 종료한다. H200 대기는 별개의 자원 상태로 남기되 무한 타일·솔버 최적화를 재개하지 않는다.

## 범위

이는 동일 full-space 연산자의 실행 자원 계약 결정이다. global map injectivity, material shear 수명, moving Pext 경계와 원래 periodic 기준선의 차이, remesh 및 조건부 W0 환원은 해결되지 않았다. 실제 t6/8/10 및 논문 최장 공전의 가짜가열·밀도·E/L 검증이 최종 목표이며, 이번 용량 관문 또는 다음 한 step 성공으로 대체할 수 없다.
