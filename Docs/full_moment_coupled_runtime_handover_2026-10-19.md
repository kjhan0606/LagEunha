# 자기일관 midpoint MPI 실행 인수인계 — 2026-10-19

## 선택한 방법

보존된 explicit fine `t=0.025` checkpoint에서 별도 적분법 분기를 만든다. 같은 공간 연산자와 fixed mass, Re=∞, 무조향, common-relative-W gauge를 사용한다. 새 방법은 endpoint `Y* = M_h(Y_base, F((Y_base+Y*)/2))`의 실제 map 잔차를 반복 평가한다. 성공 시 마지막 trial이 아니라 검증된 map 출력을 받는다. 공간 closure, 회전 모멘트 및 covariance 지도를 그대로 유지한다. 장시간 정확도를 이미 달성했다고 해석하지 않는다.

## 코드 연결

- `Exam/full_moment_coupled_runtime.inc`가 controller callbacks를 기존 `Exam/full_moment_evolution.c` native MPI 실행 경로에 연결한다. 매 반복마다 실제 midpoint state/geometry/density/velocity/gravity/face provider를 새로 평가한 뒤 native MPI/GPU factor map을 한 번 적용한다.
- 외부 잔차는 위치·공통 gauge 제거 W·centroid·속도·S·C·R·thermal·Psi remainder와 dependent Psi 조건을 포함한다. 한 step 최대10 provider, 구간별 compute wall 12h, coarse/fine 합산 최대24h로 제한한다.
- 마지막 provider state, 실제 `wn` 기반 cell transport journal, factor deck, map, area 및 raw/Psi 변화는 provider release 전에 보존한다. map 출력의 산술 midpoint는 재구성 geometry에서 마지막 provider midpoint와 지역별 density/pressure RMS≤1e-8, max≤1e-7로 비교한다. 이 일치 없이 commit하지 않는다.
- 기존 explicit midpoint가 기본값으로 남는다. `moment_integrator=self_consistent_midpoint`와 fork manifest를 함께 지정해야 새 경로가 선택된다.
- `Docs/full_moment_coupled_fork.py`는 fine `.025`의 기존 결과/q48/물리/보존 감사의 해시를 검증하고, source method와 target method가 다른 별도 fork라는 점을 기록한다. 상태 6개 payload를 serialize/import하면서 time step을 수행하지 않는다.
- native build v7 SHA-256: `7b525e6e8dcfc7fe7f53407c60f8dd800dca539a409f177f0f07cdba924d9c45`. 실제 linked symbols에는 native evolution, coupled control, midpoint decode가 포함된다. native face-stage Python/lib dependencies와 hash inventory가 build snapshot에 같이 들어 있다.

### Picard 정체 이후의 유일한 변경: safeguarded Anderson-1

- v7 coarse/fine 두 경로는 모두 10회 provider 상한에 걸려 nonlinear 상태코드 114로 종료됐다. accepted step은 0이고 시작 checkpoint는 그대로 보존됐다. fine 첫 시도(`h=.0025`, 후보 `.0275`)에서 반복 4–9 잔차가 plateau였으며 iteration 9 최대 잔차는 `W=4.07e-3`, centroid=`4.53e-3`, `S=3.63e-3`, thermal=`6.71e-3`, Psi remainder=`1.92e-2`였다. 줄어든 dt만으로 Picard 정체를 해결하지 못했다.
- 물리식, h, 질량, Re=∞, 무조향 및 수렴 임계값을 바꾸지 않고 고정점 반복만 convex depth-one Anderson map mixing으로 교체했다. map 잔차를 기존 componentwise norm과 일치하는 dimensionless vector로 기록하고, 직전 잔차와의 최소제곱 혼합을 사용한다. 과거 map의 혼합계수는 0.8 이하이며 예상 잔차 감소가 없거나 decode 검증에 실패하면 현재 Picard map으로 되돌린다.
- 새 native MPI/CUDA snapshot `Codex_full_disk_coupled_midpoint_anderson1_mpi_build_v1_20261019`의 실행 파일 SHA-256은 `a9b299fd40b75157fa6542d210cfa3f29572560e8bc88fd9b93c0e1a4d1a393e`다. Anderson 혼합계수/사용 여부/safeguard fallback과 벡터 잔차 내적을 반복별로 남긴다. 최초 패키지에서 이 실행기가 요구하는 `provenance.json:sources` 목록이 빠져 즉시 실행 실패했으나, v7 스냅샷과 일치하는 46개 파일 해시 목록을 채우고 누락된 두 원본 모듈을 복원했다. 바이너리는 바뀌지 않았다.
- 새 첫 관문 `Codex_full_disk_coupled_midpoint_anderson1_first_gate_v1_20261019`는 보존된 `.025` 상태에서 coarse `h=.005`와 fine `h=.0025`를 각각 `.03`까지만 검증한다. 최초 ID `428161` fine, `428162` coarse는 provenance `KeyError`로 1초 이내 종료됐다. Slurm 로그와 traceback을 `startup_failure_<jobid>.log`에 보존했다. 이 두 잡은 native process를 띄우지 않았고 provider 0, accepted step 0이므로 물리 관측치가 아니다. provenance를 수정한 뒤 `428213` fine, `428214` coarse로 재제출했으며 마지막 확인에서 `QOSGrpGRES`로 대기 중이었다. 각 잡은 `--gres=gpu:H200:2`, 2 MPI ranks, 8 CPU/rank, 192 GiB를 요청했고 dependency는 없다.

## 보존된 시작 상태

- 입력: `/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019/fine/checkpoints/fine_0010/accepted.npz`
- `t=.025`, 19,945 cells, source integrator explicit midpoint, target self-consistent midpoint
- checkpoint SHA-256: `465f9713ab405f19e438396abb8796510c3c51784b0c85028139248b34c28d8e`
- 기존 fine checkpoint state checks, q48 density report 및 all-region physical monitor 검증을 통과한 후 fork manifest를 만들었다. native import serialization fixture에서 IDs, mass, raw state, factors, Psi, weight-rate seed 여섯 payload가 보존됐고 RHS/accepted step은 0이었다.

## Slurm 현황 및 기록

- 첫 시도 `428093` coarse와 `428094` fine은 provider import 중 `curved_face_nodes_native`가 build snapshot에 누락되어 각각 7초/3초 후 실패했다. native runtime identity는 확인됐으며 두 잡 모두 physical provider 0, accepted step 0이다. 실패 기록은 `Codex_full_disk_coupled_midpoint_first_gate_v1_20261019`에 유지한다.
- 누락된 native face modules와 두 `.so`를 포함한 build v7로 고쳤다.
- v7 Picard 잡 `428097` coarse (`h=.005`)와 `428098` fine (`h=.0025`)는 각각 10 provider에서 정체되어 종료됐다. 둘 다 accepted step은 0이며 `.025` 시작 checkpoint가 보존됐다. fine 후보 endpoint는 `.0275`였고 coarse/fine 모두 첫 step을 수락하지 않았다. 종료 manifest는 `Codex_full_disk_coupled_midpoint_first_gate_v2_20261019/{coarse,fine}/runtime/failure/manifest.json`에 있다. 자동 dt 변경/explicit fallback은 없었다.
- Anderson-1 첫 제출 `428161/428162`는 실행기 provenance schema 오류로 시작하지 못했다. 수정 후 재제출 `428213` fine (`h=.0025`, step 12에서 `.03`)과 `428214` coarse (`h=.005`, step 6에서 `.03`)가 새 first-gate root에 대기 중이다. 최신 binary/source provenance/입력/스크립트 SHA는 `retry_submission.json`에 기록했다.
- 다음 과학 관문은 coarse와 fine 양쪽에서 `.03`까지 완료하고 nonlinear residual, factor/q48/물리 감사, 지역별 density·pressure 일치가 통과하는지 평가하는 일이다. 둘 다 통과해야만 `.05`의 나머지 구간으로 진행한다. Anderson이 실패하면 그 잔차 기록으로 다음 solver 변경을 평가하며 물리식이나 성공 기준을 완화하지 않는다.
- `428213/428214`는 `QOSGrpGRES`로 대기한다. `sacctmgr show qos normal`의 사용자당 한도는 `gres/gpu=8`이며, `sacctmgr show qos part_h200`에서 H200 파티션 QoS의 공유 한도 `GrpTRES=gres/gpu=6`을 확인했다. `h200` 파티션은 해당 `part_h200` QoS를 사용한다. 따라서 공유 파티션 한도가 대기 원인을 설명한다. 기존 잡은 취소/중복 제출하지 않는다.

## 아직 입증되지 않은 것

이 build/실행은 코드 연결과 첫 `.03` 관문을 위한 것이다. 케플러의 장시간 안정성, t6–10, 3차원 확장, 회전항이 shear error를 줄이는 인과관계는 아직 입증하지 않았다. 회전 진단에서 affine angular projection과 half-vorticity는 다른 양이며, 둘을 합산 에너지 항으로 더하지 않는다.

## 회전과 케플러 전단에 대한 물리 제약

- 케플러 전단의 속도 기울기는 대칭 변형률 성분을 포함한다. 반대칭 강체 회전만으로는 순수전단을 나타낼 수 없다. `canonical_spin_relaxation_shear_2026-10-18.json`의 정확한 아핀장 면 null test에서 강체-spin 재구성은 순수전단에 인공 접선 미끄럼과 양의 진단 열률을 만들었다 (`n=4`: `0.1056724`, `n=8`: `0.0596168`). 전체 아핀 재구성은 두 격자 모두 약 `1e-30` 이하의 열률과 `1.2e-15` 이하의 최대 미끄럼을 보였다. 따라서 강체 회전만을 기준으로 한 면 마찰/완화항을 Kepler RHS에 연결하지 않는다.
- 현재의 `C` 교차 모멘트는 반대칭 내부 각운동량과 대칭 아핀 변형을 함께 담는다. 회전 에너지 `Krot=L_internal^2/(2I)`는 이 운동에너지의 진단용 사영이며 별도 저장 에너지로 더하지 않는다. `moment_spin_ledger.h`의 도함수는 질량 변화와 완전한 moment/ALE rate를 포함한다. 단, 이 chain-rule 감사만으로 실제 moving-cell ALE 에너지 수지가 닫혔다고 보지는 않는다.
- 제조된 비균일 밀도장 압력 모멘트 감사 (`canonical_pressure_spin_work_2026-10-18.json`)는 CPU 1/4 thread와 CUDA 사이 토크 결과가 일치하고 해석적 토크 오차가 `4.14e-16` 수준임을 보인다. 이는 실제 진화가 아니며 moving-face flux 및 케플러 궤도 안정성 검증을 대신하지 않는다.
- 다음 물리 판정은 기존 Re=∞ 설정을 유지한 실제 coupled RHS에서, pressure/gravity/face/ALE transport가 orbital·affine·thermal 에너지 및 각운동량 사이에서 닫히는지 확인하는 것이다. 검증 전에는 임의의 spin 감쇠, 물리 점성, 추가 토크 또는 `Krot` 중복 가산을 넣지 않는다. 이 제약은 회전이 목표 문제를 해결했다는 주장이 아니라, 다음 긴 케플러 실행에 앞서 통과해야 할 보존 조건이다.
- 같은 해석을 현재 full-moment 입력으로 재계산했다. `G=C^T S^{-1}`, `K_aff=1/2 sum Tr(G S G^T)`, `K_rot=sum L_i^2/(2 I_i)`, `L_i=Cxy-Cyx`, `I_i=tr(S_i)`로 정의하면, t=0의 19,945셀에서 `K_aff=4.65472930943e-4`, `K_rot=4.94182826432e-5` (10.6168%)이고 t=.025에서 각각 `4.65520775610e-4`, `4.95096039864e-5` (10.6353%)이다. 따라서 이 짧은 구간에서 셀 affine KE의 약 89.4%는 최적 강체회전 사영에 포함되지 않는다. 계산 입력은 초기 `/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/initial.npz` (SHA-256 `6794d8e56c0bf8829bf95bb2ae1a633d2b690588db243ccbccdf9f8d45809c46`)와 t=.025 `/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019/fine/checkpoints/fine_0010/accepted.npz` (SHA-256 `465f9713ab405f19e438396abb8796510c3c51784b0c85028139248b34c28d8e`)다. 이는 회전만 더하는 모델보다 전체 아핀 변형 상태의 수송·응력 장부가 핵심임을 뒷받침하지만, .025까지만의 계측이므로 t6–10 안정성 증거는 아니다.
- 이미 저장된 coarse t=.025 실제 RHS 캡처의 18개 원시 모멘트 변화율로도 회전 장부를 계산했다. 재실행 없이 모든 입력 SHA를 재검증한 결과, 셀별 내부 각운동량 변화 합 `+8.84154462e-6`과 셀 중심 궤도 각운동량 변화 합 `-8.84154462e-6`이 상쇄한다. 총 원시 각운동량 변화 `-1.23e-18`은 저장된 전역 경계/중력 장부 `-2.96e-18`과 `1.73e-18` 차이다. 내부 강체회전 사영 에너지 변화율은 `-3.49e-7`; 그 안에서 토크 일은 `+2.07e-6`, 관성 변화 일은 `-2.42e-6`이다. 열에너지 변화율 `+5.23e-7`은 회전 사영 에너지 변화율과 다르다. 병진 KE `+2.19015e-5`, affine KE `-2.32974e-5`, 잔여 KE `+1.07786e-6`, 열 `+5.23055e-7`의 변화율 합은 저장된 유체 에너지 변화율과 셀별 최대 `4.24e-22`까지 일치한다. 따라서 이 순간의 내부 회전 변화는 별도 감쇠/가열항으로 해석할 수 없고, 토크·관성·병진·변형·열 사이의 실제 모멘트 교환으로 계측해야 한다. 재현 스크립트와 수치는 `Docs/audit_saved_full_moment_spin_budget.py`, `Docs/canonical_saved_full_moment_spin_budget_2026-10-19.json`에 있다. 이는 기존 단일 RHS의 instantaneous tangent 분석이며, 유한 스텝이나 장시간 회전 억제 효과를 검증하지 않는다.
- 같은 결과를 저장된 정확 곡면 지오메트리의 방사 구간별로 추가 분해했다 (`Docs/canonical_saved_full_moment_spin_budget_regions_2026-10-19.json`). inner-contact 101셀은 원시 총 각운동량 변화율 `+2.3750e-5`, 내부 스핀 변화율 `-5.5300e-7`, 강체회전 사영 에너지 변화율 `-1.9341e-6`; disk-core 11,812셀은 각각 `-2.3207e-5`, `+9.1824e-6`, `+2.5070e-6`이다. 두 영역의 큰 반대 부호 각운동량률은 내부/궤도 재분배에 대한 점검 지점을 지정하지만, 이 순간 자료만으로 face transport·압력·중력의 인과 기여를 분리하지 못한다. 다음 실제 accepted midpoint map은 같은 영역에서 `ell`, centroid orbital L, `Krot`, inertia-work, thermal 및 raw energy 변화의 midpoint-rate 대 유한차분을 함께 기록해야 한다. 보고서는 같은 저장 RHS만 읽었고 새 provider/step은 추가하지 않았다.


### 2026-10-19 scheduler dry-run 재확인

- 실제 상태 조회에서 `428213` fine과 `428214` coarse는 모두 `PENDING (QOSGrpGRES)`였고 실행시간은 0이었다.
- 기존 실행 스크립트의 임시 사본을 1 MPI rank, `--gres=gpu:H200:1`로 바꾼 `sbatch --test-only`는 두 작업 모두 10월 29일 08:55:09 시작 예상으로 응답했다. 기존 2 rank, `--gres=gpu:H200:2` 요청의 `--test-only` 예상은 10월 29일 08:55:18이었다. 테스트 전용 조회는 실제 제출이나 예약이 아니며, GPU 수를 줄여도 예측 대기시간이 개선되지 않는다는 점만 확인한다.
- 후속 queue/node 확인에서 `syn104`는 `gpu:H200:8` 중 `AllocTRES gpu=6`이었고, 두 카드는 물리적으로 비어 있었다. `part_h200` 공유 상한은 6장이므로 그 두 카드는 배정할 수 없다. `kjhan`의 `427873`, `427874`, `428252`가 각각 H200 1장을 실행 중이고, 이 첫 관문의 `428213/428214`는 각각 H200 2장 요청으로 계속 `QOSGrpGRES` 대기였다. 기존 실행은 건드리지 않았다.
- 따라서 대기 중인 2-GPU 잡은 보존하고 1-GPU 대체 잡은 제출하지 않았다. `part_h200`의 6-GPU 공유 상한에 이미 도달해 있는 동안은 같은 파티션에서 1-GPU로 줄여도 실행할 수 없다. 이 QoS 제한은 에이전트가 우회하지 않는다. 실제 재개 시점은 기존 작업 종료 후 scheduler가 정한다.
- 같은 물리 H200 노드에 등록된 `sjjlee_hgpu` 파티션도 읽기 전용 확인했다. 이 partition은 `AllowQos=sjjlee_unlim`이고 사용자 association은 `normal`만 보인다. 해당 QOS를 명시한 `sbatch --test-only`는 `allocation failure: Invalid qos specification`으로 거절됐으며 어떤 작업도 제출되지 않았다. 따라서 그 별도 파티션은 현재 계정으로 사용할 수 있는 대안으로 확인되지 않았다.


## 2026-10-19 midpoint/Anderson의 과학적 범위

현재 접근만으로 과거 `t=6–10` 벽을 넘는다고 예측할 근거는 아직 없다. `Exam/full_moment_coupled_control.c`의 Anderson-1은 같은 비선형 midpoint map 반복에서 다음 trial을 만들 뿐이고, 수렴 시 `mapped` 출력을 검증한 뒤 native evolution이 그 상태를 commit한다. 따라서 Anderson은 Picard 정체를 풀 수는 있어도 압력·중력·면/ALE 수송의 공간 closure나 그 물리 원인을 바꾸지 않는다. `Exam/full_moment_evolution.c`도 각 실행 보고서에서 `long_time_stability_verified=false`로 기록한다.

과거 t6의 entropy/pressure/work 진단은 구형 legacy/adaptive 보조기억 궤적이다. 그 보고서는 두 상태에서 순간적인 총 궤도 이탈 생성률이 양수라고 보였지만, 현재 full-moment 상태와 다른 모델·궤적이므로 새 경로의 예측값이나 개선율로 옮길 수 없다. 반대로 이 결과는 pressure/geometry가 중요할 수 있다는 경고이므로, 새 경로가 `.03`의 solver gate를 통과한 뒤 압력·밀도·위상·총에너지·각운동량을 같은 초기 조건의 시간 간격 수렴 자료로 계측해야 한다.

판정 순서는 coarse/fine 모두 `.03` 수락 및 실제 midpoint 잔차·지역 density/pressure 검증, 이어 같은 방법으로 기존 승인 endpoint까지 계속 진화, 그 뒤 동일 초기조건의 더 긴 경로를 역사적 `t=6`, `t=10`까지 비교하는 것이다. `h` 수렴, 초기 ID 기준 압력/위상/밀도 오차, orbit/affine/thermal 에너지와 각운동량 장부가 함께 닫히기 전에는 “벽을 깼다”고 결론내리지 않는다. `.03` 테스트는 이 연속 경로의 첫 실행 가능성 관문일 뿐이다.


## 2026-10-19 기존 explicit midpoint 기준선

기존 coarse (`h=.005`)와 fine (`h=.0025`) 19,945셀 궤적의 저장 checkpoint 비교는 `Docs/canonical_full_disk_completed_density_temporal_budget_2026-10-19.json`에 있다. 이 기준은 같은 모형의 coarse/fine 물리 오차 차이를 재며, 정확해 대비 시간차수나 독립 정답은 아니다. 영역별 비교 예산은 pressure/density/velocity RMS `<0.002`, phase RMS `<0.2°`; 각 checkpoint의 개별 density 예산은 RMS `<0.01`, 최대 `<0.1`이다.

### 시간 비교 결과

- `t=.01`: 8개 관심 영역 모두 coarse/fine 비교 예산 통과.
- `t=.025`: cavity density 차이 RMS `0.00200579`로 기준 `0.002`를 근소하게 초과; 나머지 보고 영역은 통과.
- `t=.05`: cavity density 차이 RMS `0.00410605`, inner-contact pressure 차이 RMS `0.00756273`으로 초과. coarse/fine 각각의 checkpoint는 개별 density RMS/max 예산을 모두 통과했다. 따라서 개별 상태의 밀도 검사는 통과했지만 시간 간격 사이의 물리 일치는 실패했다.

### 짧은 구간의 보존 장부

`t=.05` coarse/fine `physics/report.json`에서 유체에너지 변화와 `Psi` 변화가 상쇄된다. coarse는 `ΔEhyd=+9.18452489e-9`, `ΔPsi=-9.18452489e-9`, 결합 변화 `-1.94e-18`, 원점 각운동량 변화 `+1.93e-17`; fine은 각각 `+2.28859019e-9`, `-2.28859019e-9`, `+7.26e-19`, `-9.45e-18`이다. q48 `Psi` 재구성 defect는 signed `-2.59e-15` / `-2.24e-15`; P/E/L update ledger 최대 성분은 `1.94e-17` / `1.09e-17`이다.

Coarse와 fine 각각의 `fields.npz`, physics report, physical gate, density report, reference report 다섯 파일의 SHA-256을 현재 원본과 대조해 모두 일치했다. 이 전역 보존 결과는 양호하지만 지역 temporal budget 실패를 설명하거나 인과 항을 특정하지 않는다.

### 셀별 오차 위치와 열·밀도 분해

`Docs/analyze_full_disk_t05_temporal_error_concentration.py`와 `Docs/canonical_full_disk_t05_temporal_error_concentration_detail_2026-10-19.json`은 셀별 coarse/fine 오차 차이를 같은 질량으로 가중해 분해한다.

- inner-contact 101개 셀 모두의 pressure·density 오차차 제곱 기여가 `|r−0.5|<0.005` 안에 있다. ID 786의 relative pressure error는 coarse `−8.66e-3`, fine `−1.08e-4`; ID 794의 relative density error는 각각 `−1.63e-2`, `−7.95e-3`이다.
- cavity density 차이 기여의 `78.6%`는 `|r−0.5|<0.04`인 94개 셀에 있다. 반면 cavity pressure 차이는 더 넓게 퍼져 상위 118/754셀이 절반을 차지한다.
- 이상기체식 `P=(γ−1)ρu`로 상대 압력 오차차를 밀도 및 specific-thermal 기여로 항등 분해했다. 보고서는 `Docs/analyze_full_disk_t05_pressure_difference_split.py`와 `Docs/canonical_full_disk_t05_pressure_difference_split_2026-10-19.json`이다. inner-contact 기여 RMS는 밀도 `2.62e-3`, specific thermal `7.23e-3`, 두 기여 상관 `−0.976`; cavity에서는 `4.19e-3`, `4.24e-3`, 상관 `−0.996`이다. cavity net pressure 차이는 `5.55e-4`로, 큰 밀도·열 차이가 서로 상쇄된다. 두 지역의 항등식 최대 잔차는 `7e-16` 미만이다.

이 위치 분석은 접촉면 인근의 강한 시간간격 민감도를 보여주지만, 압력 복원·밀도 접촉 처리·가중치/ALE 수송 가운데 원인을 확정하지 않는다. 압력 기준만으로는 상쇄되는 열·밀도 차이를 놓칠 수 있으므로 새 accepted midpoint 계측에는 specific thermal도 함께 남긴다. 입력 hash는 analyzer가 검증했고 두 분석 모두 provider 0, accepted step 0이다.

### 새 coupled midpoint와 비교할 기준

새 coarse/fine gate는 둘 다 같은 fine `.025` checkpoint `465f9713ab405f19e438396abb8796510c3c51784b0c85028139248b34c28d8e`에서 출발한다. 실제 checkpoint 파일 hash를 확인했다. 따라서 새 gate의 차이는 서로 다른 시작 상태에서 비롯되지 않는다.

새 경로가 `.05`까지 수락되면 같은 초기 ID·지역 마스크·시간 비교 예산으로 기존 explicit 기준선과 비교한다. 차이가 남으면 pressure, density, specific thermal 및 midpoint 잔차, orbital/affine/thermal energy, angular-momentum exchange를 접촉부에서 함께 분해한다. 이 비교는 기존 결과에서 재분석했으며 새 provider나 accepted step을 추가하지 않았다.

### t=.05 가중치·셀기하와 시간간격 오차의 관계

저장된 coarse/fine accepted checkpoint를 분석하는 read-only 보고서 `Docs/analyze_full_disk_t05_weight_geometry_disagreement.py`와 `Docs/canonical_full_disk_t05_weight_geometry_disagreement_2026-10-19.json`을 추가했다. 두 checkpoint의 고정 질량 배열과 cell ID는 비트 단위로 동일하다. 각 branch에서 `actual_density = mass/area`의 최대 상대 잔차는 0이며, 따라서 이 고정질량 기준선에서 실제 밀도 차이는 셀 면적 차이와 정확히 대응한다. 가중치 비교는 additive power weight `state[2]`의 전역 공통 gauge를 제거한 뒤 수행했다.

- cavity에서 coarse/fine 가중치 차이의 질량가중 RMS는 평균 면적 대비 `2.018e-2`, 상대 면적 차이는 `4.186e-3`, 밀도 참조오차 차이는 `4.106e-3`이다. 가중치 차이와 면적 차이 상관은 `+0.394`, 가중치 차이와 밀도 참조오차 차이 상관은 `−0.395`이다.
- inner-contact에서 해당 RMS는 각각 `2.467e-3`, `1.024e-2`, `2.605e-3`이며 pressure 참조오차 차이는 `7.563e-3`이다. 가중치 차이-면적 차이 상관은 `+0.521`, 가중치 차이-밀도 참조오차 차이 상관은 `−0.945`이다. 이는 접촉부 밀도 오차와 가중치/기하의 연관을 강하게 지목하지만, pressure 차이나 오차의 발생 방향을 설명하지 않는다.
- disk-core의 가중치 차이 RMS는 `1.186e-3`, 면적 차이 `1.277e-4`, 밀도 참조오차 차이 `1.279e-4`, pressure 참조오차 차이 `1.130e-3`이다. 영역별 결과와 상위 오차 셀 ID는 JSON 보고서에 기록했다.
- coarse/fine `weight_rate_seed` 차이의 전역 평균 gauge shift는 0이지만, 공간별 차이는 남는다. 이를 평균 면적으로 나눈 질량가중 RMS는 all `4.436e-2`, cavity `7.996e-1`, inner-contact `1.127e-1`, disk-core `3.923e-2`다. 이 값은 다음 nonlinear mass-root solve의 초기 seed이지 solved/accepted `Wdot`가 아니므로 물리적 가중치 변화율로 해석하지 않는다. 분석은 provider 0, accepted step 0이다.

native 코드 경로도 확인했다. `closed_weight_root.solve`는 매 provider stage의 실제 face/Riemann 상태에 대해 모든 유체 셀의 질량 flux residual을 푸는 비선형 방정식이며, 평균 `Wdot=0` gauge와 벽 mass-flux compatibility를 검사한다. 해를 구한 뒤 `full_moment_build_candidate`가 `W_new=W_base+h Wdot`를 적용한다. 따라서 coarse/fine endpoint의 `W` 차이는 이 제약을 매 stage 풀며 시간 적분한 결과다. `W`를 단순 동결하거나 solved rate를 임의로 바꾸는 비교는 고정질량 closure를 깨므로 유효한 원인분리 실험이 아니다.

이 결과는 가중치-기하 결합을 다음 원인분리 대상으로 정하지만 인과관계는 입증하지 않는다. 유효한 후속 비교는 동일한 시작 checkpoint와 timestep에서 **가중치 시간 적분 또는 stage 평가 방식만** 바꾸되, 매 stage의 fixed-mass root, wall compatibility, 질량 잔차 및 나머지 RHS/mesh/solver 조건을 유지해야 한다. 매 accepted step의 solved `Wdot`, `ΔW`, `ΔA`, `Δrho`, specific thermal과 접촉부 pressure를 함께 기록한다. 이는 단순 coarse/fine endpoint 상관보다 강한 검증이며, 결과 없이 가중치 완화나 감쇠를 추가하지 않는다.

그 계측을 위해 `Docs/full_moment_coupled_observation.py:final_journal`에 source-level 진단을 추가했다. 다음 번 새 observer snapshot은 provider `weight_rates.bin`과 base/accepted `W`로 `ΔW−h Wdot`를 roundoff 한계와 비교하고, internal spin/orbital/total angular momentum 및 `Krot`/affine/nonspin/thermal의 유한 변화량과 실제 midpoint tangent를 함께 저장한다. 물리 RHS, weight root, 승인조건, provider 횟수는 바꾸지 않는다. Python AST 문법 검사는 통과했다. 기존 저장 t=.025 RHS로 helper 공식을 검증해 canonical spin budget과 `Krot`, `Kaffine`, internal-spin/orbital-L rates가 각각 정확히 일치했고, `Krot` rate 차이는 `1.06e-22`, 에너지 분해 최대 잔차 `4.24e-22`, 각운동량 rate 폐쇄 잔차 0이었다. 런타임 journal은 아직 없어 이 계측은 미검증이며, 현재 제출된 `428213/428214`의 이미 패키징된 실행 경로에는 들어 있지 않고 이 작업트리에서만 대기 중이다. 마지막 scheduler 조회에서도 두 잡은 `PENDING (QOSGrpGRES)`, elapsed 0이었다.

### Voronoi 시작에서 동적 Laguerre weights로의 실제 전이

초기 상태 `/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/initial.npz`의 `state[:,2]`는 모두 정확히 0이다. 즉 초기 테셀레이션은 gauge상 Voronoi다. 기존 고정질량 경로의 `t=.05` accepted coarse/fine 상태에서는 gauge-제거 weight RMS가 각각 `1.632e-5`, `1.267e-5`이며, 초기 질량 배열은 두 경로에서 bitwise 동일하게 보존됐다. 실제 stage의 closed mass-root도 마지막 coarse/fine midpoint에서 mass residual `1.07e-17`/`1.25e-17`로 수렴했다. 따라서 이 짧은 구간에서는 Voronoi 시작, 시간변화 Laguerre weights, 고정 셀 질량의 결합이 작동한다. t=.05의 공간·시간간격 오차와 기존 t=6–10 장시간 실패 때문에 안정성이나 물리 정확도까지 입증된 것은 아니다.

### 셀 내부 회전과 비강체 아핀 변형의 시간간격 민감도

`Docs/analyze_saved_full_moment_rotation_shear_timeline.py`와 `Docs/canonical_saved_full_moment_rotation_shear_timeline_2026-10-19.json`은 기존 accepted 상태에서 `Krot=ell^2/(2 tr(S))`와 전체 affine 운동에너지 내부의 비강체 잔여 `Kaffine−Krot`를 나눠 계산한다. t0의 전체 affine energy 중 강체-spin 사영은 `10.62%`; t=.025에서 coarse/fine은 `10.73%`/`10.64%`, t=.05에서는 `10.86%`/`10.66%`다. inner-contact에서는 spin 사영 비율이 t=.05 coarse/fine `52.62%`/`50.83%`이고 disk-core에서는 약 `7.5%`다. 단일 강체 회전만으로는 내부 운동을 표현하기 어렵고 전체 affine 상태가 필요하다는 기존 진단과 일치한다.

가장 관련 있는 시간간격 차이는 inner-contact다. coarse/fine의 cellwise `Krot` 차이 RMS는 fine RMS 대비 t=.025 `3.98%`에서 t=.05 `7.89%`로 증가한다. 같은 측정에서 `Kaffine−Krot` 차이는 `0.148%`와 `0.132%`로 작다. 따라서 이 구간의 접촉부 affine-energy coarse/fine 차이는 비강체 affine 잔여보다 강체-spin 사영 성분에 집중된다. 이는 회전 진단을 다음 coupled midpoint에서 따로 추적해야 할 근거이지, spin 항이 오차 원인 또는 해결책이라는 인과 증명은 아니다. 분석은 기존 저장 상태만 읽었고 새 provider/accepted step은 0이다.

같은 네 endpoint에 기존 `Docs/diagnose_full_moment_cell_rotation.py`도 실행해, 각 branch의 actual-cell steady Kepler reference에 대한 회전오차를 계산했다. 결과는 `Docs/canonical_saved_cell_rotation_{t0025_coarse,t0025_fine,t005_coarse,t005_fine}_2026-10-19/report.json`에 있다. 각 보고서는 rigid-rotation, 좌표 회전 불변성, 비등방 strain과 vorticity 구분 fixture를 통과했고 physics input hash를 확인했다. inner-contact에서 intrinsic-ell 상대 RMS 오차는 coarse/fine t=.025 `3.22%`/`0.678%`, t=.05 `6.46%`/`1.25%`다. `Krot` 오차의 reference-thermal RMS 대비 크기는 각각 `8.64%`/`1.68%`, `17.15%`/`3.09%`다. Half-affine-vorticity 상대 RMS 오차도 t=.05에서 `1.47%`/`0.341%`다. fine timestep은 접촉부 회전오차를 크게 낮추지만 t=.05까지 오차가 증가한다. 두 dt만으로 수렴 차수나 장시간 거동은 판정할 수 없다.

같은 영역에서 `ell/tr(S)`와 half-affine-vorticity 자체의 질량가중 RMS 차이는 약 `1.09–1.13`이며, 둘은 같은 회전량이 아니다. 비등방·변형 셀에서는 rigid-spin 사영을 절반 curl과 바꿔 쓰지 말고 두 값을 별도로 계측한다. 이 결과는 별도 감쇠를 추가하라는 근거가 아니라, 다음 자기일관 midpoint에서 세 번째 timestep(`h=.00125`)까지 동일 시작 상태로 비교해 시간오차 추세를 판정할 근거다. 새 진화는 실행되지 않았다.

### 최신 Slurm 확인

2026-10-19 14:36 KST의 `squeue` 조회에서 `428213` fine과 `428214` coarse는 둘 다 `PENDING (QOSGrpGRES)`, elapsed `0:00`이었다. 제출·취소·재시작은 하지 않았다. 현재 사용할 수 있는 증거는 저장 상태 분석이며, 다음 동적 midpoint 관문은 동일 job handle이 실제 실행될 때까지 기다린다.

### 회전 오차와 밀도·열·압력 오차의 위치 조건부 비교

`Docs/analyze_saved_rotation_error_conditioned.py` 및 `Docs/canonical_saved_rotation_error_conditioned_v2_2026-10-19.json`은 위의 네 accepted endpoint를 다시 읽는다. 동일 cell ID와 저장된 intrinsic-ell 오차를 대조하고, 실제 moving-cell centroid의 질량가중 저차 반지름 추세와 방위각 Fourier 모드 m=1,2를 제거한 뒤 spin 오차와 상대 밀도·specific-thermal·압력 오차의 상관을 계산했다. 각 상태에서 이상기체 압력 오차 항등식 `e_P=(1+e_rho)(1+e_u)-1`도 최대 `7.1e-16` 이내로 확인했다. 원인분석이 아니라 큰 공통 공간 패턴이 상관을 설명하는지 확인하는 진단이다.

핵심 접촉부(inner-contact, 101 cells)에서 spin–density 상관은 raw/공간 조건부 순서로 `.025 coarse +.756/+.438`, `.025 fine +.903/+.640`, `.05 coarse +.825/+.410`, `.05 fine +.905/+.600`이다. 따라서 fine 경로 두 endpoint에서는 저차 공간 패턴을 제거해도 중간 정도의 같은 부호 연관이 남지만, coarse 경로의 연관은 더 약하다. spin–specific-thermal 조건부 상관은 각각 `−.441`, `−.658`, `−.172`, `−.590`으로 밀도와 반대 부호다. spin–pressure 조건부 상관은 `−.239`, `−.653`, `+.886`, `+.287`로 시간간격과 시간에 따라 크게 변한다. 압력 오차가 밀도·열 오차의 상쇄에 민감하므로 spin 상관만으로 압력 오차의 원인을 정할 수 없다.

같은 contact region에서 density와 specific-thermal 상대오차 자체는 네 endpoint 모두 강하게 반대 부호(`−.9949`에서 `−.9998`)다. 질량가중 `RMS(e_rho+e_u)`는 `.025 coarse/fine 3.430e-4/2.601e-4`, `.05 coarse/fine 7.552e-3/1.679e-4`다. 정확한 이상기체 압력 오차는 `e_P=(1+e_rho)(1+e_u)-1`이며, 이 선형합과 작은 곱항이 `.05` fine에서의 작은 압력 오차를 만든다. 따라서 `.05` coarse의 큰 pressure 오차는 EOS 항등식 위반이 아니라 density/thermal 오차 상쇄가 fine 상태보다 덜 정확한 것이다. 같은 contact region의 spin 오차 RMS coarse/fine은 `.025`에서 `2.897e-9/6.097e-10`, `.05`에서 `5.825e-9/1.124e-9`다. finer dt가 각각 약 4.75배와 5.18배 낮추지만, fine에서도 `.025→.05` 증가가 남는다. `.05` 상대 압력 오차 RMS는 coarse `7.571e-3`, fine `1.376e-4`로 차이가 크다. 상위 10% 제곱오차 cell의 spin–density 교집합은 `.025 coarse/fine 0/8`, `.05 coarse/fine 0/9`; spin–pressure 교집합은 `0/4`, `0/2`다. 이 위치 중첩과 조건부 상관은 회전 오차와 contact density 오차를 후속 추적할 근거지만, 회전이 전단 오차를 일으키거나 억제한다는 인과 증거는 아니다. 네 endpoint 모두 t=.05 이하라 Kepler 장시간 안정성의 증거도 아니다.

### 자기일관 midpoint 첫 관문 실행 중

2026-10-19 약 14:44 KST 조회에서 fine job `428213`은 H200 노드 `syn104`에서 실행 중이고 coarse job `428214`는 `PENDING (QOSGrpGRES)`다. fine runtime에는 후보 midpoint `t=.02625`의 `passed_native_tentative_admissibility`가 저장됐다. 이 후보의 mass array는 초기와 bitwise 동일하고 P/E/L ledger 잔차는 기록된 각 bound보다 작으며 최소 열에너지도 양수다. 그러나 manifest의 `accepted`는 `false`다. 따라서 이는 nonlinear iteration 중 tentative map 상태이지 accepted step이나 `.03` gate 통과가 아니다. 다음 판단은 기존 잡의 동일 handle과 최종 residual/admissibility 기록을 확인한 뒤 한다. 추가 제출, 취소, 재시작은 하지 않았다.

약 14:47 KST의 새 로그에서는 fine job이 계속 실행 중이고 elapsed는 `6:25`; coarse `428214`는 `PENDING` 사유가 `QOSGrpGRES`에서 `Nodes required for job are DOWN, DRAINED or reserved for jobs in higher priority partitions`로 바뀌었다. coarse는 그대로 두고 재제출하지 않았다. fine의 `rhs_2.json`은 실제 moving stage `t=.02625`에서 19,945셀/4,336,408 face nodes를 대상으로 `passed_actual_full_moving_stage`를 기록했다. closed mass-root residual `5.24e-16`, generator-velocity defect 0, wall flux와 중력/전역 장부가 각 bound 이내이며 provider `Wdot` 변화는 outer iteration 0→3에서 `1.072e-3 → 3.35e-16`으로 줄었다. 별도 cell transport journal도 보존됐다. native log에는 두 MPI rank가 서로 다른 NVIDIA H200 NVL UUID를 선택했다고 기록돼 있어, `service_GPU_calls=0`은 service callback이 CUDA를 호출하지 않는다는 뜻이지 GPU 미사용이 아니다. 실제 owned CUDA/MPI stage는 두 rank 각 약 `1.3–1.5 s`, 전체 약 `4.34 s`였다. 두 rank factor audit도 `passed_owned_native_factor_audit`다.

이 자료는 실제 stage와 GPU/MPI 경로의 정상 실행 및 fixed-mass root 수렴을 증명하지만 아직 시간 상태를 수락하지 않는다. `accepted_time_steps=0`이고 `.02625` admissibility의 `accepted=false`가 마지막으로 확인된 값이다. 이 관문은 nonlinear midpoint의 map/residual·factor·물리 승인 후에야 step을 수락한다. 따라서 현 시점은 첫 accepted step이나 `.03` 통과가 아니다.

### 독립적인 다음 실험 평가

별도 계획 검토에서도 가장 판별력 높은 다음 비교는 새 spin 물리항을 도입하는 것이 아니라, 승인된 공통 fine `.025` checkpoint에서 현재 자기일관 midpoint 경로를 `.03` gate 뒤 기존 조건으로 `.05`까지 완료해 explicit fine 기준선과 비교하는 것으로 정리됐다. 비교량은 상관계수 자체가 아니라 같은 cell ID에서 `.025`부터 늘어난 spin/density/specific-thermal/pressure 오차와 실제 `Vdot`·체적-map 잔차다. 이것으로 midpoint 시간결합이 density drift와 spin error까지 줄이는지, pressure 상쇄만 바꾸는지 구별한다. 현재 fine 첫 map은 아직 accepted 전이고 coarse가 대기 중이므로 새 모델 분기나 재제출은 하지 않는다. 이 비교가 통과해도 t6–10 안정성은 별도 장시간 검증 대상이다.

fine runtime의 `coupled_0011/iteration_00.json`은 첫 fixed-point map residual을 기록한다: X `.132`, gauge W `.00635`, centroid `.132`, velocity `.0229`, S `.0176`, C `.00212`, R `1.43e-7`, thermal `.00461`, Psi remainder `.101`. componentwise 수락 한계보다 크고 `converged=false`, `accepted=false`이며 Anderson history는 아직 비어 있다. 이는 initial trial=base의 첫 residual이므로 solver 실패 판정이 아니다. 다음 residual 및 Anderson update가 수렴하는지 계속 확인한다.

약 14:49 KST에는 fine `428213`이 `RUNNING` elapsed `10:46`, coarse `428214`가 계속 `PENDING (QOSGrpGRES)`였다. fine의 두 번째 provider 후 `iteration_01.json`에서 dimensionless Anderson residual vector norm²가 `204.76 → 0.8395`로 감소했고 X `.00106`, centroid `.00266`, velocity `.000272`로 줄었다. 반면 W `.00289`, S `.01740`, Psi remainder `.05496`, thermal `.00395`는 여전히 매우 높아 componentwise 수렴 기준을 충족하지 않는다. `converged=false`, `accepted=false`; 고정밀도 map은 수렴 방향으로 움직였으나 step은 아직 수락되지 않았다. 이 관측을 성공 판정으로 대체하지 않고 같은 잡의 후속 iteration을 확인한다.

약 14:56 KST 조회에서 fine `428213`은 `RUNNING` elapsed `15:41`, coarse `428214`는 `PENDING (QOSGrpGRES)`였다. 세 번째 provider/Anderson 평가 후 residual vector norm²는 `.002619`로 더 감소했다. 각 stage의 mass-root residual은 `1.0–1.2e-17`이고, factor audit는 두 rank 모두 통과했다. 세 tentative endpoint 모두 질량 배열 bitwise 동일, 양의 thermal, ledger bound 통과를 기록했지만 수락되지는 않았다. 현재 iteration 2의 component maxima는 W `.00420`, centroid `.00442`, S `.00350`, thermal `.00615`, Psi remainder `.01460`이고 strict componentwise 기준보다 여전히 크다. norm 감소는 크지만 일부 성분은 직전 반복보다 커져, 전체 맵 수렴이나 step 수락으로 판정하지 않는다.

약 14:59 KST 조회에서 fine은 `RUNNING` elapsed `18:34`, coarse는 계속 `PENDING (QOSGrpGRES)`였다. 세 번째 provider stage `rhs_6.json`과 rank 0/1 factor audit는 통과했고, native CUDA/MPI 로그에 해당 epoch의 H200 실행이 보인다. 네 번째 provider stage가 시작됐지만 아직 `iteration_03.json`은 생성되지 않았다. iteration 2 residual은 여전히 미수렴이며 accepted step은 0이다. 현재 실행 handle을 유지하고 후속 residual을 확인한다.

2026-10-19 약 15:01 KST fine `428213`은 `RUNNING` elapsed `20:52`; coarse `428214`는 `PENDING (QOSGrpGRES)`다. 네 번째 provider와 factor audit도 통과했고 질량은 각 tentative map에서 초기와 bitwise 동일하다. 그러나 iteration 3의 residual norm²는 직전 `.002619`에서 `.003139`로 다시 증가했다. W `.00388`, centroid `.00434`, S `.00348`, thermal `.00643`, Psi remainder `.01871`로 정체/비단조가 보이며 `converged=false`, `accepted=false`다. 이는 아직 10-provider step budget의 최종 결론은 아니지만, 단순 Anderson 혼합으로 수렴이 보장되지 않는다는 경고다. 다음 map 잔차와 provider budget/terminal manifest를 확인하고, 실패하면 이 기록을 바탕으로 수학식이 아닌 nonlinear solver 전략을 재평가한다.

기존 Picard fine v2의 저장 `iteration_*.json`과 같은 시작 상태/`h=.0025`에서 비교하면, 새 Anderson-1의 iteration 3도 이미 옛 Picard 반복 4–9의 plateau에 매우 가깝다. 옛 plateau 대비 W/centroid/S/thermal/Psi 차이는 대략 3–5% 수준뿐이다. 전체 Anderson norm²는 크게 내려갔지만 이는 초기에 큰 X·velocity 성분이 줄어든 효과가 포함되며, 물리적으로 지배적인 다섯 component residual은 여전히 `1e-3–2e-2`다. 따라서 새 mixer가 solver bottleneck을 풀었다고 해석할 수 없다. 다만 새 실행은 아직 provider 4까지이며, terminal residual은 아니다. step budget 끝까지 확인한 후에만 정체 여부를 확정한다.

2026-10-19 약 15:06 KST에는 fine 잡이 elapsed `26:00`으로 계속 실행 중이고 coarse는 `QOSGrpGRES` 대기다. iteration 4에서 norm²가 `.0004624`로 내려왔으며 W `.00267`, centroid `.00237`, S `.00164`, thermal `.00333`, Psi remainder `.00571`까지 낮아졌다. iteration 3의 비단조 증가 뒤 다시 감소한 실제 진전이다. 그러나 모든 핵심 최대 residual은 여전히 해당 `1e-10`/`1e-8` 한계보다 여러 자릿수 크고 `accepted=false`; provider 6에 해당하는 epoch 12 tentative 상태도 질량 bitwise 동일/양의 열에너지이나 미수락이다. solver는 아직 진행 중이며 convergence 여부는 열려 있다.

fine provider 6 이후 iteration 5에서도 수렴 판정은 실패했다. norm²는 `.0004624 → .0009733`으로 다시 증가했고 W `.00395`, centroid `.00412`, S `.00312`, thermal `.00588`, Psi remainder `.01357`가 iteration 3 쪽으로 되돌아왔다. 특히 낮은 잔차 iteration 2/4와 높은 잔차 iteration 3/5가 교대하는 양상이 나타난다. 저장된 residual dot product에 구현된 depth-one 공식으로 mixing weight를 다시 계산하면 iteration 3/4/5에 `0.524/0.068/0.611`이며, 혼합에도 residual이 왕복한다. 이는 현재 map에서 oscillatory solver mode의 가능성을 보여 주지만 provider limit 전이므로 terminal failure는 아니다. Stage mass root residual `8.59e-15`, rank factor audits 및 물리 admissibility는 통과하고 질량 배열은 bitwise 고정이다. `accepted=false`, coarse `428214`는 여전히 `QOSGrpGRES` pending이다.

fine의 provider 7 iteration 6에서는 norm²가 다시 `.0001415`로 감소했고 W `.00163`, centroid `.00161`, S `.00115`, thermal `.00228`, Psi remainder `.00408`까지 낮아졌다. 따라서 residual envelope는 oscillatory지만 감소하고 있다. 다만 componentwise 수락 한계보다 아직 5–7 orders 크고 `converged=false`, `accepted=false`다. rank factor audit와 실제 moving-stage 질량 root는 계속 통과했다. 2026-10-19 약 15:15 KST의 `squeue`에서는 fine `428213`이 `RUNNING` elapsed `34:55`, coarse `428214`는 `PENDING (QOSGrpGRES)`였다. 이후 provider 8–10의 residual과 terminal report가 현재 step의 판정을 결정한다.

### Fine Anderson-1 관문 종료: 10회 반복으로 미수렴

2026-10-19 15:27 KST에 fine `428213`은 provider budget 10회 후 Slurm `FAILED`, exit code `114`, elapsed `47:20`으로 종료됐다. terminal `control.json`은 `iterations=10`, `provider_calls=10`, `converged=false`, `final_checks_passed=false`, `accepted=false`, compute elapsed `2835.89 s`를 기록한다. 15:27의 마지막 iteration 9에서 dimensionless residual norm²는 `2.23747e-4`였고, 최대 성분은 gauge W `1.62723e-3`, centroid `1.59947e-3`, velocity `8.12950e-5`, S `1.20817e-3`, C `6.18457e-5`, R `3.75441e-9`, thermal `2.37725e-3`, Psi remainder `5.44534e-3`다. 종속 Psi 잔차는 0이다. 모든 핵심 map 성분이 요구 기준보다 여전히 크게 초과하며 어떤 time step도 수락되지 않았다.

잔차 norm²는 iteration 7 `1.10361e-4`에서 iteration 8 `2.97948e-4`로 반등했고 iteration 9에서 `2.23747e-4`로 일부 낮아졌다. 전체 10회 중 마지막 9번에 Anderson mixing이 사용됐고 계수는 iteration 9에서 `0.4391`; safeguard fallback은 없었다. 따라서 이 첫 midpoint의 병목은 provider 실행 실패나 mass root 실패가 아니라, 현재 depth-one convex Anderson fixed-point 반복의 map 수렴 실패다. 추가 단순 반복이나 작은 timestep만으로 해결됐다고 볼 근거는 아직 없다.

모든 tentative 후보에서 질량 배열은 초기와 bitwise 동일하고 thermal은 양수이며 P/E/L ledger와 두 rank native factor audit는 통과했다. 마지막 tentative state도 `passed_native_tentative_admissibility`이나 `accepted=false`다. `runtime/failure/manifest.json`은 이전 accepted `.025` checkpoint가 보존됐고 retry와 timestep 변경이 없었음을 확인한다. native 로그는 실제 두 MPI rank가 H200 NVL을 사용했으며 provider 10을 포함한 실제 GPU/MPI stage를 실행했음을 확인한다. 즉 이번 결과는 짧은 계산이나 CPU 대체 실행이 아니라, H200에서 수행된 완전한 19,945셀 fine 경로의 nonlinear solver 실패다.

같은 첫 관문의 기존 coarse 작업 `428214`는 이 확인 직후 H200에서 시작해 실행 중이었다. 새 job 제출·취소·재시작은 하지 않는다. coarse의 terminal residual을 같은 형식으로 읽어 timestep 민감도와 fine 실패 재현성을 판별한 뒤 다음 solver 변경을 정한다. 물리 RHS, Re=∞, 질량 제약, 수락 기준은 유지한다. 이 결과는 장시간 Kepler 안정성에 대한 진전이나 반증이 아니라, 그 검증으로 가기 전 현재 self-consistent midpoint solver가 첫 `.0275` 후보조차 수락하지 못했다는 사실을 확정한다.

### Coarse Anderson-1 관문 종료 및 fine/coarse 대조

2026-10-19 16:15 KST에 coarse `428214`도 10회 provider 뒤 Slurm `FAILED`, exit code `114`, elapsed `47:30`으로 끝났다. 마지막 iteration 9의 norm²는 `2.28283e-3`; 최대 성분은 gauge W `3.74289e-3`, centroid `4.33306e-3`, velocity `2.35990e-4`, S `3.33433e-3`, C `1.87758e-4`, R `1.19925e-8`, thermal `6.38376e-3`, Psi remainder `1.09540e-2`다. `control.json`은 `converged=false`, `final_checks_passed=false`, `accepted=false`, provider 10, compute elapsed `2845.72 s`를 기록한다. 마지막 세 residual은 감소–반등–감소했지만 모든 주요 성분은 수락 기준보다 여전히 크다.

같은 `.025` 시작 checkpoint, 19,945셀, 동일 물리/수락조건에서 비교하면 coarse 마지막 norm²는 fine의 `10.20`배다. W/centroid/velocity/S/C/R/thermal/Psi 최대 residual은 coarse가 fine보다 각각 `2.30/2.71/2.90/2.76/3.04/3.19/2.69/2.01`배 크다. 따라서 `h=.005`가 `h=.0025`보다 nonlinear solve에 더 어렵지만, fine timestep도 10 provider 안에 수렴하지 못해 단순 dt 절반으로는 solver 병목이 해결되지 않았다.

두 branch 모두 actual moving-stage mass root, 두 rank native factor audit, 질량 bitwise 보존, 양의 thermal 및 P/E/L ledger를 통과했고 실제 2-rank H200 NVL 계산이었다. 그래도 한 step도 수락되지 않았으며 두 failure manifest 모두 기존 accepted `.025` 시작점을 보존하고 `retry=false`, `timestep_reduced=false`다. 그러므로 이번 이중 관문은 물리 RHS가 잘못되었다는 증거가 아니라, 현 depth-one Anderson/Picard fixed-point iteration이 첫 step의 map 잔차를 수락 수준까지 줄이지 못한다는 재현 가능한 solver 실패다.

다음 코딩 작업은 수락 기준이나 물리식을 바꾸지 않고 nonlinear solve만 개선해야 한다. 잔차 이력이 비단조이고 coarse에서는 일부 반복의 Anderson mixing이 예측감소 검사를 통과하지 못해 pure map으로 돌아갔으므로, 단순히 provider 한도를 늘리는 것보다 rank-revealing QR을 쓰는 multi-history Anderson(작은 history 깊이부터) 또는 safeguarded multisecant/Broyden을 먼저 평가하는 편이 판별력이 높다. 혼합 후보는 같은 고정질량/양의 covariance 상태공간과 positivity guard를 유지하고, 비교 관문은 이 보존된 `.025` checkpoint에서 fine/coarse를 다시 실행해 iteration별 동일 잔차 성분·step 수락 여부를 확인해야 한다. 성공해도 다음은 `.05` 비교이며, 최종 목표인 long-term Kepler와 역사적 `t=6–10` 검증은 여전히 남아 있다. 현재 지시는 코드 변경이 driver에 있다는 사용자 지시에 따라 여기서는 기존 수치결과와 다음 solver 요건을 기록하는 데 그친다. 구현용 수식·데이터구조·수락불변조건은 `Docs/next_nonlinear_solver_design_2026-10-19.md`에 적었다.

### Nonlinear residual 좌표계 감사 — 2026-10-19

`Exam/full_moment_coupled_state.c`와 `Exam/full_moment_coupled_runtime.inc`를 대조한 결과, Anderson 벡터의 Ψ remainder 성분은 각 provider의 `origin/scale/projection/sigma/dependent`로 재정규화되는데, runtime은 그 벡터와 직전 provider의 벡터를 직접 내적한다. 따라서 metric이 실제로 변했다면 Gram 자료는 공통 좌표계에 있지 않다. 기존 반복별 `potential_norms` payload는 `/dev/shm`의 provider root가 `FMP_RELEASE`에서 삭제되며, fine `428213`의 영구 산출물에도 없어서 변화량을 사후 측정할 수 없다. 이 사실은 기존 실패를 설명한다고 단정하지 않지만, 기존 norm² 곡선을 하나의 고정 목적함수로 해석하거나 현재 Type-II 제안에 그대로 재사용할 수 없게 한다.

Type-II는 accepted base의 첫 provider metric을 timestep 안에서 고정해 solver용 잔차를 만들고, 현재 provider metric으로 계산한 실제 componentwise acceptance gate는 그대로 별도로 적용해야 한다. provider dependent mask/rank가 바뀌면 기존 이력을 폐기하고 metric epoch를 재시작한다. 가속 대상은 residual `21*n` 성분에 대응하는 per-cell raw state 20개와 `Psi`만으로 한정한다. 공통 weight gauge를 candidate에 투영하고 IDs/mass는 고정하며 derived factors를 재구성한다. 전역 ledger는 외삽하지 않고 최신 map 값으로 두며 final acceptance는 실제 새 map의 ledger로 판단한다. 코드 경로 검색에서 `weight_rate_seed`는 trial midpoint에 보존된 뒤 convex mixer에서만 섞이고, native factor가 실제 `weight_rates`를 map output seed로 쓴다. audited RHS/provider는 trial seed를 읽지 않으므로 현재 fixed-point에서는 외삽 대상에서 제외하고 최신 map seed를 전달한다. 구체적인 구현 요건은 `Docs/next_nonlinear_solver_design_2026-10-19.md`에 정리했다. 코어 소스와 시뮬레이션은 이 감사에서 수정하거나 다시 실행하지 않았다.

반복 JSON의 최대 residual cell index도 accepted `.025` geometry와 대조했다. iteration 2–9의 144개 성분별 최대 위치 중 127개가 `r=.5` interface 안쪽 cavity, 10개가 interface-crossing inner-contact, 7개가 inner disk다. 가장 빈번한 ID 720과 747은 각각 radial support `[.449593,.485216]`, `[.452422,.485158]`인 cavity-side 이웃이다. 이 두 시작 상태의 density fit은 condition `319/293`, dual residual `2.78e-16/2.91e-16`으로 통과했으므로 base-fit rank 실패는 보이지 않는다. moving midpoint fit은 실행 산출물에 없어 원인은 미정이다. 따라서 다음 solver 계측은 이 cavity-side shell의 actual-provider fit/rank와 잔차를 함께 보존해야 한다.

최종 Slurm 상태를 다시 조회한 결과 `428213`은 `FAILED`, exit `114`, 종료 `2026-10-19 15:27:59`; `428214`는 `FAILED`, exit `114`, 종료 `2026-10-19 16:15:29`다. 두 작업 모두 terminal이며 지금 실행 중인 동일 job은 없다. 새 작업은 제출하지 않았다.

### 기존 수락 경로의 연장 판정 — 2026-10-19

비선형 midpoint 실험과 별도로, 기존 continuous explicit-midpoint 경로는 fine `h=.0025`로 `.05`까지 20 step을 수락했고 각 저장 checkpoint의 물리 monitor가 통과했다. fine inner-contact pressure RMS는 `t=.01/.025/.05`에서 `1.64e-4/2.52e-4/1.38e-4`, 마지막 velocity-relative RMS는 `3.10e-4`, phase RMS는 `0.0176°`다. 압력 오차는 `.05`까지 증가하지 않았지만 phase/velocity 지표는 증가했다. Coarse `h=.005`의 `.05` inner-contact pressure RMS는 `7.57e-3`으로 fine보다 약 55배 크다. 그러므로 오래된 조건부 coarse continuation은 장기 정확도 판정용으로 쓰지 않는다.

이전 coordinator manifest의 `conditional_end=1.0`은 fine 경로를 연장하지 않았다. 원본 `run_full_disk_continuous_bundle.py --mode continue`는 오직 `coarse_0010`의 `.05` checkpoint를 읽어 `h=.005`, step 10→200으로 진행한다. 이는 10월 18일 hold/cleanup 기록이 가리키는 obsolete continuation이다. 그 프로세스는 종료됐고, 새 guarded continuation으로 갈아타야 한다.

목표에 맞는 후속 경로는 fine `fine_0020/accepted.npz` (`t=.05`)에서 동일 `h=.0025`와 동일 physics/gates로 restart해 `.1/.25/.5/.75/1.0`를 저장·검증하는 것이다. 이 restart 입력의 SHA-256은 `b6570af3551041958eb5d85bb038d0c88c99e39cefe9f79c0bf2c14d271a737e`이며 `latest_fine.npz`와 byte-identical이다. manifest는 20 accepted steps/40 provider calls를 기록한다. checkpoint physical gate가 통과했고 mass는 initial과 bitwise 동일, minimum thermal은 `2.025e-8`, P/E/L residual은 각각 기록된 bounds보다 작다. fine 원본 `.0→.05`는 20 step에 15,945초 걸렸다. 같은 평균속도를 단순 적용하면 `.05→1.0`은 약 84시간이므로, 이 값은 재측정 전의 planning estimate이며 48시간 Slurm 제한 내 여러 accepted-checkpoint 구간으로 나눠야 한다. 새 restart runner는 현재 코드에 없다. Driver가 fine restart, 대상 시각별 physical/density/rotation/ledger gates, immutable input hash 및 실패 시 accepted checkpoint 보존을 구현한 뒤, H200 `--gres=gpu:H200:1` 짧은 continuation으로 backend를 확인하고 긴 구간을 제출한다. 각 구간은 실패 시 자동 retry/dt 변경 없이 멈춘다. 이 경로가 `.05` 이산화 오차와 contact growth를 통과해야만 역사적 `t=6–10` 목표로 이어간다.

### 최종 관문 상태 갱신 — 2026-10-19

과거의 위 `PENDING`/`RUNNING` 문장은 당시 진행 기록으로 보존한다. 현재 terminal 상태는 `sacct`에서 확인했다: `428213 kepler-aa1-fine`은 `FAILED 114:0`, 2026-10-19 14:40:39–15:27:59 (47:20); `428214 kepler-aa1-coarse`는 `FAILED 114:0`, 15:27:59–16:15:29 (47:30). 두 manifest 모두 provider budget 10회에서 `converged=false`, `accepted=false`, accepted step 0으로 끝났고 기존 `.025` checkpoint를 보존했다. 현재 진행 중인 이 coupled midpoint 작업은 없다.

현재 `Exam/full_moment_coupled_control.c`와 `Exam/full_moment_coupled_runtime.inc`에는 Type-II가 아직 없고 Anderson-1의 직전 residual 한 개만 저장한다. 다음 구현 요건과 기존 실패 분석은 `Docs/next_nonlinear_solver_design_2026-10-19.md`에 정리했다. 현 단계에서 같은 solver로 잡을 다시 제출하거나 단순 반복 수를 늘리는 것은 다음 판별력 있는 조치가 아니다.

한편 기존 explicit H200 경로의 `.05` 결과는 coarse/fine 비교와 restart 등가성 검증이 끝났지만 density temporal 기준을 통과하지 못했다. 이를 이유로 `.1` 또는 t1 continuation을 시작하지 않는다. 앞의 2026-10-19 계획 절에 있는 장기 continuation은 density/temporal gate 실패 및 새 Type-II 관문이 해결되기 전의 조건부 계획으로 읽어야 한다.

### 같은 accepted 시작점의 세 번째 explicit timestep 분기 — 실행 중

두 timestep만으로 contact의 density/spin 오차 수렴을 판정할 수 없어, 기존 fine accepted checkpoint `fine_0010/accepted.npz` (`t=.025`, SHA-256 `465f9713ab405f19e438396abb8796510c3c51784b0c85028139248b34c28d8e`)에서 같은 native H200 MPI/CUDA binary로 두 개의 독립 fork를 시작했다. Python fork builder는 계산을 수행하지 않고, 원본 manifest·실제-cell/q48 reference·density preparation·physical monitor·state checks의 hash를 묶었다. IDs/mass가 원래 IC와 정확히 같은지와 열 양수·P/E/L·독립 q48 조건을 확인했으며, typed native import 전후에 six payload가 보존되도록 native verifier를 통과시켰다.

- `428541`: 동일 H200 native `h=.0025` 대조 branch, `.025→.05` 10 step, 1 rank/1 GPU.
- `428542`: 세 번째 timestep `h=.00125`, `.025→.05` 20 step, 1 rank/1 GPU.

두 작업은 같은 입력을 읽는 독립 branch이며 dependency는 없다. 첫 제출 `428538/428537`은 wrapper의 경로 오기 때문에 각각 0/1초 만에 종료됐다. 그 로그와 출력은 보존했고 native process/provider/RHS는 시작되지 않았다. 경로를 v10에서 실제로 사용했던 runtime identity wrapper로 고친 새 v2 run directory에서 위의 `428541/428542`를 제출했다. 최초 조회에서는 두 잡 모두 `RUNNING` (syn104, elapsed 00:24)이었다. 완료 후에는 `.05`의 density·specific thermal·pressure·intrinsic-spin을 같은 cell ID와 같은 `.025` 시작점에서 비교한다. 새 분기 역시 `.05`까지의 시간수렴 결과일 뿐 장시간 안정성 증명은 아니다.

### Native GPU 첫 provider 확인

최신 로그에서 두 branch 모두 정확한 v10 binary SHA `d9c79ec0f868ebb3f39f0e95ec3777465d0055b4fab079e076ae59e4cf619760`로 실행 중이고 실제 `NVIDIA H200 NVL` UUID가 기록됐다. `428541`은 provider 1에서 19,945개 셀/91,749,864 gravity nodes를 처리했고 owned GPU stage 2.178 s, MPI/GPU 5.073 s, 전체 RHS 322.965 s였다. `428542`는 같은 크기에서 2.105 s, 4.972 s, 320.086 s였다. 첫 midpoint half-state 후보는 각각 `t=.02625`와 `.025625`에 생성됐지만 `accepted=false`이며 정상적인 중간 상태다. 아직 어느 branch도 full step을 commit하지 않았으므로 accepted-step 수는 0으로 유지한다.

그 뒤 각 branch의 두 번째 RHS와 full factor/endpoint gates가 통과해 첫 step을 commit했다. `428541`은 `accepted_0011` (`t=.0275`, `h=.0025`, 2 providers, step wall 828.907 s), `428542`는 `accepted_0021` (`t=.02625`, `h=.00125`, 2 providers, 826.772 s)다. 양쪽 모두 native report에서 `fixed_mass=true`, `Re_infinity=true`; admissibility report에서 mass bitwise initial, positive thermal, P/E/L ledger bounds 통과를 확인했다. 최저 thermal은 각각 `2.0246745e-8`/`2.0246758e-8`이다. 상세 hash와 값은 `/gpfs/kjhan/LagForce/Codex_full_disk_explicit_fine_timestep_refinement_fork_v2_20261019/first_accepted_steps.json`에 보존했다. 마지막 Slurm poll에서는 두 작업 모두 계속 `RUNNING`; 각 branch는 한 step만 통과했고 `.05` endpoint gate는 아직 남아 있다.

### 진행 갱신 — 세 번째 accepted step

이후 Slurm elapsed `47:07` 시점에도 두 작업은 `syn104`에서 `RUNNING`이었다. `428541`은 source `.025` 이후 세 step을 accepted하여 `step=13`, `t=.0325`, `h=.0025`에 도달했다. `428542`도 세 step을 accepted하여 `step=23`, `t=.02875`, `h=.00125`에 도달했다. 최신 step wall은 각각 약 `870.5 s`, `867.0 s`; 총 provider 수는 branch마다 6이다. 두 report 모두 19,945 cells, `fixed_mass=true`, `Re_infinity=true`, `accepted_commit=true`, `long_time_stability_verified=false`를 기록한다. 세 번째 step의 admissibility에서 질량 bitwise 초기 동일, 양의 thermal, P/E/L 및 combined-energy roundoff 장부 통과를 확인했다. 실제 GPU UUID와 v10 binary SHA는 앞 기록과 동일하다. 원본 `.025` checkpoint와 과거 출력은 그대로 보존되어 있다.

기존 coarse/fine 시간 비교의 정확한 실패 위치도 재확인했다. `canonical_full_disk_completed_density_temporal_budget_2026-10-19.json`에서 `h=.005` 대 `.0025` temporal budget은 `t=.025` cavity density RMS `0.00200579`, `t=.05` cavity density RMS `0.00410605`, inner-contact density RMS `0.00260498`, inner-contact pressure RMS `0.00756273`로 모두 해당 `0.002` 기준을 넘었다. 따라서 이번 `.025→.05`의 `.0025` 대 `.00125` 비교가 그 부족을 줄이는지 확인한다. 현재 단계는 각 스텝의 admissibility일 뿐, 아직 density/reference physical checkpoint audit이나 temporal comparison은 완료되지 않았다.

완료 시 두 분기 모두의 `.05` accepted state와 source/manifest hashes를 확인한다. `.0025` 결과는 기존 uninterrupted fine `.05` accepted state와 restart identity를 대조하고, `.00125` 결과는 actual-cell reference/q48 density 및 physical monitor를 새로 통과시킨 뒤 `.0025`와 같은 ID·질량 기준으로 비교한다. density, pressure, velocity, phase 기준을 모두 통과하기 전에는 더 긴 continuation을 시작하지 않는다. 이 gate를 통과해도 t=6–10 장기 안정성의 증명은 아니며, 기존 adaptive-QR Voronoi 장기 실행은 t≈10에서 이미 physical stability fail로 기록되어 있다.

v10 실행 코드를 확인한 결과, `checkpoint()`는 `step==r.steps`인 최종 스텝도 무조건 checkpoint로 지정한다. 그래서 manifest의 `moment_validation=0`이어도 각 `.05` endpoint에서 `configure(...,check)`를 통해 independent q48 density, actual-cell steady reference, physical monitor가 수행된다. 이 monitor가 실패하면 native code는 해당 후보를 `accepted_XXXX`로 commit하지 않는다. 따라서 endpoint 완료 시 `candidate_audit_0020`/`candidate_audit_0040`의 `density`, `reference`, `physics`, `physical_gate.json`을 먼저 보고, 통과한 경우에만 accepted state 비교를 진행한다. 현재로선 별도 대규모 geometry/density postprocess job을 중복 제출할 필요가 없다.

회전항의 즉시 RHS 장부도 보존된 같은 `t=.025` 상태에서 별도로 확인했다: intrinsic-spin rate와 centroid-orbital rate의 전역 합은 `-1.23e-18`, saved angular-momentum ledger rate는 `-2.96e-18`이고 차이는 `1.73e-18`; rotational-energy chain-rule 최대 잔차는 `1.29e-26`이다 (`canonical_saved_full_moment_spin_budget_2026-10-19.json`). 이는 한 순간의 내부 spin·토크·관성·열 분해가 대수적으로 맞는다는 증거다. finite-step 회전 에너지 전달이나 t=6–10 궤도 안정성의 증거는 아니므로, endpoint 물리 기준과 이후 장기 진화가 여전히 필요하다.

새 `h=.0025` native 분기의 source `.025`와 accepted `t=.0325` (`accepted_0013`)를 같은 ID/질량으로 대조해 finite-step spin/orbit 장부도 검사했다. `intrinsic_ell=C_xy-C_yx`, `origin_L=D_xy-D_yx`, centroid orbital term `q×p/m`로 계산하면 내부 spin 합은 `+7.80660768e-8`, centroid orbital 합은 `-7.80660767e-8` 변했고, 전체 origin angular momentum 변화는 `5.21e-18`, 두 변화의 closure 잔차는 `4.10e-17`이다. 셀별 spin 변화 절댓값 합은 `2.16838e-6`이라 내부 spin 자유도가 실제로 움직이며 궤도 각운동량과 상쇄된다. 이 결과는 유한 스텝의 보존 전달을 확인하지만, 연속 Kepler 해에 대한 spin 오차나 전단 오차 감소, 장기 안정성을 아직 입증하지 않는다.

같은 accepted state로 사용자의 Voronoi→Laguerre/고정질량 조건도 확인했다. 원 IC `t=0`에서 모든 power weight는 정확히 0이고, 이어진 같은 물리 경로의 `t=.0325`에서는 19,945개 셀 모두 weight가 nonzero (`max |w|=7.22222e-5`)다. initial과 final의 ID 및 mass 배열은 bitwise 동일하고, manifest 조건은 무조향·`Re=∞`다. 따라서 이 구현은 Voronoi 초기 상태에서 가중치를 진화시키면서 셀 질량을 보존하고 있다. 이는 원하는 허용 제약을 구현했다는 증거이지, 그 가중치 응답이 물리 오차를 줄이거나 t=6–10을 안정화한다는 증거는 아니다.

matched-time raw spin 비교도 가능해졌다: `h=.0025 accepted_0012`와 `h=.00125 accepted_0024`는 모두 `t=.03`이며 IDs/mass가 bitwise 같다. fine 분기의 intrinsic-spin mass-weighted RMS 대비 분기 차이는 `2.9495e-4`, 최대 셀 차이는 `2.16012e-10`이다. 두 분기의 global spin 합 차이 `-3.37780e-8`은 centroid-orbital 합 차이 `+3.37780e-8`과 상쇄되어 global origin-L 차이는 `-1.38e-17`이다. `Krot=.5 ell²/trace(S)` 회전 투영 에너지의 분기 RMS 차이는 fine RMS의 `1.0615e-3`; 이는 affine kinetic energy 안의 진단 성분이며 별도 저장 에너지가 아니다. specific-thermal 분기 차이는 fine RMS의 `7.5701e-5`다. weight의 두 branch 간 차이는 max `1.66196e-6`, mass-weighted RMS `3.11771e-8`이다. 초기 회전/weight 상태가 작은 h 변화에 어떻게 반응하는지 보는 진단으로는 유망하지만, 이 단계에는 independently audited density/thermal/pressure/phase가 없어 전체 timestep budget이나 수렴 차수를 판정하지 않는다.

Slurm elapsed `1:16:43`의 최신 확인에서 두 original branch 모두 `RUNNING`이며 각각 다섯 step을 accepted했다. `428541`은 `t=.0375`, `h=.0025`, cumulative 10 providers, step wall `873.79 s`; `428542`는 `t=.03125`, `h=.00125`, cumulative 10 providers, step wall `870.76 s`다. 다섯 번째 후보의 admissibility에서 fixed mass bitwise, 양의 열, P/E/L 및 combined-energy roundoff 장부 통과를 확인했다. 둘 다 `.05` endpoint full q48/reference/physical monitor에는 아직 미도달이다. 별도 2-rank `428614`는 계속 `PENDING (QOSGrpGRES)`다.

실행속도는 장기 목표에 직접 영향을 준다. 현재 첫 네 accepted step의 평균 wall은 `.0025`에서 약 `862 s/step`, `.00125`에서 약 `850 s/step`이다. 이 속도를 고정 적용하면 1 rank 기준 `.05→6`에 `h=.0025` 약 24일, `.05→10` 약 40일이며 `h=.00125`는 대략 두 배다. 이는 단순 선형 추정이며 장시간 실제 비용을 보증하지 않는다. `.05` 물리/temporal gate를 먼저 끝내고, 이후 장기 확장을 위해 같은 physics와 엄격한 gate를 유지한 2-rank H200 scaling 및 단계별 허용 `h` 증가를 검증해야 한다. `h=.005`는 이전 `.05` density budget이 실패했으므로 검증 없이 큰 step으로 바꾸지 않는다.

장기 계산을 위한 2-rank scaling 분기를 같은 `.025` native source와 `h=.0025`로 추가 준비해 Slurm `428614`로 제출했다. 새 output은 `run_h0025_r2`, target은 같은 `.05`이며 요청은 `--gres=gpu:H200:2`, 16 CPUs, 256G, 6h다. 동일 v10 binary SHA 및 기존 immutable restart manifest를 쓰지만, 새로운 runtime/output 디렉터리와 MPI run token을 사용한다. independent branch라 기존 두 계산과 dependency는 없다. 제출 직후 상태는 `PENDING (QOSGrpGRES)`여서 아직 계산이나 GPU identity 증거는 없다. 시작 후 두 rank의 서로 다른 실제 H200 UUID, 동일 `.05` q48/reference/physical gate, rank-1 및 기존 fine `.05`와의 결과 차이를 확인한다. 이 benchmark는 정확도 분기를 대체하지 않으며, 같은 physics로 병렬 실행이 장기 연장 시간을 줄이는지 판정하기 위한 것이다.

기존 MPI2 경로의 authoritative restart-equivalence 보고서도 확인했다: `/gpfs/kjhan/LagForce/Codex_full_disk_native_restart_equivalence_v10_20261019/report.json`은 `h=.005`, `t=.01→.015` 한 스텝에서 두 rank가 uninterrupted 경로와 IDs/mass/state/factors/Psi/weight seed/time/masks 전부 exact equality, 입력 payload exact copy, initial/final q48 physical gates pass를 기록한다. `runtime_identity.json`은 실제 서로 다른 두 H200 UUID와 동일 v10 binary/operator를 확인했다. 이것은 MPI2 restart 및 실행 경로의 증거이지만 현재 `h=.0025`의 scaling·cross-rank 물리오차·장기 안정성 증거는 아니므로 `428614` 비교는 여전히 유효하다.

### Slurm 진행 재확인 — 2026-10-19 18:29 KST

Slurm 로그인 호스트 `syntax`에서 `sacct`로 조회했다. `428541`과 `428542`는 `syn104`의 H200 1장씩으로 계속 `RUNNING` (elapsed `1:22:41`)이고, `428614`는 2장 요청을 유지한 채 `PENDING (QOSGrpGRES)`다. 두 원본 런타임에는 아직 새 accepted checkpoint가 없다: `h=.0025`는 `step=15,t=.0375`; `h=.00125`는 `step=25,t=.03125`다.

다만 두 런타임 모두 `rhs_11.json`이 생겼다. 이는 각 다섯 번째 수락 step 뒤 다음 step의 첫 실제 moving-stage 평가이며, 질량 root residual은 약 `1.1e-17`, 전역 ledger residual은 저장 허용치 안이다. 이 provider 성공은 full step 수락이나 `.05` endpoint density/reference/physical gate 통과가 아니다. 작업을 재시작하거나 중복 제출하지 않고 다음 accepted checkpoint와 `.05` endpoint audit을 기다린다.

### 다음 midpoint 평가의 admissibility — 2026-10-19 18:30 KST

새 `half_audit_0016/admissibility.json`과 `half_audit_0026/admissibility.json`은 각각 h=.0025의 `t=.03875` 및 h=.00125의 `t=.031875` 중간 시험 상태를 `passed_native_tentative_admissibility`로 표시한다. 양쪽 모두 `mass_bitwise_initial=true`, 열에너지 양수, P/E/L 및 결합에너지 잔차가 기록된 roundoff bound보다 작다. 하지만 `accepted=false`는 정상적인 중간 후보 상태이며 현재 committed checkpoint는 여전히 각각 `t=.0375`, `t=.03125`다. 18:30 KST `sacct` 재조회에서 두 1-GPU 잡은 `RUNNING` (elapsed `1:24:21`), 2-GPU 잡은 `PENDING (QOSGrpGRES)`였다. 다음 판단은 full-step provider 완료 뒤에만 한다.

### h=.00125 분기의 여섯 번째 수락 step — 2026-10-19 18:36 KST

`run_h00125/runtime/accepted_0026.json`은 `step=26,t=.0325,h=.00125`, 12 physical providers, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`를 기록한다. step wall은 `873.34 s`; 장기 안정성 flag는 계속 false다. 대응하는 `candidate_audit_0026/admissibility.json`은 tentative report라 `accepted=false`로 표시되지만, native accepted-step report와 checkpoint가 authoritative commit을 나타낸다. h=.0025 분기에는 아직 `accepted_0016.json`이 없고, native log의 최신 tentative endpoint `.04`도 수락 표시가 아니므로 두 분기의 물리 비교는 수행하지 않았다. 직전 `sacct` 조회에서는 두 rank-1 H200 잡이 모두 `RUNNING` (elapsed `1:29:18`)이고 2-GPU 잡은 `PENDING (QOSGrpGRES)`였다.

### 양쪽 분기의 다음 step 수락 — 2026-10-19 18:37 KST

`accepted_0016.json`도 생성되어 h=.0025 branch가 `step=16,t=.04`를 12 provider와 `881.05 s` step wall로 수락했다. 이는 직전 h=.00125 branch `step=26,t=.0325` 수락과 함께 각각 고정 질량 및 `Re=∞` 조건을 확인한다. 두 분기는 아직 같은 시각이 아니므로 시간수렴 비교를 하지 않았고, `.05` endpoint의 q48 density/reference/physical gate도 아직 남아 있다. 18:37 KST `sacct`에서는 두 rank-1 잡 모두 `RUNNING` (elapsed `1:31:21`), 2-GPU scaling 잡은 `PENDING (QOSGrpGRES)`였다.

현재 제출 스크립트의 시간 제한은 h=.0025 branch `6 h`, h=.00125 branch `10 h`다. 18:37 KST의 elapsed `1:31:21`과 최신 수락 step wall `881.05/873.34 s`를 기준으로 하면 각각 `.05`까지 4 step(약 59분), 14 step(약 3시간 24분)이 남았다. 이는 기존 속도를 선형 적용한 일정 추정이며 endpoint q48/reference/physical audit 비용과 이후 step 변동을 포함하지 않는다. 지금 측정된 여유 시간 안에서 종점 도달이 예상되지만, 실제 accepted checkpoint와 Slurm 종료 상태로 계속 확인한다.

### 수락 상태의 회전·가중치 독립 계측 — 2026-10-19 18:42 KST

두 native checkpoint directory를 원본 IC와 직접 대조한 read-only 감사는 `Docs/canonical_live_native_spin_weight_audit_2026-10-19.json`에 저장했다. t=.04 (h=.0025)와 t=.0325 (h=.00125) 모두 19,945 IDs 및 mass 배열이 initial과 정확히 같고, initial의 모든 0 weight가 현재 모두 nonzero가 됐다. 따라서 이 짧은 경로에서 Laguerre weight가 변하면서 고정 질량은 유지된다.

두 상태에서 intrinsic-spin 변화와 centroid orbital angular-momentum 변화는 각각 `+3.84118e-7/-3.84118e-7` 및 `+2.49483e-7/-2.49483e-7`로 상쇄한다. 원점 총 각운동량 변화는 `3.44e-17` 및 `-3.52e-18`, 스핀+궤도 closure 잔차는 `3.19e-21` 및 `4.83e-21`이다. rigid-spin projection은 affine kinetic energy의 약 `10.64–10.65%`; 나머지 약 89.35%에는 케플러 전단을 포함하는 비강체 affine 변형이 남는다. 이는 보존 장부와 회전/가중치 자유도의 실제 진화를 확인하지만, 전단 오차가 줄었다거나 시간간격 수렴 또는 장기 안정성을 달성했다는 증거는 아니다. 두 checkpoint는 서로 다른 시각이라 교차 분기 정확도 비교에도 사용하지 않았다.

이 계측 시점의 `sacct` 조회에서 428541/428542는 계속 `RUNNING` (elapsed `1:36:19`), 428614는 `PENDING (QOSGrpGRES)`였다. 가장 최근 accepted 상태는 h=.0025의 t=.04 및 h=.00125의 t=.0325다.

### 다음 step의 midpoint 상태 — 2026-10-19 18:43 KST

`sacct` 재조회에서 428541/428542는 `syn104`에서 각각 H200 1장으로 `RUNNING` (elapsed `1:37:41`), 428614는 여전히 `PENDING (QOSGrpGRES)`다. h=.00125 branch의 `half_audit_0027/admissibility.json`은 tentative t=`.033125`에서 고정 질량, 양의 열, P/E/L 및 결합에너지 roundoff bound 통과를 기록하지만 `accepted=false`다. h=.0025 native log에도 다음 step의 tentative t=`.04125`가 있으나 아직 accepted 파일은 없다. 두 branch의 최신 committed 시각은 각각 `.04`와 `.0325`이며 `.05` endpoint audit 전이다.

### 두 번째 중간 후보의 admissibility — 2026-10-19 18:44 KST

`half_audit_0017` (h=.0025, tentative t=`.04125`)과 `half_audit_0027` (h=.00125, tentative t=`.033125`) 모두 `passed_native_tentative_admissibility`, mass bitwise initial, 양의 thermal, P/E/L 및 combined-energy 잔차가 각 bound 이내다. 둘 다 아직 `accepted=false`이므로 committed 시각은 `.04`와 `.0325`다. `sacct`에서 두 rank-1 H200 잡은 `RUNNING` (elapsed `1:38:36`), 2-GPU 잡은 `PENDING (QOSGrpGRES)`였다. 다음 확인은 full-step map과 native commit 완료 뒤에 한다.

### h=.00125 분기의 다음 accepted step — 2026-10-19 18:50 KST

`run_h00125/runtime/accepted_0027.json`에서 `step=27,t=.03375,h=.00125`, `physical_providers=14`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, 1 rank가 확인된다. step wall은 `875.08 s`; `long_time_stability_verified=false`다. 같은 step 번호의 `candidate_audit_0027/admissibility.json`은 tentative audit로 `accepted=false`를 기록하지만, native accepted report와 `accepted_0027` checkpoint가 실제 commit 근거다. audit의 mass-bitwise flag는 true이고 최소 thermal은 양수이며 P/E/L 잔차는 약 `1e-17`로 기록 bound보다 작다.

같은 확인 시각에 h=.0025의 마지막 accepted report는 여전히 `step=16,t=.04`; native log에는 다음 후보 `.0425`까지만 있으며 새 `accepted_0017.json`은 아직 없다. 따라서 분기는 비교 가능한 동일 시각에 도달하지 않았다. `sacct` 18:51 KST 조회에서는 428541/428542가 syn104에서 계속 `RUNNING` (elapsed `1:45:23`), 428614는 `PENDING (QOSGrpGRES)`였다. `.05` endpoint density/reference/physical audit 및 두 timestep 결과 비교는 미완료다.

### h=.0025 분기의 다음 accepted step — 2026-10-19 18:52 KST

`run_h0025/runtime/accepted_0017.json`에서 `step=17,t=.0425,h=.0025`, `physical_providers=14`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`가 확인된다. step wall은 `880.33 s`; 장기 안정성 flag는 여전히 false다. `candidate_audit_0017/admissibility.json`은 `accepted=false`의 tentative audit지만 mass-bitwise 유지, 양의 최소 thermal, ledger 잔차 약 `1e-17`을 기록하며, 실제 수락은 native accepted report/checkpoint가 증명한다.

따라서 최신 committed 시각은 h=.0025에서 `.0425`, h=.00125에서 `.03375`이며 아직 같은 시각의 시간간격 비교가 불가능하다. h=.00125 provider 로그에는 다음 full-moment provider 시작이 기록됐지만 다음 step은 아직 수락되지 않았다. `sacct` 18:52 KST에서 428541/428542는 계속 `RUNNING` (elapsed `1:46:26`), 428614는 `PENDING (QOSGrpGRES)`다. `.05` endpoint density/reference/physical gate와 timestep 간 결과 비교는 여전히 남아 있다.

### 다음 timestep 후보의 진행 — 2026-10-19 18:58 KST

`run_h00125/runtime/half_audit_0028/admissibility.json`의 tentative `t=.034375`는 `passed_native_tentative_admissibility`, `mass_bitwise_initial=true`, 양의 최소 thermal, 그리고 P/E/L 및 combined-energy 잔차가 각 roundoff bound 안임을 기록한다. 이 상태는 아직 `accepted=false`이며 committed 시각은 계속 `.03375`다.

h=.0025 native log에는 provider 15 완료와 다음 tentative 후보 `t=.04375`가 있다. 그러나 해당 half audit 및 accepted checkpoint는 아직 보이지 않아 committed 시각은 `.0425`다. 18:58 KST `sacct`에서 428541/428542는 계속 `RUNNING` (elapsed `1:51:38`), 428614는 `PENDING (QOSGrpGRES)`였다. 따라서 양쪽 분기의 endpoint 비교와 `.05` density/reference/physical gate는 미완료다.

### h=.0025 다음 midpoint 후보 admissibility — 2026-10-19 18:59 KST

`run_h0025/runtime/half_audit_0018/admissibility.json`의 tentative `t=.04375`도 `passed_native_tentative_admissibility`다. mass bitwise initial, 양의 최소 thermal, P/E/L 및 combined-energy 잔차가 기록 bound 안이다. 그러나 `accepted=false`이므로 h=.0025 committed 시각은 `.0425` 그대로다. h=.00125의 `half_audit_0028`도 `t=.034375`에서 같은 tentative admissibility만 통과했고 commit은 아직 `.03375`다. `sacct` 18:59 KST에서 두 rank-1 H200 잡은 `RUNNING` (elapsed `1:53:15`), 2-GPU 잡은 `PENDING (QOSGrpGRES)`였다.

### h=.00125 accepted step 28 — 2026-10-19 19:05 KST

`run_h00125/runtime/accepted_0028.json`은 `step=28,t=.035,h=.00125`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, 16 physical providers, 19,945 cells, 1 rank를 확인한다. step wall은 `870.25 s`; `long_time_stability_verified=false`다. `candidate_audit_0028/admissibility.json`은 별도의 tentative audit (`accepted=false`)이며 mass-bitwise initial, 양의 최소 thermal, P/E/L 잔차가 각 bound 이내다. native accepted report/checkpoint가 commit의 근거다.

h=.0025 branch의 native log에는 provider 16 및 다음 tentative `t=.045`가 있지만 새 accepted checkpoint는 없다. 따라서 최신 committed 시각은 h=.0025 `t=.0425`, h=.00125 `t=.035`다. 두 시각이 다르므로 시간간격 결과 비교는 아직 불가하며 `.05` endpoint physics gate도 미실행이다. `sacct` 19:05 KST에서 rank-1 잡 428541/428542는 `RUNNING` (elapsed `1:58:57`), 2-GPU 잡 428614는 `PENDING (QOSGrpGRES)`였다.

### h=.0025 accepted step 18 — 2026-10-19 19:07 KST

`run_h0025/runtime/accepted_0018.json`은 `step=18,t=.045,h=.0025`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, 16 physical providers, 19,945 cells, 1 rank를 기록한다. step wall은 `880.38 s`; `long_time_stability_verified=false`다. candidate admissibility audit의 `accepted=false`는 tentative audit 형식이며, 실제 commit은 accepted JSON과 checkpoint가 증명한다. 질량은 bitwise initial, thermal 양수, ledger 잔차는 기록 bound 안이다.

동일 시각이 아니므로 현재 committed 상태 h=.0025 `t=.045`, h=.00125 `t=.035`를 직접 비교하지 않았다. `.05`까지 각각 2 step과 12 step이 남았다. 최근 step wall을 그대로 적용한 선형 추정은 약 29분과 2시간 54분이며, endpoint q48 density/reference/physical audit 비용 및 step 변동을 포함하지 않는다. `sacct` 19:07 KST에서 428541/428542는 `RUNNING` (elapsed `2:00:38`), 428614는 `PENDING (QOSGrpGRES)`였다. `.05` 비교 및 장기 안정성은 여전히 미검증이다.

### 다음 step midpoint 후보 — 2026-10-19 19:12 KST

h=.00125 `half_audit_0029/admissibility.json`의 tentative `t=.035625`는 `passed_native_tentative_admissibility`다. 질량은 initial과 bitwise 동일하고, thermal 양수, P/E/L 및 결합에너지 잔차는 bound 안이지만 `accepted=false`; committed 시각은 `.035`다. h=.0025 native log에는 다음 step의 provider 17 및 tentative `t=.04625`가 기록됐으나 accepted commit이나 해당 midpoint audit은 아직 없다. committed 시각은 `.045`다.

`sacct` 19:12 KST에서 428541/428542는 계속 `RUNNING` (elapsed `2:06:15`), 428614는 `PENDING (QOSGrpGRES)`였다. `.05` endpoint 및 timestep 비교는 미완료다.

### 장기 안정성 목표의 명시적 범위 — 2026-10-19 19:14 KST

사용자가 목표 시간을 `t>6`으로 명확히 갱신했다. 따라서 `.05` endpoint는 수렴성·초기 물리 게이트일 뿐이며 목표 종료 시각이 아니다. 두 시간간격 분기의 `.05` 게이트를 마친 뒤 체크포인트 연속 진화를 통해 반드시 `t>6`까지 가고, 그 전 구간에서 궤도 위상, 밀도/압력/속도 오차, 보존량 및 셀 회전·가중치 장부를 검증해야 한다. `long_time_stability_verified`는 `t>6`까지의 실제 증거가 생기기 전에는 false로 취급한다.

19:14 KST 현재 committed 시각은 h=.0025 `t=.045`, h=.00125 `t=.035`다. h=.0025의 `half_audit_0019` tentative `t=.04625`는 `passed_native_tentative_admissibility`, mass bitwise initial, 양의 minimum thermal, P/E/L 및 결합에너지 잔차 bound 내지만 `accepted=false`다. h=.00125의 최신 committed `t=.035` 뒤 tentative midpoint `t=.035625`도 이전 audit에서 같은 admissibility 통과만 했고 commit은 아니다. `sacct`에서 428541/428542는 `RUNNING` (elapsed `2:07:34`), 428614는 `PENDING (QOSGrpGRES)`였다. 어느 분기도 아직 `.05`에 도달하지 않았다.

### `t>6`까지의 현재 성능 외삽 — 2026-10-19 19:17 KST

실제 accepted report에서 h=.0025의 step wall은 `880.38 s` (현재 committed `t=.045`), h=.00125는 `870.25 s` (`t=.035`)다. 동일 비용이 계속된다고 단순 가정하면 strict `t>6`까지 각각 2,383 step / 4,773 step, 약 24.28일 / 48.08일의 계산시간이 필요하다. 이는 장기 비용을 보여주는 선형 외삽일 뿐, 장기 물리 안정성 또는 step 비용이 일정하다는 증거가 아니다. endpoint 진단 및 checkpoint/job 재시작 overhead도 제외된다.

19:17 KST `sacct`는 428541 (h=.0025) `RUNNING`, elapsed `2:10:51`, Timelimit `6:00:00`; 428542 (h=.00125) `RUNNING`, 같은 elapsed, Timelimit `10:00:00`; 428614 (2 ranks/2 H200) `PENDING`, Timelimit `6:00:00`를 보고한다. 이 allocation들은 예상대로 `.05` gate에는 충분한 headroom이 있지만 각각의 현재 job allocation 하나로 `t>6`에 도달할 수 없다. `.05` 비교가 통과하면 검증된 branch를 checkpoint에서 여러 allocation으로 연속 수행하거나, 정확성을 보존하는 성능 개선을 평가해야 한다. 긴 목표를 만족하는 물리적 안정성 검증은 계속 필수다.

### CUDA/MPI 가속 경로 감사 — 2026-10-19 19:22 KST

갱신된 목표의 성능 조건을 현재 source와 live log로 점검했다. 실제 실행은 1 MPI rank / 1 H200이다. `Exam/full_moment_evolution.c`에서 rank 수를 읽고 rank 1–2만 허용한다. `Exam/full_moment_face_provider.py`는 rank-0 Python service에서 `full_disk_stage_memory.evaluate`를 호출하고 report에 `service_MPI_calls=0`, `service_GPU_calls=0`를 명시한다. 이 service가 `build_geometry`, `prepare_states`, `records_build`, `solve_from_transport` 전역 face/root 경로를 실행한다. 반면 `Exam/full_moment_native_stage.c`만 device를 선택하고 셀 partition 뒤 `lagmoment_owned_stage`를 호출하며 MPI `Gatherv`로 결과를 모은다.

로그에서 provider RHS는 대략 `320–324 s`; 같은 호출의 native owned kernel은 약 `2.3 s`, native MPI/GPU stage 전체는 약 `5.2 s`다. accepted step 하나에 provider 두 호출(`physical_providers`가 step마다 +2)과 native stage 두 번이 들고, 최근 step wall은 `874.28 s`다. 따라서 현재 MPI/CUDA 경로에서 계측되는 전체 owned stage는 accepted step wall의 약 1.2%뿐이다. 이는 H200이 느리다는 뜻이 아니라, 전역 기하·root service가 대부분 CPU/Python 측에 남아 GPU/MPI 가속 범위가 전체 계산의 작은 부분이라는 source+runtime 증거다.

현재 accepted 상태는 h=.0025 `step=19,t=.0475` (`accepted_0019.json`), h=.00125 `step=29,t=.03625` (`accepted_0029.json`)이며 둘 다 `fixed_mass=true`, `Re_infinity=true`, `long_time_stability_verified=false`다. 428541/428542는 `RUNNING`; 두 rank/두 H200의 428614는 여전히 `PENDING`이며 아직 rank-scaling 측정값이 없다. 효율 목표를 충족했다고 말할 수 없다. 유효한 다음 성능 판단은 428614가 실제 실행된 뒤 step wall 및 provider/native stage 분해를 rank 1 대비 비교하는 것이고, 장기 계산 전에는 전역 geometry/root의 GPU/MPI 분할을 주요 최적화 대상으로 평가해야 한다. 현재 `.05` 수렴 관문은 두 branch 모두 아직 통과하지 않았다.


### 최신 진행 — 2026-10-19 19:39 KST

h=.0025 branch `428541`은 Slurm에서 `COMPLETED`로 바뀌었고 `accepted_0020.json`은 `step=20,t=.05`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, `long_time_stability_verified=false`를 기록한다. 같은 시각의 `candidate_audit_0020/physical_gate.json`은 all/inner/cavity/disk_core/inner_contact/outer_contact/dynamic_buffer/wall 8개 영역을 모두 통과했다. q48 밀도 준비 보고는 `passed_full_stage_Psi_density_preparation`, 독립 정상장 참조 보고는 `passed_steady_reference_on_moving_cells`다. 이로써 h=.0025 branch의 `.05` checkpoint/초기 물리 gate는 확인됐지만, 이것은 시간간격 수렴 비교나 장기 안정성 증거가 아니다.

h=.00125 branch `428542`는 계속 `RUNNING`이며 최신 accepted report는 `step=30,t=.0375`, step wall `886.70 s`; 같은 checkpoint는 고정질량 및 `Re=∞`를 기록한다. 따라서 동일한 `.05` 시각에서의 temporal/density 비교는 아직 수행할 수 없다. 2-rank/2-H200 scaling job `428614`는 계속 `PENDING (QOSGrpGRES)`이며 측정값은 없다. Slurm accounting은 `428541` 종료 상태 `COMPLETED`, `ExitCode=0:0`을 확인했다.

### 2-rank H200 비교 잡 실행 확인 — 2026-10-19 19:48 KST

`squeue/sacct` 조회에서 `428614`는 `syn104`에서 `RUNNING`이며 2 H200, 16 CPU, 256 GB를 할당받았다. `428542`도 계속 `RUNNING`이다. `428614/runtime/restart_import.json`은 `verified_restart_provenance_no_state_update`, `time=.025`, `requested_h=.0025`, `numerical_updates=0`, `all_source_files_verified=true`를 기록한다. 양쪽 MPI rank의 native log도 `t=.025` 재시작 import 완료를 출력했다. 아직 새 물리 timestep은 수락되지 않았다.

확인 시 native log의 최신 항목은 19:46:54의 import 완료이며, provider rank 로그는 비어 있고 GPU UUID 및 native stage timing은 아직 없다. 소스상 CUDA 장치 선택 로그는 첫 native face stage 안에서 출력되므로 이 시점에 UUID가 없는 것만으로 CUDA 미사용이나 멈춤을 판정할 수 없다. 첫 provider/native stage가 끝난 뒤 실제 rank별 GPU UUID와 timing을 확인한다. 2-rank 성능 향상은 아직 측정되지 않았다.

### h=.00125 step 31 및 비교 잡 초기화 — 2026-10-19 19:49 KST

`run_h00125/runtime/accepted_0031.json`은 `step=31,t=.03875,h=.00125`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, `long_time_stability_verified=false`, 19,945 cells를 기록한다. 이번 step wall은 `878.84 s`다. fine branch는 `.05`보다 9 accepted step 뒤에 있다. 동일시각 timestep 비교와 장기 안정성 판정은 아직 할 수 없다.

마지막 Slurm 조회에서 428542는 `RUNNING` (elapsed `2:42:30`), 428614는 `RUNNING` (elapsed `5:00`, 2 H200)였다. 2-rank `provider_rank0.log`에는 service 준비(`37.18 s`)와 약 4.34M node resource preflight가 나타났지만 첫 native face-stage timing, 실제 UUID별 CUDA 사용, accepted step은 아직 기록되지 않았다. 따라서 MPI/GPU scaling은 미측정 상태다.

### 2-rank provider 반복 계측 — 2026-10-19 19:51 KST

다음 Slurm 조회에서 428614는 계속 `RUNNING` (elapsed `6:56`, 2 H200), 428542도 `RUNNING` (elapsed `2:44:26`)이었다. 428614의 rank-0 provider는 약 4.34M face nodes에 대해 geometry/face-state/mass-root 반복 0–2까지 진행했다. 반복별 합산시간은 각각 `74.67 s`, `40.55 s`, `40.45 s`이며, 반복 2의 mass residual은 `5.24e-16`이다. 따라서 시작 지연으로 보였던 구간은 실제 CPU provider 계산 중이었다.

2-rank 잡은 아직 첫 accepted step에 도달하지 않았고 native CUDA stage의 rank별 UUID 및 MPI/GPU 시간도 아직 로그에 없다. 현재까지는 GPU 2장 할당만 검증됐으며, 처리량 향상은 미측정이다. h=.00125 branch의 최신 accepted report는 step 31, `t=.03875`이며 `.05` 동일시각 비교까지 9 step 남았다.

### 2-rank 첫 native stage — 2026-10-19 19:52 KST

428614의 두 MPI rank가 서로 다른 H200 UUID를 실제 선택했다: rank 0 `9281c2516d6db647453f180365da4dbd`, rank 1 `cc32769c2315f4a28e744ad3a43159e6`. 첫 native stage에서 rank별 owned kernel은 각각 `1.300 s`와 `1.598 s`, collective MPI/GPU stage wall은 `4.432 s`였다. 같은 t=.025 source checkpoint에서 시작한 1-rank control의 첫 stage는 kernel `2.178 s`, MPI/GPU `5.073 s`였다. 따라서 해당 native stage wall은 약 12.6% 줄었지만 아직 전체 step의 가속률은 아니다.

같은 첫 RHS provider 구간은 2-rank `299.643 s`, 1-rank control `322.965 s`로 7.2% 짧았다. 1-rank control은 `accepted_0011.json`에서 step 11 (`t=.0275`)을 `828.91 s`에 수락했다. 2-rank branch는 현재 첫 tentative midpoint(`t=.02625`)까지만 진행했고, 아직 step 11 commit/full-step wall은 없다. 실제 checkpoint·물리 결과·전체 step wall을 비교한 뒤에만 MPI scaling을 평가한다. Python geometry/face/root provider가 약 300 s/RHS로 남아 있는 것이 주요 병목이라는 점은 이 초기 측정에서도 유지된다.

### 2-rank 첫 step midpoint geometry 재구성 — 2026-10-19 19:56 KST

`428614`는 Slurm에서 계속 `RUNNING` (elapsed `11:40`, 2 H200)이며 `428542`도 `RUNNING` (elapsed `2:49:10`)이다. 2-rank 잡의 첫 midpoint audit 뒤 provider service가 새 midpoint state용 configure/resource preflight를 시작했고, 19:55:07에 약 4.33M node preflight를 남겼다. native log에는 아직 step 11 full RHS/commit이 없으며, `accepted_*.json`도 생성되지 않았다. `428542` fine branch의 최신 accepted report는 step 31, `t=.03875`다.

### rank 1 대 rank 2 동일 상태 step 11 비교 — 2026-10-19 20:01 KST

2-rank `428614`가 step 11 (`t=.0275,h=.0025`)을 실제 수락했다. authoritative `accepted_0011.json`은 `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, `long_time_stability_verified=false`, step wall `708.88199 s`다. 같은 `t=.025` source checkpoint와 동일 binary에서 시작한 1-rank control step wall은 `828.90670 s`; 2-rank는 첫 동일 조건 step에서 `14.48%` 빨랐다. 비교는 1-rank 1 H200/8 CPU와 2-rank 2 H200/16 CPU allocation의 결합 효과이며 GPU-only speedup으로 분리하지 않는다.

MPI 결과의 수치 일치도 확인했다. 두 accepted state의 `ids`, `mass`, `state`, `factors`, `weight_rate_seed` 파일은 bytewise 동일하다. `Psi`의 16-byte long-double 수치도 전 셀에서 정확히 같고, 저장 padding만 다르다. header는 run token만 다르며, 두 `candidate_audit_0011/admissibility.json`도 정확히 같다 (`passed_native_tentative_admissibility`; 이 tentative 보고의 `accepted=false`는 최종 commit을 뜻하지 않으며, 실제 수락은 별도의 accepted report가 증명한다). 따라서 rank 수를 바꿔도 이 첫 accepted state의 물리 수치가 변하지 않았음을 확인했다.

2-rank 잡은 다음 step들을 계속 수행 중이다. 이 한 step의 14.48% 향상은 실측이지만 장기 throughput은 아직 미확정이다. h=.00125 branch는 마지막으로 확인된 accepted `step=31,t=.03875`에서 `.05` 시간간격 비교를 향해 가고 있으며, 두 branch 모두 `t>6` 검증은 남아 있다.

### t=.03875 세이브드 스핀·가중치·각운동량 계측 — 2026-10-19 20:04 KST

h=.00125 `accepted_0031`의 저장 state만 읽어 t=0 Voronoi IC와 비교했다. cell spin은 독립 DOF가 아니라 raw 18 moment 안의 `ell=C_xy-C_yx`; rigid-rotation energy `Krot=ell²/(2 Tr(S))`는 affine kinetic energy 안의 projection이다. t=.03875에서 이 projection은 affine energy의 `10.6366%`이며 나머지 affine deformation/shear가 `89.3634%`다. 이는 회전/변형 분해의 bookkeeping이지, 회전이 Kepler shear 오차를 줄였다는 인과 증거는 아니다.

`fixed_mass`는 initial과 bitwise 동일하다. 초기 가중치 전부 0에서 현재 19,945개 cell 모두 gauge-removed 가중치가 nonzero로 변했다 (`RMS=9.63e-6`, min `-1.09e-5`, max `7.85e-5`), 즉 zero-weight Voronoi 극한에서 벗어난 weighted geometry다. 내부 spin 합은 `+2.7250e-7` 변했고 centroid-orbital 각운동량은 `-2.7250e-7` 변했으며, origin 기준 총 각운동량은 초기 `13.883106850908922` 그대로였다. spin–orbit 교환 장부는 보존된다. 그래도 이 검증은 `t=.03875`까지만이며, shear-error suppression 및 `t>6`는 미검증이다.

### t=.04 fine checkpoint 및 spin/weight 수렴 중간점검 — 2026-10-19 20:05 KST

h=.00125 `accepted_0032`는 `t=.04`를 `827.65 s`에 수락했다 (`fixed_mass=true`, `Re_infinity=true`, 장기 안정성 flag false). h=.0025 `accepted_0016`도 동일 `t=.04`에 있다. 두 branch는 같은 t=.025 source checkpoint와 동일 binary에서 갈라졌으며 ids와 masses가 bytewise 같다.

저장 raw moments만 이용한 동일시각 spin 비교에서 h=.0025 대비 h=.00125의 mass-weighted RMS 차이는 intrinsic spin `0.089%`, 회전율 `0.080%`, `Krot` `0.312%`, affine energy `0.0347%`였다. gauge-removed weight field 차이는 fine branch weight RMS의 `4.22%`다. origin 각운동량 총합 차이는 `-3.82e-17`이었다. 이는 spin/affine/weight 지표의 `.04` 중간 temporal comparison일 뿐, full density·velocity physical gate나 `.05` 공식 비교를 대체하지 않는다.

마지막 scheduler/log 확인에서 `428614`는 계속 실행 중이었다. rank-2 branch는 step 11 `t=.0275`를 수락했고 step 12의 half-step RHS(`t=.02875`)까지 진행했다. 같은 시각의 전체 step scaling과 `.05` endpoint는 아직 남아 있다.

### step 12 full-RHS 진행 — 2026-10-19 20:09 KST

Slurm 조회에서 `428614`는 `RUNNING` (elapsed `25:47`, 2 H200), `428542`는 `RUNNING` (elapsed `3:03:17`)이었다. 2-rank step 12는 half-step tentative 이후 full-RHS provider에 들어갔고, rank-0 log에는 약 4.33M face node에 대한 geometry 반복 0–2가 진행됐다. native log에는 아직 step 12 accepted report가 없으므로 step 11 실측속도를 sustained throughput으로 일반화하지 않는다. fine branch는 계속 accepted `step=32,t=.04`이며 `.05`까지 8 step 남았다.

### rank 2 step 12 대조군 재현성과 throughput — 2026-10-19 20:12 KST

2-rank job은 step 12 (`t=.03,h=.0025`)을 `742.50500 s`에 수락했다. 같은 시작 state의 1-rank control은 `870.59941 s`; 두 번째 연속 step에서는 `14.71%` 빨랐다. step 11–12 평균은 rank 1 `849.753 s`, rank 2 `725.693 s`, 즉 `14.60%` 단축이다. 이는 2 H200/16 CPU와 1 H200/8 CPU를 비교한 결합 scaling이다.

step 12도 numerical state 재현성이 확인됐다: `ids`, `mass`, `state`, `factors`, `weight_rate_seed`가 bytewise 같고 `Psi` long-double 숫자는 모두 같았다. header의 run token만 다르다. 두 tentative admissibility JSON도 같으며 `accepted=false`는 audit 형식의 후보 상태를 뜻한다. 별도 `accepted_0012.json`이 실제 commit을 증명한다. 장기 안정성과 `.05` temporal gate는 여전히 미검증이다.

### 2-step MPI throughput 분해 — 2026-10-19 20:13 KST

step 11–12 전체 평균 1-rank `849.75 s`, 2-rank `725.69 s` (14.60% 개선) 중 provider RHS 합산시간은 각각 step당 평균 `617.13 s`, `571.81 s` (7.34% 개선)였다. 2-rank step wall의 약 `78.8%`가 여전히 Python provider 계산이다. 개별 owned MPI/GPU stage는 1-rank `5.07–5.16 s`, 2-rank `4.43–4.47 s`로 빨라졌으나 step당 수 초 수준이라 전체 시간을 지배하지 않는다. 다음 성능 개선은 H200 추가만으로 해결되지 않으며 전역 geometry/face/root provider 구간이 핵심이다.

현재 두 step 평균을 일정하다고 놓는 단순 외삽으로는 `.05`부터 strict `t>6`까지 h=.0025 기준 약 2,381 step, 약 20.0일이다 (1-rank step wall을 쓰면 약 23.4일). 장기 step 비용과 물리 안정성은 전혀 검증되지 않았으므로 planning estimate로만 취급한다. `.05` temporal/physics gate 뒤에는 체크포인트 연속 restart 계획이 필요하고, 먼저 step wall이 여러 시점에서도 유지되는지 관측한다.

`428614`는 step 12 이후에도 `RUNNING`이며 step 13 provider의 resource preflight가 20:13 KST에 남았다 (약 4.34M face nodes). 아직 step 13 accepted report는 없고, job은 중단/재시작하지 않았다. h=.00125 branch는 마지막 accepted `step=32,t=.04`다.

### 두 분기 다음 후보 진행 확인 — 2026-10-19 20:16 KST

한 차례 `squeue/sacct` 확인에서 두 계산 모두 계속 `RUNNING`이었다: 428614 rank-2/2 H200 elapsed `32:11`, 428542 rank-1/1 H200 elapsed `3:09:41`. 잡을 취소하거나 재시작하지 않았다.

rank-2 branch의 최신 accepted 상태는 여전히 `step=12,t=.03`; native log는 다음 tentative/midpoint `t=.03125`를 보이고, provider log는 새 약 4.34M-node 계산의 geometry/mass-root 반복 0–3까지 기록했다. 아직 `accepted_0013.json`은 생성되지 않았다.

h=.00125 branch의 최신 accepted 상태는 `step=32,t=.04`. `half_audit_0033`의 `t=.040625` 후보는 질량 bitwise 유지, 양의 thermal, P/E/L 잔차 roundoff bound 통과를 기록하지만 `accepted=false`다. native log는 이어서 candidate `t=.04125`까지 진행했으며, `accepted_0033.json`은 아직 없다.

따라서 rank-2와 fine branch의 최신 committed 시각은 각각 `.03`과 `.04`이고, 둘 다 아직 `.05` endpoint에 도달하지 않았다. 한-rank h=.0025 대조군은 앞서 `.05`에서 종료됐다. `.05` 시간간격 비교 및 물리 gate 전까지 장기 안정성 판정은 계속 보류한다.

### 현재 Slurm 잡의 종료 범위 및 장기 연속 실행 준비 — 2026-10-19 20:17 KST

활성 잡의 `params.dat`와 제출 기록을 직접 확인했다. 428614 rank-2 h=.0025는 source step 10에서 `moment_steps=20`, 428542 fine h=.00125는 source step 20에서 `moment_steps=40`으로 설정되어 둘 다 설계상 `t=.05`에서 끝난다. `submission.json`의 목표 시각도 `.05`다. 따라서 이 잡들은 `t>6`까지 자동으로 이어지는 장기 잡이 아니다.

`.05` gate와 공통시각 두 timestep 비교가 통과하면, 실제 accepted `.05` checkpoint에서 state와 provenance manifest를 새 run directory에 fork하고, 각 Slurm allocation의 종료 전에 재시작 가능한 연속 구간을 설계해야 한다. 현재 입력 생성기는 `.025` checkpoint SHA/시각에 고정되어 있으므로 그대로 재사용하면 안 된다. 장기 run용 checkpoint-fork builder/runner를 준비해야 한다. 이번에는 active jobs나 solver source를 수정하지 않았다.

### h=.00125 accepted step 33 — 2026-10-19 20:18 KST

이전 점검 직후 fine branch가 `accepted_0033.json`을 생성했다: `step=33,t=.04125,h=.00125`, `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, `long_time_stability_verified=false`, step wall `883.6906 s`. 대응 admissibility 보고는 후보 형식상 `accepted=false`지만 mass bitwise flag는 true, minimum thermal은 양수, P/E/L 잔차와 combined-energy 변화는 각각 기록된 roundoff bound 안이다. commit 근거는 별도 accepted report다.

fine branch는 `.05`까지 7 accepted steps 남았다. rank-2 h=.0025 branch는 최신 확인에서 accepted `step=12,t=.03`이며 step 13 midpoint `t=.03125`와 provider RHS epoch 5까지만 보인다. 따라서 두 branch의 공통 `.05` gate는 여전히 미완료다.

### 활성 provider 단계 재확인 — 2026-10-19 20:19 KST

`squeue/sacct`에서 428614와 428542는 계속 `RUNNING`이었다 (elapsed 각각 `35:43`, `3:13:13`). rank-2 native log에는 아직 `accepted_0013.json`이 없지만, provider log는 새 RHS의 약 4.33M face-node resource preflight와 geometry 반복 0을 20:19:34까지 기록했다. Fine branch는 `accepted_0033` (`t=.04125`) 이후 새 provider call을 시작했다. 두 작업 모두 실행 중이며, 아직 다음 accepted step은 나오지 않았다.

### rank-2 accepted step 13 및 3-step 재현성 — 2026-10-19 20:23 KST

428614는 step 13 (`t=.0325,h=.0025`)을 `734.77168 s`에 수락했다. 같은 입력에서 시작한 1-rank 대조군 step 13은 `870.51120 s`다. step 11–13 평균은 rank 1 `856.67244 s`, rank 2 `728.71956 s`; 2 H200/16 CPU가 1 H200/8 CPU보다 세 연속 step에서 평균 약 `14.93%` 빨랐다. GPU-only speedup으로 분리하지 않는다.

동일 t=.0325 checkpoint의 serialized fields를 비교했다. `ids`, `mass`, `state`, `factors`, `Psi` long-double 수치, `weight_rate_seed`가 모두 numeric-equal이며, header는 run token만 다르다. step 13 후보 admissibility는 고정 질량, 양의 thermal, P/E/L 및 combined-energy roundoff bound를 통과했다. 따라서 3개 연속 step에서 rank 수 변경이 허용된 timestep state를 바꾸지 않았다는 증거가 추가됐다.

같은 시각 fine branch는 accepted `step=33,t=.04125` 뒤 step 34 후보 `t=.041875`를 계산 중이다. 마지막 scheduler 확인에서 두 잡은 계속 `RUNNING`; `.05` endpoint gate와 long-time 검증은 여전히 미완료다.

### step 14 / step 34 provider 진행 중 — 2026-10-19 20:26 KST

다음 `squeue/sacct` 확인에서도 428614와 428542는 모두 `RUNNING`이었다 (elapsed 각각 `42:51`, `3:20:21`). 428614는 `accepted_0013` 이후 step 14용 provider 계산을 시작해 geometry 반복 0–1까지 기록했으나, 아직 `accepted_0014.json`은 없다. 428542는 accepted step 33 (`t=.04125`) 상태에서 step 34 midpoint admissibility (`t=.041875`)를 통과했으며, 다음 provider call이 시작됐지만 `accepted_0034.json`은 아직 없다. 두 작업 모두 끝내거나 재시작할 이유가 없으며 계속 실행 중이다.

### h=.0025 t=.05 재시작 후보 provenance 봉인 — 2026-10-19 20:28 KST

기존 rank-1 h=.0025 accepted checkpoint `run_h0025/runtime/accepted_0020` (step 20, `t=.05`)는 전체 8개 물리 영역 gate, q48 밀도 검증, 고정 질량, 양의 thermal 및 P/E/L 장부 검증을 통과한 상태다. 이를 선택 가능한 continuation seed로 남기기 위해 새 디렉터리 `restart_seal_h0025_t0050/manifest.json`에 provenance를 봉인하고 `verify()`로 다시 확인했다. 모델/integrator, `h=.0025`, time/step, 초기 질량 및 모든 상태 파일 SHA, checkpoint audit/factor evidence가 확인됐다. seal/verify는 `physical_providers=0`, `accepted_time_steps=0`이며 state를 바꾸거나 simulation step을 수행하지 않는다.

이 checkpoint는 아직 장기 branch로 채택하지 않았다. fine h=.00125의 같은 `.05` gate와 두 timestep의 공통시각 density/velocity/phase 비교가 끝난 뒤에만 h=.0025 seed를 사용한다. 현재 마지막 scheduler 확인에서는 두 활성 잡이 여전히 `RUNNING`이었다 (428614 elapsed `44:48`, 428542 elapsed `3:22:18`); rank-2는 step 14 midpoint/provider 계산 중이고 fine branch는 step 34 candidate를 계산 중이다.

### native h=.0025 `.05`와 기존 동일-timestep endpoint 대조 — 2026-10-19 20:31 KST

기존 `continuous_two_dt_to_t1/fine/checkpoints/fine_0020`는 h=.0025, `t=.05` endpoint다. 새 native rank-1 h=.0025 `candidate_audit_0020/candidate.npz`와 read-only 비교했다. IDs와 mass는 exact equality다. flattened `state` relative RMS 차이는 `7.58e-14` (max absolute `4.58e-12`), `Psi` `4.70e-14`, factors `3.16e-9`, weight-rate seed `1.73e-7`다. 비트 단위 전체 state equality는 아니지만 매우 근접하다.

두 endpoint의 8-region physical-gate 보고서 모두 통과했고, 32개 pressure/velocity/phase 지표 중 최대 상대 차이는 `5.61e-7` (wall pressure maximum)이다. independent q48 최대 오차는 과거 `4.34e-17`, 새 native `4.70e-17`로 둘 다 작다. 따라서 같은 h=.0025 native endpoint에서 이미 수행한 spin/weight/rotation 해석을 다시 실행할 필요는 낮다. 새 h=.00125 native branch의 `.05` 결과와 matched-time 비교는 여전히 필요하다.

### fine accepted step 34 및 활성 잡 진행 — 2026-10-19 20:32 KST

428542 fine branch는 `accepted_0034.json`에서 step 34 (`t=.0425,h=.00125`)를 `859.96171 s`에 수락했다. fixed mass/Re=∞ 조건을 유지하며 장기 안정성 flag는 false다. 대응 candidate audit은 mass bitwise initial, 양의 최소 thermal, P/E/L과 combined-energy roundoff bound를 통과했다. `.05`까지 6 step 남았다.

428614 rank-2 branch는 latest accepted `step=13,t=.0325`; step 14 half audit `t=.03375` admissibility는 통과했지만 아직 `accepted_0014.json`은 없다. 마지막 scheduler 확인에서 두 작업 모두 `RUNNING` (elapsed `48:35` 및 `3:26:05`)이었다. 현재 `.05` endpoint와 native fine/coarse matched-time gate는 여전히 미완료다.

### step 14 / 35 계산 재확인 — 2026-10-19 20:33 KST

428614와 428542는 모두 `RUNNING` (elapsed `50:12`, `3:27:42`). Rank-2의 최신 commit은 여전히 step 13 (`t=.0325`); step 14 midpoint `t=.03375`는 통과했고 provider RHS epoch 7은 `296.25 s`에 끝났으며, 이어지는 provider 반복 2까지 mass-root residual은 `1.91e-16`이다. step 14 commit은 아직 없다.

Fine branch는 accepted step 34 (`t=.0425`) 뒤 다음 provider call을 시작했지만 step 35 accepted report는 아직 없다. 두 branch는 계속 진행 중이며 `.05`까지 각각 coarse 7, fine 6 accepted steps 남았다.

### step 14 full-step candidate 및 fine step 35 provider 시작 — 2026-10-19 20:34 KST

428614 rank-2 job은 `RUNNING` (elapsed `50:41`). Step 14의 provider RHS epoch 8이 약 `296.92 s`에 완료됐고 native log는 tentative candidate `t=.035`를 기록했지만 `accepted_0014.json`은 아직 없다. Provider nonlinear weight/mass-root outer iteration 3 residual은 `1.91e-16`이다.

428542 fine job도 `RUNNING` (elapsed `3:28:11`). Accepted step 34 (`t=.0425`) 뒤 다음 provider call이 새 약 4.338M-node geometry의 resource preflight를 시작했다. `accepted_0035.json`은 아직 없다. 두 branch는 후보와 accepted report를 계속 구분한다.

### rank-2 accepted step 14 및 4-step 재현성 — 2026-10-19 20:35 KST

428614는 step 14 (`t=.035,h=.0025`)를 `740.68844 s`에 수락했다. 1-rank control의 같은 step은 `877.85657 s`였다. 연속 step 11–14 평균 wall은 rank 1 `861.96847 s`, rank 2 `731.71178 s`; 2 H200/16 CPU allocation이 1 H200/8 CPU보다 약 `15.11%` 빨랐다.

동일 `t=.035`의 `ids`, `mass`, `state`, `factors`, `Psi` long-double 값, `weight_rate_seed`가 모두 numeric-equal이고 header의 차이는 run token뿐이다. 후보 admissibility는 mass bitwise, 양의 minimum thermal, P/E/L 및 combined-energy roundoff bound를 통과했다. 4개 연속 accepted step에서 rank 1/2 결과 일치가 확인됐다. h=.00125는 마지막 accepted `step=34,t=.0425`로 `.05`까지 6 step 남았으며, `t>6` 검증은 미완료다.

### step 15 / 35 계산 현황 — 2026-10-19 20:38 KST

마지막 scheduler 조회에서 428614 rank-2와 428542 fine은 계속 `RUNNING`이었다 (elapsed `54:19`, `3:31:49`). Rank-2는 step 15 provider geometry를 시작해 iteration 0을 기록했으나 아직 midpoint audit/accepted report는 없다. Fine branch는 step 35용 RHS epoch 29 (`320.87 s`) 및 tentative `t=.043125`까지 갔으며 provider outer iteration 3까지 residual `1.16e-17`로 기록했지만 `half_audit_0035`와 `accepted_0035.json`은 아직 없다.

### fine step 35 midpoint audit — 2026-10-19 20:40 KST

마지막 scheduler 조회에서 428614 rank-2와 428542 fine은 `RUNNING` (elapsed 각각 `56:23`, `3:33:53`)이었다. Rank-2 branch의 latest accepted는 step 14 (`t=.035`); step 15 provider outer iteration 3까지 로그가 왔으나 `half_audit_0015`와 `accepted_0015.json`은 아직 없다.

Fine branch는 accepted step 34 (`t=.0425`)에서 다음 step 35 midpoint `t=.043125`의 `half_audit_0035/admissibility.json`을 생성했다. mass bitwise initial, 양의 minimum thermal, P/E/L 및 combined-energy 잔차가 bound 안이지만 `accepted=false`; 따라서 step 35는 아직 commit되지 않았다. 두 endpoint의 `.05` gate 및 `t>6` 검증은 미완료다.

### 두 branch의 다음 midpoint admissibility 확인 — 2026-10-19 20:41 KST

428614/428542는 scheduler에서 계속 `RUNNING` (elapsed `57:56`, `3:35:26`). Rank-2 step 15 half audit `t=.03625`는 fixed mass, positive thermal, P/E/L 및 combined-energy bound를 통과했지만 `accepted=false`; 이어 provider outer iteration 3 residual은 `1.05e-17`, accepted step 15은 아직 없다.

Fine step 35 half audit `t=.043125`도 같은 admissibility 기준을 통과했지만 `accepted=false`; 다음 provider resource preflight는 약 4.34M nodes다. 아직 `accepted_0035.json`은 없다.

### step 15 / 35 full-stage 진행 확인 — 2026-10-19 20:45 KST

428614와 428542는 scheduler에서 계속 `RUNNING` (elapsed `1:01:48`, `3:39:18`). Rank-2 step 15 midpoint `t=.03625`는 통과했으며 이어지는 provider 반복 2 mass-root residual은 `1.10e-17`; 아직 accepted step 15는 없다.

Fine step 35 midpoint `t=.043125`도 통과했다. Native log는 full-step tentative `t=.04375`까지 갔으나 아직 `accepted_0035.json`은 없다. Provider 반복 2 residual은 `1.11e-17`. 어느 후보도 accepted로 세지 않는다.

### step 15/35 accepted, rank scaling 및 실제 병목 감사 — 2026-10-19 20:50 KST

`428614`는 step 15 (`t=.0375,h=.0025`)를 `739.51955 s`에, `428542`는 step 35 (`t=.04375,h=.00125`)를 `826.68479 s`에 수락했다. 두 report 모두 `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, `long_time_stability_verified=false`다. 각 branch 모두 `.05`까지 accepted step 5개가 남았다. 20:50 KST scheduler 조회에서 두 작업은 계속 `RUNNING` (elapsed `1:06:11`, `3:43:41`)이었고 다음 step report는 아직 생성되지 않았다.

동일 `t=.0375`의 기존 1-rank h=.0025 control과 새 2-rank 결과를 직접 대조했다. `ids`, `mass`, `state`, `factors`, `weight_rate_seed` 파일은 byte-identical이다. `Psi` raw file bytes는 다르지만 19,945개 값을 `longdouble`로 읽은 배열은 전부 numeric-equal이며 max abs 차이는 0이다. Header는 run별 값이 달라 byte-identical하지 않다. 연속 step 11–15 평균은 rank 1 `864.331 s`, rank 2 `733.273 s`, 약 `15.16%` wall-time 개선이다. 이 비교는 GPU 1→2와 CPU 8→16을 함께 바꿨으므로 GPU만의 scaling으로 해석하지 않는다.

실제 job 설정 및 호출 경로를 읽어 느린 지점을 구분했다. native MPI/CUDA executable이 state step과 `backend.stage`를 수행하지만, `full_moment_provider_service.py`의 `FACE_STAGE`가 Python `full_moment_face_provider.py`를 통해 `full_disk_stage_memory.py`의 `build_geometry`, `prepare_states`, `records_build`, `solve_from_transport`를 실행한다. 반복 geometry/face/root 계산은 Python provider 경로에 남아 있다. 기록된 반복당 geometry는 cache warm 뒤 약 `16 s`, face-state는 약 `10 s`, mass-root는 약 `14 s`; 첫 geometry는 약 `48–50 s`다. 같은 시각 native owned stage는 rank당 약 `1.3–1.6 s`, MPI/GPU 구간은 약 `4.4 s`였고 전체 provider RHS는 약 `298 s`였다. 따라서 현 하이브리드 경로는 실제 H200/MPI를 사용하지만, timestep 전체는 Python provider의 CPU geometry/face/root 작업에 지배된다.

20:52 KST 코드 출처를 대조해 정정한다. active v10 executable은 사용자의 native C/MPI/CUDA 코드 스냅샷을 컴파일했으며 주요 Kepler/update/native-stage/provider/CPU-MPI C 파일의 SHA가 저장소 `Exam/`와 일치한다. 다만 v10의 `full_moment_evolution.c`와 `full_moment_provider.h`는 현재 저장소보다 오래되어 coupled-control 변경이 없고, 실행 provider service/face-provider도 현재 `Docs/` 버전과 다르다. 특히 두 active run은 `LAGFORCE_NATIVE_FACE_STAGE` 선택 없이 old Python face path를 사용한다. 따라서 이번 실행은 사용자 코드 계열의 명시적 midpoint baseline이지만 저장소의 최신 provider 경로 전체와 동일한 빌드는 아니다.

native face/root 구현은 새로 만들 필요가 없다. 사용자의 `Docs/native_full_disk_stage.py`가 기존 물리 stage를 감싸 native straight-face 및 transport-root 함수를 선택하며, 저장소 `Docs/full_moment_face_provider.py`의 `LAGFORCE_NATIVE_FACE_STAGE=1` opt-in으로 이미 연결돼 있다. 고정 저장 상태 `t=.025`, 19,945 cells에서 실행한 기존 full-provider 등가성 보고서는 node inputs/extra hash, 전체 closed-root iteration history, `stage`, `incoming_Q`, offsets, generator, weight rates, `Psi_rates`, ledgers, mass, IDs 및 initial-state hash가 baseline과 모두 exact-value equal임을 확인했다. 그 2-RHS/0-step replay의 geometry+root 구간은 `2.4590x`, 전체 provider는 `1.5578x` 빨라졌다. 이건 저장 상태 backend equivalence이지 시간 적분이나 장기 안정성 증거는 아니다. Frozen module bundle의 `native_full_disk_stage.py`는 저장소판과 journal observer 래퍼만 다르고, explicit run에선 journal을 요청하지 않으므로 기존 등가성 범위와 충돌하지 않는다.

구 old Python 경로 accepted-step 시간을 단순 외삽하면 rank-2 h=.0025 현재 `t=.0375→6` 약 `20.2`일, h=.00125 현재 `t=.04375→6` 약 `45.6`일이다. 저장 상태 native face/root 성능비 `1.5578x`를 전체 step에 그대로 적용해도 대략 `13.0/29.3`일 규모이며, 이는 한 지점 replay를 사용한 계획 추산이다. 지금 잡들은 `moment_steps=20/40`, 종료 `t=.05`인 시간수렴 관문이다. 이 둘은 건드리지 않는다. `.05` density/physical gate 뒤의 연장은 새 native face 구현을 발명하는 대신 저장소의 검증된 opt-in 코드를 versioned provider snapshot으로 실행한다: 현재 `Docs` provider files, proven native module bundle, `LAGFORCE_NATIVE_FACE_STAGE=1`, 필요한 `PYTHONPATH`를 job identity에 묶고 exact restart checkpoint에서 provider/backend identity와 accepted gate를 다시 확인한다. coupled-control solver를 택하는 경우에는 v10 binary 대신 현재 `Exam/`의 coupled source를 대상으로 해야 한다. 더 큰 timestep은 h-refinement와 물리 오차 gate 없이 채택하지 않는다.

### CUDA ledger 및 MPI stage 데이터 이동 최적화 — 2026-10-19

사용자 CUDA bridge `Docs/moment_owned_bridge.cu`의 `lagmoment_ordered_ledgers`는 wall 이외의 약 4.3M node도 GPU 단일 thread에서 훑고, 전체 flux/owner 및 stage/hydro 배열을 다시 GPU로 복사했다. ledger에 기여하는 항은 외부 wall node, frozen cell, active cell뿐이다. CUDA 경로를 이 세 그룹의 필요한 5개 operand만 canonical index 순서대로 compact 전송하도록 바꿨다. GPU 합산은 기존 left-to-right 순서를 유지하며, 물리식·ledger 정의·CPU 경로는 바꾸지 않는다. `LAGFORCE_PROFILE_LEDGER=1`은 compact 개수, 전송 byte 및 함수 시간을 기록한다.

`Exam/full_moment_native_stage.c`는 각 rank가 계산한 disjoint cell 범위를 rank 0으로 `MPI_Gatherv` 세 번 보내고 있었다. 실행기는 모든 rank가 같은 host인지 이미 확인하고 출력 파일은 `/dev/shm`의 공유 MAP_SHARED 파일이므로, rank 0이 파일을 만들고 각 rank가 자기 canonical slice를 직접 쓰도록 바꿨다. 기존 provider-epoch barrier 뒤 ledger 계산이 시작되며 output ordering은 그대로다. 이는 같은 노드의 shared mapping을 전제로 하므로 multi-node 실행에 적용하지 않는다.

두 수정은 source snapshot `/gpfs/kjhan/LagForce/User_full_disk_mpi_compact_ledger_build_v1_20261019`에 고정했다. CUDA library SHA-256 `45ab78db7d4ff1b81ab1aadbb72635deedef6f4885870abf6f6ee21919ab2d5e`; MPI stage object SHA-256 `5acc661dcc4a680466c85a9c135df07c1f4b23e20ecf987fb61bbd9af801d346`. nvcc 및 Intel MPI C 컴파일을 거쳐 MPI/CUDA executable도 link했고, SHA-256은 `0a4be604ea5fe8f4ca800d757a6ca9c39bcc721515d510cd272d740e049f6985`다. 별도 runtime provider snapshot `User_full_moment_native_face_provider_runtime_mpi_compact_v1_20261019`의 46개 provider/module hash도 모두 일치한다. H200 저장 capture의 ordered-ledger CPU/CUDA 동등성 작업 `428717`은 `QOSGrpGRES`로 대기 중이며, 컴파일/link 성공만으로 수치 동등성이나 속도 향상을 주장하지 않는다.

최근 산출물에서 기존 기준 작업은 `428542` fine step 37, `t=.04625`, `h=.00125` 및 `428614` rank-2 step 17, `t=.0425`, `h=.0025`까지 각각 accepted다. 마지막 scheduler 조회에서도 두 작업은 `RUNNING`, component job `428717`은 `QOSGrpGRES` pending이었다. 기준 작업은 원래 v10 실행과 기존 Python face 경로를 유지하며 `.05` 물리/density gate 전이다. source 수정이나 대기 중인 CUDA component check가 이 잡들을 변경하지 않았다. 다음 적용 전 조건은 component ledger bitwise check, MPI direct-mapped output의 actual two-rank restart equality, 그리고 `.05` checkpoint physical/temporal gates다.

### 재시작 provenance의 CUDA operator 버전 고정 — 2026-10-19

장기 continuation용 `Docs/full_moment_restart_provenance.py`가 CUDA operator의 구 v4 경로/SHA만 허용해 compact-ledger candidate로 봉인·검증할 수 없는 것을 확인했다. 기본값은 구 v4 그대로 유지하고, `EUNHA_MOMENT_OPERATOR_PATH`와 `EUNHA_MOMENT_OPERATOR_SHA`가 함께 지정되면 그 버전의 실제 파일 SHA를 검사하고 manifest에 경로와 SHA를 기록하도록 보완했다. 경로 override가 있을 때도 SHA를 생략하면 구 v4 기대값으로 검증되므로 실패한다.

기존 경로/SHA manifest는 경로 필드가 없더라도 구 v4 기본 operator에서만 계속 검증되도록 하위 호환을 유지했다. 이 helper 수정은 현재 실행 중인 job의 snapshot에는 반영되지 않는다. Python AST 및 구 v4/candidate operator identity 검증만 수행하며, physics step은 실행하지 않았다. 실제 restart seal/import는 `.05` 물리·density·temporal gate가 끝난 뒤 새 versioned provider snapshot에서 수행해야 한다.

### Versioned provider snapshot v2 for candidate restarts — 2026-10-19

기존 provider snapshot은 보존하고 새 경로 `User_full_moment_native_face_provider_runtime_mpi_compact_v2_20261019`를 만들었다. 46개 source hash를 모두 확인했으며, 변경된 restart helper SHA-256은 `d75efb6328379f55c41dedb35ab5956486ff8b95a6b2f25bf116a56ccbba9929`이다. 이 provenance의 MPI/CUDA executable SHA-256 `0a4be604ea5fe8f4ca800d757a6ca9c39bcc721515d510cd272d740e049f6985`와 candidate operator SHA-256 `45ab78db7d4ff1b81ab1aadbb72635deedef6f4885870abf6f6ee21919ab2d5e`를 실제 파일에서 다시 확인했고 snapshot 안 helper의 operator identity 검사도 통과했다.

이는 continuation에 필요한 입력 snapshot을 고정한 것이다. 아직 이 binary/provider로 time step이나 MPI restart equality를 실행하지 않았고, active `.05` 기준 잡에도 반영하지 않았다. 적용은 기존 physical/density/temporal gate와 ordered-ledger component check 이후에 한다.

### 재시작 호환성 확인 및 `.05` 기준 분기 진행 — 2026-10-19 21:31 KST

수정한 helper를 기존 `h=.0025`, `accepted_0020`, `t=.05` 구형 manifest에 실제 적용해 `verify()`를 통과시켰다. manifest에 새 `operator_path` 필드가 없어도 기본 v4 operator SHA `51014e...`로만 성공하는 하위 호환을 확인했다. candidate 라이브러리에 대해서는 별도로 정확한 경로와 SHA가 새 manifest identity로 선택되는 것을 확인했다.

마지막 `sacct` 확인에서 `428542`와 `428614`는 계속 `RUNNING`, `428717`은 `QOSGrpGRES` pending이었다. 최신 승인 checkpoint는 fine `step=38,t=.0475,h=.00125` (약 `874.67 s/step`) 및 rank-2 `step=18,t=.045,h=.0025` (약 `739.54 s/step`)다. 두 branch 모두 `.05`까지 각 2 step 남았다. 로그는 다음 provider RHS 계산을 보여주며 실패 표시는 없다. 후보 CUDA ledger check 결과는 아직 생성되지 않았다.

### Runtime snapshot identity-launcher 보완 — 2026-10-19

v2 snapshot 대조 중 runtime identity launcher `run_full_moment_with_identity.py`가 빠진 것을 발견했다. v2는 보존하되 실행 입력으로 사용하지 않는다. 현재 저장소판 launcher를 추가한 완전한 v3 snapshot `User_full_moment_native_face_provider_runtime_mpi_compact_v3_20261019`을 만들고 47개 source SHA, candidate executable SHA, CUDA operator SHA를 다시 검증했다. launcher도 `EUNHA_MOMENT_OPERATOR_SHA` override를 읽으며, provider의 `LAGFORCE_NATIVE_FACE_STAGE=1` 선택 및 restart helper가 모두 포함된다. v3 runner SHA-256은 `39ef5a5d696bff84c901999ce0cb2d15b6600b49b15fe531b5f22a6f8881a425`이다.

v3는 정적 identity 검증을 통과한 완전한 실행 snapshot이지만, 실제 MPI runtime identity와 restart equivalence는 아직 미검증이다. `.05` 기준 gate와 CUDA ledger component gate 후 이 v3 경로를 사용할 것.

### 물리 체크포인트 연속 구간 준비기 — 2026-10-19

`Docs/prepare_full_moment_restart_segment.py`는 사용자가 승인한 체크포인트/manifest와 initial/physics 파일을 먼저 검증하고, v3 provider snapshot 및 binary/operator SHA를 확인한 뒤 새 run directory에 `params.dat`, `submit.sh`, `segment_plan.json`을 만든다. 생성된 H200 job은 `LAGFORCE_NATIVE_FACE_STAGE=1`로 native face/root 경로를 선택하고, endpoint physical/density gate가 포함된 최종 accepted step을 확인한 뒤 다음 restart manifest를 봉인한다. 기존 디렉터리는 거부한다. 준비기는 job 제출이나 적분을 하지 않는다.

설정 smoke는 기존 승인 `h=.0025,t=.05,step=20` checkpoint를 잠깐 봉인해 후보 operator로 검증한 뒤, 새 `t=.15,step=60,ranks=2` 구간 설정을 `/tmp` 임시 경로에서 생성했다. absolute step, 2-rank/2-H200 설정, native-face 선택, endpoint seal 명령을 확인했고 생성된 `submit.sh`에 `bash -n`이 통과했다. 이 smoke는 provider call 0, accepted step 0이며 후보 물리 정확도나 MPI runtime restart 동등성 증거가 아니다.

### `.05` 물리 관문 직전 상태 — 2026-10-19 21:38 KST

마지막 scheduler 조회에서 `428542` fine 및 `428614` rank-2는 계속 `RUNNING`, `428717` ledger component 검증은 계속 `QOSGrpGRES` pending이었다. 최신 accepted state는 두 branch 모두 `t=.0475`: fine `step=38,h=.00125` (마지막 accepted step `874.67 s`), rank-2 `step=19,h=.0025` (`738.15 s`). `.05`까지 fine 2 step, rank-2 1 step 남았다. `.05` physical/density/time-convergence 보고서는 아직 없다. 다음 Slurm 확인 전 native/provider 로그는 두 분기의 provider 작업이 진행 중임을 보였고 오류는 기록되지 않았다.

### `.05` 관문 계산 진행 — 2026-10-19 21:40 KST

428542/428614는 모두 계속 실행 중이고 428717은 `QOSGrpGRES` pending이다. fine branch 최신 accepted는 `step=38,t=.0475`; 다음 step 39의 midpoint `t=.048125`와 provider 계산이 진행 중이다. rank-2 최신 accepted는 `step=19,t=.0475`; 다음 step의 provider RHS 계산이 시작됐다. 둘 다 오류 없이 진행 로그가 있으며 `.05` accepted report와 endpoint physical/density audit은 아직 없다. 후보 component test 결과도 아직 없다.

### 재시작 manifest에 source 실행 신원과 params 봉인 — 2026-10-19

`Docs/full_moment_restart_provenance.py`의 새 seal은 source run의 `runtime_identity.json`과 `params.dat`도 검증한다. identity는 native exit 0, 런처 통과, ranks/GPU UUID, 실제 실행 binary SHA, rank별 mapped CUDA operator SHA, 실제 provider snapshot의 source hashes를 확인한다. `params.dat` SHA도 런처가 기록한 값과 비교한다. 새 manifest에는 source binary/operator/provider identity 및 runtime identity/params 파일 경로와 SHA가 들어간다. 과거 manifest는 기존 default v4 operator에서의 backward-compatible 검증을 유지한다.

v4 provider snapshot `User_full_moment_native_face_provider_runtime_mpi_compact_v4_20261019`에는 이 helper를 반영해 source map 47개 SHA를 갱신했고, launcher·executable identity도 재검증했다. 현재 restart helper SHA-256은 `0022469012f3e4b9f198b5701ee02af99206675d3bf4858dd84ed426760565c3`이다. `.05` 기존 checkpoint를 봉인할 때 source runtime identity가 binary `d9c79e...`, operator `51014e...`, 한 rank와 실제 H200 UUID를 증명하고 36개 provider source hashes와 `params.dat` SHA가 일치했다. 그 manifest로 후보 v4 snapshot을 사용한 `.15/step=60`, 2-rank continuation 설정을 다시 생성했고 generated `submit.sh`의 `bash -n`이 통과했다. 이 작업은 provider call 0, accepted step 0, Slurm submission 0이다.

### Continuation 준비기 source identity gate 및 최신 진행 — 2026-10-19 21:45 KST

`prepare_full_moment_restart_segment.py`는 이제 source manifest에 runtime identity, source params, binary/operator/provider identity 필드가 없으면 구간을 만들지 않는다. 따라서 레거시 manifest는 관찰/호환 검증은 가능하지만 새 장기 분기 입력으로 바로 사용할 수 없고, 최신 helper로 identity-bound seal을 먼저 만들어야 한다. v4 provider snapshot에서 새 source identity/params seal 및 설정 생성 smoke와 bash 문법 확인이 통과했다.

21:45 KST scheduler 조회에서 428542/428614는 여전히 `RUNNING`, 428717은 `QOSGrpGRES` pending이었다. fine 최신 accepted는 `step=39,t=.04875,h=.00125`; rank-2는 `step=19,t=.0475,h=.0025`다. 두 분기 모두 `.05`까지 accepted step 1개가 남아 있다. endpoint 물리·밀도 보고와 ledger 결과는 아직 없다.


### 2026-10-19 22:08 KST — native H200 t=.05 timestep-refinement result

The user directed the next work to use their native code. The completed refinement runs used the existing `eunha2.exe` v10 binary (`d9c79ec0...`), native MPI time control, the owned GPU stage/CUDA operator, two actual H200 ranks where requested, fixed mass and Re=∞. Python remained the provider/orchestration process for the audited geometry and physical RHS; this was not a Python-only integrator. The h=.00125 run (`428542`) and h=.0025 single-rank reference (`428541`) both started from the same accepted t=.025 state. Both completed 20 interval steps to t=.05. The independent 2-rank h=.0025 run (`428614`) also completed to t=.05.

At t=.05 both h=.0025 and h=.00125 branches passed the existing absolute density bounds (mass-weighted RMS <.01 and maximum <.1) in every region. The cavity density error was RMS/max `0.00760656/0.0303563` for h=.0025 and `0.00699214/0.0260636` for h=.00125. Inner-contact density RMS was `0.00484692` and `0.00415882`. Both checkpoints passed their physical monitor, independent q48 density preparation, fixed-mass bitwise check, positive thermal, and P/E/L plus combined-energy roundoff ledgers. The h=.00125 minimum thermal was `2.0246703e-8`.

A saved-state comparison of the two paths over the shared `.025→.05` interval passed the unchanged temporal budgets in all eight regions: relative pressure/density/velocity error-difference RMS <.002 and phase RMS <.2 degrees. Selected values: all-cell density `7.11442e-5`, cavity density `9.50670e-4`, inner-contact density `8.98871e-4`, inner-contact pressure `1.64219e-4`, and inner-contact phase `0.0057335°`. Specific-thermal error difference relative to the reference RMS was `9.04286e-4` in inner contact; intrinsic-spin (`ell`) was `4.72177e-3` under its own very small reference scale. These are same-time temporal diagnostics, not exact-solution or long-time validation. Reports are preserved under `/gpfs/kjhan/LagForce/User_full_moment_h0025_h00125_temporal_comparison_v1_20261019`.

The separate restart equivalence test from the accepted 2-rank `.05` checkpoint first attempted the native face-stage path. Job `428752` imported the state with the unmodified user v10 executable and made no accepted step; its provider failed at the first RHS because the generated script omitted `LAGFORCE_NATIVE_STRAIGHT_FACE_LIB` (and the paired transport-root library path). This is a launch configuration failure, not a physical failure. Its logs are preserved. The dependent candidate `428753` did not run and remains `DependencyNeverSatisfied`.

`Docs/prepare_full_moment_restart_segment.py` now verifies both native provider libraries and exports their exact snapshot paths. Physics code and acceptance limits were not changed. New 2-rank H200 jobs use the same accepted source state, v10 reference versus its compact-ledger MPI/CUDA candidate, identical v4 provider sources, native face-stage selection, and a two-hour wall limit: reference `428757` was RUNNING at first query and candidate `428758` waits on `afterok:428757`. The original failed attempt and never-satisfied dependent job remain preserved. The planned comparison checks IDs, mass, state, factors, Psi, weight-rate seed, physical/density audits, and step wall time.

No result here verifies stability past t=.05. After the native-face restart equivalence and runtime comparison, continue the same fine h=.00125 path through t=.1 with the same absolute and temporal gates; use accepted checkpoints and preserve them on failure. Historical t=6–10 remains the eventual Kepler validation target.


### 2026-10-19 22:42 KST — user v10 native continuation submitted from t=.05 to .1

The user explicitly directed that their existing code be used. The continuation therefore reuses the accepted h=.00125 checkpoint `Codex_full_disk_explicit_fine_timestep_refinement_fork_v2_20261019/run_h00125/runtime/accepted_0040` at t=.05, the user's native v10 executable SHA-256 `d9c79ec0f868ebb3f39f0e95ec3777465d0055b4fab079e076ae59e4cf619760`, and the same CUDA operator SHA-256 `51014e2fb7de34bb8eb1ed99c8ed0e22ecbbe5690a67327f8fec9892d5968471`. Fixed mass, Re=∞, no steering, explicit midpoint, and h=.00125 remain unchanged. The native C/MPI time controller and CUDA-owned stage perform time integration; Python remains the geometry/physical RHS provider and orchestration layer. The tested v10 reference run used the same native-face provider family, with `LAGFORCE_NATIVE_FACE_STAGE=1` and its straight-face/root shared libraries explicitly exported.

The accepted t=.05 state was sealed in the new experiment root `/gpfs/kjhan/LagForce/User_full_moment_kepler_nativeface_t005_to_t01_v1_20261019`. Seal verification bound the accepted checkpoint, original initial/physics files, the producing run's native runtime identity, and its audits. The prepared segment targets absolute step 80 (t=.1), i.e. 40 more h=.00125 steps, using one H200, one MPI rank, 8 CPU threads, 128 GiB, and an 18-hour limit. The generated submit script passed `bash -n`; no physics code or timestep was altered.

Slurm accepted job `428765`. At submission it had not yet been queried for queue state. Its run directory and plan are in `/gpfs/kjhan/LagForce/User_full_moment_kepler_nativeface_t005_to_t01_v1_20261019/run_h00125_to_t01`; `submission.json` stores the sealed source, hashes, and submitted job identity. Continue from this exact run and preserve all accepted checkpoints.

The separate v10-versus-compact-ledger one-step check is complete: jobs `428757` and `428758` both produced an accepted step to t=.0525 with matching checkpoint values and meaningful payload bytes. The only byte-level discrepancy was nonsemantic padding in 16-byte x86 long-double `Psi` storage; numeric Psi values match. Whole-step speedup was only 0.086%, so the long continuation uses the user's original v10 binary. This remains a short-horizon restart/equivalence result; `long_time_stability_verified=false`, and the historical t=6–10 wall has not been crossed.


### 2026-10-19 22:43 KST — continuation job started

The single scheduler check after submission found job `428765` `RUNNING` on `syn104`, with one H200, 8 CPUs, 128 GiB, no dependency, and the 18-hour limit. `launch.log` records native runtime rank 0/1 with h=.00125, absolute step target 80, and the expected v10 executable SHA. The runtime directory currently contains restart import and provider logs but no accepted step yet (elapsed at scheduler observation was 1:08); `runtime_identity.json` is emitted by the wrapper at completion. The run is progressing and no restart or duplicate submission was made.

The Slurm accounting check also confirms the reference/candidate one-step equivalence jobs `428757/428758` both ended `COMPLETED 0:0` in 13:34/13:31. Their exact numeric checkpoint equivalence and padding-only Psi byte difference are recorded in the preceding section and in `User_full_moment_mpi_compact_step_equivalence_v1_20261019/retry_submission.json`.


### 2026-10-19 22:57 KST — h=.00125 accepted step and matched rank-scaling run

The continued one-rank user-v10 run `428765` accepted step 41 at t=`0.051250000000000004`, h=.00125, after two physical provider evaluations. Its accepted-step wall was `641.4426 s`; `fixed_mass=true`, `Re_infinity=true`, and `long_time_stability_verified=false`. The intermediate half-state at t=.050625 passed its native tentative admissibility check: mass remained bitwise equal to initial, minimum thermal was positive, and all four P/E/L residual components plus combined energy were inside their recorded roundoff bounds. This half-state was not counted as accepted; the later native commit report is authoritative. The native log confirms H200 NVL UUID `cc32769c2315f4a28e744ad3a43159e6` and the pinned user v10 executable.

To choose the most efficient production rank count at the fine timestep, a matched one-step scaling job `428770` was submitted from the exact same accepted t=.05 checkpoint/manifest, target t=.05125, h=.00125, using the same v10 executable, CUDA operator, and native-face provider, but 2 MPI ranks / 2 H200 GPUs / 256 GiB. It is independent of the running one-rank continuation `428765`; no dependency was set. The 2-rank segment has a two-hour limit. Compare accepted physical state and endpoint gates as well as `step_wall_seconds` before changing the long continuation's rank count.

The one-rank continuation remains targeted at t=.1 (absolute step 80). No long-time stability claim is supported at these times.


### 2026-10-19 23:00 KST — step-42 progress and rank-2 queue state

At the latest scheduler observation, `428765` remained RUNNING on `syn104` with one H200. The checkpointed native state is step 41/t=.05125 (`641.4426 s`); the next step had begun provider call 3, whose fourth geometry iteration reported weight-rate change `8.6046e-16`. There was no accepted step 42 yet.

The matched h=.00125 2-rank benchmark `428770` remains `PENDING (QOSGrpGRES)` and has no allocation or runtime output. The same queue query showed the main job continues; no job was canceled or restarted. Keep the rank-2 benchmark pending until Slurm grants its requested two H200 GPUs, then compare it to the already accepted one-rank step 41 from the exact same t=.05 source state.


### 2026-10-19 23:02 KST — step-42 half-state admissibility

The latest observed scheduler state was `428765 RUNNING` on `syn104`, elapsed 19:41, and `428770 PENDING (QOSGrpGRES)` with no allocation. The native log records step 42 provider call 3 at 223.72 s and its tentative midpoint time t=.051875. `half_audit_0042/admissibility.json` passed native tentative admissibility: mass bitwise equal to initial, minimum thermal `2.02466996e-8`, all four P/E/L residuals inside their respective bounds, and combined-energy change `1.80e-17` inside its `3.48e-11` roundoff bound. It explicitly reports `accepted=false`; step 42 has not yet committed. Do not use this half-state as a restart point.

The 2-rank matched scaling benchmark remains queued for two H200 GPUs. Preserve its pending job and continue the one-rank t=.05→.1 trajectory.


### 2026-10-19 23:09 KST — fine native continuation accepted step 42

At the scheduler observation (`428765 RUNNING` on `syn104`, elapsed 27:12), the fine user-v10 continuation had accepted both step 41 at t=.05125 (`641.4426 s`) and step 42 at t=.0525 (`683.2838 s`). Each step used two physical provider evaluations; cumulative providers=4. Step 42 has `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, while `long_time_stability_verified=false`. The half candidate t=.051875 and full candidate t=.0525 each passed native tentative admissibility with bitwise fixed mass, positive thermal, and P/E/L plus combined-energy residuals within their reported bounds. Their `accepted=false` audit flags identify tentative diagnostics; `accepted_0042.json` and the native commit line confirm the actual step commit.

The matched 2-rank/2-H200 h=.00125 one-step benchmark `428770` is still `PENDING (QOSGrpGRES)` without allocated resources or runtime. The continuous one-rank branch remains active with target t=.1.


### 2026-10-20 — dependent native continuation t=.1 to .15 submitted

To continue the validated fine path without restarting from t=0, a serial dependent segment `428781` was submitted with `afterok:428765`. It requests one H200, one MPI rank, 8 CPUs, 128 GiB, and 18 hours. Only after job 428765 exits successfully will the driver verify its exact source `runtime/accepted_0080` plus `restart_seal/manifest.json`, create a provenance-verified segment, and run from t=.1/step80 to t=.15/step120 at unchanged h=.00125. It uses the same user-v10 executable, operator and v10-reference native-face provider. The generated builder is a hashed copy (`bec96fe5…`); the driver passed `bash -n` and is sealed in `chain_t010_to_t015_plan.json`. No physics source was modified and no step has run in this dependent job yet.

The independent 2-rank/2-H200 matched scaling job `428770` remains pending `QOSGrpGRES`; it is preserved while the rank-1 physical path advances.


### 2026-10-19 — current fine branch and chained segment queue check

At scheduler observation 2026-10-19 23:13 KST, `428765` remained RUNNING on `syn104`, elapsed 31:38, with the user v10 rank-1/H200 path. Accepted state remained step 42/t=.0525. Step 43 had completed provider call 5 and logged the tentative half time t=.053125; no accepted step 43 was present yet.

The rank-scaling comparator `428770` remained `PENDING (QOSGrpGRES)`, with no allocation. The continuation job `428781` (`t=.1→.15`) was `PENDING (Dependency)` on `afterok:428765`, with no allocation. Thus the t=.1→.15 segment cannot start before the t=.1 segment exits successfully and creates its sealed source checkpoint.


### 2026-10-19 23:25 KST — accepted progress, provider cost, and executable lineage audit

One scheduler check found `428765 RUNNING` on `syn104`, elapsed 43:52. The same H200 rank-1 run has accepted step 43 at t=.05375 (`664.953 s`); step 44 is still tentative at midpoint t=.054375. Its `half_audit_0044/admissibility.json` passed tentative checks: mass is bitwise unchanged, minimum thermal is positive, and P/E/L plus combined-energy residuals are inside their recorded bounds. `accepted=false` remains authoritative until the native commit file is written. The matched 2-rank job `428770` remains pending `QOSGrpGRES`; dependent `428781` remains pending on `afterok:428765`.

The last two native RHS RPCs took 202.895 s and 227.549 s. Their stage reports separately recorded 150.747 s and 176.016 s; each report covers three or four face/root outer iterations over roughly 4.34 million nodes. The provider uses seven CPU workers. Per outer iteration, `face_state` took about 9.9 s and the mass root about 4.8–5.8 s; cached geometry after the first iteration took about 1.8–1.9 s, versus 34–35 s for the first geometry pass. The owned CUDA stage took about 2.25 s of kernel time and 5.1 s including its MPI/GPU stage. This confirms that the slow portion is chiefly the repeated CPU provider work, not the CUDA kernel. These timings are diagnostics, not a reason to weaken any physics gate.

The user asked that the user-authored code remain the basis. An exact lineage audit shows the active executable uses the checkout's `eunha2.c` byte-for-byte (`b704041c…`) and a pinned earlier FullMoment native-stage snapshot. However, the active `full_moment_native_stage.c` hash (`8038edb3…`) differs from the current checkout (`8c939495…`); the current checkout's `full_moment_evolution.c` also has the optional coupled-control path absent from the v10 binary. The v10 link contains the FullMoment C objects and `liblagmoment_owned.so`, but does not link `Exam/exam_gpu.cu`/`exam_gpu.o`. Its CUDA stage is the separate owned-moment operator (`moment_owned_bridge.cu`); geometry and physical RHS still come through the Python provider. Therefore this is a pinned user-EUNHA2/FullMoment baseline, not the complete current repository GPU path. Preserve this distinction in reports.

Do not interrupt or reinterpret the already-running explicit-midpoint trajectory. At its sealed t=.1 endpoint, compare the current repository source and the pinned v10 baseline from the same accepted checkpoint before selecting a continuation binary. Keep equations, fixed-mass/Re=∞ configuration, timestep, and acceptance gates identical for that comparison; do not substitute the legacy `Exam/exam_gpu.cu` force kernel unless its FullMoment state/flux mapping is shown to be equivalent. Prioritize the CPU provider path for speed work, while preserving the current accepted checkpoint chain.


### 2026-10-19 23:31 KST — fine continuation accepted step 44

The latest scheduler observation found `428765 RUNNING` on `syn104`, elapsed 49:24. It committed step 44 at t=`0.055`, h=`0.00125`, with 2 physical RHS evaluations and `step_wall_seconds=686.206645`. `accepted_0044.json` confirms `accepted_commit=true`, `fixed_mass=true`, `Re_infinity=true`, and `long_time_stability_verified=false`. The saved tentative admissibility report has minimum thermal `2.0246691e-8`; pressure RMS/max are `4.12686e-5/4.73772e-4` globally and `6.41827e-5/2.65911e-4` in the inner region. Its P/E/L residuals and combined-energy change are within the recorded roundoff bounds. The report's `accepted=false` is the pre-commit audit flag; `accepted_0044.json` is the subsequent native commit evidence.

The run has now advanced four accepted fine steps from t=.05 through t=.055. It has not yet reached the t=.1 endpoint density/physical gates or restart seal. Job `428770` remains queued for `QOSGrpGRES`; dependent continuation `428781` remains pending on `afterok:428765`. Preserve both jobs and all accepted state files.

### 2026-10-19 23:51 KST — accepted step 46 and current-checkout replay queue

At the scheduler observation, `428765` remained `RUNNING` on `syn104` (elapsed `1:11:26`). The pinned user-v10 fine branch committed step 46 at `t=0.0575`, `h=0.00125`, after 2 physical provider evaluations; `step_wall_seconds=662.419996`, `fixed_mass=true`, `Re_infinity=true`, `accepted_commit=true`, and `long_time_stability_verified=false`. The preceding midpoint at `t=0.056875` passed tentative admissibility with positive minimum thermal `2.02466872e-8`, bitwise fixed mass, and P/E/L plus combined-energy residuals inside their recorded bounds; that midpoint itself had `accepted=false` and is not an accepted state.

The first runtime replay of the current repository build is job `428805`: current checkout `eunha2.c` SHA `b704041c…`, current-checkout executable SHA `545b4280…`, current provider snapshot, the existing owned CUDA operator, and the exact sealed accepted `t=.05/step40` checkpoint; it requests 2 H200 ranks for one step to `t=.05125`. At the same scheduler observation it was still `PENDING (QOSGrpGRES)` and had produced no runtime output. The matched v10 rank-2 benchmark `428770` was also `PENDING (QOSGrpGRES)`; dependent continuation `428781` remained `PENDING (Dependency)` after `428765`. No existing job was changed or canceled.

The accepted fine-step wall time implies roughly 36.4 days from `t=.0575` to `t=6` if held constant (`4754` steps at `662.42 s/step`); this is a simple planning estimate, not a forecast or stability result. It makes the current Python CPU geometry/physical-RHS path a major optimization target. The pending current-checkout replay has not yet tested runtime or speed, and no result here establishes long-time stability or that the derived rigid-spin projection reduces Kepler shear error.

### 2026-10-19 23:57 KST — step 47 midpoint and H200 group quota

At the 23:57 scheduler query, `428765` remained `RUNNING` on `syn104` (elapsed `1:16:32`). Step 46/t=.0575 remains the latest accepted checkpoint. The step 47 midpoint at `t=.058125` passed native tentative admissibility: minimum thermal `2.02466841e-8`, mass bitwise equal to initial, and all P/E/L plus combined-energy residuals within the saved bounds. `accepted=false`; no step 47 endpoint has been accepted.

Jobs `428805` (current-checkout 2-rank replay) and `428770` (matched v10 2-rank benchmark) remain `PENDING (QOSGrpGRES)`; `428781` remains dependency-pending. A read-only `sacctmgr show qos` query reports the `part_h200` QoS group limit as `gres/gpu=6`. The queue reason is consistent with the H200 group limit, so the replay has not started despite the requested H200 nodes; no alternate GPU partition was substituted and no jobs were canceled or altered. Recheck the queue before changing resource requests.

### 2026-10-20 00:00 KST — step 47 full candidate awaiting commit

A subsequent scheduler check found `428765` still `RUNNING` on `syn104` (elapsed `1:21:39`). Step 47's second provider call completed in `226.749748 s`; the native log then created a tentative full-step candidate at `t=0.05875`, `h=0.00125`. `factor_15_rank0.json` exists, but `candidate_audit_0047/` and `accepted_0047.json` do not yet exist, so this candidate is not accepted. The `t=0.058125` half-state passed its admissibility gate as recorded above. `428805` and `428770` remain `QOSGrpGRES` pending under the observed six-GPU `part_h200` group cap; no job or partition was changed.

### 2026-10-20 00:04 KST — accepted fine step 47

`428765` remained `RUNNING` on `syn104` (elapsed `1:22:45`). It committed step 47 at `t=0.05875`, `h=0.00125`, after 2 provider evaluations; wall time was `689.784474 s`. The accepted report confirms fixed mass, `Re_infinity=true`, and `accepted_commit=true`; `long_time_stability_verified=false`. The endpoint admissibility audit has positive minimum thermal `2.02466816e-8`, bitwise fixed mass, and P/E/L plus combined-energy residuals inside the recorded roundoff bounds. Its `accepted=false` describes the pre-commit audit; `accepted_0047.json` is the authoritative subsequent commit.

The last observed fine-step rate now implies about 37.9 days from `t=.05875` to `t=6` at constant `689.78 s/step` (4753 steps). This reinforces that runtime acceleration remains necessary; it is not a stability forecast. Current-checkout replay `428805` and matched rank-2 run `428770` were still `QOSGrpGRES` pending, with `part_h200` configured at a six-GPU group limit. Neither has runtime output, and no alternate partition was used.

### 2026-10-20 00:05 KST — measured RHS time split identifies CPU target

The accepted step 47's second RHS was at `t=0.058125`. Its saved `rhs_14.json` reports `total_seconds=174.398` and `compiled_seconds=6.037`. Across four outer iterations, geometry took `34.666 s` on the first pass and `1.79–1.93 s` on cached passes; face-state work took `9.89–9.94 s` per pass; mass-root solves took `4.79–5.68 s` per pass. The enclosing native RHS call took `226.750 s`. Its owned H200 stage used `2.264 s` kernel time and `5.156 s` for the logged MPI/GPU stage.

Thus the GPU stage kernel is about 1% of that RHS call; kernel-only tuning cannot materially change the roughly 690-second accepted timestep. The next performance work should profile and reduce the CPU geometry, face-state, mass-root, and provider-task/serialization costs while preserving the same stage equations and acceptance gates. The `52.35 s` difference between the stage report and enclosing RPC is a scope difference requiring its own timing before attributing it to one subtask. Current-checkout job `428805` remains queued, so these timings describe the active v10 baseline only.

### 2026-10-20 00:09 KST — step 48 midpoint admissibility

At the scheduler observation, `428765` remained `RUNNING` on `syn104` (elapsed `1:27:49`); step 47/t=.05875 remains the latest accepted endpoint. Step 48's midpoint at `t=0.059375` passed tentative admissibility: minimum thermal `2.02466810e-8`, mass bitwise equal to initial, and P/E/L plus combined-energy residuals inside their saved bounds. It reports `accepted=false`, so the midpoint is not an accepted checkpoint. The full step to `t=.06` is still pending. Current-checkout replay `428805` and matched 2-rank benchmark `428770` remain `PENDING (QOSGrpGRES)`; `428781` remains dependency-pending.

### 2026-10-20 00:14 KST — accepted fine step 48

`428765` remained `RUNNING` on `syn104` (elapsed `1:34:08`). It committed step 48 at `t=0.06`, `h=0.00125`, after 2 provider evaluations; accepted-step wall time was `660.791382 s`. The report confirms `fixed_mass=true`, `Re_infinity=true`, and `accepted_commit=true`; `long_time_stability_verified=false`. The midpoint at `t=0.059375` and the endpoint admissibility report both passed their gates, with bitwise fixed mass, positive minimum thermal (`2.02466810e-8` at midpoint; `2.02466785e-8` at endpoint), and P/E/L plus combined-energy residuals inside their saved bounds. The endpoint report's `accepted=false` is its pre-commit audit flag; `accepted_0048.json` confirms the subsequent commit.

The fine path has advanced 8 accepted steps from t=.05 to .06. The current-checkout replay `428805` and matched rank-2 benchmark `428770` remain `PENDING (QOSGrpGRES)`; dependent `428781` remains pending after `428765`. This is still far short of the t>6 stability objective.


### 2026-10-20 00:31 KST — accepted step 49 and scoped native-face OpenMP benchmark

The fine user-v10 path committed step 49 at t=0.06125 with h=0.00125 and `step_wall_seconds=665.259530`; fixed mass and Re=∞ are confirmed, while `long_time_stability_verified=false`. At the latest scheduler observation, job `428765` remained RUNNING on `syn104`; jobs `428770` and `428805` were still `PENDING (QOSGrpGRES)`, and dependent `428781` remained dependency-pending.

A code audit found that `native_straight_faces.cpp` has OpenMP parallel loops in `lag_line_cuts` and `lag_line_fill` and its shared library links `libgomp`, but the persistent provider sets `OMP_NUM_THREADS=1`. Therefore those native loops were executing with a one-thread team. An isolated provider overlay is prepared at `/gpfs/kjhan/LagForce/Codex_native_face_openmp_v1_20261020`: it sets the OpenMP team to the available maximum (capped at 8) only around `native_straight_faces.evaluate`, then restores the provider's one-thread setting. The physical stage, face formula, root solver, CUDA operator, and acceptance gates are unchanged. Python syntax checks and an import/API smoke check passed; no runtime speedup is claimed yet.

Matched one-step benchmark `428822` was submitted from the same sealed t=0.05/step40 checkpoint as accepted baseline step 41. It uses the same v10 executable, one H200, h=0.00125, fixed mass and Re=∞, with only the scoped native-face OpenMP setting changed. The baseline step wall is 641.442581 s. Compare step acceptance and conserved-state gates first, then runtime; the test does not establish long-time stability.


### 2026-10-20 00:34 KST — OpenMP benchmark preflight retry

The first submission, `428822`, exited after 1 second in the identity wrapper before launching MPI because the isolated overlay had no colocated `provenance.json`; `sacct` confirms `FAILED`, exit code 1, and 3.2 MiB batch MaxRSS. Its error log remains at `launch.log`. I added an overlay provenance file and verified all 48 listed source hashes. The provider import resolves to the overlay while the service still keeps `OMP_NUM_THREADS=1`; the scoped OpenMP setting is 8 and the runtime API restores the caller setting to 1.

The same sealed-checkpoint benchmark was resubmitted as `428824`; its scheduler state has not yet been polled. No simulation result or speedup is claimed.


### 2026-10-20 00:36 KST — matched OpenMP benchmark entered native runtime

Job `428824` is RUNNING on `syn104` (elapsed 1:31). The source checkpoint manifest was verified with zero state updates, and the native runtime identity line reports the expected v10 binary hash, one rank, and h=0.00125. The provider has not started its first request yet at this observation; there is no step or speed result so far.


### 2026-10-20 00:38 KST — accepted step 50 and benchmark restart audit

The baseline job `428765` committed step 50 at t=0.0625, h=0.00125, in 683.809031 s; fixed mass and Re=∞ remain confirmed, and long-time stability remains unverified.

The matched OpenMP benchmark `428824` is still RUNNING on `syn104`. Its restart admissibility report passed at t=0.05: bitwise-equal mass, positive minimum thermal, and all P/E/L plus combined-energy residuals within bounds. The independent reference audit had passed 17,500 of 19,945 cells in 52.1 s at the latest file observation. The physical provider had not started yet, so no timing comparison is available.


### 2026-10-20 00:41 KST — OpenMP benchmark first RHS profile

Benchmark `428824` has started its first full-disk provider RHS. The first two outer iterations report geometry 34.497 s then 1.697 s, face-state 10.166 s then 9.864 s, and mass-root 5.970 s then 4.860 s. The root is still converging, so this is not a completed provider or step. The first geometry timing is close to the original one-rank profile (~34.7 s), showing the current scoped OpenMP change has not yet reduced the dominant first-pass geometry preparation; a final speed claim awaits the completed matched step and the emitted provenance.


### 2026-10-20 00:44 KST — matched first-RHS OpenMP result

At the identical t=0.05 state, the 8-thread face-kernel run and the one-thread baseline produced bitwise-identical `node_inputs` and `extra` arrays. The provider stage took 153.612157 s versus 153.470680 s (0.092% slower); the enclosing provider RPC took 205.402417 s versus 204.848240 s (0.271% slower). This is a useful equivalence result but no measurable speedup. The targeted OpenMP loops are therefore too small a share of runtime to materially advance t>6. Keep the ongoing one-step job through its endpoint gate, then prioritize the serial Python line-sample preparation and per-cell face-state work.


### 2026-10-20 00:48 KST — batched per-cell line-sample candidate prepared

The next isolated optimization batches the same five straight-face sample points for each cell before calling the existing `velocity_samples` formula. The original `native_straight_faces.prepare` still builds/caches the final line constants and the same native cut/fill and physical gates remain in place. A small 4-cell/6-face parity check produced bitwise-identical constants and zero difference against the scalar path; all 49 provider source hashes and Python syntax passed. A matched t=0.05/step40→step41 H200 benchmark is prepared but not submitted while `428824` is still running. This small synthetic parity check does not substitute for the full-state endpoint validation.


### 2026-10-20 00:53 KST — OpenMP one-step endpoint and exact state comparison

Slurm accounting confirms `428824` completed successfully in 16:54. It accepted step 41 at t=0.05125 (`step_wall_seconds=778.964336`) and passed the native runtime identity check. Its `state.bin`, `factors.bin`, `weight_rate_seed.bin`, `ids.bin`, `mass.bin`, and `Psi.bin` are each bitwise identical to baseline job `428765` step 41; only `header.bin` differs because restart run tokens/epochs differ. Thus this scoped OpenMP path preserves the accepted physical state exactly.

The measured stage had no speedup: first-RHS stage time was 153.612157 s versus baseline 153.470680 s. The whole accepted step was 778.964336 s versus baseline 641.442581 s. Both jobs ran concurrently on `syn104`, so the whole-step increase is not attributable to the code change; the measured face cut/fill optimization itself is too small to show an improvement. The long-time objective remains unverified.


### 2026-10-20 00:55 KST — user-code lineage reconfirmed; batched sample benchmark submitted

The user explicitly reiterated that their code must remain the basis. I compared the active checkout `/tmp/lageunha-push/eunha2.c` and the recorded build input `/gpfs/kjhan/LagForce/Codex_full_disk_mpi_bridge_build_v2_20261019/eunha2.c`; both SHA-256 hashes are `b704041c2c6e4754c24de4827f976bde23ff10c4042e2df94bf41bb96db5b556`. The active v10 executable hash is `d9c79ec0…`; its provenance points to the original MPI bridge object compiled from that matching `eunha2.c`. The FullMoment provider/native-stage pieces remain pinned companion code, so this must not be described as the complete current-checkout CUDA path.

Submitted H200 job `428831` for one matched step from the already sealed t=0.05/step40 checkpoint to t=0.05125. It uses the same executable, operator, physical parameters, and acceptance gates; only the per-cell batched straight-face sample preparation overlay changes. This is a performance-only candidate and does not alter equations or the continuing `428765` trajectory. The candidate passed small synthetic bitwise parity, but full-state parity and any speedup remain unverified. `sbatch` accepted the job; its queue state has not yet been polled. Existing jobs were left untouched.


#### Benchmark configuration correction for job 428831

A final script audit found that the prepared batch overlay is based on the earlier OpenMP overlay: job `428831` enables both per-cell sample batching and an 8-thread scoped native face team. The prior OpenMP-only control `428824` had bitwise-identical accepted state and no measurable first-RHS speedup, but the new run is still a combined-change comparison versus the baseline. Do not attribute any observed timing change to batching alone. If the combined run is promising, compare batching with the OpenMP team set to one before adopting the performance result. The running job files were left untouched.


### 2026-10-20 01:06 KST — batched candidate first RHS and midpoint gates

Job `428831` passed the no-update restart import at t=0.05 and entered the first full RHS. Its stage report shows `faces_batched=59364`, `cells_batched=19945`, and `scalar_velocity_calls_avoided=593640`; 32 sampled faces were bitwise equal to the scalar reference with maximum difference zero. Batch preparation took 15.421943 s. At the same t=0.05 state, enclosing provider time was 188.240289 s versus baseline 204.848240 s (−8.11%); first geometry was 18.251422 s versus 34.496922 s. The candidate midpoint t=0.050625 passed tentative admissibility and the native factor audit.

This is preliminary: the candidate has not accepted step 41, and its complete numerical arrays have not been compared with baseline. Both jobs are on `syn104`; the candidate includes both the batch overlay and scoped OpenMP=8, so timing cannot yet be assigned to batching alone. Preserve the run through its accepted endpoint gate.


### 2026-10-20 01:09 KST — matched batch candidate second RHS and resource progress

The batched candidate completed its second moving RHS at midpoint time `0.050625`: provider time `197.396460 s` versus baseline `204.474126 s` (−3.46%); stage total was `144.980284 s` versus `152.663037 s` (−5.03%). First-pass geometry was `20.235029 s` versus `34.495178 s`. Batched preparation took `17.054025 s`; 32 sampled faces again matched the scalar formula bitwise. The endpoint candidate at t=`0.05125` and its factor audit are present, but the final acceptance/global audit and accepted-state comparison are not yet complete.

At the 12:49 scheduler observation, job `428831` remained RUNNING on `syn104`; its MPI step reported AveCPU `00:41:27`, MaxRSS `23,962,724 KiB`, AveRSS `12,497,120 KiB`, so the process is making progress and is within its 128 GiB memory request. This remains a combined batch+OpenMP comparison sharing `syn104` with baseline `428765`; use provider/stage timings only as preliminary evidence until endpoint parity and a batch-only control are complete.


### 2026-10-20 01:12 KST — batch candidate endpoint pre-commit gates passed

The one-step candidate's t=`0.05125` endpoint passed tentative admissibility: minimum thermal `2.0246700214e-8`, bitwise fixed mass, and P/E/L plus combined-energy residuals inside saved bounds. The independent density reference audit also completed successfully. `admissibility.json` has `accepted=false` because it is the pre-commit audit; there is still no `accepted_0041.json` or committed checkpoint. The native commit and full numerical-array comparison against baseline remain pending.


### 2026-10-20 01:15 KST — job 428831 accepted exact step 41

Slurm accounting reports job `428831` COMPLETED, exit `0:0`, elapsed `00:16:27`, one H200 on `syn104`; live identity passed with the expected `eunha2.exe` hash (`d9c79ec0…`), owned operator hash (`51014e2f…`), and provider entrypoint. The step accepted at `t=0.05125`, `h=0.00125`, rank 1, fixed mass and Re=∞. Its physical monitor passed in every reported region. `long_time_stability_verified=false`.

All numerical checkpoint arrays (`ids.bin`, `mass.bin`, `state.bin`, `factors.bin`, `Psi.bin`, `weight_rate_seed.bin`) are bitwise identical to baseline `428765` step 41. `header.bin` differs only in restart run token/epoch. Thus the batch+OMP overlay preserves the accepted physical result exactly at this checkpoint.

The two RHS stage totals were `281.531178 s` candidate versus `306.133717 s` baseline (−8.04%); enclosing provider times were `385.636749 s` versus `409.322366 s` (−5.79%). However accepted-step wall was `754.705 s` versus `641.443 s` (+17.66%), with both jobs sharing `syn104`; that full-step comparison is not clean. The candidate also combined batching and OpenMP=8. Treat the reduced stage timings as promising, not a proven end-to-end speedup. Next isolate batch preparation with OpenMP=1 when H200 queue capacity permits.


### 2026-10-20 01:17 KST — baseline progress and batch-only queue hold

The continuing user-source baseline `428765` has accepted step 53 at t=`0.06625` (`step_wall_seconds=686.157075`); fixed mass/Re=∞ remain active and t>6 stability remains unverified. The existing rank-2 benchmark `428770` and current-checkout replay `428805` still report `PENDING (QOSGrpGRES)`. The prepared batch-only OMP=1 comparison therefore remains unsubmitted to avoid adding another H200 request ahead of the existing queue; its scripts and hashes are sealed in `/gpfs/kjhan/LagForce/User_full_moment_native_face_batch_omp1_benchmark_t005_step41_v1_20261020/`.


### 2026-10-20 01:19 KST — batch-only OMP=1 control submitted

Submitted H200 job `428845` for the same sealed t=0.05/step40 → step41 comparison. The user-source `eunha2.c`, executable, owned operator, provider overlay, checkpoint, h, and physics are hash-verified unchanged; both `OMP_NUM_THREADS` and `LAGFORCE_NATIVE_FACE_THREADS` are set to 1, isolating batching. The job requests one H200, 8 CPUs and 128 GiB; script syntax and runtime-path preflight passed. It has not yet been polled. Existing jobs `428770`, `428805`, `428781`, and continuing baseline `428765` were left untouched.


### 2026-10-20 01:24 KST — OMP=1 startup failure diagnosed; isolated retry 428846 submitted

Job `428845` exited `114:0` before a physical provider call. `provider_rank0.log` shows `Service.__init__` failed its `workers <= len(os.sched_getaffinity(0))` assertion: the script had `OMP_NUM_THREADS=1` together with `I_MPI_PIN_DOMAIN=omp`, which sized the MPI process domain to one OpenMP thread while the provider needs seven workers. The preserved failure manifest confirms `provider_attempts=0`, `source_geometry` was not initialized, and the accepted input checkpoint was copied unchanged. No physical step was attempted. Intel MPI documents `I_MPI_PIN_DOMAIN=omp` as using the OpenMP thread count to size the process domain and supports `I_MPI_PIN=0` to disable process pinning.

The retry `428846` keeps OMP and native face threads at 1, keeps seven provider workers, removes `I_MPI_PIN_DOMAIN`, sets `I_MPI_PIN=0`, and first launches a tiny MPI affinity assertion requiring at least seven visible CPUs. It uses a new output directory and the same hash-verified user-derived executable, operator, sealed t=.05 checkpoint, and physical parameters. Its queue state has not yet been polled; `428845` logs remain preserved.


### 2026-10-20 01:25 KST — fixed OMP=1 affinity preflight passed

Job `428846` is RUNNING on `syn104`. The added MPI-launched affinity preflight passed with eight visible CPUs `[2, 3, 20, 21, 38, 39, 44, 45]`, and the native runtime reported the expected executable SHA `d9c79ec0…`. At the 1:19 observation, it was still in restart audit and had not started the physics provider. Failed attempt `428845` is accounted as `FAILED 114:0` in four seconds; it made zero provider calls, and its failure manifest confirms the accepted source state was preserved.


### 2026-10-20 01:28 KST — batch-only OMP=1 restart audit passed

Job `428846` passed restart tentative admissibility and the checkpoint physical monitor at t=`0.05`; minimum thermal is `2.0246703258e-8`, mass matches bitwise, and the import receipt records zero numerical updates. The MPI rank now sees eight CPUs. At elapsed `00:03:56`, it remained RUNNING on `syn104` with `AveCPU=00:21:06` and MaxRSS about 4.5 GiB; the first physical provider RHS had not yet been logged and no accepted step exists.


### 2026-10-20 01:33 KST — isolated batch-only first RHS and midpoint pass

With scoped native-face OpenMP at one thread, job `428846` measured first RHS `188.084692 s` versus baseline `204.848240 s` (−8.18%); stage report was `136.450184 s` versus `153.470680 s` (−11.09%). First geometry was `18.265901 s` versus `34.530403 s` (−47.1%), essentially matching the batch+OMP8 run. The per-cell path batched all 59,364 faces / 593,640 scalar sample calls; 32 verification faces matched bitwise, maximum difference zero, batch preparation `15.403841 s`. This isolates batching as the source of the first-RHS geometry/runtime reduction.

The midpoint at t=`0.050625` passed its admissibility and factor gates with positive thermal, bitwise mass and roundoff-bounded ledgers. Job `428846` remains in progress; second RHS and accepted endpoint parity are not yet available.


### 2026-10-20 01:37 KST — isolated batch-only second RHS and endpoint candidate

Job `428846` completed RHS 2 at t=`0.050625` with OMP=1: provider time `189.610109 s` versus baseline `204.474126 s` (−7.27%); stage total `137.135587 s` versus `152.663037 s` (−10.17%). First-pass geometry was `18.695557 s` versus `34.495178 s`; batched preparation took `15.674038 s`, and 32 sampled faces again matched bitwise. The full step endpoint candidate at t=`0.05125` was generated, and its native factor audit passed with positive thermal.

Across both RHS calls so far, provider time totals `377.694801 s` versus `409.322366 s` (−7.73%), and stage totals `273.585771 s` versus `306.133717 s` (−10.63%). The candidate pre-commit global admissibility/physical gates and accepted checkpoint comparison have not yet completed; no step is accepted yet.


### 2026-10-20 01:40 KST — batch-only endpoint admissibility and density preparation passed

The t=`0.05125` candidate from `428846` now has a passed native tentative admissibility report and `passed_full_stage_Psi_density_preparation`. Minimum thermal is `2.0246700214e-8`, mass is bitwise initial, and P/E/L plus combined-energy residuals lie within bounds. Its physical-monitor report, native commit and accepted-array comparison remain pending; no accepted endpoint is claimed yet.


### 2026-10-20 01:44 KST — batch-only OMP=1 accepted endpoint and user-source parity

Slurm accounting confirms job `428846` COMPLETED with exit `0:0`, elapsed `00:16:24`, one H200 on `syn104`, and MaxRSS `17,139,208 KiB`. The native runtime identity check passed with no identity errors and the expected executable SHA `d9c79ec0…`. The executable was built from the user's `eunha2.c` source SHA `b704041c…`, matching the active checkout. The FullMoment provider/stage remains the hash-pinned v10 companion implementation; this result is not a claim that the entire current checkout's CUDA/provider path was replayed.

The batch-only candidate accepted step 41 at `t=0.05125`, `h=0.00125`, with fixed mass and Re=∞. Restart, midpoint, endpoint admissibility, density preparation, factor, and endpoint physical-monitor gates passed. Its numerical arrays `ids.bin`, `mass.bin`, `state.bin`, `factors.bin`, `Psi.bin`, and `weight_rate_seed.bin` match the baseline step-41 checkpoint bitwise; `header.bin` differs only in restart run token/epoch. `long_time_stability_verified=false`.

With OMP and the scoped native-face team both at one thread, the two RHS provider times total `377.694801 s` versus `409.322366 s` baseline (−7.73%); stage totals are `273.585771 s` versus `306.133717 s` (−10.63%). Each RHS batched 59,364 faces / 593,640 scalar sample calls; the 32 sampled-face checks were bitwise equal with zero maximum difference. This isolates the per-RHS improvement to the batch preparation. The full accepted-step wall was `748.684 s` versus `641.443 s` baseline, but both jobs shared `syn104`, so this is not a clean end-to-end comparison and no whole-step speedup is claimed.

The benchmark JSON now records final acceptance and accounting. The long-time trajectory remains the governing test: job `428765` is still running, while rank-2 benchmark `428770` and current-checkout replay `428805` remain queued for `QOSGrpGRES`; dependent Kepler job `428781` remains pending. These jobs were left untouched.


### 2026-10-20 01:44 KST — continuing user-source baseline progress

The continuing baseline job `428765` is RUNNING on `syn104` with one H200, 8 CPUs, and 128 GiB. Its last accepted checkpoint is step 55 at `t=0.06875` (`h=0.00125`), recorded at 01:33:50 KST; `step_wall_seconds=685.206 s`. The latest Slurm observation was elapsed `03:03:07`. The two-GPU matched test `428770` and the current-checkout replay `428805` remain pending on `QOSGrpGRES`; dependent job `428781` remains pending. No existing job was cancelled or altered, and no new H200 work was submitted ahead of those queued jobs.


### 2026-10-20 01:55 KST — source-lineage correction, GPU bottleneck, and spin audit

A source audit refined the meaning of “user-source baseline.” Running job `428765` uses binary `d9c79ec0…` from `Codex_full_disk_native_evolution_build_v10_20261019`. It contains the current `eunha2.c` (`b704041c…`) and `full_moment_kepler.c` (`c564f196…`), but its native evolution and stage files are earlier snapshots: `full_moment_evolution.c=dd6e0dde…` and `full_moment_native_stage.c=8038edb3…`, versus current checkout `f9b2e74f…` and `8c939495…`. Therefore `428765` is a useful continuing v10 reference trajectory, not proof of the whole current native source path. Do not interrupt it.

The distinct pending job `428805` is the current-checkout two-rank/two-H200 replay of step 41. Its recorded compile commands compile `Exam/full_moment_evolution.c`, `Exam/full_moment_native_stage.c`, and the other native C/MPI units directly from `/tmp/lageunha-push/Exam`; the benchmark selects explicit midpoint. It still uses the version-pinned v10ref Python geometry/provider and the hash-pinned CUDA operator, so it isolates the newer native C/MPI layer rather than replacing every companion component. `428770` is the separate v10 two-rank comparison. At the 01:55 scheduler observation, both remained `PENDING (QOSGrpGRES)` and `428765` remained RUNNING.

The runtime log shows the current throughput limit is not the H200 kernel. At step 56 (`t=0.07`), each provider RHS took about `227–230 s`; the owned CUDA kernel processed about `91.6 million` gravity nodes in `2.26–2.36 s` per RHS, and the enclosing MPI/GPU stage took about `5.2 s`. The accepted step took `687.869 s`. Thus the GPU kernel accounts for under 1% of that step; CPU/provider geometry, roots, and associated work dominate. The earlier batch-only line-sample optimization is promising, but its one-step `t=0.05` timing does not establish speedup later in the trajectory. Prioritize CPU/provider profiling and then measure the current two-GPU native path before proposing more GPU count.

Added `Docs/audit_kepler_spin_state.py`, a read-only checkpoint diagnostic, and ran it on accepted steps 41, 55, and 56. It decodes the existing centered moments: `S=Q-qq^T/m`, `C=D-qp^T/m`, `I=tr(S)`, `L_int=C_xy-C_yx`, and the rigid-rotation kinetic projection `K_rot=L_int^2/(2I)`. From `t=0.05125` to `0.07`, total `Lz=sum(D_xy-D_yx)` stayed `13.883106850908922`; total energy changed by `1.63e-9` (about `3.4e-10` relative). `K_rot` increased from `4.95129136e-5` to `4.95215787e-5` and is about `10.64%` of the internal kinetic moment. The kinetic partition identity closes to about `1e-20`. This is only an early-time v10 diagnostic; `long_time_stability_verified` remains false.

Rotation is already represented implicitly in the evolving off-diagonal position/velocity raw moments; the helper `moment_spin_ledger.h` explicitly defines `K_rot` as a diagnostic projection, and `msl_rate` has no production call site in `Exam` or `Voro/Laguerre`. No extra spin energy should be added on top of the transported kinetic moments. At Re=∞, do not add an arbitrary damping/relaxation term; use this ledger to measure how the existing model transfers energy among rotation, affine shear, residual motion, and thermal energy. The pending current-checkout replay remains the next authoritative native-source comparison.


### 2026-10-20 02:00 KST — Anderson-1 coupled midpoint jobs terminated without accepted steps

Final Slurm accounting for jobs `428213` (fine `h=0.0025`) and `428214` (coarse `h=0.005`) is `FAILED 114:0`, elapsed `00:47:20` and `00:47:30`. Their native failure manifests show ten provider attempts and `failed_native_transaction_preserved`; fine stopped at step 11, coarse at step 6. Neither branch has an accepted checkpoint in this gate, and the `.025` source checkpoint remains preserved. This supersedes the earlier handover entries that still described these jobs as pending.

The saved last iterations show why merely raising the provider cap is costly: fine Anderson vector norm-squared fell from `204.76` to `0.8395`, then `0.00262`, but oscillated and ended at `2.24e-4`; its final maximum residuals include `W=1.63e-3`, centroid `1.60e-3`, `S=1.21e-3`, thermal `2.38e-3`, and Psi remainder `5.45e-3`. Coarse ended at norm-squared `2.28e-3`, with `W=3.74e-3`, centroid `4.33e-3`, `S=3.33e-3`, thermal `6.38e-3`, and Psi remainder `1.10e-2`. Both remain far above the unchanged componentwise convergence limits; no step was accepted.

This rejects Anderson-1 at the tested h values under a ten-provider budget; it does not reject the physical midpoint equations. The next coupled-method attempt must improve the safeguarded nonlinear solver using the stored residual evidence (for example, a bounded multi-history Anderson/Broyden method) while retaining the same physics, fixed-mass constraints, and acceptance thresholds. Do not treat a larger budget or loosened residual limit as a solution. Until such a method passes, the explicit fine branch remains the only continuing accepted trajectory, but is far too slow for `t>6` without major acceleration.


### 2026-10-20 02:12 KST — interim convex-Anderson draft rejected before simulation

The user's instruction to use their code was checked against the active checkout. `/tmp/lageunha-push/eunha2.c` remained at SHA-256 `b704041c2c6e4754c24de4827f976bde23ff10c4042e2df94bf41bb96db5b556`, matching the input recorded for the running v10 H200 reference. No edits were made to `eunha2.c` or its CUDA kernels.

An interim two-history convex mixer was implemented and locally fixture-tested, but a review against `Docs/next_nonlinear_solver_design_2026-10-19.md` showed that it violated the approved solver coordinates and Type-II design. It was never built into a run binary or submitted to H200. That interim result is superseded by the Type-II implementation below and must not be used as a physical result.


### 2026-10-20 02:27 KST — Type-II solver implemented to match the saved design

The native controller now uses a maximum depth-three Type-II Anderson fit. Rank 0 forms `ΔR`, runs reorthogonalized rank-revealing modified Gram-Schmidt QR, drops the oldest dependent columns, rejects coefficients with magnitude above 10, and broadcasts the selected depth/coefficients/rank/condition. The candidate is `M_k - ΔM γ`; powers-of-two backtracking checks finite state plus centered S/R positive-definiteness and thermal positivity, with a pure-map fallback. This remains a solver-only change: the timestep, physical RHS, Re=∞, fixed-mass condition, actual-provider component gates, and ten-provider step limit are unchanged.

The Anderson residual coordinate uses the accepted timestep's iteration-0 provider metric for every fit in that timestep. Actual-provider fit metrics still determine the existing componentwise acceptance test. The runtime saves each actual `potential_norms.bin` payload, saves the fixed metric for each metric epoch, records a canonical 64-bit metric fingerprint, dependent-cell count, actual `sigma` range/rank margin, raw Ψ-remainder numerator and fixed-metric residual norm. A dependent-mask change starts a new metric epoch and clears history. Two consecutive actual residual norms exceeding four times their predictions also clear the history.

Only `U=state[0:20*n]` and `Psi` are accelerated. The candidate is projected back to the accepted base's common W gauge; IDs/mass remain the base arrays, the four global ledger entries and `weight_rate_seed` come from the newest verified map, and Cholesky factors are rebuilt from candidate moments. The control report stores per-iteration actual residual maxima, metric fingerprints, QR depth/rank/condition, coefficients, chosen λ, predicted norm and fallback/miss reason.

Local verification passed 18 native fake-callback fixtures, including depth-two QR, dependent-column removal, coefficient rejection, Type-II controller integration, positivity backtracking, metric-epoch reset, prediction-miss reset, fixed IDs/mass/global ledger/seed preservation, factor rebuild, and provider/transaction invariants. They passed AddressSanitizer and UndefinedBehaviorSanitizer with leak detection. `Exam/full_moment_evolution.c` compiled as an MPI object with the cluster Intel MPI compiler; no MPI executable was run on the login node. Whitespace scan found no trailing whitespace. The 18 fixtures exercise the controller with synthetic callbacks only; they are not a physics test.

At the 2026-10-20 02:12 KST scheduler observation, the continuing H200 baseline `428765` was RUNNING and current native replay `428805` plus rank-2 comparison `428770` were `PENDING (QOSGrpGRES)`; dependent `428781` remained pending. No job was changed or added. The Type-II implementation has not yet been built into a simulation executable or run on H200, so no physical convergence or accepted step is claimed. Next, build with the unchanged user `eunha2.c` and run fine `h=.0025` and coarse `h=.005` from the preserved `.025` checkpoint when H200 capacity permits. Compare fixed-metric predictions, actual-provider maxima, accepted state and physical monitors; preserve all existing checkpoints and jobs.


### 2026-10-20 02:38 KST — user-source Type-II build prepared; first gate queued for capacity

The current user `eunha2.c` remains SHA-256 `b704041c2c6e4754c24de4827f976bde23ff10c4042e2df94bf41bb96db5b556`. A separate GPFS build was produced at `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_build_v1_20261020/`; its binary SHA-256 is `3b85ebea99e4192209cd4b2b4277c10fa4ed3719f34560a73ce243ee424dd2ae`. It recompiles the staged user main source and links the current Type-II controller/runtime. Static symbol checks found `full_moment_evolution_run`, `full_moment_coupled_control`, and `full_moment_coupled_type2_fit`; `ldd` reports no missing libraries. This is a build check only, not a physics execution.

The existing provider service v1 already implements commands 15–18 (`COUPLED_METADATA`, `VERIFY_FORK`, `MIDPOINT_OBSERVE`, `COUPLED_JOURNAL`). The initial compile attempt exposed that the new runtime was accidentally resolving a stale provider header; the current matching header was staged and the build then completed. The Python provider, native geometry libraries and CUDA operator remain pinned to their recorded hashes. No old build or provider directory was modified.

Prepared, but not submitted, two separate H200 branches under `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_first_gate_v1_20261020/`. Both fork from the sealed `t=0.025` checkpoint (SHA `465f9713…`), use `self_consistent_midpoint` with the Type-II depth-three controller, and target `t=0.03`: fine `h=0.0025` (global step 12) and coarse `h=0.005` (global step 6). Each requests 2 MPI ranks, `--gres=gpu:H200:2`, 192 GiB, and 12 hours. Actual provider metrics and existing component acceptance limits remain authoritative; the runs make no long-time stability claim.

At the 02:38 KST scheduler observation, reference job `428765` was RUNNING at one H200; step 60 at `t=0.075` was accepted and step 61's first RHS was in progress. Current-checkout replay `428805` and rank-2 reference comparison `428770` remained PENDING on `QOSGrpGRES`; dependent job `428781` remained PENDING. The Type-II branches were not submitted while these queued H200 requests were blocking capacity. The governing `t=6–10` Kepler stability objective is still unverified; the current `t=0.075` reference trajectory is far short of it.


### 2026-10-20 02:45 KST — intrinsic-spin diagnostic extended to t=0.075

Read-only `audit_kepler_spin_state.py` analysis of the accepted v10-reference states at `t=0.05125, 0.06875, 0.07, 0.075` found the kinetic partition identity closes to `1e-20`. From `t=0.05125` to `0.075`, total `Lz` is unchanged at `13.883106850908922`; total energy drifts by `2.14e-9` (`4.53e-10` relative). The rigid-rotation projection `K_rot` rises from `4.95129136e-5` to `4.95244436e-5` (0.0233% relative); its share of internal kinetic energy moves from `10.6353%` to `10.6396%`. Maximum absolute intrinsic-spin `omega` stays near `9.11`, and minimum thermal remains positive near `2.02466e-8`. This supports the early-time bookkeeping of existing cell rotation; it does not test long-time spin stability, and the checkpoint remains on the older v10 native stage.

At the 02:44 KST scheduler observation, `428765` was still RUNNING, while `428770` and `428805` remained PENDING on `QOSGrpGRES`. Type-II fine/coarse jobs remain prepared and unsubmitted.


### 2026-10-20 02:47 KST — throughput implication and scheduler projection

At the observed accepted-step wall of about 686 s with `h=0.00125`, a constant-step extrapolation from `t=0.075` to `t=6` is 4,740 steps or about 37.6 days on one H200. This is only a fixed-step arithmetic estimate; it excludes later step rejection/growth and does not predict Type-II throughput. It confirms the explicit baseline cannot be the practical route to the requested long-time gate. The Type-II gate and rank-scaling evidence are needed before selecting a larger-step long-run plan.

The 02:47 `scontrol` projection put rank-2 job `428770` at an estimated start of 14:37 KST; it remains `PENDING (QOSGrpGRES)`. Existing jobs remain untouched. No new H200 job was submitted.


### 2026-10-20 02:51 KST — one-H200 Type-II fine gate submitted

To use an available backfill opportunity without disturbing existing jobs, submitted supplementary single-rank fine Type-II gate `428882` at `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_rank1_fine_v1_20261020/`. It uses the exact staged user `eunha2.c` SHA `b704041c…`, Type-II binary SHA `3b85ebea…`, sealed `t=0.025` checkpoint, `h=0.0025`, and targets `t=0.03`/global step 12. It requests one H200, one MPI rank, 8 CPUs, 128 GiB, and four hours. This is a supplementary solver-convergence run; the prepared two-rank fine/coarse branches remain necessary for MPI validation and later advancement.

Initial scheduler state is `PENDING (Resources)`, so it has not started. `scontrol show node syn104` reports all 8 H200 GPUs allocated (`AllocTRES gpu=8`); existing jobs, including the two-rank pending jobs `428770` and `428805`, were not changed. The single-rank run will not be described as an MPI scaling result or long-time stability evidence.


### 2026-10-20 03:00 KST — stale provider pin caught; corrected user-checkout Type-II package prepared

Job `428882` ended `FAILED`, elapsed one second, exit `1:0`. Its runner rejected the Type-II binary before starting `mpirun`: the submitted provider v1 provenance pinned binary `545b4280…`, while the executable was `3b85ebea…`. The traceback ends at the runner's binary-SHA assertion. No provider service, MPI rank, CUDA operator, or physical step ran; this is a launch-preflight failure, not a solver result.

Created a separate provider snapshot at `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_provider_v1_20261020/`, preserving the old provider directory. It copies the existing 47 hash-verified provider files, adds the three coupled Type-II helper modules plus `curved_disk_cell.py`, and changes the coupled observation helper to import that local pinned geometry module. Its 51 source hashes, binary SHA, build source manifest SHA, and current user source SHA are recorded in `provenance.json`. The build source manifest confirms the executable was compiled from `/tmp/lageunha-push/eunha2.c` SHA `b704041c…`; binary SHA remains `3b85ebea…`.

Prepared corrected two-rank fine/coarse branches in the new directory `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_first_gate_v2_20261020/`. Fine targets global step 12 at `h=.0025`; coarse targets global step 6 at `h=.005`; both fork from the unchanged `.025` checkpoint and request `--gres=gpu:H200:2`. Static preflight passed: all 51 provider files match their hashes, Python source parses, provider wire command IDs match the user's C header, the run parameters/scripts match the preparation manifest, and the build manifest's source/binary hashes match the current user checkout. No physical simulation was run in this preflight.

At 03:00 KST, `428765` and `428592` were RUNNING on `syn104`; `428770` was PENDING `(Resources)`, `428805` PENDING `(Priority)`, and `428781` PENDING on its dependency. The corrected fine/coarse branches remain prepared but unsubmitted so they do not jump ahead of the existing GPU queue. No existing job or checkpoint was changed. The `.03` Type-II acceptance gate has therefore not run, and nothing here establishes stability at `t>6`.


### 2026-10-20 03:02 KST — corrected current-checkout fine/coarse branches submitted

Submitted the corrected branches as Slurm jobs `428884` (fine, global step 12, `h=.0025`) and `428885` (coarse, global step 6, `h=.005`). Both use the current-checkout Type-II binary SHA `3b85ebea…`, the new 51-source provider snapshot, the preserved `.025` checkpoint, two MPI ranks, `--gres=gpu:H200:2`, 192 GiB, and 12 hours. Their run directories and immutable inputs are separate; neither consumes the other's output.

The first scheduler observation after submission reports both `PENDING (Priority)`. They have not started physics. Existing jobs and checkpoints remain untouched. The required next gate is still a fully accepted `t=.03` fine/coarse comparison; this is not evidence for `t>6` stability.


### 2026-10-20 03:13 KST — coupled checkpoint continuation path built and fixture-verified

The original Type-II executable rejected coupled restart manifests and hard-capped coupled targets at `t=.03`; that would have made a successful first gate impossible to extend. Added a separately versioned continuation path in `Exam/full_moment_evolution.c` and `Exam/full_moment_coupled_restart_provenance.py`. It accepts exactly one of an initial explicit-to-coupled fork or a sealed coupled restart, verifies the same `h`, the next staged endpoint, full source/runtime/provider/operator identity, paired fine/coarse gate, accepted control/audit reports, and all checkpoint payload hashes. It serializes the six state payloads into a fresh native run token without advancing time. Coupled restart bootstraps also rerun the independent q48 audit. The next endpoint is `.05` after `.03`, then every `.025` through `t=10`. This interval fits the current 100-provider fine / 50-provider coarse branch budgets at the existing 10-provider-per-step limit. A branch cannot be sealed for continuation without a paired fine/coarse comparison artifact.

Built this path into `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_restart_build_v1_20261020/eunha2.exe`, SHA `272cc1ce…`. Its source manifest confirms the unchanged user `eunha2.c` SHA `b704041c…`; the new runtime provider snapshot is `/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_restart_provider_v1_20261020/` with 52 hash-pinned sources. The build completed with the toolchain warning that `ifx` ignores `-fno-strict-aliasing`; linking succeeded, required native symbols are present, and `ldd` reports no missing libraries.

A synthetic accepted-checkpoint fixture passed seal→verify→serialize: it bound the exact checkpoint to a fine/coarse gate record and hashed comparison artifact/script, preserved all six numerical payloads, exercised the new-token serialization header, rejected a changed `h`, and rejected a target beyond the manifest cap. The native C import path is compile-verified but has not yet run under MPI. The fixture performed zero physics providers and accepted zero time steps; it is software-path evidence only. No `.03` physical checkpoint exists yet, so this continuation binary has not run on GPU.

At 03:13 KST, `428884` (fine) remained `PENDING (Priority)` with estimated start 18:38; `428885` (coarse) remained `PENDING (Priority)` with no estimated start. Existing H200 scale/replay jobs still have earlier estimated starts. The older user-source baseline `428765` accepted step 63 at `t=.07875` in `682.81 s`; step 64 had entered its first RHS at the latest log observation. This slow explicit trajectory is not the current-checkout Type-II branch and remains far from the long-time target.


### 2026-10-20 03:23 KST — runtime-language audit and performance evidence corrected

The current Type-II path is a hybrid runtime, not a Python-free solve. User
`eunha2.c` dispatches into the native C/MPI evolution driver (staged user
source SHA `b704041c…`). Each MPI rank starts the configured Python provider
with `posix_spawnp`; the provider uses a seven-process `ProcessPoolExecutor`
per rank for owned CPU preparation. Its global face/root/audit path is Python.
With `LAGFORCE_NATIVE_FACE_STAGE=1`, it hands prepared arrays to
`native_full_disk_stage`; the state update and owned GPU stage remain under
native MPI/CUDA control. The submitted branches request two ranks, eight CPUs
per rank, and two H200s. Describe this as native MPI/CUDA evolution with a
substantial Python CPU provider, not as Python used only for analysis.

The available instrumented profile at
`/gpfs/kjhan/LagForce/Codex_provider_profile_once_20261018/report.json` is a
different one-provider replay: 1,590 cells, 272,336 nodes, `t=5e-5`, with
concurrent load. It measured 363.04 s provider wall time. Exclusive CPU hot
spots were density face traces 91.49 s, velocity face samples 66.52 s,
moment-flux accumulation 64.90 s, and exact Euler-star pressure 39.05 s.
These numbers prioritize likely optimization work but are not throughput
measurements for the 19,945-cell Type-II branches.

At 03:21 KST, `squeue --start` still showed `428884` and `428885` pending;
the fine branch's estimated start remained 18:38 and the coarse branch had no
estimate. No physical step has run in these branches. Do not submit a second
profile job onto the same H200 queue before the paired gate. Once these exact
branches start, retain their native logs and provider timing reports as the
first full-size Type-II baseline. Target per-face CPU sampling/flux assembly
for later acceleration, and compare any native-batched replacement against
the current provider outputs before using it in a physics run.


### 2026-10-20 03:39 KST — native face speedup and moving batch overlay separated

The saved-state native-provider equivalence artifact at
`/gpfs/kjhan/LagForce/User_full_moment_current_checkout_provider_v1_20261019/native_face_provider_equivalence.json`
records `t=0.025`, 19,945 cells, two physical provider calls, and zero accepted
time steps. Native versus baseline outputs are exactly equal for stage,
incoming flux, offsets, generator, weight rates, Psi rates, ledgers, mass, IDs,
and node input/output hashes. On this one saved state, geometry plus transport
root fell from 165.505 s to 67.305 s (2.459x); whole-provider speedup was
1.558x. This meets the predeclared 2x geometry/root objective for that saved
state, but is not a time-evolution or long-time stability result.

The separate moving replay `428846` at `t=0.050625` compares a batched line
sample overlay against its already-native parent provider. The baseline and
candidate use identical native face-module hashes. Across the two replayed
RHS calls, measured geometry plus mass-root time fell from 108.345 s to
75.946 s (1.427x); total provider time fell 7.73% and stage time 10.63%.
The accepted step at `t=0.05125` preserved numerical arrays bitwise. Its
whole-step wall comparison is not clean because both jobs shared `syn104`.
These numbers evaluate the extra batching overlay, not the original native
geometry/root substitution; they do not replace the 2.459x native-provider
comparison above.

The `.03` fine/coarse Type-II jobs `428884` and `428885` were still pending at
the 03:37 KST scheduler observation. No physical job was changed or added.
`Docs/create_full_moment_coupled_pair_gate.py` now rejects malformed or
nonfinite four-entry P/E/L ledger checks, malformed energy bounds, and any
continuation forecast outside the approved 15-step/150-provider and combined
24-hour budget. Syntax and CLI checks passed; it has not yet processed a
physical `.03` pair because neither branch has run.


### 2026-10-20 03:40 KST — rank-2 comparison entered RUNNING

The 03:39 scheduler observation shows `428770` (rank-2 comparison) entered
RUNNING, while reference `428765` remained RUNNING and current replay `428805`,
Type-II fine/coarse `428884`/`428885`, and dependent `428781` remained pending.
The `.03` branches still have no accepted provider calls. No job was submitted,
canceled, or modified in this observation.


### 2026-10-20 03:42 KST — rank-2 restart reached its first provider call

Job `428770` is running on `syn104` with the requested two H200s, two MPI ranks,
16 CPUs, and 256 GiB. Its restart-import report verifies the unchanged source
checkpoint at `t=0.05`/step 40 and records zero numerical state updates. The
native log confirms both ranks started; rank 0's provider log records the first
call for 19,945 cells. No accepted step-41 or completed-batch record exists yet,
so this is not a rank-scaling result.

The live reference `428765` has now accepted through step 66 at `t=0.0825`;
step 66 wall time is 682.64 s. At that constant cost and `h=0.00125`, the
remaining 4,734 steps to `t=6` project to about 37.4 compute days, excluding
queue time and future cost changes. This confirms the explicit reference is a
monitor/comparator, not a practical long-time path. The `.03` Type-II pair and
current replay remain pending. No job was changed or submitted.


### 2026-10-20 03:45 KST — rank-2 first-RHS timing and MPI serialization point

Rank-2 job `428770` completed its first RHS at `t=0.050625` with
`FULL_MOMENT_RHS=187.016 s`. Both H200s were active and divided the owned stage
as 9,972/9,973 cells; the distributed CUDA stage reported 4.512 s. The matching
one-rank baseline first RHS is 204.848 s, only a 1.095x provider-time ratio
(8.71% reduction). Its face event sums geometry+mass-root to 55.769 s versus
53.775 s for the one-rank baseline. The rank-2 step has not yet produced an
accepted or completed record, and both jobs share `syn104`; do not treat these
first-RHS timings as a final scaling result.

Source-path inspection explains the weak early scaling: bounds/density/velocity
preparation is divided by stable-ID cell ranges across ranks, but `GEOMETRY`
and `FACE_STAGE` requests are issued only by rank 0 and each needs the global
face input. All ranks then participate in the native owned GPU stage. Thus MPI
currently halves the cell-local CUDA work while leaving large global geometry,
face-state, and transport-root preparation on the root rank. Use the completed
rank-2 step to quantify the whole-step effect; any distributed global-face
redesign is a separate exact-parity performance change, not part of the
currently running comparison.


### 2026-10-20 03:54 KST — rank-2 parity gate completed

Slurm job `428770` completed with exit `0:0` in 784 s. It accepted step 41 at
`t=0.05125`, `h=0.00125`, with `step_wall_seconds=619.555`; the completed-batch
report records 780.058 s. Fixed mass, Re=∞, the eight-region physical monitor,
independent q48 density audit, runtime identity, and restart seal all passed.
The build's staged `eunha2.c` has the same SHA-256 as the current user source:
`b704041c2c6e4754c24de4827f976bde23ff10c4042e2df94bf41bb96db5b556`.

Compared with the one-rank accepted step from the same `t=0.05` source,
`ids`, `mass`, `state`, `factors`, `Psi`, and `weight_rate_seed` are all exactly
equal as numerical arrays. `ids`, `mass`, `state`, `factors`, and
`weight_rate_seed` files are byte-identical. `Psi.bin` has different raw file
hashes only because this host stores each 80-bit `long double` in a 16-byte
slot: all ten value bytes match for every cell while six padding bytes vary;
the loaded `Psi` arrays compare exactly. The checkpoint header differs only
in the expected run token.

The two rank-2 RHS calls total 371.560 s versus 409.322 s for the one-rank
baseline (9.23% lower); accepted-step wall is 619.555 s versus 641.443 s
(3.41% lower). Both ran on `syn104` while sharing the host, so this is a weak,
non-isolated scaling observation rather than a clean speedup claim. It verifies
same-step numerical parity across 1/2 ranks, but does not explain or resolve
the long-time Kepler drift.

Summing the two provider event groups separates the serial face/root work from
the remaining RHS path. Geometry + face-state + mass-root total 172.078 s for
rank 2 and 167.835 s for rank 1 (+2.53%); the remaining RHS time is 199.482 s
versus 241.488 s (-17.39%). The two distributed GPU stages total 9.003 s, or
2.42% of rank-2 RHS time. The remainder includes owned CPU preparation,
gather/serialization, and stage transfer, so it must not be labeled CPU kernel
time alone. This confirms that adding GPUs has little leverage until the
rank-0 face path and provider-side CPU costs are reduced.

At 03:54 KST, fine/coarse Type-II jobs `428884`/`428885` and replay `428805`
were still `PENDING (QOSGrpGRES)`; Slurm projected their start to 16:41 KST.
Dependent job `428781` remains pending on its prerequisite. None has started a
physical step. No existing job was modified.


### 2026-10-20 03:55 KST — spin and energy ledger extended through t=0.0825

Read-only `Docs/audit_kepler_spin_state.py` analysis on rank-1 reference
checkpoints at steps 41, 55, 56, 60, and 66 (`t=0.05125` through `0.0825`)
finds total `Lz=13.883106850908922` at every saved point. Total energy changes
from `4.737666648431901` to `4.7376666514205485`, a relative drift of
`6.31e-10`. The rigid-spin projection `K_rot` rises from `4.95129136e-5` to
`4.95291361e-5` (`+0.0328%`), while affine kinetic energy changes by `-0.0231%`;
the rigid-spin share of internal kinetic energy moves from `10.6353%` to
`10.6411%`. Maximum intrinsic-spin `|omega|` stays near `9.10998` and minimum
thermal remains positive (`2.02466e-8`). The kinetic-partition closure residual
stays below `4.7e-20`.

This extends the early spin bookkeeping only. `K_rot` remains a projection of
the existing affine kinetic energy, not a new stored energy; these `t<0.1`
states do not establish the cause of density drift or Kepler stability at
`t>6`.


### 2026-10-20 03:59 KST — reference spin ledger advanced one accepted step

The still-running reference `428765` accepted step 67 at `t=0.08375`
(`668.96 s`). Extending the read-only spin audit to that checkpoint gives
`E=4.737666651569337`, `Lz=13.883106850908922`, and minimum thermal
`2.02466161e-8`. Relative to step 41 (`t=0.05125`), energy change is
`6.62e-10`, `K_rot` is higher by `0.0345%`, and affine kinetic energy is lower
by `0.0240%`; kinetic-partition residual is `1.59e-20`. This remains only an
early-time conservation/rotation bookkeeping check.


### 2026-10-20 04:08 KST — preserve the proven face-sampling batch in checkout

The latest queue observation in this turn found the `.03` Type-II jobs
`428884`/`428885` and current-checkout replay `428805` still pending
`QOSGrpGRES`, with projected start `2026-10-20T16:41:27`; dependent job
`428781` remained pending on its dependency. Since `squeue --start` lists
pending jobs only, a follow-up Slurm accounting check was made: `428765` is
still `RUNNING`, elapsed `5:28:56`. These sealed jobs were not modified. Their
provider snapshot v1 does not contain the later batched line-sampling overlay.

Copied the exact tested `parallel_line_prepare.py` from
`Codex_native_face_batch_v1_20261020` into `Docs/parallel_line_prepare.py` and
connected it from `Docs/native_full_disk_stage.py`. It groups the same five
straight-face velocity samples by owner cell, then leaves the original native
cut/fill, transport root, GPU stage, and physics audits in place. The full
one-step batch-only benchmark `428846` had already accepted the same t=.05→.05125
state with all six numerical checkpoint arrays bitwise identical to baseline;
its two provider RHS calls were 7.73% faster and stage totals 10.63% lower.
Those timings are useful provider evidence, not an isolated whole-step speedup.

`Docs/full_moment_face_provider.py` now resolves a provider-local native stage
before falling back to the pinned shared module directory. This lets a future
hash-pinned provider bundle use the checked-in batch stage without changing
the active immutable job snapshots. Python compilation and an import smoke
check verified that the local stage loads and its `prepare` callback is the
batched implementation. The new helper SHA matches the tested module exactly.
No H200 job, `eunha2.c`, CUDA kernel, physical operator, or accepted state was
changed in this source integration. The updated checkout source has not yet
been rebuilt or run on H200; the existing Type-II gate remains the next
physics result, and t>6 Kepler stability remains unverified.


### 2026-10-20 04:12 KST — immutable batched Type-II provider v2 staged

Created a separate provider snapshot at
`/gpfs/kjhan/LagForce/User_full_moment_current_checkout_type2_provider_v2_batch_20261020/`.
It copies the 51-source Type-II v1 snapshot, overlays the checked-out
`full_moment_face_provider.py` and `native_full_disk_stage.py`, and adds the
previously tested `parallel_line_prepare.py`. The new manifest verifies all
52 source hashes, records the v1 provenance SHA
`8d2017e2c1be3e1b3642ebfa461863d8f307a163ff4e6aadda5b8a7229c3edf5`, and keeps
the Type-II binary unchanged at SHA
`3b85ebea99e4192209cd4b2b4277c10fa4ed3719f34560a73ce243ee424dd2ae`.

Provider-local import smoke verification selected the new snapshot's native
stage and confirmed that its straight-face `prepare` callback is the batched
implementation. No physics provider, accepted step, or Slurm job was run or
submitted for this packaging check. Existing fine/coarse jobs `428884`/`428885`
remain on their immutable v1 provider and are still pending `QOSGrpGRES`;
`428805` is pending for the same reason, `428781` is dependency-pending, and
`428765` remains RUNNING on `syn104` at elapsed `5:30:59`. Use provider v2 only
for a subsequent run after evaluating the current `.03` pair; it has no H200
performance result yet and it proves no long-term Kepler stability.

The 04:13 KST `squeue` recheck confirms no queue transition: `428884`,
`428885`, and `428805` remain pending `QOSGrpGRES`; `428781` is still pending
on its dependency, and `428765` remains RUNNING on `syn104` at elapsed
`5:32:14`.

The 04:16 KST queue check again found `428884`, `428885`, and `428805` pending
`QOSGrpGRES`, with the scheduler's earlier 16:41:27 KST start estimate still
in effect; `428781` remains dependency-pending. `428765` is RUNNING on `syn104`
at elapsed `5:34:41`. The `.03` physics branches have not started.

At 04:17 KST, `squeue` again showed `428884`/`428885`/`428805` pending
`QOSGrpGRES` and `428781` dependency-pending; `428765` remained RUNNING on
`syn104` at elapsed `5:34:41`. The fine/coarse runtime directories still
contain only their sealed `params.dat` and `submit.sh`; no provider call has
started in the `.03` pair.


### 2026-10-20 04:22 KST — verify user-source lineage and correct execution-path description

The current checkout is `/tmp/lageunha-push`, HEAD
`fbc5dfff6d167f108f76d329903f6498af9cacc3`. Its `eunha2.c` SHA-256 is
`b704041c2c6e4754c24de4827f976bde23ff10c4042e2df94bf41bb96db5b556`, matching
the current Type-II build manifest's `user_main_sha256`; the resulting binary
SHA-256 is `3b85ebea99e4192209cd4b2b4277c10fa4ed3719f34560a73ce243ee424dd2ae`.
The sealed fine/coarse scripts invoke this binary and pin its hash. Their
parameters require fixed mass and `Re=∞`.

The source trace narrows what “using the user's GPU code” means here:
`eunha2.c` dispatches `FullMomentKepler` configurations to
`full_moment_evolution_run`; the evolution path calls
`full_moment_native_stage`, which calls `lagmoment_owned_stage` from the pinned
`liblagmoment_owned.so`. The face/provider service is Python. Thus these runs
use the user's exact main source and an MPI/CUDA FullMoment extension, but they
do not run through the pre-existing `Exam/exam_gpu.cu` solver path. Earlier
wording that implied the latter was too broad. The prepared batched provider
v2 is also not used by these already submitted jobs; their scripts pin
provider v1.

At 04:22 KST, `428884`/`428885`/`428805` remained pending `QOSGrpGRES`,
`428781` remained dependency-pending, and reference `428765` remained RUNNING
on `syn104` (elapsed `5:40:39`). No queued job was modified and no new
simulation was submitted during this provenance audit.


### 2026-10-20 04:25 KST — distinguish active branch from the open science tree

A separate read-only comparison found that the open workspace's source repo is
`/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/src`, branch `master` at
`75aaeedc52a2183077ad45f73f4a1c1471b11a6d`; its `eunha2.c` SHA-256 is
`2a80b30d5c71c0908c0075865451f309ec3d35170ee9f71cf69e0a3bf09ccb57`.
The H200 jobs listed above were built instead from the separate Git worktree
`/tmp/lageunha-push`, branch `codex/physics-validation-20261017` at `fbc5dff`,
with `eunha2.c` SHA-256 `b704041c...`. Both report the same GitHub remote, but
they are distinct code states and the full-moment branch is not present in the
open workspace's `eunha2.c`. The interactive source tree's local `master` also
lags its GitHub `master` ref observed during this audit (`8e25dfe`, 2026-09-29).

Do not describe the submitted jobs as built from the open workspace's exact
source. They remain valuable branch-specific FullMoment evidence; any result
must retain that provenance. The earlier 04:22 section's phrase “current
checkout” refers specifically to `/tmp/lageunha-push`, not the open workspace.
