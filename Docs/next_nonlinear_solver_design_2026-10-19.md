# 다음 비선형 solver 설계 — 2026-10-19

이 메모는 첫 self-consistent midpoint gate에서 fine/coarse 양쪽이 10회 provider 제한에 도달해 실패한 뒤, 물리식과 acceptance criteria를 바꾸지 않고 solver만 개선하기 위한 구현 요건이다. 코드는 driver가 작성한다.

## 확인된 결과

- fine `h=0.0025` (`428213`)와 coarse `h=0.005` (`428214`)는 같은 `t=0.025` checkpoint, 같은 19,945셀 물리상태, 같은 H200 2-rank 실행 바이너리를 사용했다.
- 둘 다 provider 10회 후 `FM_COUPLED_NOT_CONVERGED` (`control.json` status 4, Slurm exit 114)로 종료했다. `converged=false`, `final_checks_passed=false`, `accepted=false`; accepted time step은 0이다.
- fine의 마지막 norm²는 `2.23747e-4`; W/centroid/S/thermal/Psi 최대 map 잔차는 `1.63e-3/1.60e-3/1.21e-3/2.38e-3/5.45e-3`다.
- coarse의 마지막 norm²는 `2.28283e-3`; 같은 잔차 성분은 `3.74e-3/4.33e-3/3.33e-3/6.38e-3/1.10e-2`다. coarse norm²는 fine의 `10.20`배이며 각 주요 성분도 약 `2–3.2`배 크다.
- 두 경로 모두 mass root, fixed mass, positive thermal, P/E/L ledger, native factor audits 및 실제 H200 GPU/MPI stage를 통과했다. 문제는 provider 또는 mass closure의 실패가 아니라 outer midpoint fixed-point residual의 미수렴이다.
- 입력 checkpoint 파일은 현재도 SHA-256 `465f9713ab405f19e438396abb8796510c3c51784b0c85028139248b34c28d8e`다. 양쪽 failure manifest도 이전 accepted 상태 보존, retry 없음, timestep 변경 없음을 기록한다.

## 왜 다음에는 더 긴 단순 반복이 아니라 여러 이력을 써야 하는가

현재 `Exam/full_moment_coupled_control.c`는 직전 map 하나만 보유하는 convex Anderson-1이다. `FullMomentCoupledResidual`은 current/previous residual의 norm²와 dot product만 저장한다. residual은 큰 초기 감소 후 여러 번 반등했고, coarse에서는 일부 반복의 혼합계수가 0으로 걸러져 pure Picard map이 사용됐다. 따라서 단순히 provider 횟수만 늘리면 같은 진동/plateau를 더 비싸게 따라갈 수 있다.

실제 dimensionless residual norm²의 전체 이력은 다음과 같다.

| provider | fine `h=.0025` | coarse `h=.005` |
|---:|---:|---:|
| 1 | `204.764` | `826.251` |
| 2 | `0.839507` | `14.0409` |
| 3 | `0.00261864` | `0.0707845` |
| 4 | `0.00313948` | `0.107931` |
| 5 | `0.000462418` | `0.0263809` |
| 6 | `0.000973338` | `0.0180627` |
| 7 | `0.000141502` | `0.00662139` |
| 8 | `0.000110361` | `0.00219973` |
| 9 | `0.000297948` | `0.00695746` |
| 10 | `0.000223747` | `0.00228283` |

초기 transient는 크게 감소하지만, 이후 fine/coarse 모두 반복마다 반등하고 terminal residual norm은 수락 기준을 만족하지 않는다. 특히 norm² 감소만으로는 충분하지 않다. 각 componentwise 최대 residual을 함께 봐야 하며, 두 branch 모두 W/S/thermal/Psi 성분이 끝까지 `1e-3–1e-2` 수준에 남는다.

## 구현 전에 고쳐야 할 좌표계 문제: Ψ 잔차 metric이 반복마다 바뀐다

소스 `Exam/full_moment_coupled_state.c`를 다시 추적하면 21성분 벡터 중 Ψ remainder 성분은 현재 provider의 `FullMomentPotentialNorm`을 사용한다. 구체적으로 매 반복의 fit에서 얻은 `origin`, `scale`, `projection`, `sigma`, `dependent`로 모멘트 변화의 다항식 성분을 빼고 `mass*sigma`로 나눈다. 이 값들은 midpoint geometry와 density fit에 따라 바뀔 수 있다. 그런데 runtime은 iteration `k`에서 만든 벡터와 iteration `k-1`에서 다른 provider metric으로 만든 벡터의 dot product를 Anderson 계수에 사용한다. 좌표변환이 실제로 달랐다면 저장된 norm²·dot·mixing coefficient는 하나의 공통 최소제곱 문제의 Gram 자료가 아니다. 이전의 비단조 norm² 이력을 같은 유클리드 목적함수의 수렴 곡선으로 해석할 수도 없다. 현재 실행에서 metric 변화가 얼마나 컸는지는 payload가 저장되지 않아 정량 판정할 수 없다.

산출물도 이 변환을 사후 재구성하기에 부족하다. fine `428213`의 `runtime/coupled_0011/base_metadata`에는 `areas.bin`, `gamma.bin`, `ids.bin`만 있고, 반복 폴더 `coupled_midpoint_0011_00..09`에는 `admissibility.json`만 있다. iteration JSON에는 잔차 최대값과 Gram scalar만 있으며 `potential_norms.bin`이나 그 hash는 없다. provider가 만든 metadata는 `/dev/shm` 아래 service root에 생성되고 `FMP_RELEASE` 때 지워진다. 따라서 현재 보존된 norm² 변화만으로 이 metric 문제를 사후 해소할 수 없다.

Type-II multi-history를 구현할 때는 두 용도를 분리한다.

1. **실제 수렴/수락 판정:** 기존대로 매 반복의 실제 provider fit metric과 dependent-Ψ 제약을 사용한다. 현재 componentwise gate와 물리 acceptance 조건은 변경하지 않는다. 반복 기록에는 실제 fit의 dependent cell 수, `sigma` 범위, raw Ψ-remainder 분자 최대값 및 provider metadata hash를 함께 남긴다.
2. **Anderson 최소제곱 좌표:** 해당 timestep의 iteration 0 provider, 즉 accepted base 상태에서 만든 per-cell `origin/scale/projection/sigma/dependent`를 고정해 모든 `r_k`에 동일하게 적용한다. Type-II의 모든 `ΔR`·dot product·예측 norm은 이 `r_AA`만 사용한다. 좌표 변환 payload를 timestep 산출물에 보존하고 각 iteration에서 같은 metric fingerprint를 기록한다.

현재 provider의 dependent mask가 고정 mask와 달라지거나 고정 변환이 유효하지 않으면 서로 다른 좌표계의 이력을 섞지 않는다. 이 경우 기존 이력을 모두 비우고 새 metric epoch에서 다시 시작하며, 전환 횟수와 영향 cell을 기록한다. 반복적으로 mask/rank가 바뀌면 Anderson을 강행하지 말고 그 step을 실패시켜 provider rank 전이를 따로 조사한다. 물리 수렴 판정은 언제나 actual-provider metric으로 다시 계산하고, `r_AA`의 감소만으로 step을 수락하지 않는다.

Type-II candidate의 상태 규약도 명시해야 한다. residual vector는 per-cell `state[0:20*n]`와 `Psi[0:n]`에 대응해 길이 `21*n`이지만, `FullMomentState`에는 그 밖에 전역 장부 `state[20*n:20*n+4]`, `weight_rate_seed`, derived factors, 고정 IDs/mass가 있다. 이 부가 배열을 전부 같은 Anderson 식으로 외삽하면 residual least-squares가 제약하지 않은 양을 바꾸게 된다. 가속 대상은 오직 `U=(state[0:20*n],Psi[0:n])`로 한정한다. `W`는 잔차에서 공통 gauge를 제거하므로 외삽 candidate도 accepted base와 같은 gauge로 투영한다. IDs/mass는 base의 bitwise 동일 배열을 유지하고, factors는 candidate의 raw moments에서 다시 재구성한다. 전역 장부는 현재 map 값 그대로 보존하며 최종 ledger gate는 항상 검증된 새 map에 적용한다.

`weight_rate_seed` 사용처를 coupled native path에서 검색한 결과, midpoint는 trial seed를 보존하고 현재 convex mixer는 이를 이전 map seed와 섞지만, native factor는 매 map의 실제 `weight_rates`로 output seed를 덮어쓴다. audited RHS/provider 경로는 trial seed를 읽지 않는다. 따라서 현재 구현에서는 이 값이 fixed-point 물리 unknown이 아니다. Type-II candidate에는 seed를 외삽하지 말고 최신 map의 seed를 그대로 전달한다. 추후 seed가 RHS에 실제 영향을 주도록 바뀌면 그때 고정 metric residual과 수렴/감사 기준을 추가해야 한다.

기존 Anderson-1 기록(`428213/428214`)의 `anderson_vector_norm2`, `anderson_dot_previous`는 이 감사가 반영되기 전의 값이다. 그러므로 Type-II 설계의 효과를 주장하는 기준선으로 재사용하지 말고, metric 변화가 미계측이라는 한계를 붙여 solver 진동의 정성적 단서로만 취급한다. 새 비교에서는 fixed-metric norm²와 actual-provider component maxima를 나란히 기록한다.

## 저장 residual의 공간 위치 감사

Fine/coarse 각 반복 JSON의 성분별 최대 cell index를 accepted `t=.025` checkpoint의 ID·정확 곡면 support와 대조했다. ID 배열은 mesh 및 initial checkpoint와 일치한다. iteration 2–9에서 유효한 성분별 최대 위치 144개 중 127개가 초기 반경 `r<.4`인 inner disk가 아니라 `r=.5` 원의 안쪽에 놓인 cavity cell이다. 나머지는 10개 inner-contact cell과 7개 inner-disk cell이며 disk-core 최대는 없다. index/ID 720은 두 branch 합계 최대 위치로 75회, ID 747은 35회 나타난다. 두 cell 모두 초기 site 반경 `.46875`, support가 각각 `[.449593,.485216]`, `[.452422,.485158]`이어서 `r=.5` 경계를 가로지르지 않는 접촉부 바로 안쪽 cell이다. 따라서 현재 기록의 최대 map residual은 interface crossing cell보다 contact 안쪽 이웃에 더 자주 집중된다.

이 위치가 시작 시 potential-density fit rank 실패 때문이라는 증거는 없다. `.025` accepted 상태의 ID 720/747 fit은 각각 dual residual `2.78e-16/2.91e-16`, condition `319/293`, `sigma/rank_error` 약 `6.5e9/1.2e10`으로 통과했다. 비교용 inner-contact ID 795는 condition `1414`, `sigma/rank_error≈646`으로 통과했다. 다만 이 값은 시작 fit만 설명하며 각 nonlinear midpoint provider의 fit 변화는 저장되지 않았다. 새 구현은 cavity-side 최댓값 cell들의 반복별 actual-fit `sigma`, rank margin, fit condition 및 component residual을 기록해 solver-coordinate 문제와 moving-fit conditioning을 분리해야 한다. 공간 상관만으로 물리 RHS의 원인이나 t=6–10 붕괴를 설명하지 않는다.

runtime callback `Exam/full_moment_coupled_runtime.inc`는 이미 물리 스케일에 맞춘 길이 `21*n`의 signed residual vector를 계산하지만, 현재는 rank 0에서 직전 벡터 하나만 저장한다. 19,945셀에서 long-double residual vector 하나는 약 6.7 MB이므로 3–4개 이력을 보관하는 비용은 이 문제 규모에서 감당할 수 있다.

## 제안: positivity safeguard가 있는 Type-II multi-history Anderson

현재 trial `x_k`, 실제 map 출력 `m_k=M_h(x_k)`, 실제 map residual `r_k=m_k-x_k`를 쓴다. 최근 residual difference 열 `Δr_j=r_{k-j+1}-r_{k-j}`에 대해 rank-revealing QR로 다음 least-squares 문제를 푼다.

\[
\gamma_k=\arg\min_\gamma\|r_k-\Delta R_k\gamma\|_2^2,
\qquad x^{AA}_{k+1}=m_k-\Delta M_k\gamma_k,
\]

여기서 `ΔM_j=m_{k-j+1}-m_{k-j}`다. 이는 계수 합이 1인 Anderson Type-II affine extrapolation이다. 현재의 비음수 convex 계수 제한을 그대로 늘리기보다, 작은 다중 이력이 주는 외삽력을 사용하고 실제 유효 상태 조건으로 제한한다.

외삽 결과가 양의 covariance/열 조건을 벗어나면 검증된 map 출력 `m_k`에서 시작하는 선분 감쇠를 사용한다:

\[
x_{k+1}(\lambda)=m_k+\lambda(x^{AA}_{k+1}-m_k),\qquad \lambda\in(0,1].
\]

`full_moment_centered_decode`가 모든 셀에서 S/R 양의 정부호와 양의 열을 확인하고, 상태·Psi가 finite인 가장 큰 `λ`를 선택한다. powers-of-two backtracking에서 유효한 `λ`가 없으면 `λ=0`인 pure map을 쓴다. `m_k`가 물리적으로 유효하고 두 map 상태가 같은 외부 mass/ID 배열을 쓰므로 이 감쇠로 fixed mass를 바꾸지 않는다. 다음 반복은 언제나 새 실제 midpoint provider로 평가한다.

Anderson 예측 norm은 candidate 선택에만 쓰며 acceptance를 대체하지 않는다. 수락은 기존 componentwise 실제 map residual, factor, 실제 midpoint agreement 및 보존 장부 조건이 모두 통과할 때만 한다.

### 구현 조건

1. 직전 이력 하나 대신 길이 `21*n` residual history와 대응하는 physical `U_k`/`M(U_k)` ring buffer를 보관한다. `U`만 외삽하고 IDs/mass, 4개 global ledger, derived factors 및 `weight_rate_seed`는 위 규약대로 처리한다. history는 **현재 timestep 안에서만** 유지하고, 다음 timestep에서는 비운다.
2. rank 0에서 `ΔR` 열을 만들고 고정밀도 rank-revealing QR로 least-squares 계수를 구한다. 거의 종속인 가장 오래된 열부터 제거해 depth를 낮춘다. 계수와 선택된 depth만 MPI ranks에 broadcast한다.
3. 먼저 history depth `m≤3`으로 구현한다. `m=4`는 QR rank, coefficient size, 실제 residual trend를 보고 선택한다. history vector 하나는 약 6.7 MB이므로 각 rank의 state/map ring buffer 메모리와 함께 측정해도 충분한 headroom이 있어야 한다.
4. 계수가 nonfinite/과대이거나 QR rank가 불충분하면 depth를 낮춘다. extrapolated state가 decode를 통과하지 못하면 `λ=1,1/2,1/4,...` 순으로 줄이고, 모든 backtrack이 실패하면 `λ=0` pure map으로 돌아간다. 실제 다음 residual이 반복적으로 예측보다 커지면 history를 초기화한다.
5. 출력 iteration record에 history depth, QR rank/condition, 전체 mixing 계수, 선택된 `λ`, fixed-metric predicted norm², 다음 iteration의 fixed-metric norm², actual-provider component maxima, fallback 사유, metric epoch/fingerprint를 저장한다. 실제 Ψ `sigma` 요약과 dependent-mask 변경도 기록한다. 실패한 timestep에는 재현에 필요한 실제 provider metric payload를 보존한다.
6. 기존 positivity/SPD decode, exact fixed-mass bitwise check, factor audit, actual-midpoint agreement, P/E/L checks 및 10-provider limit을 그대로 유지한다. 실패 시 시작 accepted checkpoint만 보존하며 retry/dt 변경/legacy fallback을 금지한다.

## 검증 순서

새 구현은 보존된 동일 `t=.025` checkpoint에서 fine과 coarse 두 branch로 다시 검증한다. 두 branch 모두 실제 `.03` endpoint까지 accepted step을 완료하고 모든 solver·물리 gate를 통과해야 다음 `.05` 비교로 이동한다. 둘 중 하나라도 10 provider 안에 수렴하지 않으면 그 결과로 다음 solver 변경을 판단하며 수렴 기준은 완화하지 않는다. `.03` 통과도 t6–10 안정성의 증거는 아니다.

## Bulk 물리 평가 반영 — spin·밀도·압력과 다음 관문

2026-10-19의 bulk Astra 평가는 접촉 영역에서 저차 공간 추세를 제거한 뒤에도 spin error와 density error의 부분상관이 남는다고 보고했다: `.025` coarse/fine `0.438/0.640`, `.05` coarse/fine `0.410/0.600`. fine `.05`에서는 상위 spin·density error 11셀 중 9셀이 겹쳤지만 coarse에서는 겹치는 셀이 없었다. 따라서 이는 접촉부의 국소 연관 신호이지, 전체 영역에서 보편적인 관계나 인과 증거는 아니다. 셀 간 공간 의존성과 고차 공간 구조도 제거되지 않았다.

추적 가능한 원자료는 `Docs/canonical_saved_rotation_error_conditioned_v2_2026-10-19.json`과 `Docs/canonical_saved_rotation_error_cospatiality_2026-10-19.json`이다. 압력·밀도 예산과 tangent 대조는 각각 `Docs/canonical_full_disk_completed_density_temporal_budget_2026-10-19.json`, `Docs/canonical_full_disk_saved_t0p025_tangent_growth_2026-10-19.json`에 있다.

같은 평가에서 `.05` inner-contact pressure 상대 RMS는 coarse `7.571e-3`, fine `1.376e-4`로 약 55배 차이가 났다. 반면 intrinsic-spin **오차** RMS는 `5.825e-9` 대 `1.124e-9`(약 5.18배), density RMS는 `7.040e-3` 대 `4.847e-3`(약 1.45배) 차이다. fine의 specific thermal과 density 오차는 각각 약 `4.97e-3`, `4.85e-3`인데 pressure 오차는 훨씬 작다. `e_P=e_rho+e_u+e_rho*e_u` 관계의 상쇄가 pressure를 작게 만들 수 있으므로 압력만으로 상태의 정확성을 판정하지 않는다. 회전항이 압력 개선을 일으켰다고도 결론내리지 않는다.

별도 저장 tangent 분석은 coarse inner-contact `.025→.05` pressure error 변화를 RMS `7.3877e-3`로 측정했고, `.025` 순간 tangent의 total prediction은 RMS `1.8065e-4`, signed mean은 관측 `-7.3707e-3`과 반대인 `+1.5320e-4`였다. 따라서 그 순간 tangent는 뒤의 접촉 압력 변화를 예측하지 못한다. 이는 순간 volume/thermal 계측을 유한시간 원인으로 확대해석하지 말아야 함을 뒷받침한다.

Astra가 권한 다음 물리 비교는 동일 `.025` accepted checkpoint에서 self-consistent midpoint fine/coarse 경로를 기존 explicit 경로와 비교하고, `.03` 및 이후 `.05`까지 ID별 density·specific thermal·pressure·spin 오차의 증가분과 실제 체적 변화 계측을 나란히 보는 것이다. 그러나 Astra 응답의 Slurm 현황은 최신이 아니었다. 현재 `sacct`상 `428213` fine과 `428214` coarse는 모두 `FAILED`, exit `114`, 각각 elapsed `47:20/47:30`이며 두 작업 모두 accepted step 0이다. 따라서 해당 비교는 아직 수행되지 않았다.

### 다음 실행 관문

1. driver는 위 Type-II 설계(고정된 timestep 내부 Ψ metric, rank-revealing QR, positivity safeguard 및 실제 provider 기준 수락)를 구현한다. 기존 failed outputs와 `.025` 시작 checkpoint를 보존한다. 여기서 수렴 임계값·물리식·`h`를 바꾸지 않는다.
2. 같은 보존 `.025` 상태에서 fine `h=.0025`, coarse `h=.005` 두 branch를 H200 native MPI/CUDA 경로로 각각 실행한다. Re=`∞`, 무조향, 고정 mass/ID bitwise 불변을 유지하고 최대 10 provider 안에 실제 `.03` endpoint step이 accepted되어야 한다.
3. 성공한 branch는 기존 explicit 결과와 같은 ID·영역 기준으로 density, specific thermal, pressure, intrinsic-spin 오차 증가분과 실제 체적 변화 journal을 비교한다. pressure만 개선되고 density/thermal 또는 spin이 악화되면 성공으로 판정하지 않는다. fine/coarse 간 차이와 모든 기존 conservation·q48·실제-cell gate를 함께 통과해야 `.05` 비교로 간다.
4. `.05`의 temporal density·pressure·velocity·phase와 회전/열 budget가 모두 기존 기준을 통과하기 전에는 `.1`, `1`, 또는 역사적 `t=6–10` 연장을 시작하지 않는다. Type-II solver가 통과해도 최종 장시간 목표 달성으로 간주하지 않는다.

이 단계는 회전 모델을 새로 바꾸는 실험이 아니라, 이미 포함된 회전/모멘트 상태에서 시간결합 방식이 밀도·열·spin drift를 줄이는지 판별하는 비교다. 여기서 분리 가능한 개선이 확인된 뒤에야 응력에 의한 내부 회전 완화나 가변 weight 법칙을 독립 변수로 추가한다.
