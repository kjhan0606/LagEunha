# Kepler t=8 진단·수정·Astra 검토 — 2026-10-20

## 프로젝트와 목표

사용자 CUDA+MPI 코드: `/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`, origin `kjhan0606/LagEunha`. 새 Python 유체 solver로 대체하지 않는다. 이 디렉터리의 Python은 native 실행·스냅샷·감사 결과 분석 도구다.

목표는 과거 t6–10 실패를 넘어 GIZMO 논문의 최장 Kepler 진화를 검증하는 것이다. Fig7 t600은 그림 비교 관문이고, 논문 본문은 약1000회전까지의 적분도 언급한다. 각 회전수의 기준 반경과 연화 중력을 명시해야 한다. 이번 원반은 Rin2, Rout8, GM1, eps0.047이므로 내측 주기17.77889이며 t8은0.45회전에 불과하다. 현재 결과로 장시간 성공을 선언할 수 없다.

조건: 생성점 조향 없음, Kepler의 Re∞, W=0일 때 Voronoi+개별 셀 고정질량. 동적 가중치·회전은 공통 힘/토크/에너지/질량 계약을 갖출 때 연구 가능하다.

## t8에서 입증된 것

근거: `baseline_t8_true_area.json`, `audit_snapshot.py`, `snapshot_offsets.json`.

스냅샷 `.../Kepler_AdaptiveStencil_CleanLong_H200_Halo8_InitBounds_20261020/production/khout.013378.dat`: t8.00013686,65536셀, native stride576.

독립적인 전역 Voronoi 재구성에서 주기 영상의 영향을 받지 않음을 확인한64240셀의 저장 부피는 float 정밀도(약6e-8)로 일치한다. 저장 rho=m/V도 같은 수준에서 일치한다. 이 시각의 그림/셀 부피 오류를 중심 가열의 설명으로 삼을 근거는 없다. 과거 RK stage의 부피까지 검증한 것은 아니다.

| 반경 | 평균 접선속도 상대오차 | rms(vr/vc) |
|---|---:|---:|
| R<0.5 |0.54479|0.30922|
| 0.5–1 |0.28181|0.21926|
| 1–2 |0.05076|0.09335|
| 2–4 |1.8621e-5|8.5778e-5|
| 4–8 |3.2159e-7|2.7487e-6|

원반은 이 시각 정확하지만 중심 배경은 이미 크게 손상됐다. 최대 P=0.4150934는 초기1e-6의 약42만배, 위치 R0.027957. 전체 U=0.00531127(초기0.000864), R<2의 U≈0.00436414. 작은 전역 총에너지 오차~9.3e-5는 이 국소 궤도에너지→열 전환을 배제하지 않는다.

공간 간격 h0.09375에 eps0.047(0.50h). 중심 최소 주기0.064로 원반보다 훨씬 빠르다. 이것은 강한 해상도 의혹이며, 원인으로 확정한 것은 아니다. 같은 점질량 힘의 중력-only 대조가 필요하다.

GPU stage의 volume floor=1e-3*(4h)^2=1.40625e-4. 최대 P셀의 실제 V=7.64473e-5로 floor의0.544배. t8 스냅샷에서 이 floor 아래인 R<0.5 셀은 하나이고0.5<R<8에는 없다. stage/final 부피 불일치는 검사해야 하지만 전체 가열을 이 하나로 귀속할 수 없다.

`baseline_event_history.json`: pair-pressure 보정 최초 t0.395491, 음의 U floor 최초 t0.769161, 전체 U가 초기 두 배를 넘은 최초 t2.17403. 가열은 t8에서 처음 시작된 현상이 아니다.

## 현재 코드에서 수정한 확정 결함

`voronoi_tessellate_kernel`의 이웃 density 입력과 자신의 density 출력이 같은 device 배열이었다. 블록 실행 순서에 따라 서로 다른 시간의 density를 읽으며 restrict 조건도 위반했다. 영구 버퍼에 입력 density를 복사한 뒤 같은 CUDA stream의 kernel에서 읽도록 수정했다. 식·질량·조향·점성은 바꾸지 않았다. 업로드된 density의 시각 자체를 고차화한 변경은 아니다.

원본의 같은 t6 스냅샷 한 스텝 두 번:392739개 face-end의 기하는 비트 단위 동일하지만 pstar 최대 차이2.4884e-7, 속도 최대 차이7.8244e-7. 수정 후 MPI1에서65536셀의 저장 상태/기울기와392739개 face-end의 pstar/work/기하가 모두 정확 일치한다. MPI4 반복도 동일한 검사에서 모두 정확 일치한다.

근거: `compare_late_repeat.json`, `compare_fixed_repeat.json`, `compare_fixed_r4_repeat.json`.

MPI1↔MPI4의 처음 차이(pstar 최대2.5855e-6, 속도 성분 최대2.5425e-5)는 R>10 외곽으로 국한됐다(`rank_difference_localization.json`). 원인은 MPI1에서 periodic padding의 새 기울기를 갱신하지 않던 NID>1 조건이었다. 이 조건만 제거한 뒤 MPI1 반복 및 MPI1↔MPI4의65536셀 상태/기울기와392739면 기록이 전역에서 비트 단위로 일치했다(`compare_periodic_repeat.json`, `compare_periodic_r1_r4.json`). 이 t6 한 스텝의 원인 개입 검증은 완료했다. 장시간 MPI 불변성까지 입증한 것은 아니다. `GFS_GPU_SYNC_HALO_STATE` unset은 기본 ON이며 MPI1의 주기 영상에도 필요하다. 기준선 로그에도4랭크 HALO_SYNC가 확인된다.

Astra가 제기한 pressure 동시 쓰기 의혹은 함수 경계까지 대조한 결과 같은 결함으로 확인되지 않았다. A의 out_pressure 인자는 미사용이고 실제 pressure 쓰기는 A 이후 같은 stream의 별도 `pressure_stress_kernel`(B)에 있다. 이 정정은 Astra에게 전달했다.

## 기하/일/물질부피를 구분하는 새 근거

근거: `t8_face_work_audit.json`, `audit_face_work.py`, 수정 binary의 t8 한 스텝 재생 `fixed_t8_audit/`. 이미 가열된 상태의 순간 연산자 검사이며 정상상태의 장기 검증은 아니다.

392809개 face-end에서 일 속도의 기하학적 속도와의 차이가1e-12보다 큰 항은188개이며 모두 R<2에 있다. R2–8에서는 per-cell work/geometry volume-rate 차이/V rms가 약1e-15이다. 따라서 현재 원반의 일반 면이 전부 Riemann 속도로 일한다는 설명은 틀리다. 실제 source의 phase1/use_muscl/extreme/geom_fv 분기를 구분해야 한다.

그렇지만 같은 t8 기하 위에서 **해석적인 원운동 속도장**을 대입해 계산한 기하 Vdot/V는 0이 아니다(실제 진화를 바꾸지 않는 사후 계산): rms R<0.5에서27.70, R2–4에서0.07051, R4–8에서0.0004591. 이 원운동장의 연속 div(v)는 정확히0이다. 이는 정확한 Voronoi 부피가 고정질량 유체의 물질부피와 일치하지 않을 수 있음을 보여준다. 이미 왜곡된 t8 기하의 결과이므로 초기 성장의 인과 증명으로 확대하지 않는다.

앞으로 두 잔차를 따로 측정한다.

1. R_work=Vdot_geom−Vdot_work: 실제 압력 일과 기하 부피의 호환성.
2. R_material=Vdot_geom−V div(v_fluid): 조향 없는 생성점의 Voronoi 셀이 물질 셀인가의 문제.

과거 인수인계의 셀 미끄러짐→부피/압력 요동→가역 불안정 기제는 두 번째와 관련된다. 첫 번째 잔차가0이어도 두 번째와 그 압력 응답은 남을 수 있다. 실제 비균일 P, Riemann 소산, floor가 있는 상태의 열 변화는 전체 Udot+P Vdot 장부로 분해해야 한다.

## 운전자 판단과 Astra의 권고를 반영한 다음 묶음

1. **입력/halo 계약:** density race와 local periodic padding 결함을 수정했고 같은 t6 한 스텝에서 반복 및 MPI1/4 전역 정확 일치를 통과했다. 이후 연산자 변경도 이 관문으로 회귀 검증한다.
2. **실제 CUDA 장부:** 기존 t6/t8 체크포인트의 짧은 재생에서 raw/floored V, R_work, R_material, 압력력/압력 일, U·KE·potential·floor, 각운동량 이동을 같은 셀 ID로 잇는다. 균일 P 극한의 entropy 식을 실제 비균일 상태 전체에 대입해 원인을 단정하지 않는다.
3. **중력/해상도 분리:** 같은 Plummer 힘과 초기 원운동의 중력-only RK 대조로 orbital E/L/위상/반경을 검사한다. 필요한 경우 같은 연속 IC와 같은 eps에서 중심부 초기 해상도 대조 한 쌍을 추가한다. 이는 초기 이산화를 바꾸며 기존 개별 질량이 동일한 실험은 아니다. 실제 셀 적분으로 초기 질량을 정한다. eps를 키워 같은 문제를 해결했다고 주장하지 않는다.
4. **잔차에 맞는 물리 연산자:** 입력 결함과 공간/시간 오차를 분리한 뒤, 물질부피/질량/압력 일 계약 또는 보존적 고차 모멘트 경로를 선택한다. 과거 동적 Laguerre 물질부피 추종의 가중치 유효성 한계를 그대로 반복하지 않는다. 단순 geometric-work 스위치만으로 고정질량 ALE 계약이 닫히지는 않는다.

회전은 누락된 각운동량/운동에너지를 담을 수 있다. 공유 면 작용점 torque, orbital+intrinsic L, E_spin=S²/(2I)의 I-dot 항, 변형과 모멘트 수송을 함께 정식화해야 한다. Re∞에서 임의 회전 감쇠를 넣지 않는다. 비연화 Kepler의 spin Ω/4·strain3Ω/4 및 등방성 moment 에너지 비9를 Plummer 중심에 적용하면 안 된다. 일반적으로 spin=Ω+rΩ'/2, strain=|rΩ'|/2이고 Plummer 중심은 강체 회전에 가까워 spin의 상대적 역할이 커진다. 실제 구현의 필요성은 앞 장부가 판별한다. 기존 full-moment 연구 자료는 활용하되 현재 CUDA 본체의 장시간 검증과 혼동하지 않는다.

## 가속 수정의 검증 범위

512개 후보 배열을 없애고 같은 linked-cell 순서로 즉시 clipping한다. t6의392739개 face-end 기하가 원본과 비트 단위 일치했다. thread stack20384→1952바이트. 공유 GPU 조건의 전체 실행241s→259s, median tess20.38→20.52ms였으므로 실제 속도 향상은 입증되지 않았다. 독립 GPU 측정 전에는 속도 향상 수치를 주장하지 않는다.

native RK4 끝의 CPU volume 갱신을 완전 이웃 탐색의 GPU raw geometry 경로로 이식했다. `GFS_GPU_FINAL_VOLUME=1`의 opt-in이며 기본은 아직 OFF다. 아래 정확도·성능 관문을 통과했다. force-stage의 face merging/volume floor는 별도 미해결 항목이다.

## 후속 판별 실험 결과와 해결 방향

### 중력의 시간적분만으로 중심 붕괴를 설명하기 어렵다

`gravity_replay/`: 현재 native kepler_accel/kepler_phi 함수 본문을 그대로 추출하고 hash를 기록한 C RK4 대조. 생산 로그의15345개 timestep(인쇄 정밀도)을 t8.99965까지 적용했다. 원래 IC의 중심88개 격자점에서 최대 상대 orbital E오차1.578e-7, L오차1.582e-7, 반경오차2.522e-7, 최종 위상오차5.10e-5rad였다. dt/2에서는 E/L오차가 약32배 작아졌다. 전체 유체 RK4를 그대로 replay한 것은 아니므로 cell-average 중력/source 결합까지 배제하지 않는다. 그러나 점 중력 시간적분 단독은 coupled 유체의54% 궤도오차 설명으로 매우 약하다.

### 초기 첫 공간 연산자부터 잘못된 압력·열이 생긴다

`initial_face_work_audit.json`: 진화를 시작하기 전 t0 첫 힘 평가에서 모든262144면의 work와 geometric 속도는 일치한다. R<.5에서는 pstar/P0 최대392.87, 양의 추가 열률2.351e-4, 초기 Ucore1.160e-6이다. 즉 재시작한 손상 상태만의 효과가 아니고, 첫 공간 flux부터 문제가 존재한다. 이 순간 기하 Vdot/V rms1.440은 연속 원운동의 정확한 div0과 다르다. work/geometry 불일치만으로 가열을 설명할 수 없다.

기존 `gradient_method=1` LSF로 바꾸면 R2–4의 양의 추가 열률은1.779e-5→1.228e-8(약1449배 감소), pstar/P 최대12.60→1.036이다. R4–8 열률은약3912배 감소했다. 반면 중심 pstar/P는316.43, 추가 열률1.797e-4로 여전히 크다. LSF는 필요성이 입증된 개선이나 전체 해결이 아니다.

### 기존 Green–Gauss 설정은 선형 속도장을 재현하지 못한다

`affine_probe/results.json`: 균일 rho/P, 중력 OFF, v=(.1x,-.1y)의 native 1스텝 patch test에서 GG+OoA1/3의 기울기는 CPU와 GPU 모두 정확값보다22.222% 작았다. GG+OoA0는상대오차약5e-8, LSF는5.2e-14로 통과했다. 정사각형4이웃에서 기존 보정항 kp+km−i−j가 상쇄되지 않아 약7/9의 기울기를 만드는 구조와 일치한다. 기존 LSF를 새 비교의 기준으로 삼는 것이 타당하다. CPU/GPU 공통 수치식 결함이며 CUDA만의 문제가 아니다.

### 강체회전도 왜곡 셀의 제한자에서 손상된다

`irregular_affine_probe/results.json`: 독립 검증한 실제 t6 셀의 위치/부피를 사용하되, 균일 rho1/P1e-6과 정확한 v=.1J(x-center)로 초기화한 native 진단. 중력 없이 첫 순간 연산자를 검사하며 Kepler 장시간 안정성 검증이 아니다. 정상 격자의 LSF 강체회전은 기울기오차약5e-14, pstar/P≈1을 통과했다.

왜곡 격자에서는 같은 LSF라도 R.5–1의 pstar/P 최대30.60, R2–4 최대28.60이 생겼다. 중심 기하 Vdot/V는~1e-16으로0인데도 가열된다. 일부 기울기의 상대오차는약.64다. 정확 affine strain에서도 같은 수준의 오차를 보였다. **Cartesian 속도 성분을 이웃 중심의 min/max로 제한하는 TVD가 선형/회전장을 깎는 것**이 유력하다. 왜곡된 face centroid가 이웃 중심의 값 범위를 벗어날 수 있으므로 선형장조차 new extremum으로 잘못 처리할 수 있다.

후속 velocity-only/merge-only 개입을 완료했다. 주 근거는 `irregular_affine_probe_v2/results.json`이다. 운전자가 v1 진단 입력의 KH trailer에 이전시간간격 대신 Gamma를 기록한 오류를 발견해, 기존 파일은 보존하고 동일 particle payload에 dt/dtold1e-7을 넣은 별도 입력으로 재검증했다. 기하/LSF/TVD 수식은 이 필드와 독립이지만 pair-pressure dt cap은 다르므로 v2 수치를 사용한다. 생산 baseline과 restart reader의 결함은 아니었다.

`affine_limiter_intervention/first_gradient_audit.json`은 실제 native 첫 면 기하에 현재 LSF/TVD 수식을 독립 적용했다. R<4의 제한 전 gradient 최대상대오차는1.14e-14 이하이고, 제한 후에는최대.637이다. R.5–1의14/96셀에서 clipping이 일어나며 최소alpha.203, 정상속도 jump가~8e-17에서.00955까지 커진다. 순수 강체회전에 가짜 대칭 strain도 생긴다. 따라서 제한자가 affine 재현성을 깨는 인과 연결이 확인됐다.

v2 native R2–4 첫 연산자 대조:

| 속도 clipping | 작은 면 병합 | pstar/P0 최대 | 양의 추가 열률 |
|---|---|---:|---:|
| ON | ON |28.599625|5.83357e-8|
| OFF | ON |1.000615|1.72402e-12|
| ON | OFF |28.599625|5.82975e-8|
| OFF | OFF |1.000000111|2.11933e-15|

이는 **진단 전용 개입**이다. `GFS_AFFINE_UNLIMITED_VELOCITY` 및 `GFS_AFFINE_UNMERGED_FACES`는 EUNHA_IC=affine_probe에서만 허용하고 다른 IC에서는 abort한다. 실제 충격파/생산 Kepler에서 모든 제한자를 제거한 해결로 간주하지 않는다. 해결 후보는 강체회전/affine 재현성을 갖는 재구성·제한 규칙이며, 해석적 Kepler 속도로 강제하거나 Re를 낮추지 않는다.

중심 R<2에는 두 개입 후에도 일부 압력이 남았다. `affine_limiter_intervention/pair_residual.json`은 실제 native `gfs_pair.h`를 C wrapper로 그대로 호출해 잔차를 분해했다. 이78개 face-end의 추가압력은 근접 셀의 pair restoring pressure이며, Pstar−mean(P)−Ppair의 최대값은 초기P 대비8.96e-14이다. 정상 생성점 속도차는2.47e-17 이하인데도, 현재 보정은 정지/강체회전 근접쌍에서도 열적 반발압력을 준다. 최대 추가압력15.694P0이다. 이 항은 알려진 근접쌍 보호장치의 작동이며 새 점성 원인을 발명한 것이 아니다. 제한자 수정만으로 중심부 모든 힘이 물리적이 되는 것은 아니다. 충격에서의 안전성을 보존하면서 smooth 회전/전단을 침범하지 않는 근접쌍·기하·보존 계약이 필요하다.

### 마지막 부피 갱신의 정확도와 속도

`final_volume_gpu/periodic_area_audit.json`: 주기 영상을 포함한 독립 전역 Voronoi 재구성으로 모든65536셀의 이웃 완전성을 기하학적으로 증명했다. totalarea576. GPU raw volume은 t6/t8 모두 최대상대오차약6e-8, 1e-6 초과셀0. 기존 CPU end-volume은 고정3x3 후보 탐색으로 t6의24셀, t8의36셀을 잘못 계산하며 최대상대오차9.16/21.37이었다. GPU force의 adaptive complete search를 CPU final path까지 사용하지 않던 실제 결함이다. 앞의64240개 안전한 내부셀 t8 감사와 범위가 다르며, 기존 CPU가 전역에서 정확했다는 뜻이 아니었다.

동일 binary에서 flag0/1로 t6→6.01009의15스텝을 비교한 벽시계는29.521초→16.168초(관측1.826배), CPU final timing은약1.0→0.1초/step이었다. 실행 중인 기준선과 H200을 공유한 순차 측정이므로 전용 GPU 성능/강한 scaling의 입증은 아니다. CPU의 잘못된 부피를 고치므로 여러 스텝 후 그 셀까지 정확히 같아야 한다고 요구하지 않는다. pressure/stress/U는 GPU geometry 함수에서 갱신하지 않고 원래 caller의 EOS/에너지 장부에 맡겼다. 비양수/비유한/halo누락 geometry는 abort한다.

### 작은 면 병합과 공유면의 상호성

`face_reciprocity_audit.json`: t6의392739면끝 중13개, t8의392809개 중5개가 반대쪽 기록을 갖지 않는다. 나머지 일부 면도 양쪽 면적벡터가 정확히 반대가 아니며 최대상대차이는.0140/.00634다. t6 R<.5에서는 양쪽 면이~1e-12로 일치해, 이것을 초기 중심 가열 전체의 원인으로 삼을 수는 없다. 전체 hydro work 잔차는t6 2.96e-8/t8 5.85e-7이다.

독립 개입 `affine_limiter_intervention/face_reciprocity.json`: 동일한 실제 왜곡 geometry에서 작은 면 병합만 끄면 orphan이12→0, 최대 면벡터 절대차이가1.669e-4→2.336e-14, 최대상대차이가.00407→1.775e-10으로 줄었다. 병합이 공유면 불일치를 만드는 원인임을 확인했다. 단순 양쪽 벡터 평균보다 원래 닫힌 셀과 공통 면 support를 보존하는 방식으로 tiny-face 처리/CFL을 수정해야 한다. 병합 제거 진단의 좌표/속도는 변경하지 않았다.

**운전자 수정 우선순위:** 입력/halo → 완전한 부피·공유면 → LSF 및 강체회전/affine 보존 재구성 → smooth 회전을 침범하는 근접쌍 보정의 역할 분리 → 고정 eps에서 중심 해상도/중력 공간 결합 → 남는 물질부피 결함에 대한 가중치/모멘트 보존 폐쇄. 회전 모멘트 연구는 계속 활용하되 위의 확정된 공간 연산자 결함을 spin 감쇠로 가리지 않는다. 새로운 증거와 이 해결 후보를 한 번 묶어 Astra에 후속 평가를 요청했고 동의 및 보완을 받았다(`ASTRA_T8_FOLLOWUP_REVIEW.md`). Astra는 affine 정확성, 좌표 회전/균일 boost 불변성, rho/P 양성 및 충격 처리를 함께 요구하며, 공통 velocity 제한계수 하나나 전역 clipping 해제로 해결했다고 간주하지 말 것을 권고했다.

### 사용자 질문: 실제 Kepler 안쪽은 강체회전인가?

실행 코드 `Exam/KH/util.c:493` 이하의 IC와 `Exam/exam.c:4454`의 중력을 대조했다. Rin2는 밀도 profile의 경계이고, R<Rin에 강체회전을 별도로 강제하는 분기는 없다. 전 영역에 같은 Plummer 원궤도 속도 vphi=R(R²+eps²)^(-3/4)를 쓴다. 따라서 Ω=(R²+eps²)^(-3/4), q=-dlnΩ/dlnR=(3/2)R²/(R²+eps²). 강체회전은 R≪eps 극한(q→0)이고 R≫eps는 q→1.5의 Kepler 차등회전이다.

현재 eps.047, h.09375, 초기 최소 생성점 반경h/√2=.066291>eps이므로 강체회전 극중심을 처음부터 분해하지 못한다. 가장 안쪽 생성점도 q=.998223으로 유의한 전단을 겪는다. t8 압력 peak R≈.028에서는 q≈.393으로 회전이 더 우세하지만 정확한 강체회전은 아니다. R<.5 영역 전체를 강체회전으로 부르지 않는다. 별도의 affine rigid 진단은 native 연산자의 회전 재현성을 분리한 대조이며 실제 Kepler의 중심 곡률·물질부피 문제를 대체하지 않는다.

## 실행/파일 보존

- 기존 baseline429482는 기존 binary로 계속 실행한다. 새 build가 실행 중인 파일을 바꾸지 않는다.
- 초기 빌드430017은 QOS 대기에서 소유 확인 후 취소하고,429482의 할당 내 CPU1개로 빌드를 완료했다. MPI direct-PMI 시험 하나가 초기화 전에 실패했고 원본 log를 보존했다. 이후 Slurm step 내부 local Hydra fork로 native MPI를 실행했다.
- 구형 binary 성능 대기430019/430020, direct-PMI 대신 재사용한430021은 소유 확인 후 대기 상태만 취소했다. 다른 작업은 건드리지 않았다.
- 모든 추가 재생/분석은 Slurm 할당 내에서 끝나는 명시적 step이다. 백그라운드 감시 loop는 없다.
- 원본 binary와 소스 before/after·patch·SHA를 보존했다. 신규 fixed binary: `density_snapshot/eunha2_stable_input.exe`.
- 사용자 요청에 따라 Astra에 한 번 묶어 검토를 요청했고 답변을 받았다: `t8_review_request.md`. 우선순위는 입력→장부→중력/해상도 분리로 합의했다. 검증 결과로 운전자 판단을 보완했다. 커밋·푸시는 하지 않았다.

최종 상태 확인: 2026-10-20 19:14 KST, 기준선429482 RUNNING, t9.85185/step17136, dEtot/|Etot0|1.165e-4, 누적 floor9.045e-4. 기존 계산은 보존했다. 목표 미완료; 장시간 안정성은 아직 입증되지 않았다.

## 2026-10-20 후속 갱신: 기존 기준선 완료와 실패 원인 분리

이 절이 위19:14 상태를 갱신한다. 기준선429482는 t10.000582/step17446에서 정상 종료했다. 그러나 중심 가열·궤도 오류가 남으므로 물리 안정성 통과가 아니다. 실제 Voronoi 밀도장 그림은 `linear_velocity_limiter/baseline_t10_density.png`에 보존했다.

- 새 공통 vector velocity limiter는81개 helper 검사와 왜곡 강체회전 patch/MPI1·4 일치를 통과했으나, 같은 binary의 cold Kepler2×2에서는 중심 entropy L1이34.44→110.8로 악화됐다. 작은 면 병합 여부와 거의 무관하다. 기본 채택 보류. 상세: `linear_velocity_limiter/RESULTS.md`.
- 첫 실제 cold RHS의 C/native-helper trace에서 중심 제한계수가~1e-13으로 내려가 전체 속도 tensor를 거의 지운다. 최대 가열 면의 압축 정상속도차는 old~.77, vector~1.26, unlimited 진단~.54다. 작은 면이 관여하는 사례도 있지만 보통 길이의 면에서도 발생한다. Unlimited는 진화에 사용하지 않았다.
- raw Sedov는 old limiter에서도 양성 실패. 엄격한 RK 단계 검사로 최종 ledger보다 이른 step256/k4_input, t_start.0245808의 음의 U2셀을 포착하고 정상 step-start를 보존했다. 같은 한 구간의 full step이 원래 실패 상태를 정확 재현했고, 두 half step에서도 두 번째k4_input의 한 셀이 음수였다. 최종 U가 양수여도 중간 단계가 비물리적일 수 있다. 상세: `rk4_diagnostics/RESULTS.md`.
- 고정eps.047에서 중심88개 coarse 셀만8×8로 세분한71080셀 IC를 만들었다. 실제 주기 Voronoi 다각형의 같은 연속밀도를 적분하여 coarse/fine 총질량은 모두75.97420601059443, 총열도 반올림 범위에서 일치한다. fine최소r.008286, h.011719로 연화 중심을 분해한다.
- fine 첫 검사에서 dense seed의 중간 polygon capacity를 확인했다. 최근접8각 sector의 실제 반평면을 먼저 적용한 뒤 원래 완전 탐색을 모두 실행하도록 순서를 보완했다.48꼭짓점 제한을 늘리거나 이웃을 버리지 않았다. 이어 같은 support를 한쪽은1면, 다른 쪽은2면으로 Riemann 계산하던 결함을 발견해 연속 동일 support를1면으로 합쳤다. 상세: `raw_vertex_dedup/RESULTS.md`.
- frozen v4 fine preflight:71080셀 부피 최대상대오차9.095e-15, 공유면벡터 오차2.853e-16, surfaceclosure1.552e-17. netforce~2e-9→~1e-18. 남는 한쪽 FP미세면의 최대길이1.086e-16과 힘8.04e-19는 별도 공개했다.
- v4 MPI4 fine1스텝은 보존한 MPI1 결과와 모든 기록 필드가 비트 단위 일치했다. 새 기하로 Sedov의 동일 실패 구간을 재생해도 full 및 half 단계 음의 열은 남았다. 기하 수정이 RK4 양성 문제까지 해결했다는 뜻이 아니다.

430452 및 의존 분석430454는 완료됐다. 같은 frozen v4 binary의 coarse+old만 t.2에 도달했다. fine+old는 t_start=.10196264987779348/step164/k2_input, fine+vector는 t_start=.11947779279039178/step192/k3_input에서 음의 U로 거부됐다. 두 fine 분기 모두 ID67935,68001,68526,68592 네 셀이 로그에 기록돼 있다. 실패 분기에는 t.2 결과가 없으므로 같은 시각의 최종 정확도 비교라고 제시하지 않는다. 원본과 정상 step-start/실패 stage는 보존했다. 모든 새 GPU 요청은H2001장과 명시적 메모리를 사용했다.

동일 연속 IC/eps.047의 첫 공간 RHS에서 R<.5의 양의 Udot 합은 coarse1.8144767438e-4→fine1.8646217113e-7로 약973배 감소했다. 이는 압력 일도 포함하는 순간 지표이며 순수 비가역 열이나 장시간 수렴차수를 뜻하지 않는다. 전체 입자 수65536→71080(+8.46%), 중심 간격.09375→.01171875. 국소 세분화는 유효하지만 단계 양성 문제를 닫지는 못했다.

430482에서 fine old의 정상 step-start로 같은 한 구간(h=.0006255366384684754)을 native v5로 재생했다. full의 U/V/m은 원래 실패 상태와 비트 단위 일치했고 위치 최대차이3.55e-15/속도6.66e-16이었다. full은k2, half는첫k4에서 같은 네 셀이 음수였다. quarter4회는 동일종료시각.10258818651626195까지 floor 없이 통과했다. 이는 시간간격 민감성만 입증하며 장시간 재개나 자동 timestep 재시도 도입의 근거로 삼지 않는다. v5의 물리 연산자는v4와 같으며 audit ID필터만 추가됐다. 상세: `rk4_diagnostics/fine_cold_interval/manifest.json`.

최근 Astra 권고(`rk4_diagnostics/ASTRA_POST_REFINEMENT_REVIEW.md`)에 따라 실패 셀/이웃의 signed 열 장부와 실제 Voronoi 체적 Jacobian—압력력의 일치성을 검사했다.430502는QOSGrpGRES 대기로 소유/작업경로 확인 후 취소했고, 사용자 허용lageunha에서 CPU1스레드 읽기 전용 분석을 마쳤다. 실패 셀의 추가냉각~2.366e-8은(Pi-Pstar)s_geom으로 재구성되고, pair=0/work속도 차이~1e-21이다. adjoint 가역 기준은power1.84e-20 잔차를 통과했고 독립체적 미분은delta반감 시4배 수렴해L2relative8.73e-9였다. 이는 연산자 기준이며 장시간 안정성 증명이 아니다.

결과와10월9일 가역 sliding-cell 불안정 이력을 묶어 Astra에게 평가를 요청했다. 다음 한 묶음은 기존quarter4구간의 압력비100분기/경계U→0 판별이었다(`rk4_diagnostics/ASTRA_BRANCH_REVIEW.md`). 물리 연산자/임계값을 바꾸지 않고 선택 셀의16회force call 입력·실제flux를 기록하는v6를 작성했다. Nativeforcekernel 산술은 수정하지 않았다. 사용자가 허용한lageunha의 유휴Ada에서 동일quarter4를1회 실행했고, H200의 기존quarter4와 체크포인트별 모든 기록필드 및metadata가 비트 단위 일치했다. Native실행5.54초(성능벤치마크가 아님). 수동런처의Slurm환경변수 조기평가 오류는 MPI실행 전 발견/수정했으며 실패준비디렉터리를 보존했다.

저장 입력의 양분기 평가도 완료했다(`rk4_diagnostics/branch_trace/RESULTS.md`). Native face 재구성 최대오차는 P7.30e-18/work속도8.15e-16/Udot9.71e-21이다. ID68526의 첫 quarter는 압력비96.4/129.6/91.8/195.3에 따라 extreme면0/2/0/4개가 바뀌며 Udot가 음/양/음/양으로 전환된다. 실제 전환면의 동일입력 양분기 차가 관측 jump를 재현하여 chatter 기여를 확인했다. 그러나 U→0+의 cell-HLL 한쪽 극한은 네 셀 모두+2.985–2.988e-9로 안쪽을 향한다. 이 상태에서 유한 음의 경계 열률이 남는다는 가설은 지지되지 않는다. RK4 단계 양성 문제와 누적 열 고갈/수송 계약을 구별해야 한다. 임계값·분기·시간적분법·가중치·회전 물리는 바꾸지 않았다. 최종 장시간 목표는 미완료다.

Astra 검토는 결과 묶음 단위로 요청했다(`linear_velocity_limiter/ASTRA_REVIEW.md`). eigenaxis limiter나 임의 spin 감쇠를 추가하지 않았다. 장시간 최종 요구는 `GIZMO_COMPLETION_CONTRACT.md`에 고정했다. Fig7 t600뿐 아니라 본문의~1000orbit 언급까지 최종 목표에 남아 있다. 새 짧은 실험으로 기존 t6–10 실패나 최종 장시간 검증을 대체하지 않는다. 커밋·푸시 및 사용자 공통 context 변경은 하지 않았다.

## SSPRK3 묶음 결과 갱신

Astra의 묶음 계획에 따라 운전자가 opt-in native SSPRK3와 collective full rollback을 구현했다. 기본 RK4와 frozen reference 일치, 주기 경계/MPI 소유 이동 뒤 강제 rollback의 MPI1/4 전체 기록필드·metadata 비트 일치를 확인했다. 저장된 fine/Sedov 실패구간은 모든 중간 단계 양성으로 통과했다. Fine nominal half ceiling은 실제 dt를 바꾸지 않아 검증으로 세지 않았고, 한 번의 실효 half 대조(3→5step)를 추가했다. 국소 열 고갈은 그대로였다. Sedov는 stage 양성과 시간 민감도까지만 확인했으며 충격 정확도 관문 통과로 분류하지 않았다.

새 경로는 현재 geometry→EOS→halo→gradient/force를 준비한다. 기존 GPU gradient가 직전stage rho/P cache를 받을 수 있던 입력 준비와 차이가 있으므로 순수하게 동일RHS의 시간적분법 대조라고 주장하지 않는다. 압력비100분기·force/work·limiter 공식은 바꾸지 않았다. 고정질량/W0/무조향/Re∞ 후보이며 dynamic Laguerre/spin 완성본은 아니다.

Fresh71080셀은 사용자 허용lageunha/Ada에서 한 경로만 수행했다. t=.13509103681003889,515accepted/3072RHS/584rollback/873.38초에서 미리 정한 예산으로 종료했다. t=.2 미도달이다. 정상 시작점 `rk4_diagnostics/ssprk3_fresh_v2/fine_fresh/khdiag.step000516.ssp_rejected_start.dat`를 보존했고 재시작/새 IC 진화는 추가로 시작하지 않았다. 개별질량/W0 비트고정, floor0, minU3.43356e-12이나 네 셀 U/Uad=.01947-.02108로 심한 열 고갈이 남았다. 현재 기하의 독립 부피 감사는~5.75e-8로 float 저장 정밀도 안이다. SSP는 시간 안전장치로 유용하나 이번 장시간 해법 후보는 통과하지 못했다. 자세한 결과/한계는 `rk4_diagnostics/SSPRK3_RESULTS.md`와 `ssprk3_final_analysis_v2.json`에 있다. 최종 자료를 Astra에게 한 번 묶어 전달하여 다음 공간물리 묶음의 평가를 요청했다.

사용자의 안쪽 분해능 증가 제안은 초기 RHS에서 이미 유효한 증거(약973배 감소)가 있으나, 물리 전단률 qΩ 자체는 h를 줄여도 남는다. 현 Plummer전단률은 r=2eps/sqrt3에서 최대이며 중심에서는0이다. h(r)∝sqrt(r²+eps²)를 초기 배치 후보로 제안했다(아직 새 IC/진화 없음). 초기 셀 적분질량을 다르게 배정하고 각 질량을 고정하면 밀도장/환원조건을 보존할 수 있다. 이 제안도 기존 sliding-cell 이력과 함께 평가 중이며 장시간 안정성 증명으로 삼지 않는다.

## Native ALE shadow 묶음 완료

Astra의 다음 묶음 권고에 따라 사용자의 CUDA/MPI 본체에 actual-face conservative Euler flux를 읽기 전용으로 추가했다. 기본은 off이며 상태에 적용하지 않는다. lower-global-ID 소유 rank가 공유pair를 한 번 CUDA 계산하고 단일 packet으로 정확한 ±수송 장부를 만든다. 전체 mass/momentum/energy flux와 pressure/advection/근사Riemann 잔차를 분리한다. frozen binary SHA `fbdb5d2bea0354b30c34d82ee60198c101f51ece1d5b2d22a0888edb53217d48`.

왜곡64² uniform/rigid/shear/boost/MPI4의 첫RHS를 무진화로 기록했다. Shear의 mdot-rhoVdot max1.14e-17, 실제Vdot max1.04e-4였으며 고정mass라면 rho_dot max.4107/time이다. 이는 실제polygon을 exact Euler control volume으로 쓰는 유한격자 계약의 차이이며 모든fixedmass particle법의 불가능성/비수렴 정리는 아니다. 2point native Gauss와 독립3point polynomial 적분은~1.5e-15 natural-scale 일치, boost오차≤5.07e-15 relative, MPI1/4 fullflux 비트 일치. 실제native uniform P의float EOS오차 때문에 residual3.54e-12를 machinezero로 분류하지 않았다.

fine의 보존한 두 정상상태도 첫RHS만 각각~3초에 기록했다(진화 없음). 기존3step fine구간을 계측 observer검증 목적으로만 한 번 재생하여 frozenSSPv2와 전체기록필드/metadata 비트 일치를 확인했다. 실제fullALE 열률은 한상태에서가열, 다른상태에서큰냉각을 보이므로 곧바로 대체하지 않았다. 개별mass/W/thermalU는 그대로다.

추가실행 없이 저장된면/현재속도로 exactpolygon moments와 독립unlimited affine fit를 계산했다. uniform-density affine-cell 모델의 Kvariance 추정은 대표셀에서 현재U의1697/3683배이며, 실제native limited gradient/실제spin/알 수 없는subcell rho의mass moments라고 부르지 않는다. 첫상태 nativeUdot -2.369e-8 중 fitted rotationpressurepower -2.750e-8 +strain +4.520e-9로 남는 일이 -7.135e-10이다. 하지만 생성점/기하중심 기준 mean/rotation 일배분이 크게 달라 단순열보정은 하지 않는다. Fulladvective momentum의 orbital torque와 필요한 intrinsic torque 합도~2.12e-22로 닫혔으나 spin-energy 진화식은 아직 아니다.

전체 구현/가정/출력/한계는 `rk4_diagnostics/ALE_SHADOW_RESULTS.md`; 원자료/분석은 `ale_shadow_patches_v1/`, `ale_shadow_real_v1/`, `ale_shadow_observer_v1/`에 보존했다. 이번 완성자료를 Astra에게 한 번 묶어 평가 요청했다. 다음 모델은 질량-부피-열 및 angular/affine moment 계약을 함께 닫아야 하며, 현재 장시간검증은 미완료다. 새 graded진화·DAE 재활성·임의감쇠·커밋/푸시는 하지 않았다.
