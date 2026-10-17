## 결론

**B(영역 중심) 정식화는 차가운 Kepler 원반의 단기 셀 규모 모드 억제에서 유망합니다.** 다만 현재 증거는 5–8궤도 수준이고, 뜨거운 배경의 성장은 남아 있으며, 재배치 코드에는 실제 질량·내부에너지 생성 경로가 있습니다. 따라서 C 포트는 기존 TCF 경로와 분리한 CPU 연구 분기로 설계하고, **무원천 보존 재배치·영역 인공점성·재시작 계약**을 통과한 뒤 승격하는 것이 맞습니다. 문서와 소스는 읽기만 했습니다.

## 1. 정식화와 재배치

### (a) 고정 관성과 재배치 질량

연결이 고정된 동안 \(U=\sum_T K_Tm_T^\gamma A_T^{1-\gamma}/(\gamma-1)\), \(F_a=\sum_TP_T\,\partial A_T/\partial x_a\)는 일관된 변분 계입니다. 각 삼각형 넓이는 평행이동과 회전에 불변이므로 내부 압력력의 합력과 토크 합은 0입니다. 이산 운동량 \(\sum_a m_a v_a\)와 각운동량 \(\sum_a m_a x_a\times v_a\)는 이 힘에 대해 보존됩니다. \(m_a\)가 초기 질량으로 고정됐다는 사실은 이 수학적 보존을 깨지 않습니다. 정식화와 힘은 [staggered_verification_2026-10-09.md:3](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:3), [stagP.py:63](/tmp/lageunha-push/Docs/polybnd/stagP.py:63), [stagP.py:67](/tmp/lageunha-push/Docs/polybnd/stagP.py:67)에 있습니다.

다만 이는 **고정 질량 생성점에 대한 보존**입니다. 입자 관성 \(m_a\)는 처음의 인접 영역 질량 ⅓ 합으로 정하고 계속 고정하지만, 재배치되는 영역 질량 \(m_T\)에는 운동량 변수가 없습니다. 코드의 중력 퍼텐셜 에너지도 고정 \(m_a\)를 씁니다([stagP.py:64](/tmp/lageunha-push/Docs/polybnd/stagP.py:64), [stagP.py:71](/tmp/lageunha-push/Docs/polybnd/stagP.py:71)). 그러므로 힘에 대한 이산 운동량·각운동량 보존은 성립해도, 재배치되는 유체 질량의 운동량 및 중력 질량까지 일관되게 표현한다는 주장은 아직 성립하지 않습니다. 자가중력을 붙일 때는 영역 질량을 입자 힘에 배분하는 규칙과 관성·중력 질량을 함께 다시 정의해야 합니다.

### (b) 엔트로피와 규모

옛 영역의 겹침으로 새 영역의 질량과 내부에너지를 적분하면, 전체 겹침이 정확하고 누락·중복이 없을 때 두 양은 보존됩니다. 서로 다른 상태를 새 영역 안에서 평균내는 것은 엔트로피를 증가시키므로, 이 재배치는 **명시된 비가역 투영 소산**으로 정당화할 수 있습니다. 이상기체에서 \(s=c_v\ln K+\mathrm{const}\)이므로 정확한 장부는

\[
\Delta S_{\rm remap}
=c_v\left(\sum_{T'}m_{T'}\ln K_{T'}-\sum_Tm_T\ln K_T\right)
\]

로 누적해야 합니다. 현재 출력의 질량 가중 평균 \(K\)는 엔트로피 그 자체가 아닙니다([stagP.py:90](/tmp/lageunha-push/Docs/polybnd/stagP.py:90), [stagP.py:98](/tmp/lageunha-push/Docs/polybnd/stagP.py:98), [stagP.py:117](/tmp/lageunha-push/Docs/polybnd/stagP.py:117)).

중요한 단서가 있습니다. `m_new <= 0`이면 코드가 `Sig`와 `P0`로 새 질량·에너지를 채웁니다. 이는 보존 재배치가 아니며, 검증 문서도 무질서 상태의 이 경로를 질량 생성 결함으로 기록합니다. 이 fallback은 금지하고, 겹침 실패는 재시도 또는 단계 실패로 처리해야 합니다([stagP.py:94](/tmp/lageunha-push/Docs/polybnd/stagP.py:94), [stagP.py:97](/tmp/lageunha-push/Docs/polybnd/stagP.py:97), [staggered_verification_2026-10-09.md:15](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:15)).

요청서의 \(2\text{–}4\times10^{-5}\) / 10 코드 시간과 영역당 약 0.9회/시간 재배치를 그대로 쓰면, 평균적으로 재배치 한 번당 \(\Delta\langle K\rangle/\langle K\rangle\sim2.2\text{–}4.4\times10^{-6}\)입니다. 단순한 동일 밀도·작은 차이 혼합 모형에서는 이는 이웃 상태의 유효 대비가 대략 \(0.2\text{–}0.3\%\)인 크기와 맞습니다. 정확한 대비는 실제 겹침 질량·에너지 분포에서 계산해야 합니다.

반지름 \(R=0.84\)의 비연화 Kepler 주기는 약 4.84 코드 시간입니다. 위 드리프트율을 100궤도까지 선형 외삽하면 평균 \(K\) 변화가 약 \(0.10\text{–}0.20\%\)입니다. 이는 규모 추정이지 장기 예측은 아닙니다. 검증 문서의 96² 실행은 \(t=26\)에서 \(\langle K\rangle/K_0-1=2\times10^{-5}\)를 보고하므로, 이 수치만 선형 외삽하면 약 \(0.037\%\)입니다([staggered_verification_2026-10-09.md:9](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:9), [staggered_verification_2026-10-09.md:14](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:14)). 관측 창과 정규화를 통일하고 100궤도 누적을 직접 측정해야 합니다.

### (c) 단계 내 연결 고정과 flip 시점

**현재처럼 RK4 한 스텝의 네 단계에서 연결을 고정하고, 승인된 스텝 끝에서만 재삼각분할·재배치하는 방식이 타당합니다.** 그래야 한 스텝 동안 같은 \(U(A)\)의 매끄러운 힘을 적분합니다. 코드도 네 RK4 단계 뒤에 연결을 재구성하고 겹침 재배치를 합니다([stagP.py:78](/tmp/lageunha-push/Docs/polybnd/stagP.py:78), [stagP.py:83](/tmp/lageunha-push/Docs/polybnd/stagP.py:83), [stagP.py:86](/tmp/lageunha-push/Docs/polybnd/stagP.py:86)).

단계마다 재삼각분할하면 RK4의 네 단계가 서로 다른 연결과 다른 상태 투영을 보게 됩니다. 이 경우 시간 4차 정확도를 기대할 수 없고, 재배치 횟수와 혼합도 늘어납니다. flip을 스텝 중간에 반드시 처리해야 한다면 이벤트 시각에서 적분을 나누고, 그 지점에서 보존 투영한 뒤 남은 구간을 새 연결로 적분하는 방식이 필요합니다.

### (d) 수렴 차수와 100궤도

코드 주석이 명시하듯 겹침 재배치는 **보존적 1차 공간 재매핑**입니다. 연결 고정 구간의 RK4 4차 정확도가 전체 해법의 차수가 되지는 않습니다. 반복 재배치가 누적하는 오차를 공간 해상도와 시간 간격으로 분리해 측정해야 합니다([stagP.py:5](/tmp/lageunha-push/Docs/polybnd/stagP.py:5), [stagP.py:106](/tmp/lageunha-push/Docs/polybnd/stagP.py:106)).

현재 64² \(t=40\)과 96² \(t=26\)의 결과는 각각 약 8궤도와 5궤도입니다. 이 결과로 100궤도를 보증할 수 없습니다([staggered_verification_2026-10-09.md:7](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:7), [staggered_verification_2026-10-09.md:9](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:9)). 64²/96²/더 높은 해상도와 \(dt\), \(dt/2\)를 비교하고 모드 진폭·정확한 엔트로피·최소 영역 넓이·에너지 잔차의 누적 추세가 수렴하는지 봐야 합니다.

## 2. 뜨거운 배경의 잔여 성장

영역 넓이는 연결이 유지되는 동안 그 **이산 삼각형**의 정확한 기하 넓이입니다. 이것만으로 실제 유체 물질 부피를 모든 비아핀 흐름에서 보장하는 것은 아닙니다. 특히 연결이 바뀌면 서로 다른 옛 상태를 새 삼각형 안에서 섞고, 그 과정에서 새 압력 분포가 생깁니다. 따라서 잔여 성장의 우선 가설은 **flip 재배치가 만든 압력 잡음이 씨앗을 제공하고, 뜨거운 배경의 압력 응답이 이를 키운다**는 것입니다. 다만 삼각형 압력력의 분산 오차와 RK4 스텝 끝에서만 처리하는 flip의 위상 오차도 경쟁 가설입니다.

자료는 차가운 원반에서 64² 심은 모드가 \(0.0300\to0.0313\)으로 8궤도 유지되고, Nyquist 모드도 안정한 반면, \(c_s/d=0.72\)인 안쪽 배경에서 성장률 \(0.15\text{–}0.2\)를 보입니다. \(c_s/d=0.23\)으로 낮춘 경우에는 \(0.02\text{–}0.03\)입니다([staggered_verification_2026-10-09.md:7](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:7), [staggered_verification_2026-10-09.md:10](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:10), [staggered_verification_2026-10-09.md:12](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:12)). 이전 면 중심 시험에서도 셀 규모 선형 불안정이 \(c_s/d\)와 연결되고 전단률에는 둔감하다고 보고됐지만, 이는 같은 영역 중심 정식화의 원인 증명은 아닙니다([HANDOVER_2026-10-04.md:222](/tmp/lageunha-push/Docs/HANDOVER_2026-10-04.md:222)).

판별은 다음 로그를 flip 이벤트마다 기록하는 것으로 시작해야 합니다.

- flip 직전·직후의 모드 진폭, \(\sum m\ln K\), 압력 잡음, 에너지·질량 잔차
- flip 사이에서 성장하는 모드와 flip 시각의 점프 크기
- 동일 초기 조건의 짧은 고정 연결 대조 구간과 flip 시각을 국소화한 구간
- \(dt\) 절반 및 해상도 증가에 따른 성장률 변화
- 주기 아핀 전단·균일 압력 시험: 연결 고정 상태에서 힘과 압력 잡음이 생기는지

\(c_s/d\approx0.07\) 원반 조건은 요청된 실험에서 이미 짧은 시간 직접 시험됐으므로 단순 외삽만은 아닙니다. 다만 시험은 비주기 점 집합과 정지된 바깥 영역을 사용합니다. 코드 주석도 주기 영상점 설명과 달리 실제 초기화는 비주기라고 적고 있습니다([stagP.py:36](/tmp/lageunha-push/Docs/polybnd/stagP.py:36), [stagP.py:43](/tmp/lageunha-push/Docs/polybnd/stagP.py:43)). 그러므로 이 결과를 주기계나 100궤도 결과로 옮겨 말할 수는 없습니다.

## 3. 얇은 영역과 겹침 실패 규칙

**최소 넓이로 압력을 자르거나 빈 겹침에 질량을 넣는 규칙은 채택하지 않는 것이 맞습니다.** 현 코드는 스텝 내부에서 면적이 0 이하인지 검사하지만, 작은 양수 면적은 막지 않습니다([stagP.py:79](/tmp/lageunha-push/Docs/polybnd/stagP.py:79)). 검증에서는 \(0.005d^2\) 미만 영역이 \(P/P_0=10^2\text{–}10^5\) 폭주를 일으켰습니다([staggered_verification_2026-10-09.md:15](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:15)).

권고 규칙은 다음과 같습니다.

1. 겹침 면적 적분의 누락 또는 무겹침 새 영역을 만나면 전체 스텝을 거부하고, 더 작은 \(dt\) 또는 이벤트 국소화로 재시도합니다. 반복 실패는 진단과 함께 중단합니다.
2. 얇은 영역은 우선 보존적인 국소 재삼각분할로 피합니다. 유효한 대체 연결이 없을 때는 인접 영역을 **기하 영역 자체로 합친 다각형 거시영역**을 후보로 둡니다. 이때 상태는 \(m\)과 내부에너지 \(U\)를 합산하고, 합집합 넓이에서 \(K,P\)를 다시 구해야 합니다. \(K\)를 이웃에 단순 복사하거나 작은 영역을 유지한 채 변수를 흡수하면 질량·부피 장부가 모호해집니다.
3. 현재 Python 코드의 `Sig/P0` fallback은 삭제합니다. 질량과 내부에너지는 재배치 전후 전체 합을 검사하고, 엔트로피는 별도로 기록합니다.
4. 거시영역 병합이 아직 구현되지 않은 초기 C 단계에서는 품질 기준에 걸린 스텝을 실패 처리합니다. 압력·넓이 바닥값으로 숨기지 않습니다.

## 4. 충격에서의 영역 \(Q\)

압력만으로 충격을 포착한다고 가정해서는 안 됩니다. Wilkins는 다차원 von Neumann–Richtmyer 인공점성으로 충격을 다루고, Caramana 등은 다차원 점성이 압축에서만 소산하고 확장·강체 회전·충격 접선 방향에서 불필요하게 작동하지 않도록 하는 조건을 제시합니다([Wilkins 1980](https://www.sciencedirect.com/science/article/pii/0021-9991(80)90161-8), [Caramana et al. 1998 전문](https://cnls.lanl.gov/~shashkov/papers/ed_vis.pdf)).

이 정식화의 최소 보존형 영역 점성은 삼각형별 \(Q_T\ge0\)를 두고 다음처럼 같은 기하 연산자에 넣는 것입니다.

\[
F_a=\sum_T(P_T+Q_T)\frac{\partial A_T}{\partial x_a},
\qquad
\dot U_T=-(P_T+Q_T)\dot A_T .
\]

그러면 노드 운동에너지 변화율과 영역 내부에너지 변화율이 반대칭으로 닫히고, 압축 \(\dot A_T<0\)에서 점성 가열은 \(-Q_T\dot A_T\ge0\)입니다. \(Q_T\)는 팽창·강체회전에서 0이어야 하며, 부드러운 아핀 흐름의 가열을 억제하는 압축 센서가 필요합니다.

TCF의 \(Q_f=Z_{\rm eff}\,\mathrm{minmod}(\Delta u_{\rm rec},\Delta u_{\rm cell})\)와 QBP의 재구성 잔차는 **조건과 진단 아이디어**로 활용하되, 현재 Voronoi 면의 \(Q_f\) 공식을 영역 \(Q_T\)에 그대로 복사하면 안 됩니다. QBP는 현재 선형·이차 속도장을 0으로 만드는 Voronoi 면 입력입니다([plan_thermo_compatible_flux.md:55](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:55), [exam.c:679](/tmp/lageunha-push/Exam/exam.c:679)). 새 영역에서는 삼각형의 속도 기울기와 인접 삼각형 경계에 맞춰 다시 유도해야 합니다.

검증 순서는 **1D Sod → 2D 축정렬·회전 Sod → Gresho**를 권합니다. Sod에서 충격 위치·후충격 상태·엔트로피 증가와 총에너지 장부를 확인하고, Gresho에서 \(Q\)가 매끄러운 회전에 가열을 만들지 않는지 확인한 뒤 Kepler로 진행합니다. 프로젝트 TCF 계획도 Sod 등 충격 검증을 별도 관문으로 두고 있으며, 영역 검증 문서는 Sod·Gresho가 아직 미검증이라고 기록합니다([plan_thermo_compatible_flux.md:86](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:86), [staggered_verification_2026-10-09.md:17](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:17)).

## 5. 프로젝트 원칙 및 Laguerre와의 관계

이 방식은 **입자를 생성점으로 유지하고, 조향 없이 압력 에너지의 변분력으로 움직인다**는 점에서 입자=생성점·조향 금지 원칙과 맞습니다. 열역학 부피는 Voronoi 셀이 아니라 생성점으로 만든 Delaunay 영역의 넓이입니다. 따라서 “부피는 기하가 정한다”는 원칙에는 부합하지만, 현재 Voronoi 부피가 열역학 부피라는 TCF 구조와는 다른 이산 모델입니다. 현재 코드의 centroid shift 옵션은 실제 생성점 이동이므로 이 분기에서는 꺼야 합니다([plan_thermo_compatible_flux.md:9](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:9), [exam.c:14687](/tmp/lageunha-push/Exam/exam.c:14687)).

총에너지는 연결 고정 구간에서 반이산적으로 보존되고, RK4에서는 시간 이산 오차가 남습니다. 재배치는 내부에너지 합을 보존할 수 있지만, fallback이나 부정확한 겹침은 그 성질을 깨뜨립니다. 따라서 “정확 보존”은 **영역 상태와 재배치까지 포함한 조건부 목표**로 남겨야 합니다.

Laguerre 가중치는 이 모델의 삼각형 넓이를 직접 조절하지 않으므로 열역학 부피 제약 변수로서는 역할이 없습니다. 가중 정규 삼각분할로 연결 선택에 개입할 수는 있지만, 넓이 기반 힘과의 변분 일관성을 다시 유도해야 합니다. 가중치 없는 영역 정식화는 25차 답신의 권고와 이어지지만, 현 Laguerre 물질부피 추종 경로를 재개하는 것은 아닙니다([audit_impl25_astra.md:21](/tmp/lageunha-push/Docs/audit_impl25_astra.md:21), [audit_impl25_astra.md:36](/tmp/lageunha-push/Docs/audit_impl25_astra.md:36), [HANDOVER_2026-10-04.md:252](/tmp/lageunha-push/Docs/HANDOVER_2026-10-04.md:252)).

자가중력에서는 특히 질량 배분을 다시 정해야 합니다. Voronoi 밀도를 진단용으로 계산하는 것과 그것을 실제 중력원으로 쓰는 것은 다릅니다. 고정 입자 관성 \(m_a\)와 재배치 영역 질량 \(m_T\) 중 어느 쪽을 중력원으로 쓸지, 힘·에너지와 함께 명시해야 합니다.

## 6. C 구현 계획

| 경로 | 재사용할 부분 | 추가·분리할 부분 |
|---|---|---|
| Voronoi 면 루프 / TCF | RK4 호출 구조, 보존 장부와 진단 패턴, QBP의 매끄러운 흐름 판별 아이디어 | TCF 면 힘·에너지 식을 영역 압력력에 섞지 않습니다. 현 면 루프는 Voronoi 면의 \((P_i,P_j)\), TCF \(Q\), 면 일 교환식을 계산하는 별도 모델입니다([exam.c:8527](/tmp/lageunha-push/Exam/exam.c:8527), [exam.c:8638](/tmp/lageunha-push/Exam/exam.c:8638)). |
| RK4 blend | 입자 위치·속도 단계, 외력 평가, 스텝 제한, 재시도 골격 | 별도 기본 꺼짐 분기를 둡니다. 각 단계에서 고정 연결의 \(A_T,\partial A_T/\partial x,P_T\)를 계산하고 힘을 노드에 합산합니다. 상태 변수 \(m_T,K_T\)는 스텝 중 고정하고, 종료 후 한 번 재삼각분할·재배치합니다. |
| GPU 밀도 | 기존 Voronoi 밀도를 출력·비교 진단으로 유지 가능 | GPU의 현재 압력 EOS는 입자별 Voronoi 부피·밀도·`stress.K`를 사용합니다. GPU Voronoi 넓이 바닥값도 있습니다. 영역법의 열역학 부피에는 재사용하지 말고, CPU 검증 뒤 삼각형 상태와 영역→노드 산포를 별도 구현합니다([exam_gpu_extract.c:134](/tmp/lageunha-push/Exam/exam_gpu_extract.c:134), [exam_gpu.cu:2075](/tmp/lageunha-push/Exam/exam_gpu.cu:2075), [exam_gpu_extract.c:330](/tmp/lageunha-push/Exam/exam_gpu_extract.c:330)). |
| 재시작 | 현재 입자 스냅샷의 전체 구조체 저장·복원 방식 | 영역 상태는 입자 구조체에 없으므로 별도 버전 스냅샷에 연결 키(생성점 ID), \(m_T,K_T\) 또는 \(m_T,U_T\), 병합 상태와 누적 장부를 저장합니다. 현재 KH 스냅샷은 전체 입자 stride와 stress 필드를 보존하지만 영역 연결은 저장하지 않습니다([kh.c:20](/tmp/lageunha-push/Exam/KH/kh.c:20), [kh.c:124](/tmp/lageunha-push/Exam/KH/kh.c:124)). |

특히 현재 재시도 wrapper는 입자 버퍼와 TCF 누적 장부를 복사·복원합니다([exam.c:14133](/tmp/lageunha-push/Exam/exam.c:14133), [exam.c:14186](/tmp/lageunha-push/Exam/exam.c:14186)). 새 영역 배열과 재배치 진단도 재시도 대상에 포함해야 합니다. 현 RK4는 각 단계에서 입자 `ie`와 `stress.K`를 갱신하는데, 이를 영역 상태의 위치에 억지로 저장하지 말고 영역 고유 배열로 분리합니다([exam.c:14435](/tmp/lageunha-push/Exam/exam.c:14435), [exam.c:14665](/tmp/lageunha-push/Exam/exam.c:14665)).

최소 진행 순서는 다음과 같습니다.

1. **1D 주기 Lagrangian 구간:** 균일 압력 정지, 음향파 수렴, Sod와 영역 \(Q\)를 확인합니다. 질량·운동량·총에너지 장부가 시간 간격 축소에 따라 수렴하고 점성 엔트로피가 음수가 아니어야 합니다.
2. **2D 주기 Delaunay:** 영상점을 포함한 실제 주기 연결, 균일 압력·선형 압력 힘 항등식, 순힘·순토크, 아핀 전단, flip 겹침 장부를 확인합니다. 모든 flip에서 질량·내부에너지 잔차를 기록하고 fallback 발생 0, 전역 \(\Delta S\ge0\)를 요구합니다. 이어 Sod 회전 대조와 Gresho를 통과시키고, 심은 셀 규모 모드의 시간·공간 수렴을 봅니다.
3. **비진공 Kepler:** 문서의 \(0.5\le r\le2\) 원반, 배경 0.01, 균일 압력, 스플라인 연화 조건을 재현합니다. 현재 Python 결과의 64² 8궤도 모드 성장률 약 \(10^{-3}\), 96² 약 0, 압력 잡음 0.8–1.2%가 비교 기준입니다. 그 뒤 해상도·\(dt\) 수렴으로 장시간 추세를 확인하고 100궤도 시험을 별도 통과시킵니다([staggered_verification_2026-10-09.md:7](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:7), [staggered_verification_2026-10-09.md:9](/tmp/lageunha-push/Docs/staggered_verification_2026-10-09.md:9)).

따라서 첫 C 구현의 범위는 **CPU 1D→2D 주기 영역 상태·힘·보존 재배치·진단**으로 고정하는 것이 좋습니다. TCF는 현재 기준선과 비교 경로로 유지하고, GPU·Laguerre 연결·Kepler 장기 계산은 앞 단계의 보존 및 견고성 관문을 통과한 뒤 진행하십시오.