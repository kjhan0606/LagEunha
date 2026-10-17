# 결론

현재 코드만으로 셀 규모 가장자리 모드의 원인을 하나로 확정할 수는 없습니다. 우선순위는 **① Kepler 중심·softening 설정이 두 코드에서 같은지 확인, ② 같은 상태에서 C와 Python의 힘·열역학 업데이트가 실제로 같은 연산자인지 확인, ③ 단계별 기하 재생성과 위상 flip이 만드는 작은 교란을 분리**하는 순서가 타당합니다.

RK4의 네 단계와 최종 가중치는 수식상 표준 RK4와 같습니다. 하지만 C는 각 단계 사이에 주기 경계 처리·마이그레이션을 하고, 그 위치에서 Voronoi 기하와 힘을 다시 계산합니다. 또 MUSCL 재구성, 상태별 압력 계산, TCF/인공 점성, 제한자와 단계 바닥값 등이 조건에 따라 Python의 “고정 엔트로피·정확한 에너지 기울기” 기준과 달라집니다. 반면 Python도 매 RK 단계마다 기하를 다시 만들므로 **위상 flip 자체가 C에만 있는 차이는 아닙니다.**

가장 먼저 맞춰야 할 구체적 차이는 Kepler 외력입니다. Python 기준은 원점 기준이고 `EPS=0.01`입니다([dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:6), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:15), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:21)). C는 `EUNHA_IC=kepler`일 때만 그 힘을 켜고, 중심 기본값은 `(2,2)`, softening 기본값은 `0`입니다([exam.c](/tmp/lageunha-push/Exam/exam.c:6384), [exam.c](/tmp/lageunha-push/Exam/exam.c:6386), [exam.c](/tmp/lageunha-push/Exam/exam.c:6391), [exam.c](/tmp/lageunha-push/Exam/exam.c:6411)). 따라서 실행 환경에서 중심과 softening을 명시적으로 맞췄다는 확인이 없으면, 같은 Kepler 문제를 적분하고 있다고 볼 수 없습니다.

## 기준 Python과 C의 차이

### 1. RK4 산술은 같지만, C의 각 단계는 순수한 상태 평가만 하지 않습니다

Python `dyn.py`는 네 번의 힘 평가와 최종 조합을 명시적인 표준 RK4 식으로 수행합니다([dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:27), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:29), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:33)). C도 `k1`부터 `k4`를 누적 갱신하고 마지막에 같은 \(1,2,2,1\) 가중치로 결합합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14333), [exam.c](/tmp/lageunha-push/Exam/exam.c:14404), [exam.c](/tmp/lageunha-push/Exam/exam.c:14475), [exam.c](/tmp/lageunha-push/Exam/exam.c:14563)).

차이는 중간 처리입니다. 각 중간 위치·속도 갱신 후 C는 `postStage`를 부르고([exam.c](/tmp/lageunha-push/Exam/exam.c:14335), [exam.c](/tmp/lageunha-push/Exam/exam.c:14343), [exam.c](/tmp/lageunha-push/Exam/exam.c:14406), [exam.c](/tmp/lageunha-push/Exam/exam.c:14414)), 그 뒤 기하와 힘을 다시 계산합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14356), [exam.c](/tmp/lageunha-push/Exam/exam.c:14427), [exam.c](/tmp/lageunha-push/Exam/exam.c:14498)). KH 경로에서 쓰는 `periodic_postStage_blend`는 좌표를 상자 안으로 감고 입자를 마이그레이션합니다([kh.c](/tmp/lageunha-push/Exam/KH/kh.c:277), [kh.c](/tmp/lageunha-push/Exam/KH/kh.c:292), [exam.c](/tmp/lageunha-push/Exam/exam.c:13586)). 즉 KH의 기본 콜백 자체가 압력·속도 제한을 적용하는 것은 아니지만, 단계별 경계 래핑과 재배치가 있습니다. 벽 경계 콜백은 반사·클램프까지 수행하지만 KH 기본 경로에서는 선택되지 않습니다([exam.c](/tmp/lageunha-push/Exam/exam.c:13644)).

경계에서 떨어진 한 MPI rank의 KH 계산이면 이 처리는 대체로 좌표 표현 변경과 이동 관리입니다. 반대로 입자가 경계를 넘거나 rank가 달라지면 단계 상태에 실제 불연속이 끼어들 수 있습니다. 이 효과는 기록된 진화 설정과 함께 평가해야 합니다.

### 2. 기하를 단계마다 다시 만드는 것은 양쪽 모두의 동작입니다

C의 상태 갱신은 각 단계에서 `updateDenW2Pressure2DBlend`를 부르고, 그 안에서 Voronoi 셀을 다시 구성합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14241), [exam.c](/tmp/lageunha-push/Exam/exam.c:14356), [exam.c](/tmp/lageunha-push/Exam/exam.c:4605)). 원판 경계가 켜져 있으면 셀을 원판에 제한합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:4637)). 힘 계산기 역시 새 연결 리스트를 만든 뒤 셀 힘 루프를 수행합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:10192), [exam.c](/tmp/lageunha-push/Exam/exam.c:10222)).

Python `dyn.py`도 `acc(X)`를 각 RK 단계에서 호출하며 매번 `pg.Geom(X, ...)`를 만듭니다([dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:16), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:29)). `dynpar.py` 역시 각 호출에서 Delaunay 이웃과 원판 절단 셀을 구성합니다([dynpar.py](/tmp/lageunha-push/Docs/polybnd/dynpar.py:47), [dynpar.py](/tmp/lageunha-push/Docs/polybnd/dynpar.py:66), [dynpar.py](/tmp/lageunha-push/Docs/polybnd/dynpar.py:79)).

그러므로 단계 사이에 Delaunay/Voronoi 이웃이 바뀌는 것은 RK4 내부의 같은 근사 다각형을 유지하지 않기 때문이라는 **C 전용 오류**로 볼 수 없습니다. 다만 힘이 위상 변화에서 매끄럽지 않거나 미분 불연속이면 RK4의 고차 수렴이 낮아지고, 가장자리처럼 셀이 작고 변형이 큰 위치에서 작은 모드 교란이 생길 가능성은 있습니다. C의 셀 부피가 정확한 원판 기하와 맞는다는 기존 검증은 부피 검증입니다. 그것만으로 매 단계의 이웃 위상, 면 기울기와 압력 힘까지 Python과 동일함을 보장하지는 않습니다.

### 3. C의 힘·압력 경로는 Python의 정확한 \(-\nabla U\)와 조건부로 다릅니다

Python 기준은 \(P_a=K_a(m_a/V_a)^\gamma\)를 계산하고, 셀 부피 미분을 모아 \(F=\sum_a P_a\,\partial V_a/\partial x\)를 만듭니다([dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:16), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:17), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:18), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:20)). `onestep.py`는 면의 정확한 기하 기울기와 고정 엔트로피를 사용한다고 명시하고 있습니다([onestep.py](/tmp/lageunha-push/Docs/polybnd/onestep.py:1), [onestep.py](/tmp/lageunha-push/Docs/polybnd/onestep.py:13), [onestep.py](/tmp/lageunha-push/Docs/polybnd/onestep.py:33)).

C에는 그와 별도로 MUSCL 재구성 경로가 있습니다. `use_muscl` 설정에 따라 면 상태를 재구성하고, 힘 계산 단계에서 그 상태를 사용합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:7189), [exam.c](/tmp/lageunha-push/Exam/exam.c:7706), [exam.c](/tmp/lageunha-push/Exam/exam.c:8130), [exam.c](/tmp/lageunha-push/Exam/exam.c:8178)). 기하·밀도·압력 계산은 기본적으로 두 패스로 동작해 새 부피에서 압력·음속을 갱신한 뒤 이웃 복사본을 다시 만듭니다([exam.c](/tmp/lageunha-push/Exam/exam.c:4566), [exam.c](/tmp/lageunha-push/Exam/exam.c:5484), [exam.c](/tmp/lageunha-push/Exam/exam.c:5515)). C의 열역학 변수는 엔트로피 모드에 따라 \(K\rho^\gamma\) 또는 내부에너지와 부피로 압력을 계산하며, 비양수 압력에는 바닥값 경로가 있습니다([exam.c](/tmp/lageunha-push/Exam/exam.c:5523), [exam.c](/tmp/lageunha-push/Exam/exam.c:5536), [exam.c](/tmp/lageunha-push/Exam/exam.c:5544)).

TCF가 켜진 경우에는 추가 차이가 분명합니다. 코드상 TCF는 단일 면점, `SEDOV_PHASE1` 총에너지 경로, 특정 다른 옵션 비활성을 요구하고([exam.c](/tmp/lageunha-push/Exam/exam.c:7883)), 재구성 상태에서 \(Z\)를 구합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:7893)). 면 힘은 평균압·인공 압력 \(Q\)의 법선 성분에 접선 압력차 항을 더하고, 에너지 플럭스도 별도로 계산합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:8453), [exam.c](/tmp/lageunha-push/Exam/exam.c:8527), [exam.c](/tmp/lageunha-push/Exam/exam.c:8539)). 따라서 이 TCF 경로 전체를 Python의 정확한 보존력과 동등하다고 놓아서는 안 됩니다. 다만 기존 감사에서 “가역 TCF 힘”을 별도로 구성해 비교했다는 사실과, 현재 기본 코드의 모든 TCF 옵션이 그 비교 구성에서 활성이라는 사실은 구분해야 합니다.

### 4. `postStage`, Kepler 가속도, 내부에너지 갱신의 구체 동작

`kepler_accel`과 Hill 항은 각 RK 단계의 현재 위치(코리올리 항은 현재 속도)에서 계산해 속도 미분에 더합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:6384), [exam.c](/tmp/lageunha-push/Exam/exam.c:6394), [exam.c](/tmp/lageunha-push/Exam/exam.c:14314), [exam.c](/tmp/lageunha-push/Exam/exam.c:14317), [exam.c](/tmp/lageunha-push/Exam/exam.c:14319)). 네 단계 모두 같은 패턴입니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14385), [exam.c](/tmp/lageunha-push/Exam/exam.c:14456), [exam.c](/tmp/lageunha-push/Exam/exam.c:14527)). 따라서 **단계 위치와 Kepler 가속도 위치가 맞지 않는 문제는 소스에서 보이지 않습니다.** 다만 C와 Python의 중심 좌표·softening 불일치는 별도로 남습니다.

`phase1_ie_stage`는 일반 모드에서는 `die * dt`를 반환하지만 `SEDOV_PHASE1`에서는 \(m\,v\cdot a_\mathrm{hydro}\)를 `die`에서 빼 내부에너지 증가량을 만듭니다([exam.c](/tmp/lageunha-push/Exam/exam.c:6635)). Kepler 중력과 Hill 항은 속도에는 포함되지만 이 에너지 보정에는 직접 들어가지 않습니다. Python 기준은 내부에너지를 적분하지 않고 \(K\)를 고정하므로 이 경로를 쓰는 실행과는 상태변수가 다릅니다.

### 5. `w2`, `dtold`, 재구성·제한자·HG는 설정에 따라 남을 수 있습니다

RK4 시작에서 C는 `w2old`와 `w2backup`을 저장하고([exam.c](/tmp/lageunha-push/Exam/exam.c:14136)), 각 단계에서 저장값을 복원한 뒤 `w2lin` 또는 동적 `w2`가 켜졌을 때만 단계별 가중치를 설정합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14348), [exam.c](/tmp/lageunha-push/Exam/exam.c:14349), [exam.c](/tmp/lageunha-push/Exam/exam.c:14419)). 활성 조건은 별도 설정과 `av_mode`, `kappa` 등에 달려 있습니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14205), [exam.c](/tmp/lageunha-push/Exam/exam.c:14216)). 따라서 `w2backup` 자체를 모드의 원인으로 볼 수는 없고, 해당 실행에서 관련 모드가 실제 활성화됐는지 확인해야 합니다.

힘 계산기는 `GAS_dtold(simpar)`를 읽으며([exam.c](/tmp/lageunha-push/Exam/exam.c:7200)), 면 속도 계산에서도 `get2dUpqradRk4(..., dtold)`를 사용합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:6611), [exam.c](/tmp/lageunha-push/Exam/exam.c:7685), [exam.c](/tmp/lageunha-push/Exam/exam.c:8260)). KH는 매 스텝 뒤 실제 선택된 `dt`를 `GAS_dtold`에 저장합니다([kh.c](/tmp/lageunha-push/Exam/KH/kh.c:327), [kh.c](/tmp/lageunha-push/Exam/KH/kh.c:330)). 일부 면 재구성·제한 경로가 활성화된 실행에서는 이력이 Python의 고정 `dt`와 달라질 수 있습니다.

HG 분할 적분은 래퍼가 관련 옵션을 감지했을 때 `HG(dt/2) → RK4(dt) → HG(dt/2)` 순으로 실행하며([exam.c](/tmp/lageunha-push/Exam/exam.c:14021), [exam.c](/tmp/lageunha-push/Exam/exam.c:14044), [exam.c](/tmp/lageunha-push/Exam/exam.c:14055), [exam.c](/tmp/lageunha-push/Exam/exam.c:14062)), `hg_apply` 뒤에는 내부에너지에서 압력과 음속도 갱신합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:13983), [exam.c](/tmp/lageunha-push/Exam/exam.c:13987)). HG가 꺼져 있으면 이 분할 경로는 원인 후보에서 제외해야 하지만, 두 패스 압력 갱신과 일반 단계의 재구성은 별개 경로입니다.

### 6. 시간 간격: C는 기본적으로 적응, Python은 고정

C는 면 음속·상대속도에 따른 CFL, 추가 면 기하 제한, 속도차 제한을 계산합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:8755), [exam.c](/tmp/lageunha-push/Exam/exam.c:8757), [exam.c](/tmp/lageunha-push/Exam/exam.c:8765), [exam.c](/tmp/lageunha-push/Exam/exam.c:8785)). Kepler 제한은 `EUNHA_KEPLER_DTETA > 0`일 때만 최소 간격에 반영됩니다([exam.c](/tmp/lageunha-push/Exam/exam.c:6541), [exam.c](/tmp/lageunha-push/Exam/exam.c:6544), [exam.c](/tmp/lageunha-push/Exam/exam.c:6562)). 고정 간격 환경변수도 안정성 제한을 넘을 수 없게 검사합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:14291), [exam.c](/tmp/lageunha-push/Exam/exam.c:14298)). Python의 기준 실행은 단일 `dt`를 모든 단계에 씁니다([dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:5), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:27), [dyn.py](/tmp/lageunha-push/Docs/polybnd/dyn.py:30)).

래퍼의 재시도·간격 축소도 기본 RK4와 별도 조건부 기능입니다. HG 분할이나 `GFS_STEP_RETRY`가 켜진 경우에만 실패 상태를 복원하고 간격을 줄여 재계산합니다([exam.c](/tmp/lageunha-push/Exam/exam.c:13994), [exam.c](/tmp/lageunha-push/Exam/exam.c:14021), [exam.c](/tmp/lageunha-push/Exam/exam.c:14082)). 따라서 “C가 적응 간격을 쓴다”는 사실과 “매 단계 제한 위반을 자동 재시도한다”는 사실은 같지 않습니다.

## 가장자리 모드 원인 후보 순위

순위는 소스상 가능한 영향과 현재 기록의 실험 결과를 함께 반영한 **검증 우선순위**입니다. 어떤 환경변수가 실제 실행에서 설정됐는지 확인하지 못한 항목은 조건부입니다.

1. **Kepler 외력 불일치: 중심 또는 softening 설정.** C의 기본 `(2,2), eps=0`과 Python의 원점·`eps=0.01`은 다릅니다. 실제 실행에서 환경변수로 일치시켰는지 확인되지 않으면 가장 먼저 배제해야 합니다.
2. **힘 및 열역학 연산자 차이.** Python은 고정 엔트로피와 기하 에너지 기울기를 사용하지만, C는 활성 설정에 따라 재구성·제한·TCF 면 플럭스·내부에너지 진화를 씁니다. 원판의 여러 셀 폭을 가로질러 생기는 모드를 씨앗하거나 증폭률을 바꿀 수 있습니다.
3. **가장자리의 단계별 위상 flip 및 기하 수치 차이.** 두 코드 모두 단계마다 기하를 재구성하므로 flip 자체는 공통입니다. 그러나 flip 근처에서는 힘의 매끄러움과 수치 구현 차이가 작은 교란의 위치와 크기를 바꿀 수 있습니다. 부피 일치만으로 이 항목은 닫히지 않습니다.
4. **적응 시간 간격과 이력 의존 경로.** 가변 `dt`와 `dtold` 기반 면 속도·재구성은 Python 고정 `dt`와 다릅니다. 공간적 증폭률이 달라지는지 같은 간격 실험으로 분리해야 합니다.
5. **`w2`, 제한자, HG 압력 재계산.** 관련 플래그가 비활성이면 제외할 수 있습니다. Q/HG를 끈 실험에서 성장이 그대로였다는 기록은 점성 경로의 우선순위를 낮추지만, 그것만으로 MUSCL·힘·외력·간격 차이를 배제하지는 않습니다.

## 한 단계 오차와 장기 성장의 판정

`onestep.py`가 보고한 속도 차이 \(3\text{–}7\times10^{-10}\)은 단일 단계에서 매우 작고, Python 속도 증가분의 약 \(10^{-4}\)라는 기존 기록과 부합합니다. 이 정도 차이도 불안정 모드의 초기 씨앗이 되어 이후 커질 수는 있습니다. 그러나 **초기 씨앗 차이만으로 서로 다른 지수 성장률을 설명할 수는 없습니다.** 같은 선형 성장률의 모드라면 씨앗 차이는 주로 모드가 눈에 띄는 시점과 진폭의 초기 오프셋을 바꿉니다.

요청 문서의 수치를 그대로 보면 C의 \(m=18\)은 \(t=16\)의 0.002에서 \(t=28\)의 0.038로 증가해 구간 평균 성장률이 약 \(0.25/\text{시간}\)입니다. Python의 \(m=18\)은 \(t=22\)의 0.0002에서 \(t=42\)의 0.0012로 약 \(0.09/\text{시간}\)입니다. 이 비교는 시작 시점·실행 설정이 완전히 일치했다는 증거가 아니므로 확정적인 성장률 대조는 아닙니다. 그래도 그대로라면 C/Python의 차이는 단순한 시작 씨앗만이 아니라 **힘의 선형화된 응답, 외력·상태 설정, 또는 기하/시간 적분 차이**를 포함할 가능성이 있습니다.

## 원인 판별 최소 시험 순서

1. **실행 계약을 먼저 맞춥니다.** C와 Python의 Kepler 중심, `eps`, 질량·초기 위치·속도, \(K\), 원판 경계, `use_muscl`, 엔트로피 모드, `av_mode`, `kappa`, TCF·Q·HG·QBP·χ 및 `GFS_FIXED_DT`를 기록해 비교합니다. 특히 C의 실제 `EUNHA_KEPLER_CX/CY/EPS` 값과 Python의 원점 기준 좌표를 확인합니다.

2. **동일 스냅숏에서 네 단계별 RHS를 대조합니다.** 최종 스텝 차이뿐 아니라 각 \(k_i\)에서 위치·속도·부피·압력·밀도·면 연결·힘·dt를 비교합니다. 차이가 Kepler 가속도에서 시작하는지, 기하·힘 평가에서 시작하는지, `postStage` 뒤 처음 생기는지 분리합니다. 기존 `onestep.py`는 3스텝 최종 비교가 중심이므로 이 항목을 단계별 로그로 확장해야 합니다([onestep.py](/tmp/lageunha-push/Docs/polybnd/onestep.py:44), [onestep.py](/tmp/lageunha-push/Docs/polybnd/onestep.py:47)).

3. **같은 상태에서 고정 간격 가역 실행을 양쪽으로 시작합니다.** 문서가 지정한 \(dt=0.02\)를 두 코드에서 똑같이 쓰고, Python이 사용하는 중심과 softening을 C에도 명시합니다. Q/HG·추가 제한·`w2`를 비활성화한 동일 조건에서 모드 진폭의 \(\log A_m(t)\) 기울기를 비교합니다. 이 기울기가 같고 최초 진폭만 다르면 씨앗 차이가 주 설명입니다. 기울기 자체가 다르면 RHS 또는 시간 적분 구현 차이입니다.

4. **동일 상태 쌍을 교차 시작합니다.** C의 \(t=20\) 상태를 C와 Python 양쪽에서, Python 상태도 양쪽에서 시작해 성장 기울기를 비교합니다. 어느 코드에서 시작했는지가 아니라 어느 적분기가 성장률을 결정하는지 판별할 수 있습니다.

5. **위상 flip과 가변 `dt`를 마지막에 분리합니다.** 먼저 각 단계의 이웃 목록·flip 시점과 해당 셀의 힘 점프를 대조합니다. 그 후 두 코드 모두 가능한 동일 고정 `dt` 대 적응 `dt` 비교를 수행합니다. 최종적으로 하나의 옵션군씩 원래 설정에 추가해 어떤 항목이 성장률을 바꾸는지 확인합니다.

이 감사는 요청하신 소스만 읽은 정적 검토입니다. 파일은 수정하지 않았고, 원인을 확정하는 실행 시험도 수행하지 않았습니다.