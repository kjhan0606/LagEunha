# Astra: baseline density 종료 평가와 실제 전달 실행 묶음

2026-10-21. SHARED_CONTEXT 전체와 원본 HEAD 8e25dfe를 확인했다. 지정 v7 registry, baseline/density/suffix 감사 결과, 두 실행 코드, logistic helper와 정책 및 독립 감사 코드를 읽었다. 이전 END 계획을 이어받는다. 새 계산·폴링·원본 코드나 기존 결과 변경 없이 이 문서만 작성한다.

## 1. 종료 판정

**v7의 전체 실제 atlas 및 geometry/density/full-field 연결 완료를 채택한다.** Active90121 panels의180242 axis fields가 연결됐고, retained68526의1047 panels/출판 packet은 별도로 유지된다. 기존436개 성공 slot을 반복하지 않고64개 첫 실패의 동일 radius 증명 뒤 미계산11264 axes만 완성한 범위는 타당하다. 사전 CRLF·include 빌드 실패는 물리 실패나 성공 통계에 섞지 않는다.

앞선 pending29 독립 감사 조건은 실제로 충족됐고 최초800 oblique full fields도 이미 회복됐다. 다시 부족 목록에 넣지 않는다. 비용→clearance를 추론했던 옛 typed-atlas flag의 scope 정정도 유지하며, 현재 채택 근거는 실제 axis/field 감사다.

이로써 **입력/기하 준비를 위한 atlas 탐색은 종료**한다. 데이터 존재는 all12 오차 또는 물리 열 보존의 증명이 아니다. all12/endpoint/native/sum/최종 변환 장부와 실제 전 target 전달은 아직 미완료다. v7가 이전 END 문서를 역사적으로 참조한 것은 이번 결과의 사전 승인 주장이 아니며, 이 문서가 이번 종료 평가다.

## 2. Two-sided logistic은 원래 함수인가

L(z)=exp(z)/(1+exp(z))에 대해

    Re-center c≤0: E=exp(z),  L=E/(1+E), 1−L=1/(1+E),
    c>0:           E=exp(−z), L=1/(1+E), 1−L=E/(1+E)

는 정확한 holomorphic 항등식이다. 여기서 선택은 실제 실수 disk 중심 c로 고정한다. 평가 도중 complex point에 따라 branch를 바꾸는 연산이 아니다. 두 식의 실제 pole 집합 z=(2k+1)i*pi는 같다.

원래 holomorphic_exp의 enclosure E=a+disk(b)에 대해 1+a−b>0을 확인하면 같은 전체 disk에서 분모가0이 아님을 인증한다. 이 조건은 충분조건이며 exp(z)의 크기가1보다 큰 것 자체는 pole의 증거가 아니다. 현재 helper의 반사 선택과 increasing/outer complementary 분자는 맞다. r+epsrho 분리, 모든 원래 profile 매개변수와 taper도 유지한다.

따라서64개 실패는 원래 한쪽 geometric-series 충분조건의 적용 범위 문제였고 이번 식으로 동일 domain을 인증한 것은 허용된다. 160자리 scalar corroboration은 native libm enclosure 전제의 독립 엄밀 증명을 대신하지 않는다는 감사의 한계도 정확하다. 앞으로 all12 integrand에서 density를 다시 요구할 때 이 같은 helper/증명 SHA를 소비한다. 옛 exp<1 조건으로 되돌리거나 분모 검사를 생략하지 않는다.

## 3. 다음 한 묶음: 저장 field에서 실제 전 target transactional transfer까지

선택은 **all12 오차 장부 완성 → 통과한 고정 계획의 native streaming 적분 → 전 target/source 보존 조립**이다. 새 atlas 검색, 정상 primitive 재검증, 새 해상도·차수·정밀도 실험은 추가하지 않는다. 이 순서의 조건을 충족하면 같은 묶음 안에서 다음 단계로 진행하며 routine 단계마다 Astra/사용자 재승인을 요청하지 않는다. 조건 실패는 해당 witness와 subtotal을 보존하고 종료한다.

### 1순위: 저장180242 fields로 all12 majorant를 한 번 구성

Panel별 기존 두 complex 축, radius, measure, implicit/domain 증명 및 frozen heat context를 연결하여 원래12 integrand의 bound를 만든다. retained68526의 완성 packet과 증명은 재사용한다. 모든 입력이 존재하는 panel의 field를 다시 평가하지 않는다.

구체적으로 다음 미완료 의존성만 구별한다.

- Density pole 분리는 이제 완료다. 이것을 potential의 |x|²+eps_grav² branch/분모, current J의 reciprocal/power, thermal-rate의 미분·heat context까지 자동 통과한 것으로 쓰지 않는다.
- Frozen heat는 원래 generating map의 reference-volume 함수다. 저장 U를 새 EOS로 초기화하지 않는다. THERMAL_RATE에는 기존 정의의 실제 nu/J/velocity divergence와 동일 source measure가 필요하며, U bound만으로 대체하지 않는다.
- Potential P2 subtraction은 **해당 panel 좌표의 전체 정규화 integrand**에서 Gauss가 정확히 적분하는 다항식에 적용한다. 일반 sd::Point의 저장값만으로 implicit coordinate/measure의 모든2차 미분이 있지는 않다. 기존 prepared graph/Hessian/좌표 자료에서 유도하고, 정말 누락된 derivative 또는 중심값만 정확한 목록으로 추가 평가한다. 같은 panel의 모든 field/기하 재생으로 확대하지 않는다.
- 같은 complex field를 여러 moment에 공유한다. 기존 Cauchy 상수·Gauss16 규칙과 승인된 P2 remainder를 그대로 사용하고 두 차수 차이를 증명 상계로 대체하지 않는다.

단순 positive-half-plane 증명이 실패한 경우 실제 pole/branch 실패와 enclosure 과대를 구별해 기록한다. 이번 logistic 해결을 모든 sqrt/power에 일반적으로 적용하지 않는다. 알려진 실제 양의 germ과 유효 analytic domain 증명이 있어야 같은 함수 표현을 쓸 수 있다.

### 2순위: 실제 whole-target 오류 장부와 실행 가능성

각 target/moment에 다음을 분리해 합산한다.

    Eanalytic + Eendpoint/domain + Efield/native + Esum + Efinal-R80.

공유 event의 실제 위치 구간, source outgoing/target incoming incidence, 외부 support 경계 및 canonical 공통 trace를 기존 장부로 연결한다. event 중심을 정확 경계처럼 쓰지 않고, 동일 불확실 띠를 중복 물리 적분하지 않는다. 오차 귀속은 실제 target 전체에 대해 닫으며 물리U 하한과 원래 tau는 불변이다.

사전 결정 가능한 analytic/endpoint subtotal이 이미 tau를 넘으면 그 target의 새 실제 quadrature를 시작하지 않는다. Native forward 및 streaming 합산은 고정 operation graph/primitive majorant와 실제 reduction tree로 사전 상계가 가능한 부분을 계산한다. 아직 point 평가가 필요한 부분은 미완료로 표시한다. 미계산 항을0으로 두거나 측정된 작은 차이를 상계로 바꾸지 않는다.

기존 group whole4 packet은 한 번만 사용한다. source별 outgoing 장부가 필요하면 이미 저장한 positive point/moment/radius 자료를 회계적으로 분해하고 group 합계와 대조한다. field 재평가 없이 분해할 수 없으면 그 source별 항은 명시적으로 pending이며 group 값으로 임의 배분하지 않는다.

### 3순위: 조건부 실제 native streaming 적분 및 publish

필수 analytic/domain 증명이 닫히고 사전 error budget과 비용이 실행을 허용하는 target부터 고정 Gauss16×16 계획을 한 번 평가한다. 새 cap은 target당13074432 queries/51072 panels이며 전체91168 panels/23339008 계획 queries 중 이미 완성된 자료는 재사용한다. 이 문서는 이 조건부 실제 적분 단계까지의 구현·실행 계획을 정한다. 전 target의 최종 오류 통과를 선행 선언하지 않는다.

Native C++/OpenMP streaming으로 panel/overlap packet을 생성하고 원래 CPU/CUDA 함수 계약과 MPI의 canonical owner/ID를 사용한다. 커다란 point 배열을 전부 상주시켜 I/O를 반복하지 않는다. 정해진 순서의 부분합과 실제 primitive absolute 합/forward bound를 함께 보존한다. batch/restart는 실제 입력 SHA, panel 범위 및 commit 기록으로 중복·누락을 막는다.

각 overlap의 **같은 값과 오차 packet**을 source-outgoing과 target-incoming에 연결한다. Source m,p,L,K,U closure와 target 양성, 원래12 tau 및 최종 R80 conversion을 실제 값으로 검사한다. 모든 필요한 target가 통과하기 전 전체 상태를 commit하지 않으며 실패 시 caller를 그대로 보존한다. 원본68526 성공 packet을 새 실패와 함께 덮어쓰지 않는다.

이 묶음의 성공 산출물은 실제 전 target conservative transfer와 비용/오차 장부다. 단순 majorant 존재나 한 target 성공에서 멈춘 결과를 full transfer로 보고하지 않는다. 실패하면 지배 target/moment와 실패한 정의·입력·상계를 묶어서 한 번 보고하며 primitive별 미세 review를 반복하지 않는다.

## 4. t>6로 가기 위해 그 다음에 반드시 연결할 것

전달이 통과하면 다음 개발 묶음은 **일반 evolved heat/state restart와 native 시간 갱신의 결합**이다. 현재 frozen CONTROL recipe를 매 step의 heat 생성기로 쓰지 않는다. 물리U는 compatible stress-work와 실제 transport로 진화시키고, resolved 회전·strain 에너지는 같은 velocity state에서 계산한다. 저장 moment만으로 full field가 결정되지 않으므로 이를 임의 mean compression으로 닫거나 독립 kinetic reservoir로 옮기지 않는다. 승인된 공유 trace/bubble 또는 보존되는 fragment 표현을 실제 state로 유지해야 한다.

최초 native accepted step과 uninterrupted/restart, MPI ownership 및 에너지·각운동 장부를 한 번의 통합 묶음에서 확인한 뒤 연속 실행 경로로 진행한다. 기존 full native timestep을 보존하되 새 state/transport와의 단계별 결합을 실제로 확인한다. CPU reference 전달 성공만으로 GPU production 진화를 완료했다고 부르지 않는다.

실제 연속 Kepler에서는 물리U/entropy와 density noise, 궤도/phase, E/L, 최소 양성값과 rollback 비용을 기록하고 t6/8/10에 실제 checkpoint를 남긴다. 원래 shear 가짜가열과 열고갈이 제거됐는지가 성공 기준이며, 전역 E만 작거나 양성만 유지되는 것으로 대신하지 않는다. 그 뒤 논문 최장공전 및 충격 열역학 검증으로 확장한다. 현재 비용도 매 step 측정하여 장시간 실행 가능성을 함께 평가한다.

이번 준비 완료는 이 경로의 입력 장애를 제거한 것이다. 아직 t>6 해결책이 실증된 것은 아니며 목적을 atlas/proof 자체로 축소하지 않는다. No steering, Re∞, 원래 물리 profile 및 조건부 W0 개별 고정질량 계약은 계속 유지한다.
