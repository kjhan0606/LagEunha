# 다음 구현 후보: 단계 양성을 검사하는 native SSPRK3

2026-10-20. 완성된 압력비 분기·경계 판별을 Astra에게 한 번 묶어 평가받았다. 운전자가 코드와 실행을 담당한다.

## 선택과 범위

기존 공간 pressure/force/work 식을 유지하며 opt-in SSPRK(3,3)와 전체 시도 복원을 구현한다. 가열·냉각 분기를 이동/평활화하거나 U를 보충하지 않는다. 경계가 안쪽을 향한 네 셀의 증거를 전체 상태에 대한 양성 정리로 확대하지 않는다. 불연속 전환면의 수치 chatter 및 인공 열평형 가능성도 남는다. 이 후보는 누적 열 고갈이나 과거 sliding-cell 불안정의 장시간 해법으로 채택된 것이 아니다.

Y1=Yn+hL(Yn)

Y2=3/4 Yn+1/4(Y1+hL(Y1))

Ynew=1/3 Yn+2/3(Y2+hL(Y2)).

실제 내부에너지는 기존 Udot=die-m v·a_hydro를 사용한다. 중력은 v에 넣으며 Udot에서 중력의 평동 일을 다시 빼지 않는다. stage 시간은t,t+h,t+h/2. 주기 좌표는 원래 위치와 unwrapped 이동량으로 조합하고 질량/W0은 혼합하지 않는다.

모든 FE 블록에서 h<=min(physicalCFL,0.9 min[U/(-Udot)])을 검사한다. 첫 단계만 검사하지 않는다. FE endpoint/혼합 상태의 열·유한성과 실제 기하를 별도 검사한다. 다음 시도 ceiling=min(h/2,실패stage의bound). 그 bound는 다음 시도의 보장값이 아니므로 다시 검사한다. 한 step당 거절최대8, 두 구간의 프로세스당 RHS최대256/accepted최대64. 한도·시간전진 소실이면 정상 상태를 보존하고 종료한다.

## 복원 계약

rank별 원래 full payload/np/전역개수,dtold,단계 floor/pair 장부를 보존한다. 실패 시 원래 rank의 데이터와 소유를 복원하고 halo/모든 GPU geometry·gradient·CSR seed를 폐기한다. 전역 particle payload allgather는 하지 않는다. domain pivot은기존pmigrate에서읽기전용임을 소스로 확인한다. 시간/출력시계 및 물리누계는commit 뒤에만 전진한다. 진단용시도/RHS횟수는복원하지않고 예산에 포함한다. 모든 rank가 같은시도를채택/거절한다.

전체에너지의 정확보존은SSPRK3의볼록구조만으로보장되지않는다. 반이산hydrowork와finite-step잔차,중력일,질량/선형운동량을따로기록한다.

## 실행 묶음과 종료

1. 코드 기본경로/RK4 보존, 주기좌표·MPI이동 뒤 강제거절/복원 비교.
2. 보존된 fine-cold 정상시작점과 Sedov 정상시작점의 같은 실패구간. 각구간에기본ceiling과그절반1회만비교.
3. 양성·기하·보존·시간민감도가허용되면 fresh fine-old IC→t.2 한경로. 고갈된checkpoint를통과한것으로열역학회복을판정하지않는다.

단계floor통과·예산초과·h붕괴·분기주기에민감한endpoint·에너지관문실패면중단한다. fresh에서고갈/궤도오류가남으면시간안전장치로만보존하고장시간해법후보로기각한다. 다음공간모델은고정질량/면운동의에너지·엔트로피수송계약과sliding/shear응답을함께검증한다. t6–10,약1000공전목표는그대로다.

## 구현 중 드러난 입력 준비 차이

GPU gradient는업로드된rho/P를사용하고그뒤EOS kernel이새V의P/cs를갱신한다. 기존RK4의후속stage는직전stage cache를gradient에쓸수있다. SSP후보는현재좌표의rawgeometry→EOS→halo→기존gradient/force순으로준비한다. Riemann·limiter·force/work공식은유지하지만,이준비순서차이를포함하므로 **엄밀한동일RHS 시간적분법만의대조라고주장하지않는다**. 기본RK4는별도로보존검증한다.
