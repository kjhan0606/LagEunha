# 모멘트 생산 경로 연결 계약

계획·평가 담당 Astra의 지시에 따른 첫 상태 계약이다. 기존 particle struct, wm17 파일 바이트, 생산 기본 solver는 변경하지 않는다. 이 문서는 전체 생산 모멘트 solver의 구현 완료를 뜻하지 않는다.

## 상태와 단위

연산 layout은 Python raw18과 같은 m,q[2],Q[4 row-major],p[2],D[4 row-major],Pi[4 row-major],Ehyd다. Dxy와 Dyx는 독립이며 전체 각운동량은 Dxy−Dyx다. 내부 각운동량에는 q×p/m를 빼야 한다. U=Ehyd−tr(Pi)/2로 회전·잔류 운동을 열에 중복 계상하지 않는다.

Psi는 raw의 추가 성분으로 숨기지 않는다. 별도 scalar와 absent / diagnostic_reconstruction / transported_moment 의미를 둔다. 진단 적분값이 수송 모멘트로 자동 승격되지 않는다. generator x[2], additive power W(길이 제곱), covariance factor L[4], stable global ID, frozen flag, stage/time/geometry epoch도 raw와 분리한다. 기존 생산 weight의 의미는 실제 기하 수식을 확인한 뒤 명시적으로 연결한다. 이름이나 주석만으로 제곱/제곱근을 취하거나 음의 additive W를 거부하지 않는다.

moment_state_layout.h의 checked adapter는 finite/layout만 검사하며 signed rate에도 사용한다. Q/Pi 대칭 성분은 정확한 equality를 요구하고 평균·clipping하지 않는다. D는 대칭화하지 않는다. wm17에는 Psi가 필요하므로 absent 입력을 0으로 채우지 않고 실패한다. 검사를 모두 끝내기 전에는 출력 배열·Psi metadata를 변경하지 않는다.

## 공유면과 ghost

face key는 min(stableID),max(stableID),periodic image/subface ID다. owner는 minID 셀의 소유 rank다. 한 owner가 공통 좌표·법선·면속도·적분점과 raw18 유량 및 사용 시 Psi 유량을 계산하고 두 셀에 반대 부호로 전달한다. 이웃 rank가 독립적으로 계산한 서로 다른 반올림 유량을 각자 적용하지 않는다.

ghost에는 raw18,L,x,W,ID,frozen,stage/time/geometry epoch가 필요하다. halo는 밀도 2-hop, quadratic velocity의 현재 반경 3 sqrt(area), 기하 구성의 필요 범위의 합집합을 충족한다. rank 분할에 따라 local prior stencil ID가 달라지면 안 된다. active massroot는 전역적으로 같은 제약 문제이며 MPI 경계를 frozen Dirichlet 경계로 취급하지 않는다.

## 재시작

C struct 메모리 덤프를 serialization으로 사용하지 않는다. schema version, field order, dtype/endianness, ID,raw18,L,x,W,frozen,Psi mode/value(존재 시), closure/version, integrator,time,dt, weightseed, boundary/gravity/energy ledgers를 명시한다. seed는 재수렴용 초기 추정이며 물리 상태의 대체물이 아니다. 기하/prior cache는 재생성하고 stage epoch를 확인한다. point checkpoint에서 모멘트 checkpoint로의 암묵적 승격은 허용하지 않는다.

## 첫 검증 범위

동일 adapter harness로 실제 캡처 raw18와 wm_atom fixture를 사용해 정확 roundtrip, 비대칭 D 보존, Q/Pi 비대칭·absent Psi·NaN 거부와 실패 시 출력 불변을 CPU/OpenMP/CUDA에서 비교한다. 이 단계는 adapter 검사이며 MPI 통신이나 전체 물리 kernel의 정합성을 증명하지 않는다. 후속 통합은 Astra 평가 후 진행한다.
