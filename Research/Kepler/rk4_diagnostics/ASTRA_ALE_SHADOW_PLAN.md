# 다음 묶음: actual-face conservative ALE shadow

2026-10-20. SSP fresh 최종 정상 상태까지 한 묶음으로 Astra 검토를 받았다. SSP는 stage 양성/rollback 도구로 유지하며 현재 공간 연산자의 장시간 해법으로 채택하지 않는다. t=.135091 예산 종료를 t=.2 완료로 보고하지 않는다.

다음 구현은 실제 native GPU 공유면에서 하나의 보존형 Euler Riemann 문제로 mass/momentum/total-energy flux를 함께 산출하는 읽기 전용 shadow다. W0 고정질량 본체에는 이 flux를 적용하지 않는다. 새 IC 진화, threshold/dt scan, 추가 장시간 재생은 이 묶음에 없다.

f_m=rho(v_n-w_n), f_p=rho*v(v_n-w_n)+P*n, f_E=Eden(v_n-w_n)+P*v_n. 같은 면 상태와 같은 실제 기하 속도를 사용한다. 열률은 Udot=Edot-v·pdot+.5*v²*mdot으로 환산한다. 기존 pressure-only와의 차이를 relative advection/pressure work/pressure closure 차이로 나누고, 근사 Riemann flux와 그 중간 상태의 물리 flux가 다른 잔차도 숨기지 않는다.

구현은 canonical lower-ID endpoint 소유 rank가 각 pair를 한 번만 CUDA 계산하고 하나의 패킷을 쓴다. 분석에서 이 packet의 정확한 ±값으로 두 셀의 장부를 합산한다. 실제 상태/질량/W/누계/출력시계는 갱신하지 않는다. 무진화 RHS-only 모드는 current-stage geometry->EOS->halo->gradient/force까지만 준비한 뒤 모든 rank가 종료한다.

제한된 자료는 균일 상태, 왜곡 mesh의 rigid/trace-free affine patch, 보존된 fine 실패 직전 정상 상태 하나, SSP fresh 최종 정상 상태 하나다. 균일/affine 적분 항등식, boost 변환, canonical pair 유일성/반대 flux 및 기존 진화 observer 일치를 확인한다. 제조된 정확 affine 재구성과 실제 native limiter 재구성은 별도 기록한다. 관련 논문 https://arxiv.org/html/2205.07916 은 완전 ALE에 대한 것으로, pressure-only에 Gauss 점만 추가해 같은 보장을 주장하지 않는다. 논문 limiter의 기준은 generator midpoint가 아닌 cell mass-centroid midpoint다.

균일밀도 trace-free affine 유동에서 mdot=rho*Vdot_geom을 검증하며 비영 mdot을 제거해서 고정질량 ALE라고 부르지 않는다. Actual 기하 부피의 GCL과 개별 고정질량을 동시에 닫을 수 없는 항이 확인되면 진화부터 시작하지 않고 다음 공간 설계로 넘긴다. 과거 cs/h sliding instability 및 material-volume weight-DAE 실패는 반복하지 않는다.

안쪽 graded 초기간격 h=eta*sqrt(r²+eps²), 완만한 외곽 cap은 후속 해상도 설계로 보관한다. 현재 새 graded 진화는 추가하지 않는다. 최종 목표는 t6-10과 GIZMO 최장 약1000공전으로 유지한다.
