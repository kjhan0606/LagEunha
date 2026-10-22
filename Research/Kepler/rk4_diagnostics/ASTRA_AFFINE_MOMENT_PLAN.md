# 다음 구현: affine moments + weight와 연결한 질량 완화

**Historical/superseded:** native 초기화로 최대10.02% cold thermal bias를
확인했다. 초기화/밀도/MPI 상태 관문은 구현됐지만 이 closure의 진화는 하지
않는다. 현재 결정은 `ASTRA_EULER_SPATIAL_PLAN.md`, 결과는
`MOMENT_INITIAL_RESULTS.md`이다. 독립 kinetic tensor 진화 제안도 Euler shock
열역학 항등식 검토 후 철회됐다.

2026-10-20. 완료된 native ALE shadow 및 실제 두 상태의 work/moment 자료를
Astra에게 한 번 묶어 평가받았다. Driver가 코드와 실행을 담당한다.

## 판정과 범위

Shadow의 mass/momentum/energy 수송 관문은 통과했다. Full ALE를 기존 점속도와
thermal U에 바로 적용하는 방식은 채택하지 않는다. Affine 추정은 후보를
구현할 근거이며 현재 실제 누락 에너지/원인 확정이 아니다. Scalar spin만
추가하거나 torque power를 기존 Udot에서 빼지 않는다.

현재 native 본체에 opt-in affine-moment ALE와 weight–mass 완화 경로를 하나의
수송 계약으로 구현한다. 기존 W0·개별 고정질량 경로는 보존한다. No steering,
Re=infinity, heat/mass/weight clipping 금지를 유지한다. 실패한 hard fixed-mass
weight DAE를 재활성하는 구현으로 대체하지 않는다.

## 상태와 열역학

각 셀에 독립적인 적분 상태를 둔다:

    m, q=int rho*x, Q=int rho*x*x^T,
    p=int rho*v, D=int rho*x*v^T, E.

    c=q/m, S=Q-q*q^T/m, C=D-q*p^T/m,
    vbar=p/m, B=C^T S^-1,
    K=|p|^2/(2m)+.5 tr(C^T S^-1 C), U=E-K.

D convention은 위치 행/속도 열이다. 원점 기준 Lz=Dxy-Dyx, 내부 Lz=Cxy-Cyx.
Density는 m,q,Q를 만족하는 양의 재구성이어야 한다. 첫 thermal closure는
specific u=U/m를 셀 안에서 일정하게 두고 P=(gamma-1)*rho*u로 명시한다.
이는 일반적인 subcell thermal profile의 정확한 표현이 아니다. Generator,
geometric centroid, density-weighted c 및 vbar를 혼동하지 않는다.

## 같은 flux와 source에서 moment 진화

    mdot = -int_face fm
    qdot = -int_face x fm + p
    Qdot = -int_face xx^T fm + D+D^T
    pdot = -int_face fp + int_cell rho g
    Ddot = -int_face x fp^T + int_cell(rho vv^T+P I)
           +int_cell x rho g^T
    Edot = -int_face fE + int_cell rho v dot g

Canonical 공유면의 같은 quadrature point/flux로 모두 계산한다. Angular flux를
face-midpoint torque로 대체하지 않는다. Gravity 힘·torque·일은 같은 volume
재구성으로 적분한다.

## Gauge-fixed weight–mass 관계

    m_i + mu_i * Wtilde_i = m_i0,
    mu_i = m_i0/V_i0 > 0,
    mdot_i = -mu_i * zdot_i, zdot_i = dWtilde_i/dt.

W는 2D length^2 단위다. 각 stage의 density/velocity 재구성을 고정한 내부 solve에서
실제 ALE flux로 얻은 질량률 F_i(z)에 대해

    F_i(z)+mu_i*z_i=0

을 푼다. Mesh-speed에 대한 massflux의 단조성/Jacobian을 직접 검사하고, 그 전에
SPD를 가정하지 않는다. mu는 첫 후보에서 고정하며 parameter scan을 하지 않는다.
Wtilde=0이면 m=m0/Voronoi이며 새 동적 weight 모드에서 질량이 함께 변한다.
W를0으로 고정한 채 비영 ALE massflux를 버리는 fallback은 금지한다.
세속 weight 성장/기하 실현성이 해결됐다는 보장은 없다.

## 구현 및 검증 순서

1. Native moment state/flux/source 및 gauge-fixed weight–mass solve.
2. Exact affine patch의 moment 적분, boost/원점 이동, angular/E 장부.
3. 양의 m,S,U, density-moment 실현 가능성, mass–weight 관계와 실제 mesh.
4. 기존 연속 Kepler 초기장을 적분한 fresh moment IC. 고갈 checkpoint에 추정 K를
   추가하지 않는다.
5. 초기 RHS에서 thermal/moment 실현가능 경계의 방향 확인 후 진화.

SSPRK rollback을 재사용하되 moment 실현가능성도 모든 stage에 검사한다. Cold
subtraction 정밀도를 별도 기록한다. 위 관문이 통과하면 기존 refined71080배치로
한 fresh 경로 t=.2만 먼저 수행한다. 새 graded배치/limiter/threshold를 동시에
바꾸지 않는다. 단기 성공으로 최종 목표를 축소하지 않는다.

관계/보존 실패, heat/mass/weight clipping 필요, 기하 상실, 이전과 같은 weight
세속증가, dt붕괴/과도한rollback, 중심 열 고갈이면 후보를 중단한다. 통과한
후보만 t6-10 및 최종 GIZMO 최장 약1000공전으로 연장할 비용을 산정한다.

## Driver implementation notes: 아직 구현되지 않은 항목

- Raw x^2 moment의 cancellation을 줄이기 위해 generator 근처 local origin에서
  q,Q,D를 저장하고 원점 이동 항을 정확히 변환하는 방식을 우선 검토한다.
- m,q,Q 양성 제약은 단순 PSD S보다 강하다. Positive density reconstruction과
  움직이는 polygon support의 실현가능성을 함께 구현해야 한다. Uniform-density
  fit을 강제로 덮어써서 q,Q를 바꾸지 않는다.
- Native weighted geometry에는 기존 wfrac clamp/float W 저장/음수 weight gauge 등
  이력이 있다. 새 모드의 실제 power diagram과 weight–mass 관계를 깨는 clamp를
  숨겨 사용하지 않는다. 기존 legacy 경로와 분리해 확인해야 한다.
- Generator advection 법칙은 reconstructed physical fluid velocity와 연결하여
  명시하고, centroid 방향의 임의 steering을 도입하지 않는다.
- MPI migration/halo와 restart는 full moment/weight 상태를 함께 보존해야 한다.
  Python 진화나 별개 toy integrator로 native 구현을 대신하지 않는다.

현 상태: 위 모델의 초기화와 상태/밀도/MPI 관문은 구현·평가됐다. 이 evolving
thermal closure는 보류·대체됐으며 실제 진화를 하지 않았다. 현재 계획은
`ASTRA_EULER_SPATIAL_PLAN.md`이다. 기존 완료 자료는 `ALE_SHADOW_RESULTS.md`와
`MOMENT_INITIAL_RESULTS.md`에 있다. 목표는 계속 active다.
