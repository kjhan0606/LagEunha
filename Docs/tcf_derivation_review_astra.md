**Astra의 판정은 “변분형 TCF의 반이산 원리는 승인, 현재 문서·계획에 대한 G1 서명은 보류”입니다.** Fable이 정정한 기하 미분, 접선 힘 부호, Laguerre 평균압력은 맞습니다. 그러나 원호 힘의 계산 방식이 모호하고, §8의 VPH 오차 진단은 틀렸으며, \(Q_f\)의 충격 포착·소산 수렴 주장은 아직 증명되지 않았습니다. 무엇보다 계획서의 구현식이 유도서의 변분식과 다릅니다.

검토 대상은 `/tmp/lageunha-push`의 LagEunha 저장소입니다. 파일 수정이나 시뮬레이션 실행 없이 문서·코드·문헌을 대조하고 아래 부호와 인자를 직접 유도했습니다.

**1. 기하 미분과 쌍 힘: 정정된 식은 맞다. 원호는 정리 전·후 표현을 구분해야 한다.**

\(e=(x_j-x_i)/d\), \(u_i=v_i\cdot e\), \(\Delta u=u_i-u_j\)로 고정하겠습니다. Laguerre 면 위 점 \(r\)는

\[
|r-x_i|^2-W_i=|r-x_j|^2-W_j
\]

를 만족합니다. 시간 미분하면 면의 법선 속도는

\[
w(r)=
\frac{(r-x_i)\cdot v_i-(r-x_j)\cdot v_j}{d}
+\frac{\dot W_i-\dot W_j}{2d}.
\]

면 중심을 \(c=x_i+\phi de+t\), \(t\cdot e=0\)로 쓰면

\[
\boxed{
w_f=\phi u_i+(1-\phi)u_j
+\frac{(v_i-v_j)\cdot t}{d}
+\frac{\dot W_i-\dot W_j}{2d}
},
\qquad
\phi=\frac12+\frac{W_i-W_j}{2d^2}.
\]

따라서 \(W=0\)에서

\[
\frac{\partial V_i}{\partial x_i}\Big|_f
=A\left(\frac e2+\frac td\right),\qquad
\frac{\partial V_i}{\partial x_j}\Big|_f
=A\left(\frac e2-\frac td\right),
\]

\[
\frac{\partial V_j}{\partial x_i}\Big|_f
=-A\left(\frac e2+\frac td\right).
\]

[유도서 §2](/tmp/lageunha-push/Docs/tcf_derivation.md:21)의 세 식은 모두 맞습니다.

원호가 없는 내부 셀에서는 변분 힘을 먼저

\[
F_i^{\rm rev}
=\sum_j A(P_i-P_j)\left(\frac e2+\frac td\right)
\]

로 얻습니다. \(\sum_j Ae=0\)을 쓰면

\[
\boxed{
F_{ij}^{\rm rev}
=-A\left[\bar P e+(P_j-P_i)\frac td\right]
},\qquad \bar P=\frac{P_i+P_j}{2}.
\]

즉 접선항은 **\(-A(P_j-P_i)t/d\)**입니다. [코드의 `hq_ftx`, `hq_fty`](/tmp/lageunha-push/Exam/exam.c:7201)와 일치합니다.

원판 제한에서는 주의가 필요합니다. 반지름이 고정이고 생성점과 함께 움직이는 원호에 대해

\[
N_i=\int_{\Gamma_i^{\rm arc}}n\,ds,\qquad
\sum_j Ae+N_i=0
\]

이고, 원호의 원래 기하 미분 기여는 \(N_i\), 원래 변분 힘 기여는 \(P_iN_i\)입니다. 따라서

\[
F_i^{\rm rev}
=\sum_j A(P_i-P_j)\left(\frac e2+\frac td\right)+P_iN_i
=\sum_j F_{ij}^{\rm rev}.
\]

**정리된 \(-A[\bar Pe+(P_j-P_i)t/d]\)를 사용하면 원호 힘 \(P_iN_i\)는 이미 포함되어 있습니다. 여기에 다시 더하면 이중 계산입니다.** 반면 부피 변화에는 반드시

\[
\dot V_i=\sum_j Aw_f+v_i\cdot N_i
\]

를 포함해야 합니다. 원호의 법선 적분이 현의 법선 벡터와 같다는 Fable의 설명은 맞습니다.

따라서 §2의 “원호 힘은 자기 압력 × 현 법선 적분”은 **정리 전 변분식에 대한 설명으로만 승인**합니다. 코드의 [(H) 원호 처리](/tmp/lageunha-push/Exam/exam.c:7421)가 별도 힘 없이 \(-P_iv_i\cdot N_i\)의 에너지 일을 넣는 구조는 이 구분과 맞습니다.

**2. 에너지·엔트로피·각운동량: 적절한 조건에서 반이산 증명이 성립한다.**

\(W\)와 원판 반지름을 고정하고, 외력이 없다고 하겠습니다. 힘과 부피 미분에 같은 기하·셀 압력을 사용하면

\[
F_i^{\rm rev}=\sum_kP_k\frac{\partial V_k}{\partial x_i},
\qquad
\dot U_i^{\rm rev}=-P_i\dot V_i
\]

이므로

\[
\dot E
=\sum_i v_i\cdot F_i^{\rm rev}-\sum_kP_k\dot V_k=0.
\]

또한 \(S_i\)를 셀의 총엔트로피로 정의하면

\[
T_i\dot S_i=\dot U_i+P_i\dot V_i=0.
\]

가역 힘의 압력을 독립적으로 재구성한 면 압력으로 바꾸면 이 증명은 일반적으로 깨집니다. Fable의 해당 제한은 맞습니다.

소산항은

\[
F_{ij}^{Q}=-AQe,\qquad
H_{i,f}=H_{j,f}=\frac12AQ\Delta u
\]

로 두면

\[
(v_i-v_j)\cdot F_{ij}^{Q}+H_{i,f}+H_{j,f}
=-AQ\Delta u+AQ\Delta u=0,
\]

\[
\boxed{T_i\dot S_i=\frac12\sum_{f\ni i}A_fQ_f\Delta u_f\ge0}
\]

입니다. \(T_i>0\), \(Q_f\Delta u_f\ge0\)가 조건입니다. **열의 절반 배분 인자도 정확합니다.** 열을 절반씩 나누므로 엔트로피 증가량 자체는 온도가 다르면 서로 다릅니다.

Fable의 면별 잔차도 맞습니다. 원호 없는 내부 면에서

\[
(v_i-v_j)\cdot F_{ij}^{\rm rev}
-P_iAw_f+P_jAw_f
=-A(P_iu_i-P_ju_j).
\]

이는 셀의 면들을 합해 사라집니다. 원호 셀에서는 원호 일을 포함해야 사라집니다. 다만 “면별로 상쇄되지 않는다”는 것은 **이 특정 장부 배분에 대한 설명**입니다. 폐합 항등식을 이용하면 \(W=0\)에서 다음과 같은 반대칭 총에너지 교환식으로도 정리할 수 있습니다.

\[
\boxed{
\dot E_{i,f}=-A\left[
\frac{P_ju_i+P_iu_j}{2}
+\frac{(P_jv_i-P_iv_j)\cdot t}{d}
+Q\frac{u_i+u_j}{2}
\right],\qquad
\dot E_{j,f}=-\dot E_{i,f}.
}
\]

이는 현재 계획의 \(-A\Pi w\)와 다른 식입니다.

각운동량은 회전 불변성으로 증명됩니다.

\[
\sum_i x_i\times F_i^{\rm rev}
=\sum_kP_k\frac{\partial V_k}{\partial\theta}=0.
\]

고립계의 진정한 Voronoi/Laguerre 기하와 원형 제한에서는 성립하고, \(Q\) 힘은 중심력이므로 쌍별 토크가 0입니다. 고정 주기 상자나 회전 대칭을 깨는 제한 경계에서는 이 증명을 적용할 수 없습니다.

따라서 Fable의 각운동량 정정은 조건부로 맞지만, **VPH 가역 힘 자체가 중심력이라는 뜻은 아닙니다.** 계획서 첫 문장의 “중심력” 목표는 현재 유도와 양립하지 않습니다.

완전 이산에서는 별도 구분이 필요합니다. 총에너지 \(E_i\)를 위와 같이 합계가 0인 우변으로 RK 적분하면 합계 \(E\)의 선형 불변량은 보존할 수 있습니다. 그러나 \(U_i\)를 직접 RK4로 적분한다고 이차식인 운동에너지까지 정확히 보존되지는 않습니다. \(U=E-KE\) 복원도 엔트로피·내부에너지 양성을 자동 보장하지 않습니다. 이 구분은 [compatible hydrodynamics 문헌](https://cnls.lanl.gov/~shashkov/papers/comp.pdf)의 핵심과도 맞습니다.

**3. \(Q_f\): 엔트로피 부호 조건은 만들 수 있지만, 충격 정확도와 잡음 소멸은 아직 승인할 수 없다.**

고정된 \((\rho_L,p_L,\rho_R,p_R)\)에 대해 정확 Riemann 해는 비진공 영역에서

\[
f_L(p_*)+f_R(p_*)=\Delta u_{\rm rec},\qquad f'_L,f'_R>0
\]

를 만족하므로

\[
\frac{dp_*}{d\Delta u_{\rm rec}}>0.
\]

따라서

\[
P_c=p_*(0),\qquad Q_0=p_*(\Delta u_{\rm rec})-p_*(0)
\]

의 부호는 재구성 점프와 같습니다. **Fable의 \(P_c\) 정정은 정확해에 대해서 타당합니다.** 단, \(P_c\)는 밀도와 EOS에도 의존합니다.

“HLLC도 항상 단조”라는 문구는 범위가 너무 넓습니다. 현재 코드의 주 Davis 분기만 보면, \(c=\max(c_L,c_R)\), \(\delta=\Delta u_{\rm rec}\)에 대해 직접

\[
p_*(\delta)=
\frac{\rho_Rp_L+\rho_Lp_R}{\rho_L+\rho_R}
+\frac{\rho_L\rho_R}{\rho_L+\rho_R}
\begin{cases}
c\delta,&\delta<0,\\
(c+\delta)\delta,&\delta\ge0
\end{cases}
\]

를 얻으므로 단조성이 있습니다. 그러나 [실제 solver](/tmp/lageunha-push/Exam/exam.c:5530)는 압력비 분기, 이른 진공 판정, Roe fallback, 압력 바닥을 포함합니다. 이 전체 함수의 단조성·연속성은 별도로 확인해야 합니다.

엔트로피 증명에 필요한 것은

\[
Q\,\Delta u_{\rm cell}\ge0
\]

이지 \(Q\,\Delta u_{\rm rec}\ge0\)가 아닙니다. 문서의 연속형 식은 다음처럼 정의하면 명확합니다.

\[
Z_{\rm eff}=\frac{Q_0}{\Delta u_{\rm rec}}\ge0,\qquad
\boxed{Q=Z_{\rm eff}\operatorname{minmod}
(\Delta u_{\rm rec},\Delta u_{\rm cell})}.
\]

0에서는 적절한 극한값을 사용합니다. 이는 문서의 clip 식과 동등하며, \(Z_{\rm eff}\)가 유한하고 연속이면 0 부근도 연속입니다. 부호 불일치 시 \(Q=0\)은 에너지·엔트로피 관점에서 정합합니다. 다만 충격에서 필요한 감쇠까지 없어질 수 있으므로 충격 포착 보장은 아닙니다.

또한 일반적으로

\[
\bar P+Q_0=p_*+(\bar P-P_c)\ne p_*.
\]

따라서 “HLLC 수준의 소산”에서 곧바로 “HLLC와 같은 충격 해”를 결론내릴 수 없습니다. 양의 엔트로피 생성은 필요한 성질이지만, 진동 억제·양성·올바른 충격 수렴까지 보장하지 않습니다. [Tadmor의 엔트로피 안정성 이론](https://math.umd.edu/~tadmor/pubs/files/Tadmor_Acta2003.pdf)에서도 보존적 기본식과 소산의 정합성을 함께 요구합니다.

매끄러운 해에서 선형 정확 재구성과 규칙적인 해상도 계열을 가정하면 \(\Delta u_{\rm rec}=O(h^2)\), \(Q=O(h^2)\)라는 논리는 맞습니다. 그러나 **해상도 규모 속도 잡음에는 이 가정을 적용할 수 없습니다.** 잡음 진폭 \(\delta v\)가 줄지 않으면 대략

\[
\frac{\dot U_{\rm diss}}{V}
\sim\frac{\rho c(\delta v)^2}{h}
\]

까지 가능하므로, 가열이 해상도와 함께 줄기는커녕 증가할 수도 있습니다. §3(b)는 조건부 주장으로 낮춰야 합니다.

**4. VPH 위험: §8의 “0차 힘 오차”는 잘못된 진단이다.**

원호 없는 내부 셀에서 \(P_j=P_i\)이면 즉시 \(F_i=0\)입니다. 더 강하게, 생성점에 선형 압력

\[
P_i=P_0+b\cdot x_i
\]

를 부여하면

\[
\begin{aligned}
F_i
&=-\sum_j A\left(\frac d2e+t\right)(e\cdot b)\\
&=-\left[\int_{\partial V_i}(r-x_i)\otimes n\,ds\right]b\\
&=\boxed{-V_i b}.
\end{aligned}
\]

**정확한 면 기하와 접선항을 포함하면 불규칙 내부 셀에서도 선형 압력 힘은 정확합니다.** [Heß–Springel 부록 A](https://arxiv.org/pdf/0912.0629)의 식 A8–A11 및 선형장 시험과 일치합니다. Fable의 선형장 시험 제안은 유용하지만, 기대값은 “0차 오차의 크기 측정”이 아니라 기하·반올림 오차 수준입니다.

원호 자유 경계에서는 균일한 내부 압력과 외부 진공 사이의 압력 차이 때문에 일반적으로 힘이 0이 아닙니다. 따라서 균일 압력 무힘 시험은 내부 셀 또는 주기계로 한정해야 합니다.

실제 위험은 다음과 같습니다.

- 생성점 이동에 따른 부피 잡음이 EOS를 통해 실제 이산 압력 잡음으로 바뀝니다.
- 부피를 거의 바꾸지 않는 입자 변위 모드에는 가역 복원력이 약하거나 없습니다.
- 특히 규칙 격자의 Nyquist 모드에서 VPH 복원력이 사라지는 문제가 알려져 있습니다. Heß–Springel의 shape correction은 이 문제를 다룹니다.

선형 압력 힘의 정확성과 기하 부피의 전단 잡음은 양립합니다. 예를 들어 \(v(x)=v_i+L(x-x_i)\)이면 내부 Voronoi 셀에서

\[
\dot V_i
=V_i\operatorname{tr}L
-2\sum_f A_f\,t_f\cdot Le_f.
\]

따라서 \(\operatorname{tr}L=0\)인 선형 전단도 불규칙 셀의 부피를 바꿀 수 있습니다. **힘의 선형 정확성만으로 이 잡음이 제거되지는 않습니다.**

조향 금지 조건에서는 shape correction을 그대로 도입할 수 없습니다. 우선 균일·선형 압력 항등식, 선형 전단의 부피 변화, Nyquist 모드, 재구성 점프와 가열률의 해상도 수렴을 각각 구분해 검증해야 합니다. 소산을 줄이는 것만으로 복원력 부재까지 해결된다고 약속해서는 안 됩니다.

**5. Laguerre 확장: 평균압력과 직교 조건은 맞지만, 본문에 모순과 경계 조건 누락이 있다.**

첫 번째 유도에서 고정 \(W\)이면

\[
\frac{\partial V_i}{\partial x_i}\Big|_f
=A\left(\phi e+\frac td\right).
\]

폐합식을 이용하면

\[
\boxed{
F_{ij}^{\rm rev}
=-A\left[\widetilde P e+(P_j-P_i)\frac td\right],
\quad
\widetilde P=(1-\phi)P_i+\phi P_j.
}
\]

Fable의 인자는 정확하며 [코드의 `Ptil`](/tmp/lageunha-push/Exam/exam.c:7214)과 일치합니다.

가중치 변화에 대해서는

\[
\dot V_i^{W}\Big|_f=A\frac{\dot W_i-\dot W_j}{2d},
\]

\[
\sum_iP_i\dot V_i^{W}
=\sum_{f=(i,j)}
\frac{A(P_i-P_j)(\dot W_i-\dot W_j)}{2d}
=g\cdot\dot W.
\]

현재의 위치 힘과 \(\dot U_i=-P_i\dot V_i\)를 유지하면

\[
\boxed{\dot E=-g\cdot\dot W}.
\]

그러므로 실제 선택한 가중치 운동에 대한 조건은 **단일 스칼라 조건 \(g\cdot\dot W=0\)**입니다. \(g\ne0\)에서

\[
\dot W=q-g\frac{g\cdot q}{g\cdot g}
\]

로 투영할 수 있고, \(g=0\)이면 투영할 필요가 없습니다.

하지만 [§5 후반](/tmp/lageunha-push/Docs/tcf_derivation.md:64)의 “각 \(l\)에 대해 \(\sum_kP_k\partial V_k/\partial W_l=0\)이어야 한다”는 문장은 다른 조건입니다. 이는 \(g=0\)이라는 더 강한 조건이며, 자유로운 무관성 \(W\) 좌표의 정류 조건으로는 나올 수 있지만 **선택된 \(\dot W\)의 에너지 보존에 필요한 조건은 아닙니다.** 두 주장을 구분해야 합니다.

원판 제한까지 포함한 확장은 추가 설계가 필요합니다. \(W=0\), 공통 반지름에서는 공유 면 위 점의 두 생성점까지 거리가 같아 양쪽 원판 절단이 일치합니다. \(W_i\ne W_j\)이면 그렇지 않아, 각자의 동일 반지름 원판으로 자른 면의 길이·중심이 달라질 수 있습니다. 그러면 하나의 공유 \(A,t\)를 전제한 쌍 힘과 \(g\) 식을 그대로 적용할 수 없습니다. [현재 코드도 원판 경로를 W=0 단계로 명시](/tmp/lageunha-push/Exam/exam.c:3168)합니다.

따라서 “정적 Laguerre면 그대로 성립”은 **공유 면이 일치하는 유효한 셀 분할에서** 승인할 수 있습니다.

**6. 독창성: 기본 변분·호환 구조는 선행 연구이며, 제약 아래의 구체적 성과에 두어야 한다.**

Voronoi 기하를 이용한 열역학 정합 입자계는 [Serrano–Español 2001](https://doi.org/10.1103/PhysRevE.64.046115), 고정 질량의 변분 VPH는 Heß–Springel 2010에 이미 있습니다. 힘의 운동에너지 손실과 내부에너지 증가를 맞추는 구조도 compatible hydro의 기존 원리입니다.

문헌 구분도 필요합니다. 사용자가 언급한 **Caramana–Shashkov–Whalen 1998**은 [다차원 인공점성 논문](https://cnls.lanl.gov/~shashkov/papers/ed_vis.pdf)이고, 총에너지 compatible 알고리즘 논문에는 **Burton**이 추가됩니다. 이들을 “현재와 동일한 VPH+Riemann 조합의 직접 선행 사례”라고 한꺼번에 묶는 것은 부정확합니다.

또한 [Duque 2023](https://arxiv.org/abs/2212.12323)은 Voronoi·power diagram·Lagrangian 역학을 통합해 다룹니다. 따라서 “Laguerre로 확장” 자체만으로 독창성을 주장하기도 어렵습니다.

독창성 후보는 **고정 질량·생성점 일치·조향 금지라는 제약 아래**, 공유 자유 경계와 동적 \(W\)를 포함한 구체적인 호환 법칙, 잡음과 충격을 함께 다루는 검증된 \(Q\), 그리고 Kepler 진공 고리에서의 수렴 성과입니다. \(R_f\)가 변분식에서 사라지는 사실이나 \(g\cdot\dot W=0\)이라는 연쇄법칙 자체는 새로운 원리로 내세우기 어렵습니다.

마지막으로, G1 서명 전에 반드시 정리할 문서 간 불일치는 두 가지입니다.

첫째, [계획서 §1.4](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:47)는 직접 전개하면

\[
\boxed{
\sigma_f=A\left[
(\Pi-\bar P)(u_i-u_j)
+(P_i-P_j)(w-\bar u)
\right]
}
\]

입니다. **본문의 두 번째 항 앞 ‘−’는 ‘+’가 맞습니다.** 바로 아래 \(\Pi_{\rm ec}\) 식은 이 정정된 부호와 일치합니다. 또 \(\sigma_f\)는 \(T_i\dot S_i+T_j\dot S_j\)이지, 온도가 다른 경우의 총엔트로피 생성률 자체가 아닙니다.

둘째, [계획서 구현식](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:72)의

\[
F=-A\Pi n,\qquad \dot E=-A\Pi w
\]

에는 변분 접선 힘과 그에 맞는 에너지 교환이 없습니다. \(R_f\)를 합이 0인 열 교환으로 재배분하는 것만으로 합계 결함을 제거할 수도 없습니다. 현재 유도는 그 계획의 단순 구현이 아니라 **운동량·에너지 법칙을 함께 바꾸는 설계**입니다.

**G1 기록 제안: “Astra — 반이산 변분 구조 승인. 원호 장부, VPH 정확성 진단, \(Q_f\) 보장 범위, 동적 \(W\) 조건 및 구현 계획의 정합성을 정정한 뒤 재서명.”**