## 결론

**엔트로피 생성 \(D_f\)와 가역 PdV 일을 같은 면 장부에서 분리해야 합니다.** 셀 압력 \(P_i,P_j\)를 각각 쓴 정의는 총에너지 장부와 맞지만 셀별 부호가 정해지지 않습니다. 대칭 면 정의는 음향 근사에서 비음수지만, 주어진 기하 잔차를 남깁니다. 공통 \(P_{\rm ref}\)로 잔차를 흡수하는 식은 유도할 수 있으나, 일반 HLLC 상태에서 그 선택이 \(D_f\ge0\)까지 보장하지는 않습니다. 일반 상태의 목표는 보존 flux와 엔트로피 안정 ALE flux를 함께 쓰는 것입니다.

특히 현재 코드의 `av_mode=5`는 HLLC \(P^*\)를 힘·일에 쓰지만 이를 `dK_diss`로 분해하지 않습니다. 이 모드에서 `entropy_mode=1`만 켜면 HLLC 충격 소산이 \(K\) 증가로 전달되지 않을 수 있습니다. [HLLC 압력과 소산 분해 경로](</tmp/lageunha-push/Exam/exam.c:6883>), [가역 PdV와 \(dK_{\rm diss}\) 분리](</tmp/lageunha-push/Exam/exam.c:7052>), [앞선 자문](</tmp/lageunha-push/Docs/consult_thermal_runaway_reply_astra.md:48>)

### 1. 두 정의와 \(P_{\rm ref}\)

면 법선을 \(L\to R\), \(\Delta u=u_L-u_R\), \(\bar u=(u_L+u_R)/2\), \(\bar P=(P_L+P_R)/2\)로 두겠습니다. 오른쪽 셀의 outward 방향은 반대이므로, 주어진 셀 정의를 면 양쪽에 적용하면

\[
\begin{aligned}
q_L&=-A(P^*-P_L)(w-u_L),\\
q_R&=+A(P^*-P_R)(w-u_R).
\end{aligned}
\]

두 셀의 합은

\[
q_L+q_R
=A(P^*-\bar P)\Delta u
+A(P_L-P_R)(w-\bar u)
=D_{\bar P}-R_f,
\]

여기서

\[
D_{\bar P}=A(P^*-\bar P)\Delta u,\qquad
R_f=A(P_L-P_R)(\bar u-w)
\]

입니다. 즉, 대칭 정의 \(D_{\bar P}\)만 절반씩 배분하면 질문에 적힌 \(R_f\)가 남습니다.

한 가지 공통 기준압력 \(P_{\rm ref}\)를 써서 양쪽 열원 합을 \(A(P^*-P_{\rm ref})\Delta u\)로 만들면, 위 셀 정의와 같은 합을 얻는 선택은

\[
\boxed{P_{\rm ref}
=\bar P+\frac{(P_L-P_R)(\bar u-w)}{\Delta u}}
\qquad(\Delta u\ne0).
\]

이는 주어진 부호 규약에서 잔차를 대수적으로 흡수합니다. 더 일반적으로, 에너지 장부에 맞는 목표 \(D_f\)를 먼저 정했다면

\[
P_{\rm ref}=P^*-\frac{D_f}{A\Delta u}
\]

가 그 목표와 짝을 이루는 압력입니다. **하지만 이 식만으로 \(D_f\ge0\)가 보장되지는 않습니다.** \(\Delta u\to0\)에서는 특이해지고, 압력 기준이 셀 압력 범위를 벗어날 수도 있습니다. 따라서 음수 \(dK\)를 셀 합산 뒤 잘라내면 에너지 장부가 다시 어긋납니다.

동일한 임피던스 \(Z=\rho c\)의 선형 음향 근사에서는 \(P^*=\bar P+\tfrac12Z\Delta u\)라서

\[
D_f=\tfrac12 A Z(\Delta u)^2\ge0.
\]

좌우 임피던스가 다르면 선형화 Riemann 해의 기준압력은 산술평균이 아니라

\[
P_{\rm ref,ac}
=\frac{Z_RP_L+Z_LP_R}{Z_L+Z_R},\qquad
Z_{\rm eff}=\frac{Z_LZ_R}{Z_L+Z_R},
\]

이고 \(P^*-P_{\rm ref,ac}=Z_{\rm eff}\Delta u\), 따라서 \(D_f=A Z_{\rm eff}(\Delta u)^2\ge0\)입니다. 이는 **선형 음향 근사** 결과입니다. 비선형 HLLC에서 이 부호 관계가 자동 성립한다고 볼 수는 없습니다.

일반적인 해법은 ALE 보존 flux를 양쪽 셀에 같은 값·반대 부호로 적용하고, 엔트로피 변수 점프에 대한 양의 준정부호 소산을 더하는 것입니다. 그러면 총에너지 보존은 공유 flux에서, 엔트로피 부등식은 소산항에서 나옵니다. 이동 격자에서는 mesh-motion 항과 기하 보존 법칙도 함께 맞아야 합니다. [이동 격자 엔트로피 조건](https://arxiv.org/abs/1812.09093), [엔트로피 안정 유한체적 flux](https://arxiv.org/abs/1209.4994)

코드에는 `GFS_SCHEME_H` 조건에서 \(D_f=A[(P_i-P_j)\delta+(P^*-P_{\rm til})(u_i-u_j)]\)를 만들고 `dte`에 분배하는 별도 경로가 있습니다([Exam/exam.c:7138](</tmp/lageunha-push/Exam/exam.c:7138>)). 이는 기하 항을 포함한 에너지 장부의 유사 사례입니다. 다만 그 분기 자체가 모든 경로의 엔트로피 안정성을 보장하거나 일반 `stress.dK`와 연결되는 것은 아닙니다.

### 2. \(D_f\)에 쓸 상태

\(D_f\)는 **현재 RK 단계, 해당 면의 재구성 상태**로 계산해야 합니다. 즉 \(P_{L,R},\rho_{L,R},u_{L,R},c_{L,R}\)는 HLLC에 실제로 넣는 제한된 면 상태여야 합니다. 엔트로피 열원을 \(K_i\)로 바꿀 때에는 별도로 셀의 단계 밀도 \(\rho_i=m_i/V_i\)를 사용합니다.

비대칭 임피던스 근사를 쓴다면 \(Z_{L,R}=\rho_{L,R}c_{L,R}\)도 그 면 상태로 계산합니다. 현재 `av_mode=5` 경로는 면에서 재구성된 \(P,\rho\)를 HLLC에 넘기면서 음속에는 셀 중심 `ibp_csound`, `jbp_csound`를 넘깁니다([Exam/exam.c:6551](</tmp/lageunha-push/Exam/exam.c:6551>)). 음향 근사 진단을 한다면 재구성 상태와 EOS상 일관된 음속을 써야 합니다.

2점 Gauss 면 적분에서는 먼저 평균을 낸 뒤 곱하지 말고,

\[
D_f=A\sum_{q=1}^{2}\frac12
(P_q^*-P_{{\rm ref},q})(u_{L,q}-u_{R,q})
\]

를 계산하십시오. 코드도 각 Gauss 점에서 별도의 재구성과 HLLC 상태를 구한 뒤 힘에는 평균 \(P^*\), 일에는 공분산 보정을 사용합니다([Exam/exam.c:6362](</tmp/lageunha-push/Exam/exam.c:6362>)). 엔트로피 열원도 이 점별 적분과 같은 가중치를 써야 합니다. 면 속도 \(w_q\), 가역 일, 힘에도 동일한 기하·시간 처리를 써야 합니다. 전단 유동의 움직이는 격자에서 고차 면 적분이 격자 잡음을 줄인다는 결과도 이 접근을 뒷받침합니다. [Zier & Springel 2022](https://arxiv.org/abs/2205.07916)

진공 분기는 별도 취급이 필요합니다. 코드의 `vnR-vnL > cmax` 분기는 \(P^*=10^{-30}\)을 반환합니다([Exam/exam.c:5490](</tmp/lageunha-push/Exam/exam.c:5490>)). 이 반환값만으로 실제 진공 발생이나 물리적 충격 소산을 판정할 수 없습니다. 해당 분기 횟수와 상태를 기록하고, 팽창 해에서 \(K\) 증가가 해상도와 시간간격을 줄일수록 사라지는지 확인하십시오. 진공에서 양의 압력을 유지하는 fallback과 엔트로피 장부도 서로 일치해야 합니다.

### 3. RK4, 양성, 바닥값, 충격

이상기체에서 \(K=P/\rho^\gamma\)라 두면 비열적 PdV는 \(K\)를 바꾸지 않습니다. 비가역 열원만 있을 때

\[
u=\frac{K\rho^{\gamma-1}}{\gamma-1}
\quad\Longrightarrow\quad
\frac{dK}{dt}
=\frac{\gamma-1}{\rho^{\gamma-1}}\dot u_{\rm diss}
=\frac{\gamma-1}{m\rho^{\gamma-1}}\dot Q_i,
\]

여기서 \(\dot Q_i\)는 셀에 들어가는 열률입니다. 현재 코드의 `stress.dK` 변환은 이 식과 일치합니다([Exam/exam.c:7402](</tmp/lageunha-push/Exam/exam.c:7402>)). Springel & Hernquist의 엔트로피 방정식도 단열 구간에서 \(A=P/\rho^\gamma\)를 일정하게 두고, 충격의 인공점성 가열을 엔트로피원으로 넣어 운동에너지 손실과 열 증가를 맞춥니다. 그 유도는 SPH의 쌍별 힘 구조에 맞춘 것이므로 HLLC 면 flux에 그대로 대입하기보다 이 코드의 면 에너지 장부와 짝을 맞춰야 합니다. [Springel & Hernquist 2002](https://arxiv.org/html/astro-ph/0111016v3)

코드는 RK4의 네 단계에서 \(K\)도 적분하고, 각 단계의 기하 부피에서 \(\rho=m/V\), \(P=K\rho^\gamma\)를 다시 계산합니다([RK 단계](</tmp/lageunha-push/Exam/exam.c:12323>), [최종 RK4 조합](</tmp/lageunha-push/Exam/exam.c:12544>), [압력 갱신](</tmp/lageunha-push/Exam/exam.c:12707>)). 세 가지 점을 관리해야 합니다.

- **RK4 단계 양성:** 최종 RK4 가중합이 양수여도 중간 단계는 음수가 될 수 있습니다. 단계 바닥값 처리는 유효한 \(K\)를 보장할 수 있지만, 바닥값으로 늘어난 내부에너지를 누적 장부에 기록해야 합니다. 코드는 엔트로피 모드에서 `GAS_K_FLOOR`를 적용합니다.
- **\(10^{-6}\)의 의미:** 코드의 에너지 모드 압력 바닥 \(P_{\min}=10^{-6}\)과 엔트로피 모드의 상수 \(K_{\min}\)은 같은 바닥이 아닙니다. 압력 바닥을 유지하려면 \(K_{\min}(\rho)=10^{-6}\rho^{-\gamma}\)처럼 밀도 의존 변환이 필요합니다. 각 바닥 적용 전후의 \(I=PV/(\gamma-1)\) 차이를 에너지 주입으로 기록하십시오([압력·\(K\) 바닥 분기](</tmp/lageunha-push/Exam/exam.c:13290>)).
- **충격의 엔트로피 생성:** 충격에서 \(K\)는 올라가야 합니다. `av_mode=5`의 HLLC 소산을 \(dK\)에 추가하지 않으면 충격 후 압력이 올바르게 회복되지 않을 수 있습니다. \(P^*\)로 인한 운동에너지 손실과 \(K\)의 내부에너지 증가는 같은 RK 단계에서 짝을 이루어야 합니다.

Hopkins 2015의 MFM/MFV는 질량·운동량·총에너지의 보존 flux를 직접 갱신합니다. 고마하수 문제에는 내부에너지의 정확도를 높이려고 내부에너지를 별도로 진화시키는 dual-energy 접근도 사용하며, 이 경우 총에너지 보존은 기계 정밀도가 아니라 시간 적분 오차 수준이 됩니다. 따라서 여기서도 보존 총에너지와 \(K\)로 복원한 열에너지를 별도 비교하는 것이 유용합니다. [Hopkins 2015, 보존 flux와 dual-energy 절](https://arxiv.org/html/1409.7395v2)

### 4. 검증과 총에너지 장부

각 RK 단계에서 정렬된 입자 ID 쌍별로 \(A_q,w_q,P_q^*,P_{{\rm ref},q},u_{L,q},u_{R,q}\), HLLC 분기, \(D_{f,q}\)의 부호를 한 번씩 기록하십시오. 같은 점별 상태와 RK 가중치를 힘·가역 일·엔트로피 열원에 적용합니다. 총에너지는

\[
\mathcal E
=\sum_i\left[
\frac12m_i|\mathbf v_i|^2+
\frac{K_i m_i\rho_i^{\gamma-1}}{\gamma-1}
\right]
\]

로 계산하고,

\[
R_E=\mathcal E_{n+1}-\mathcal E_n
-W_{\rm external}-E_{\rm boundary}-E_{\rm floor}
\]

를 봅니다. \(D_f\)는 운동에너지가 열로 바뀌는 **내부 변환**이므로 외부 주입항으로 빼지 않습니다. 면별로는 \(R_f\)를, 기하별로는

\[
R_{\rm GCL}
=\Delta V_{\rm geom}
-\Delta t\sum_s b_s\sum_q A_qw_{q,s}
\]

를 따로 기록하십시오. 사용자가 제시한 회전 격자의 큰 `GFS_VOL_AUDIT` 잔차는 엔트로피 증가량만으로 총에너지 정합성을 판단하기 어렵다는 근거입니다([요청 문서](</tmp/lageunha-push/Docs/consult_entropy_dissipation.md>), [기존 기하·RK 자문](</tmp/lageunha-push/Docs/consult_rk_topology_reply_astra.md:16>)).

검증 순서는 다음이 좋습니다.

1. **면 단위 대수 검사:** 동일 \(Z\)의 압축·팽창, 비대칭 \(Z\), \(\Delta u\simeq0\), 좌우 상태 교환에서 \(D_f\) 부호와 쌍별 총에너지 잔차를 확인합니다.
2. **회전 Hill 띠:** 기존 37° 격자에서 장기 \(K\), \(E_{\rm int}\), 속도 잡음, 바닥 주입, \(R_E\), \(R_{\rm GCL}\)을 추적합니다. \(\Delta t\) 절반과 1점·2점 Gauss 적분 비교를 분리해 봅니다.
3. **Sod 충격관:** 밀도·속도·압력 \(L_1\), 충격 위치, 충격 후 \(K\) 증가를 해석해와 비교합니다. 충격 앞뒤 매끈한 영역의 \(K\) 보존도 확인합니다.
4. **진공 팽창:** 해당 분기 빈도, 최소 압력·밀도, \(K\) 증가와 총에너지 잔차를 기록하고 해상도·시간간격 수렴을 확인합니다.
5. **Kepler \(64^2\), \(256^2\):** 장기 궤도에서 \(E_{\rm int}\), 총에너지, 각운동량, \(K\) 바닥 활성 셀 수, 기하 잔차를 비교합니다.

요청대로 문서와 코드를 읽기만 했고, 파일 수정이나 검증 계산은 실행하지 않았습니다.