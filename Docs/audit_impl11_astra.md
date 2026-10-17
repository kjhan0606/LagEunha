## 결론

제공된 진단은 코드와 측정 패턴에 잘 부합합니다. `av_mode==5`의 기본 속도 제한자는 \(v_x\), \(v_y\)를 따로 제한해 속도 기울기 텐서의 두 행에 서로 다른 계수를 곱합니다. 이 처리는 좌표 회전에 불변하지 않습니다. 진공 경계의 한쪽 stencil에서 성분별 극값이 생기면 해당 성분 기울기만 과도하게 깎일 수 있습니다. 그 오차가 면 재구성 속도 점프를 남기고, 보고된 면들에서는 그 점프가 TCF 열을 직접 결정합니다. [기울기 제한자](</tmp/lageunha-push/Exam/exam.c:4836>) [TCF 재구성](</tmp/lageunha-push/Exam/exam.c:7338>) [Q 계산](</tmp/lageunha-push/Exam/exam.c:8072>)

다만 원인을 이 제한자에 확정 귀속하려면 실행 설정도 확인해야 합니다. `EUNHA_LIM_PAIR=1`이면 해당 셀 단위 TVD 제한은 건너뛰지만, `EUNHA_REC_TCOMMON=1` 경로의 후속 제한이 켜져 있으면 \(v_x,v_y\)를 각각 클램프합니다. 따라서 두 설정과 `GFS_ORACLE_EDGE` 값을 기록해야 실제 적용 경로를 특정할 수 있습니다. [분기 조건](</tmp/lageunha-push/Exam/exam.c:4953>) [쌍 제한 경로](</tmp/lageunha-push/Exam/exam.c:7357>)

## 판정 근거와 TCF 경로

속도 구배를 \(G_i=\nabla\mathbf v_i\)라 하면, 현재 제한자는 셀 주변의 \(v_x\), \(v_y\) 범위를 각각 구하고 독립 계수 \(\alpha_{vx},\alpha_{vy}\)를 적용합니다.

\[
G_i^{\mathrm{lim}}
=
\begin{pmatrix}
\alpha_{vx}&0\\
0&\alpha_{vy}
\end{pmatrix}G_i.
\]

이는 좌표축에 의존합니다. 좌표를 회전하면 속도 성분의 극값과 두 제한 계수가 달라질 수 있습니다. 진공 호에 닿은 셀에서는 이웃 속도 범위가 실제 이웃만으로 정해지므로, 회전 흐름의 한 성분이 경계에서 극값이 되면 보고된 종류의 편향이 생길 수 있습니다. 곡률로 범위를 넓히는 선택지는 `EUNHA_LIM_SMOOTH`가 켜지고 이차 적합이 매끄럽다고 판정된 셀에만 적용됩니다. [범위 확장 조건](</tmp/lageunha-push/Exam/exam.c:4890>)

TCF에서 \(d u_{\rm rec}\) 자체는 `minmod`의 입력으로 쓰입니다. 하지만 재구성의 영향이 그 입력에만 한정되지는 않습니다. MUSCL로 만든 \(\rho_L,\rho_R,p_L,p_R,u_L,u_R\)는 HLLC 상태 계산에 들어가며, `tcf_Z`도 그 재구성 상태의 압력 차로 계산됩니다. 따라서 재구성은 \(Q\)의 제한 점프 \(m\)과 임피던스 \(Z\)를 통해 열에 영향을 줍니다. 반면 TCF의 가역 힘·일 식은 셀 중심 압력과 속도를 사용합니다. [HLLC 재구성 상태](</tmp/lageunha-push/Exam/exam.c:7467>) [재구성 상태에서 Z 계산](</tmp/lageunha-push/Exam/exam.c:7503>) [TCF 가역·Q 항](</tmp/lageunha-push/Exam/exam.c:8103>)

코드의 법선 속도 차와 minmod는

\[
d u_{\rm rec}=u_L-u_R,\qquad
d u_{\rm cell}=(\mathbf v_i-\mathbf v_j)\cdot\mathbf e,
\]

\[
m=\operatorname{minmod}(d u_{\rm rec},d u_{\rm cell})
=
\begin{cases}
\operatorname{sgn}(d u_{\rm cell})
 \min(|d u_{\rm rec}|,|d u_{\rm cell}|),&
d u_{\rm rec}d u_{\rm cell}>0,\\
0,&\text{그 밖의 경우}.
\end{cases}
\]

기본 법칙은 \(Q=Zm\)이고, `GFS_TCF_QLAW=2`이면 \(Q=Zm^2/d u_{\rm cell}\)입니다. 한쪽 열 기여는 \(H^Q_f=\tfrac12 A_f Q\,d u_{\rm cell}\)입니다. 보고된 상위 면에서는 \(|d u_{\rm rec}|\ll|d u_{\rm cell}|\)이고 부호가 일치하므로 \(m=d u_{\rm rec}\)입니다. 기본 법칙에서는 잘못 남은 재구성 점프가 \(Q\)에 직접 비례합니다. [minmod와 Q 법칙](</tmp/lageunha-push/Exam/exam.c:8072>) [열 누적](</tmp/lageunha-push/Exam/exam.c:8136>)

## 권고하는 최소 변경

`av_mode==5 && gfs_tcf_on()`일 때만 속도 제한을 회전 불변 방식으로 바꾸고, 밀도·압력 제한과 TCF의 Q 법칙은 그대로 둡니다. 각 실재 면 \(f=(i,j)\)에 대해 \(\mathbf e_f=(\mathbf x_j-\mathbf x_i)/d\), \(r_{if}=c_f-\mathbf x_i\)라 둡니다. 면의 투영 속도 범위는 해당 면의 실재 국소 생성점, 예를 들어 두 인접 Delaunay 삼각형의 \(\{i,j,k_+,k_-\}\), 에서 구합니다.

\[
U_f^-=\min_{k\in S_f}\mathbf v_k\cdot\mathbf e_f,\qquad
U_f^+=\max_{k\in S_f}\mathbf v_k\cdot\mathbf e_f,
\]

\[
u_{if}=\mathbf v_i\cdot\mathbf e_f,\qquad
\delta u_{if}=\mathbf e_f^\mathsf T G_i r_{if}.
\]

각 셀 면의 제한 계수를 다음과 같이 계산합니다.

\[
\alpha_{if}=
\begin{cases}
\min\!\left(1,\dfrac{U_f^+-u_{if}}{\delta u_{if}}\right),
&\delta u_{if}>0,\\[6pt]
\min\!\left(1,\dfrac{U_f^--u_{if}}{\delta u_{if}}\right),
&\delta u_{if}<0,\\
1,&\delta u_{if}=0,
\end{cases}
\qquad
\alpha_i=\min_{f\in F_i}\max(0,\alpha_{if}),
\qquad
G_i^{\rm lim}=\alpha_iG_i.
\]

이렇게 하면 하나의 스칼라 \(\alpha_i\)가 속도 기울기 텐서 전체에 적용됩니다. 회전하면 \(\mathbf v,G_i,r,\mathbf e\)가 함께 변환되지만 모든 투영값과 \(\alpha_i\)는 같습니다. 호 경계의 빠진 쪽에는 가상 속도 0을 넣지 않고, 실제로 존재하는 이웃만 범위에 포함합니다. `EUNHA_LIM_PAIR=1` 경로를 TCF에서 사용할 경우에도 뒤의 성분별 \(v_x,v_y\) 클램프를 이 투영 제한으로 대체해야 회전 불변성이 유지됩니다.

## 시험 순서

1. 같은 비퇴화 IC를 대칭 각도가 아닌 각도(예: \(\theta=\pi/7\))만큼 위치·속도·경계 중심과 함께 회전합니다. 초기 RHS에서 면별 \(d u_{\rm rec},Z,Q,H_Q\)를 입자 대응시켜 비교하고, 짧은 진화에서도 같은 비교를 합니다. 회전 후 결과를 원래 좌표로 되돌렸을 때 차이가 수치 오차 수준인지 확인합니다.
2. Kepler 가장자리에서 \(k=0,1\)의 \(H_Q/U\), 상위 5% 면의 열 분율, 정확 기울기 대비 오차를 전후 비교합니다. 내부와 바깥 가장자리도 함께 기록합니다.
3. Sod를 실행하되 관문 C의 정의와 임계값은 그대로 유지합니다.
4. 마지막으로 Gresho를 확인합니다. 모든 비교에서 TCF 면-Z 설정과 Q 법칙을 고정해 제한자 변경 효과를 분리합니다.

요청대로 읽기 전용 검토만 했으며 파일은 수정하지 않았습니다.