## 결론

1. **남은 가장자리 오차는 순위 실패보다 한쪽 스텐실의 절단오차가 주원인입니다.** `GFS_Q2_RANK`의 순위 확장으로 2차 적합이 유지되고 확장 셀이 31개뿐인데도 결과 곡선이 같다는 관측이 이를 뒷받침합니다. 원호에 닿는 셀과 그 면 이웃에 회전 공변 3차 적합을 적용하는 것이 첫 시도로 적합합니다.
2. **3차 적합만으로 Q의 고유 잔차가 없어지지는 않습니다.** 정확한 3차 속도장에서도 \(\delta_{b'}=d^3u'''/12\)가 남습니다. 국소 3차 미분으로 이를 보정하는 방법은 해석해를 쓰지 않는 정합적 고차 보정이지만, 충격파에서 Q를 보존한다는 보장은 없습니다. 부피 압축을 판정하는 회전 불변 게이트를 별도로 두고 Sod·쌍 시험을 통과시켜야 합니다.
3. **시험은 다항식·회전 불변성 → \(t=40\) 재시작 RHS 비교 → Sod와 쌍 시험 → Gresho → 짧은 Kepler 진화 → 장기 진화 순서**를 권고합니다.

## 1. 가장자리 기울기 적합

현재 `GFS_Q2_RANK`는 반경 \(2.5h\) 안의 후보부터 완전한 거리 껍질을 추가해 rank 5를 확인합니다. 구현은 SVD가 아니라 2회 modified Gram–Schmidt와 기계정밀도 기준을 사용합니다. 따라서 순위 확장은 적합 가능성을 개선하지만, 이미 rank가 확보된 한쪽 스텐실의 기울기 절단오차를 줄이지는 않습니다. [2차 적합과 순위 확장](/tmp/lageunha-push/Exam/exam.c:4764)

속도 성분마다 다음 회전 공변 기저로 가중 최소제곱을 풉니다. \(h=\sqrt{V_i}\), \(X=r_x/h\), \(Y=r_y/h\)입니다.

\[
\phi_3(X,Y)=
\left(
X,\ Y,\ \frac{X^2}{2},\ \frac{XY}{\sqrt2},\ \frac{Y^2}{2},\
\frac{X^3}{6},\ \frac{X^2Y}{2\sqrt3},\
\frac{XY^2}{2\sqrt3},\ \frac{Y^3}{6}
\right),
\]

\[
\min_c \sum_k \frac{h^2}{|\mathbf r_k|^2}
\left[\Delta v_k-\phi_3(\mathbf r_k/h)c\right]^2,
\qquad
G_i=\nabla v(\mathbf x_i).
\]

이는 9개 계수의 2차원 3차 다항식 적합입니다. rank 판정은 가중 설계행렬 \(A\)의 특이값으로 하고,

\[
\operatorname{rank}(A)=9
\quad\Longleftrightarrow\quad
\sigma_9>\max(n,9)\epsilon_{\rm mach}\sigma_1
\]

로 두는 편이 회전 불변 조건에 더 잘 맞습니다. 이는 물리적 잔차 문턱이 아니라 수치 rank 판정입니다. 초기 반경은 \(2.5h\)로 두고 rank가 부족하면 같은 거리의 후보를 껍질 단위로 추가합니다. 3차 rank가 안 되면 rank 5인 2차 적합으로 돌아가고, 그것도 안 되면 해당 비교를 유효한 적합으로 간주하지 않아야 합니다. 선형 fallback은 아핀 정확성은 유지하지만 2차 정확성은 보장하지 않습니다.

적용 대상은 `VORO_ARC_EDGE`에 닿는 셀과 그 면 이웃으로 시작하는 것을 권합니다. 가장자리 셀만 3차로 바꾸면 경계면 반대쪽의 2차 기울기 오차가 QBP 면 잔차에 남을 수 있습니다. 원호 표식으로 셀을 고르는 것은 기하 기반 판정이므로 물리적 잔차 문턱을 추가하지 않습니다.

원판 절단은 원호 기하를 표시하고, 진공 원호 보정은 현재 밀도와 압력 적합에만 들어갑니다. 속도 기울기 적합에는 원호 경계값이 들어가지 않습니다. [원판 원호 구성](/tmp/lageunha-push/Exam/exam.c:3954) [밀도·압력의 진공 원호 보정](/tmp/lageunha-push/Exam/exam.c:4700) 진공 경계에 임의로 \(v=0\)을 넣어서는 안 됩니다. 그 조건은 움직이는 자유표면의 속도 조건이 아니라 고정 벽의 무관통 조건에 해당합니다. 원호에 실제 속도 조건을 넣으려면 경계의 법선속도 같은 물리적 조건이 별도로 있어야 합니다.

이웃 Hessian을 쓰는 Hermite 결합도 가능하지만, 단순히 \(H_i=H_j\)를 강제하면 Hessian이 공간적으로 변하는 3차장을 잘못 제약합니다. 이를 쓰려면 공유 면에서 값과 기울기를 일치시키고 Hessian은 면까지 운반하는 다항식 패치를 구성해야 합니다. 우선순위는 국소 3차 적합보다 낮습니다.

## 2. 고유 잔차와 Q의 압축 조건

면 방향 속도를 \(u(s)=\mathbf v(\mathbf x_i+s\mathbf e)\cdot\mathbf e\)라 두고,

\[
\Delta u_{\rm cell}=u(0)-u(d),\qquad
\delta_{b'}=\Delta u_{\rm cell}
+\frac d2\mathbf e^{\mathsf T}(G_i+G_j)\mathbf e
\]

라 하면, 정확한 기울기에서도

\[
\delta_{b'}=
\frac{d^3}{12}u'''(d/2)+O(d^5).
\]

따라서 3차 적합은 **기울기 적합 오차**를 줄이지만, 이 고유 잔차 자체를 없애지는 않습니다. 실제 TCF 경로도 QBP에서 `qin`만 바꾸고 `mm=minmod(qin, ducell)`을 적용합니다. 기본 \(tcf\_Z\)는 별도로 재구성 속도차 `tcf_durec`에서 구하므로 QBP 입력과 임피던스 기준이 다릅니다. [QBP의 제한 전 기울기 캐시](/tmp/lageunha-push/Exam/exam.c:673) [TCF 임피던스 계산](/tmp/lageunha-push/Exam/exam.c:7657) [TCF Q 입력과 열 계산](/tmp/lageunha-push/Exam/exam.c:8227) 이는 이전 감사의 일관성 지적과도 같습니다. [12차 감사](/tmp/lageunha-push/Docs/audit_impl12_astra.md)

국소 3차 텐서 \(T_i=\nabla^3\mathbf v_i\)를 적합해 다음 보정을 시험할 수 있습니다.

\[
\delta_5 =
\delta_{b'}-\frac{d^3}{24}
\mathbf e\cdot(T_i+T_j)[\mathbf e,\mathbf e,\mathbf e].
\]

이 보정은 국소 자료에서 얻은 미분만 사용하며, 정확한 3차 속도장에서는 \(\delta_5=0\)이 됩니다. 한쪽 스텐실의 3차 미분 오차까지 고려하면 일반적으로 잔여 차수는 \(O(d^4)\)이고, 더 대칭적이고 정확한 중점 추정에서는 \(O(d^5)\)까지 기대할 수 있습니다. 다만 충격파는 매끄러운 3차 전개가 성립하지 않으므로, 이 보정을 Q 입력 전체에 무조건 적용해도 Sod 점성을 보존한다고 단정할 수 없습니다.

또한 현재 `minmod`가 거르는 것은 **면 쌍의 접근**이지 부피 압축만은 아닙니다. 전단 흐름에서도 일부 면 쌍은 접근할 수 있습니다. 압축만 허용하는 회전 불변 게이트의 시험 후보는 제한 전 기울기 \(G_i\)로

\[
S_i=\frac{G_i+G_i^{\mathsf T}}2,\qquad
\theta_i=\operatorname{tr}G_i,
\]

\[
\chi_i=
\begin{cases}
\min\!\left(1,\dfrac{[-\theta_i]_+}{\|S_i\|_F}\right),&\|S_i\|_F>0,\\
0,&\|S_i\|_F=0,
\end{cases}
\qquad
\chi_f=\max(\chi_i,\chi_j)
\]

를 두는 것입니다. 순수한 발산 없는 전단은 \(\chi=0\), 평면 1차원 압축은 \(\chi=1\)입니다. \( [x]_+=\max(x,0) \)에는 조정 가능한 스무딩 문턱이 없습니다. 다만 수치적으로 추정한 \(\theta\)가 오차를 가지므로 비다항 Kepler 흐름에서 정확히 0인 게이트를 보장하지는 않습니다. 임의의 유한 스텐실만으로 작은 실제 압축과 적합 오차를 문턱 없이 완벽히 구분할 수도 없습니다.

TCF와의 일관성을 위한 시험식은 다음과 같습니다.

\[
m_f=\operatorname{minmod}(\delta_5,\Delta u_{\rm cell}),
\]

\[
Q_f=
\begin{cases}
\chi_f Z_{\delta}m_f,&\text{기본 Q 법칙},\\[2mm]
\chi_f Z_{\delta}\dfrac{m_f^2}{\Delta u_{\rm cell}},
&\text{Q 법칙 2},\quad\Delta u_{\rm cell}\ne0,
\end{cases}
\]

\[
Z_\delta=
\frac{\Phi(\delta_5)-\Phi(0)}{\delta_5},\qquad
\Phi(\delta)=P^*(\bar u+\delta/2,\bar u-\delta/2).
\]

\(\delta_5=0\)에서는 \(Z_\delta=\Phi'(0)\)의 극한을 사용합니다. `GFS_TCF_FACE_Z`나 `GFS_TCF_ZH`를 쓸 때는 그 분기에서의 임피던스도 동일한 입력 기준으로 맞춰야 합니다. \(Z_\delta\ge0\)이면 두 법칙 모두

\[
H_f^Q=\frac12 A_f Q_f\Delta u_{\rm cell}\ge0
\]

를 유지합니다. 이 게이트는 후보이지 충격 보존의 증명은 아닙니다. 특히 단일 쌍만 있는 시험에서 부피 압축을 추정할 스텐실이 없다면 Q가 꺼질 수 있으므로, 쌍 시험이 요구하는 불변량이 에너지 상쇄·비음수 열인지, 접근 쌍에서 반드시 비영 Q인지 먼저 구분해 판정해야 합니다. [13차 감사](/tmp/lageunha-push/Docs/audit_impl13_astra.md)

## 3. 시험 순서

1. **다항식 정확도와 회전 검사:** 원호 셀을 포함해 아핀·2차 속도장은 기울기 정확도와 \(\delta_{b'}\simeq0\)을 확인합니다. 3차장은 \(\delta_{b'}\simeq d^3u'''/12\), 보정식은 \(\delta_5\simeq0\)을 확인합니다. 좌표를 \(\pi/7\) 회전해 기울기와 RHS가 회전 공변인지 비교합니다.
2. **동일 \(t\simeq40\) 재시작 RHS 비교:** 한 스레드·같은 상태에서 기준 QBP, 3차 기울기만, \(\delta_5\) 보정, 압축 게이트를 각각 분리해 비교합니다. 면·고리별 \(\delta_{b'},\delta_5,\Delta u_{\rm cell},\theta,\chi,Z,Q,H_Q\), rank와 조건수를 기록하고 \(H_Q\ge0\), QBP 무효 면 수 0을 확인합니다.
3. **Sod Gate C와 쌍 시험:** 압축면 Q가 유지되는지, 에너지 교환이 쌍별로 상쇄되는지, 음의 열이 생기지 않는지 확인합니다.
4. **Gresho:** 전단 가열과 소용돌이 감쇠를 확인합니다.
5. **짧은 Kepler 재시작 진화:** \(t=40\)에서 \(t=45\) 또는 \(50\)까지 비교하고 가장자리 \(H_Q/U\), \(v_r/v_K\), 압력 상승, 누출 질량을 봅니다. 통과한 조합만 \(t=60\)–72 장기 진화로 넘깁니다.

이번 작업은 지정 문서와 코드의 읽기 전용 감사이며 파일 수정이나 실행은 하지 않았습니다.