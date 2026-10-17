## 결론

원판 경계를 가진 셀에는 **내부 면 HG와 원호 HG를 함께** 둡니다. 원호는 진공 쪽 이웃을 만들지 않고, 생성점에 붙은 실제 원호 기하에서 \(\xi\)와 그 기울기를 계산하는 셀 단방향 모드로 취급합니다. 3차 QBP에는 회전 불변 압축 게이트를 결합해, 순수 압축에서 `du_cell`을 Q 입력으로 복원하는 방안을 권고합니다.

현재 [sub_build](/tmp/lageunha-push/Exam/exam.c:477)은 `VORO_ARC_EDGE`에서 즉시 실패하고, [8635행의 capfree 검사](/tmp/lageunha-push/Exam/exam.c:8633)는 모든 음수 경계 ID를 배제합니다. `hg_cell`은 모든 경계를 내부 이웃 면으로 가정하므로 이 검사를 단순히 제거하면 안 됩니다. 원판 클리퍼는 원호를 짧은 조각들로 표시하고, [Area2DPolygon](/tmp/lageunha-push/Voro/voro.c:1134)은 원호의 면적만 보정합니다. 원호 HG에는 면적뿐 아니라 중심과 기하 미분도 정확히 포함해야 합니다.

### 1. 원호 아영역과 기울기

셀 \(i\)의 연결된 원호 성분 \(\Gamma_f\)에 대해 \(x_i\)를 원의 중심, \(R\)을 반경, \(\ell_f\)를 원호 길이, \(n\)을 바깥 법선이라 둡니다.

\[
B_f(p)=\frac12\int_{\Gamma_f}(r-p)\cdot n\,ds,\quad
W_f=B_f(x_i)=\frac12R\ell_f,\quad
N_f=\int_{\Gamma_f}n\,ds.
\]

셀 도심을 \(c_i\)라 하면 원호 중심이 생성점과 함께 이동하므로

\[
B_f(c_i)=W_f-\frac12(c_i-x_i)\cdot N_f,\qquad
\eta_f=\frac{B_f(c_i)}{W_f},\qquad
\xi_{if}=\log\eta_f.
\]

아영역 상태와 길이 척도는 설계 문서의 내부 면 정의를 그대로 확장합니다.

\[
\rho_{if}=\rho_i\eta_f,\quad
P_{if}=K_i\rho_{if}^{\gamma},\quad
c_{if}^2=\gamma P_{if}/\rho_{if},\quad
h_{if}=\frac{2W_f}{\ell_f}=R.
\]

기울기는 \(B=B_f(c_i)\), \(\delta=c_i-x_i\)로 두고 다음과 같이 계산합니다.

\[
g_{if,a}=\partial_{x_a}\xi_{if}
=\left(\frac1B-\frac1{W_f}\right)\partial_{x_a}W_f
-\frac{N_f\cdot\partial_{x_a}\delta+\delta\cdot\partial_{x_a}N_f}{2B}.
\]

원호를 \(r(\theta)=x_i+R n(\theta)\)로 매개화합니다. 생성점의 미분에는 원호의 강체 이동 항 \(I\)를 포함하고, 끝점 각도는 해당 끝점의 이웃 생성점과 원판의 교점 조건을 미분해 구합니다. 끝점이 \(i,j\)의 이등분선에 있다면 \(\Delta=x_j-x_i\)에 대해

\[
2R\,n(\theta)\cdot\Delta=|\Delta|^2,\qquad
d\theta=\frac{(\Delta-Rn)\cdot d\Delta}{R(Jn)\cdot\Delta},
\]

여기서 \(J(a,b)=(-b,a)\)입니다. 이 끝점 미분으로 \(W_f,N_f\)를 미분하고, \(c_i\)의 미분에는 원호의 정확한 1차 모멘트를 사용합니다. 예를 들어 원호각이 \(\Delta\theta\)인 원형분절의 면적 보정은
\[
\Delta V=\tfrac12R^2(\Delta\theta-\sin\Delta\theta)
\]
이며 도심 계산에는 그에 대응하는 1차 모멘트 보정도 필요합니다. 접선 퇴화로 분모가 0에 가까워지는 경우는 수치를 임의로 자르지 말고 기하 단계 실패로 처리합니다.

### 2. 아핀 필터, 힘·열, 분리 적분

2-hop stencil은 **내부 유체 면 그래프만** 사용합니다. 직접 이웃에는 \(w_{ij}=A_{ij}/d_{ij}^2\), 2-hop에는 \(w_{ik}^{(2)}=h_i\sum_jw_{ij}w_{jk}\)를 적용하고, 진공 표본이나 \(v=0\) 고스트는 넣지 않습니다. 가중 아핀 적합의 잔차 연산자를 \(R_i=I-P_i\)라 하면 원호 모드의 행은

\[
L_{if}=g_{if}R_i,\qquad r_{if}=L_{if}v.
\]

아핀 rank가 부족하면 실제 유체 이웃의 거리 껍질을 확장합니다. 끝까지 rank 3을 얻지 못한 모드는 유효하다고 가장하지 말고 실패로 기록합니다. 힘에는 반드시 \(L_{if}^{T}\)를 사용합니다.

원호 자체를 하나의 HG 모드로 포함합니다. `voro2d_restrict_disc`가 표시한 짧은 원호 조각마다 별도 모드를 만들면 임의의 원호 분할 수에 따라 비선형 소산이 달라지므로, **연속된 원호 조각을 하나의 실제 원호 성분으로 합쳐** 모드를 구성합니다. \(q=|r_{if}|\), \(b=(\gamma+1)h_{if}/4\)일 때 계수·힘·열은

\[
C_{if}(q)=\rho_{if}h_{if}W_f
\left[bq+\sqrt{b^2q^2+c_{if}^2}\right],
\]

\[
F_a^{\rm hg}=-C_{if}r_{if}L_{if,a}^{T},\qquad
H_{i,if}^{\rm hg}=C_{if}r_{if}^{2}\ge0.
\]

모드의 열은 그 모드를 만든 셀 \(i\)에 누적합니다. 그러면 같은 단계에서 \(\sum_a v_a\cdot F_a^{\rm hg}+H_{i,if}^{\rm hg}=0\)입니다. \(L\)이 상수·아핀 속도장을 소거하므로 병진·회전 및 아핀 속도장에서 추가 소산도 0입니다.

`GFS_HG_SPLIT`에서는 이 원호 모드를 기존 `kind=0` 셀 모드처럼 저장할 수 있습니다. \(A_0=\rho_{if}h_{if}W_f\), \(\kappa=\gamma(\gamma-1)\eta_f^{\gamma-1}/m_i\), \(\alpha=\sum_a|L_a|^2/m_a\)라 두면

\[
C=A_0\left[bq+\sqrt{b^2q^2+\kappa U_i}\right],\quad
\dot r=-\alpha Cr,\quad \dot U_i=Cr^2.
\]

따라서 모드 한 번의 갱신은
\[
v_a\leftarrow v_a-\frac{r_0-r_1}{\alpha m_a}L_a,\qquad
U_i\leftarrow U_i+\frac{r_0^2-r_1^2}{2\alpha}
\]
로 에너지를 보존합니다. 기존 로그 진폭 \(s=\log(q_0/q)\) 적분을 재사용하고, 원호 모드를 결정적 순서로 정렬해 재시작 재현성을 유지합니다.

코드 변경은 `sub_build`의 경계 분류를 `PAIR/ARC/기타 cap`으로 나누고, `hg_cell`에서 원호 모드와 내부 면 모드를 별도로 만드는 방식이어야 합니다. capfree 조건은 **일반 cap은 계속 제외하되 정확 기하를 지원하는 원호는 허용**하도록 바꿉니다. `hg_cell`에는 기존 `GFS_BULK_AFFINE` 및 b′ 모드도 있으므로 capfree 해제만으로 이것들까지 켜지 않게 별도 분기해야 합니다. 설계의 \(\Phi=0\)과 맞춰 기존 중심 기준 `sub_eval` 퍼텐셜 경로는 이 확장에 포함하지 않습니다. 관련 현재 경로는 [hg_cell](/tmp/lageunha-push/Exam/exam.c:814), [모드 적분기](/tmp/lageunha-push/Exam/exam.c:13616), [설계 문서](/tmp/lageunha-push/Docs/design_astra.md)에서 확인했습니다.

### 3. 3차 QBP와 충격 처리

현재 `GFS_Q3`는 전 셀에 3차 다항 적합의 기울기를 사용하고 rank 부족 시 2차로 fallback합니다. QBP는 그 기울기로

\[
\delta_f^{QBP}=\Delta u_{\rm cell}
+\frac d2\,e^{T}(G_i+G_j)e,\qquad
m_f=\operatorname{minmod}(\delta_f^{QBP},\Delta u_{\rm cell})
\]

을 만들고, 현재 QLAW=2 분기는 \(Q_f=Z_fm_f^2/\Delta u_{\rm cell}\)을 적용합니다. 3차 fit은 기울기 적합 오차를 줄이지만, 코드가 3차 미분 텐서 보정까지 하는 것은 아닙니다. 이 경로에는 별도 부피 압축 게이트도 없습니다. [Q3 적합](/tmp/lageunha-push/Exam/exam.c:4793), [QBP 입력](/tmp/lageunha-push/Exam/exam.c:8247), [QLAW=2](/tmp/lageunha-push/Exam/exam.c:8297)

문턱 없는 해결안으로 감사 14차의 회전 불변 압축량을 써서 Q 입력을 연속적으로 혼합합니다.

\[
S_i=\tfrac12(G_i+G_i^T),\quad\theta_i=\operatorname{tr}G_i,\quad
\chi_i=
\begin{cases}
\min(1,[-\theta_i]_+/\|S_i\|_F),&\|S_i\|_F>0,\\
0,&\|S_i\|_F=0,
\end{cases}
\qquad \chi_f=\max(\chi_i,\chi_j),
\]

\[
q_{{\rm in},f}=(1-\chi_f)\delta_f^{QBP}
+\chi_f\Delta u_{\rm cell},\quad
m_f=\operatorname{minmod}(q_{{\rm in},f},\Delta u_{\rm cell}),
\]

\[
Q_f=
\begin{cases}
Z_fm_f^2/\Delta u_{\rm cell},&\Delta u_{\rm cell}\ne0,\\
0,&\Delta u_{\rm cell}=0.
\end{cases}
\]

평면 1차원 압축에서는 \(\chi_f=1\)이 되어 Q 입력이 `du_cell`로 돌아오고, 발산 없는 전단에서는 \(\chi_f=0\)입니다. \(Z_f\ge0\)이면 면 열은 계속
\[
H_{i,f}=H_{j,f}=\tfrac12A_fQ_f\Delta u_{\rm cell}
=\tfrac12A_fZ_fm_f^2\ge0
\]
입니다. 이 방식은 조정 가능한 압축 문턱 없이 QLAW=2와 비음수 열을 보존합니다.

요청에 기록된 Sod 128² 결과는 3차판의 최소 간격이 \(t=0.1\)에 \(0.064h\)로 Gate C 기준 \(0.1h\)를 밑돌고, 2차 QBP는 통과했습니다. 따라서 압축 혼합은 결과를 대신하는 승인 근거가 아니라 다음 후보입니다. [audit_impl14_astra.md](/tmp/lageunha-push/Docs/audit_impl14_astra.md)의 순서를 따라 검증해야 합니다.

### 4. 시험 순서

1. 원호 면적·1차 모멘트 및 \(\nabla\xi\)를 유한차분과 비교하고, 원호 분할 수를 바꿔도 결과가 같은지 확인합니다. 아핀 속도 소거, 회전 공변성, \(H\ge0\), 힘–열 항등식을 검사합니다.
2. 각 원호 모드의 분리 적분에서 운동에너지와 셀 내부에너지 합을 확인합니다. 1스레드 \(t\simeq40\) Kepler 재시작에서 같은 상태의 RHS와 모드 목록이 비트 단위로 일치하는지 비교합니다.
3. 동일 \(t\simeq40\) 상태에서 Q2, Q3, 압축 혼합을 분리 비교합니다. \(Q,H_Q,\chi_f\), invalid QBP 면 수와 최소 간격을 기록합니다.
4. Sod를 \(64^2,128^2,256^2\) 순으로 Gate C 기준(최소 간격 \(>0.1h\), 충돌·보정 0, 에너지 오차 기준 포함)으로 판정합니다. 실패하면 다음 장기 시험으로 진행하지 않습니다.
5. 통과한 구성만 Gresho와 Kepler \(t=40\to45\) 짧은 재시작 시험에 적용하고, 마지막에 지정된 장기 고리 시험으로 확장합니다.

요청대로 파일은 수정하지 않았고 계산이나 시험도 실행하지 않았습니다.