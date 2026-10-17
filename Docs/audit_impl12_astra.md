## 결론

**후보 (2)를 권고합니다.** 다만 현재 `bp_assemble`을 그대로 TCF에 연결하는 것으로는 부족합니다. TCF의 `du_rec` 대신 **2차 정합 기울기로 만든 b′형 면 잔차**를 Q의 입력으로 쓰고, 기본 `tcf_Z`도 같은 잔차를 기준으로 계산해야 합니다. 그러면 Q 입력은 회전 불변이고, 매끄러운 선형·2차 속도장에서는 0이 되며, `minmod`가 Q의 부호와 크기를 `du_cell`에 맞춰 제한합니다.

이는 **전체 TCF 해법의 엄밀한 TVD 보장**과는 다릅니다. Sod의 기존 단조 재구성 제한자는 유지하고, Q의 비음수 열 생성과 Gate C 통과를 별도로 확인해야 합니다.

## 코드 판정

- `GFS_LIM_PROJ`는 면 법선 투영에 공통 스칼라 \(\alpha_v\)를 적용하므로 회전 불변이지만, 투영값의 이웃 범위를 벗어나면 가장자리에서 \(\alpha_v=0\)이 될 수 있습니다. 요청 문서의 가장자리 결과와 일치합니다. [투영 제한자](/tmp/lageunha-push/Exam/exam.c:4912)
- 기존 `EUNHA_LIM_SMOOTH`는 이 문제를 해결하지 않습니다. 2차 적합의 잔차를 `0.1` 비율로 판정하고, 곡률로 \(v_x,v_y\) 각각의 범위를 넓힙니다. 그런데 `GFS_LIM_PROJ`는 이 성분별 범위를 사용하지 않고 \(v\cdot e\)의 투영 범위를 직접 제한합니다. 따라서 현재 경로를 켜도 투영 제한은 그대로입니다. 성분별 매끄러움 판정 자체도 회전 불변이 아닙니다. [2차 적합·매끄러움 판정](/tmp/lageunha-push/Exam/exam.c:4766) [범위 확장 경로](/tmp/lageunha-push/Exam/exam.c:260)
- TCF는 재구성으로 \(d u_{\rm rec}\)와 기본 `tcf_Z`를 계산합니다. 이후 \(m=\operatorname{minmod}(d u_{\rm rec},d u_{\rm cell})\)를 Q에 씁니다. 그러므로 Q 입력만 바꾸고 `tcf_Z`가 기존 재구성 속도 점프를 계속 쓰게 두면, 제한자 오차의 영향이 남습니다. [재구성 상태와 `tcf_Z`](/tmp/lageunha-push/Exam/exam.c:7532) [minmod와 Q](/tmp/lageunha-push/Exam/exam.c:8098)
- 현재 b′ 경로는 affine-exact한 2-hop 선형 QR 기울기를 모아 공유 면 연산자를 만듭니다. `bp_assemble`은 TCF 면 루프에서 직접 쓰는 값이 아니라 별도 분할 감쇠 모드이며, affine 입력에 대한 검사는 있지만 가장자리의 2차 정확도 검사는 없습니다. 따라서 이를 TCF에 그대로 켜기보다, **제한 전 2차 적합 기울기**를 TCF에서 쓸 수 있게 보관해야 합니다. [b′ 경로 설명](/tmp/lageunha-push/Exam/exam.c:760) [면 조립](/tmp/lageunha-push/Exam/exam.c:13539)

## 권고 수식

면 \(i\to j\)에서 \(d=|\mathbf x_j-\mathbf x_i|\), \(\mathbf e=(\mathbf x_j-\mathbf x_i)/d\)라 두고,

\[
\Delta u_{\rm cell}=(\mathbf v_i-\mathbf v_j)\cdot\mathbf e,
\qquad
\delta_{b'}=
\Delta u_{\rm cell}
+\frac d2\,\mathbf e^\mathsf T(\widehat G_i+\widehat G_j)\mathbf e .
\]

여기서 \(\widehat G_i,\widehat G_j\)는 **제한 전 2차 정합 적합의 중심 기울기**입니다. 속도를 각 성분별로

\[
v_q(\mathbf x_i+\mathbf r)\approx v_{q,i}
+\nabla v_{q,i}\cdot\mathbf r
+\tfrac12\mathbf r^\mathsf T H_{q,i}\mathbf r
\]

로 적합하고, 선형·이차 항을 함께 풀어 중심 기울기를 얻습니다. 원 코드의 b′ 기울기처럼 선형 적합만 쓰면 한쪽 stencil에서 곡률 편향이 남을 수 있어, `EUNHA_LSF_ORDER=2`와 같은 2차 기저가 필요합니다.

면 방향의 매끄러운 속도를 \(u(s)=\mathbf v(\mathbf x_i+s\mathbf e)\cdot\mathbf e\)라 하면, 정확한 중심 기울기에서

\[
\delta_{b'}=u(0)-u(d)+\frac d2\bigl(u'(0)+u'(d)\bigr)
=\frac{d^3}{12}u'''(0)+O(d^4).
\]

따라서 선형·이차 변화는 상쇄됩니다. 회전해도 \(\Delta u_{\rm cell}\)과 \(\mathbf e^\mathsf T G\mathbf e\)가 불변 스칼라로 유지됩니다.

Q 입력은

\[
m_{b'}=\operatorname{minmod}(\delta_{b'},\Delta u_{\rm cell}).
\]

기본 법칙이면 \(Q=Z_{b'}m_{b'}\), `GFS_TCF_QLAW=2`이면

\[
Q=Z_{b'}\frac{m_{b'}^2}{\Delta u_{\rm cell}}
\quad(\Delta u_{\rm cell}\ne0).
\]

\(Z_{b'}\ge0\)일 때 면 열 \(H_f^Q=\tfrac12 A_f Q\Delta u_{\rm cell}\)은 두 법칙 모두 음수가 아니며, \(|m_{b'}|\le|\Delta u_{\rm cell}|\)입니다.

기본 `tcf_Z` 경로도 같은 \(\delta_{b'}\)를 사용해 일관되게 계산해야 합니다. \(\Phi\)를 코드의 `tcf_pstar`라 하고 \(\bar u\)를 두 재구성 법선속도의 평균으로 두면,

\[
Z_{b'}=
\frac{\Phi(\rho_L,p_L,\bar u+\delta_{b'}/2,\rho_R,p_R,\bar u-\delta_{b'}/2)
-\Phi(\rho_L,p_L,\bar u,\rho_R,p_R,\bar u)}
{\delta_{b'}}.
\]

\(\delta_{b'}=0\)에서는 \(\partial_{\delta}\Phi|_0\)를 씁니다. `GFS_TCF_FACE_Z` 또는 `GFS_TCF_ZH`가 활성화된 경우에는 현재 코드가 `mm`에서 별도 \(Z_q\)를 만들므로, 해당 설정을 고정하고 새 \(m_{b'}\) 경로를 비교해야 합니다.

구현상 제한된 `stress.gU`를 재사용하면 안 됩니다. 현재 문제처럼 가장자리 기울기가 이미 0으로 제한된 뒤에는 b′도 `du_cell`로 퇴화합니다. 2차 적합의 원시 기울기를 제한 전에 캐시하고, 면 루프에서 사용해야 합니다. 필요한 2차 적합이 성립하지 않는 가장자리 셀이 있으면 그 셀의 처리를 정하기 전까지 설계를 채택하지 않는 편이 안전합니다.

## 사례별 예상과 시험 순서

- **Kepler 가장자리:** \(|\delta_{b'}|\ll|\Delta u_{\rm cell}|\)이면 \(m_{b'}=\delta_{b'}\)입니다. 같은 \(Z\)를 가정할 때 기본 Q 법칙의 열은 기존 대비 대략 \(|\delta_{b'}|/|d u_{\rm rec}|\) 비율로 줄고, QLAW 2에서는 그 비율의 제곱으로 줄어듭니다. 기존 가장자리 \(H_Q/U\)인 약 1.02와 고리 2의 약 0.85가 얼마나 내려갈지는 2차 기울기 적합 오차와 실제 \(Z\)에 달려 있어 자료만으로 절대값을 예측할 수 없습니다. 진단에 \(\delta_{b'}/\Delta u_{\rm cell}\), \(Z_{b'}\), 면별 \(H_Q\)를 함께 기록해야 합니다.
- **Sod:** 초기 불연속에서 기울기가 0인 셀은 \(\delta_{b'}=\Delta u_{\rm cell}\)이므로 처음부터 Q가 사라지지는 않습니다. 이후 충격부의 적합 기울기가 Q를 줄일 수 있으므로 Gate C의 과도·미달과 압력·속도 단조성을 확인해야 합니다. 기존 MUSCL 면 상태 제한자는 그대로 둡니다.
- **Gresho:** 매끄러운 영역에서는 국소 선형·이차 성분의 Q가 줄 것으로 예상합니다. 프로파일의 전이 영역에서는 잔차가 남을 수 있으므로 소용돌이 감쇠와 열 발생을 함께 비교해야 합니다.

순서는 요청대로 **회전 불변 RHS 비교(\(\pi/7\)) → Kepler 가장자리 RHS → Sod Gate C → Gresho → Kepler 장기 진화**입니다. 모든 비교에서 Q 법칙과 face-Z 설정을 고정하고, 입자 대응 기준으로 \(d u_{\rm rec},\delta_{b'},Z,Q,H_Q\)를 기록하십시오.

검토는 읽기 전용으로 했으며 파일은 수정하지 않았습니다.