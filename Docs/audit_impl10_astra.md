## 결론

요청서의 병진 시험은 **갈릴레이 불변성 위반이나 가역 엔트로피 생성의 증거가 아닙니다.** 코드의 `dte`는 총에너지율 \(\dot E_i\)이고, 총에너지율은 균일 속도 부스트 \(\mathbf V\)에서 \(\mathbf V\cdot\mathbf F_i\)만큼 변하는 것이 맞습니다. 내부에너지율은 `dte − v·F`로 계산해야 합니다. TCF의 이 구분과 부호는 [힘·에너지·체적률 계산부](/tmp/lageunha-push/Exam/exam.c:8062), [RK4의 총에너지율 처리](/tmp/lageunha-push/Exam/exam.c:6243)에 나타납니다.

따라서 닫힌 주기 셀에서 `dte`를 \(-P_i\dot V_i\)로 바꾸면 현재 총에너지 변수의 의미가 깨집니다. 물리식으로 맞는 내부에너지율은 \(-P_i\dot V_i+H_i^Q\)이며, 총에너지율은 운동에너지 일을 더한 \(v_i\cdot F_i-P_i\dot V_i+H_i^Q\)입니다. 이는 [design_astra.md의 에너지 장부](/tmp/lageunha-push/Docs/design_astra.md:320)와 같습니다.

## 판정에 필요한 장부

면 \(f=(i,j)\)에서 \(e=(x_j-x_i)/d\), \(t\)는 면 중심에서 입자 중점까지의 접선 오프셋, \(u_i=v_i\cdot e\)라 두면 코드는

\[
\dot V_{if}
=A_f\left[\frac{u_i+u_j}{2}
+\frac{(v_i-v_j)\cdot t}{d}\right],
\]

\[
F_{if}^{\rm rev}
=-A_f\left[\frac{P_i+P_j}{2}e
+(P_j-P_i)\frac{t}{d}\right],
\]

\[
\dot E_{if}^{\rm rev}
=-A_f\left[
\frac{P_j u_i+P_i u_j}{2}
+\frac{(P_jv_i-P_iv_j)\cdot t}{d}
\right].
\]

닫힌 셀에서는 면 법선 폐합 \(\sum_f A_f e_f=0\)을 쓰면 면들을 합산한 뒤

\[
\dot E_i^{\rm rev}=v_i\cdot F_i^{\rm rev}-P_i\dot V_i,
\qquad
\dot U_i^{\rm rev}=\dot E_i^{\rm rev}-v_i\cdot F_i^{\rm rev}
=-P_i\dot V_i
\]

가 됩니다. 요청서의 측정량 \((dte+P_i\dot V_i)/U_i\)는 이 조건에서 \(v_i\cdot F_i/U_i\)입니다. UX에 비례하는 값은 바로 이 운동에너지 일 항입니다. **엔트로피율을 검사하려면** \(dte-v_i\cdot F_i\)를 계산해야 합니다.

실제 `SEDOV_PHASE1` 경로는 `die=dte`를 저장하고 RK4 단계에서 \(m_i v_i\cdot a_i\)를 빼 내부에너지를 갱신합니다([저장부](/tmp/lageunha-push/Exam/exam.c:8589), [단계 계산](/tmp/lageunha-push/Exam/exam.c:6243)). 그러므로 그 경로에 `die=-P_i\dot V_i`를 넣으면 운동에너지 일을 잘못 처리합니다.

부호도 명확히 해야 합니다. 코드의 힘 부호에서는

\[
\sum_i v_i\cdot F_i^{\rm rev}
=+\sum_i P_i\dot V_i,\qquad
\sum_i\dot U_i^{\rm rev}
=-\sum_i P_i\dot V_i.
\]

따라서 총에너지 보존식은 두 항의 합이 0입니다. 요청서의 \(\sum_iP_i\dot V_i=-\sum_iF_i\cdot v_i\)는 실제 코드의 힘 부호 기준으로는 반대 부호입니다.

Q는 별도 소산 장부로 유지하면 됩니다. 면 Q 힘과 열은

\[
F_{i,f}^Q=-A_fQ_fe,\qquad
H_{i,f}^Q=\frac12 A_fQ_f(u_i-u_j).
\]

QLAW=2이면 \(Q_f=Z_fm_f^2/(u_i-u_j)\)이므로 \(H_{i,f}^Q=\frac12A_fZ_fm_f^2\ge0\)입니다. 총에너지율의 Q 항은 \(v_i\cdot F_{i,f}^Q+H_{i,f}^Q=-\frac12A_fQ_f(u_i+u_j)\)로 코드의 `dte` Q 항과 일치합니다([Q·열 계산](/tmp/lageunha-push/Exam/exam.c:8103)). 가장자리 Q 가열은 이 열 \(H_Q\)로 따로 판정해야 합니다.

## 원호와 cap

원호 셀의 체적률은 내부 면 항에 경계 항을 더해 계산합니다.

\[
\dot V_i
=\sum_{f\in{\rm faces}}\dot V_{if}
+\sum_{\Gamma_i}v_\Gamma\cdot N_\Gamma,
\qquad
N_\Gamma=\int_\Gamma n\,ds.
\]

고정 반경 원판에 붙어 움직이는 원호의 경우

\[
N_{\rm arc}
=R(\sin\theta_1-\sin\theta_0,\ \cos\theta_0-\cos\theta_1),
\qquad
\dot V_{\rm arc}=v_i\cdot N_{\rm arc}.
\]

코드는 원호를 짧은 현으로 저장하고, `Area2DPolygon`에서 원형 분절 면적을 보정합니다([원판 제한](/tmp/lageunha-push/Exam/exam.c:3938), [면적 보정](/tmp/lageunha-push/Voro/voro.c:1134)). 현 길이와 바깥 법선의 곱은 원호의 \(\int n\,ds\)와 같으므로, 원호가 생성점과 강체 이동할 때 `tcf_dV`의 해당 항은 기하 폐합이 맞으면 정확한 체적률이 됩니다([cap·arc 체적률](/tmp/lageunha-push/Exam/exam.c:8409)). cap도 생성점과 함께 이동하는 기하 경계라면 같은 방식으로 \(v_i\cdot N_{\rm cap}\)를 더합니다.

TCF의 쌍 힘은 셀 경계 폐합을 이용해 원호·cap의 압력 효과를 흡수할 수 있습니다. 따라서 TCF 경로에 경계 `−P_i \dot V` 항을 총에너지율로 또 더하면 안 됩니다. 현재 scheme-H cap 에너지 항은 TCF와 함께 쓸 수 없도록 코드가 제한합니다([TCF 제약](/tmp/lageunha-push/Exam/exam.c:7491), [scheme-H cap 항](/tmp/lageunha-push/Exam/exam.c:8415)).

필수 검증은 면 폐합과 체적률을 함께 보는 것입니다.

\[
R_i^{\rm rev}
=\dot E_i^{\rm rev}-v_i\cdot F_i^{\rm rev}
+P_i\dot V_i=0.
\]

원호·벽 기하 미분은 유한차분과 상대 \(10^{-6}\) 이내로 비교하고, 주기 영역은 활성 cap 0을 확인하십시오([기하 관문](/tmp/lageunha-push/Docs/design_astra.md:557)). `upperrelated<0`는 원호와 일반 cap을 모두 포함하므로 덤프에서 두 종류를 분리해야 합니다. `voro2d_restrict_disc`는 입력 꼭짓점 512개와 출력 `mp`를 넘으면 조용히 잘라내므로, 이 경우는 기하 관문 실패로 처리해야 합니다([클리핑 경계](/tmp/lageunha-push/Exam/exam.c:3991)).

## RK4와 시험 순서

각 RK4 단계에서 같은 단계 상태로 기하, \(F_i\), \(\dot V_i\), \(H_i^Q\)를 계산하십시오. 내부에너지를 직접 적분하는 경로라면

\[
\dot U_i=-P_i\dot V_i+H_i^Q,\qquad
\dot E_i=v_i\cdot F_i+\dot U_i
\]

를 사용하고, 단계 코드가 이를 다시 \(v_i\cdot F_i\)만큼 빼지 않게 해야 합니다. 현재처럼 `dte`를 총에너지율로 넘기는 경로는 `dte`를 유지하십시오. 모든 기하·에너지 누적은 double 정밀도로 하고, 단계 끝에서 한 번만 `postype`으로 변환하십시오. direct-U 경로의 RK4 전역 오차는 고정 구간에서 4차여야 하며, dt를 반으로 줄였을 때 오차가 약 16분의 1로 감소하는지 확인해야 합니다. 총에너지를 상태변수로 직접 적분하면 전역 총에너지 합은 단계별 상쇄로 반올림 수준까지 보존할 수 있지만, Kepler처럼 운동에너지와 열에너지의 크기 차가 큰 경우 \(E-K\) 복원 정밀도를 따로 확인해야 합니다.

권고 시험 순서는 다음과 같습니다.

1. Q=0인 닫힌 주기 RHS에서 \(R_i^{\rm rev}\), \(\sum_i dte_i\), 면 폐합을 검사합니다.
2. UX=0/10/100 시험을 새 갈릴레이 관문으로 추가합니다. 힘과 \(\dot V_i\)는 부스트 전후 같고, `dte` 변화는 \(V\cdot F_i\), 내부에너지율과 \(\dot{\ln K_i}\) 변화는 0이어야 합니다. `dte+P\dot V`가 0이어야 한다는 기존 기준은 사용하지 마십시오.
3. 원판 원호와 cap의 체적률을 유한차분으로 검증한 뒤 RK4 시간 간격 수렴을 확인합니다.
4. 동일한 Q 설정으로 Sod를 재실행하고 그다음 Q=0 Gresho를 비교합니다. 에너지 장부를 바로잡는 것만으로는 반주기 결과가 달라질 이유가 없습니다. 인수인계의 Sod 면-Z/minmod-Q 구성은 통과했지만 QLAW=2 속도 \(L_1=0.0536\)은 실패했고, Gresho에서는 충격과 전단을 함께 통과한 구성이 없었습니다([Sod 결과](/tmp/lageunha-push/Docs/HANDOVER_2026-10-04.md:67), [Gresho 결과](/tmp/lageunha-push/Docs/HANDOVER_2026-10-04.md:92)). Q 법칙 변경과 에너지 장부 변경을 한 번에 묶지 마십시오.
5. 실제 Kepler 내부에너지/엔트로피 변화를 `dte−v·F` 및 \(H_Q\)와 대조합니다. Gresho의 기존 압력·부피 오차 문제는 이 에너지율 교체로 해결되지 않습니다([원인 진단](/tmp/lageunha-push/Docs/HANDOVER_2026-10-04.md:108)).

따라서 가장자리의 큰 Q 열은 별도 소산 문제입니다. 제시된 병진 수치와 `dte+P\dot V`만으로는 가장자리의 가역 엔트로피 생성이 입증되지 않았습니다. 읽기 전용으로 검토했으며 파일 변경이나 시뮬레이션은 하지 않았습니다.