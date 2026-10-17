**결론: (b′)는 다음 시험 후보로 채택하되, 현 면 Q의 대체안으로는 아직 승인하지 않습니다.** 아핀장을 정확히 소거하면서 힘과 열의 장부를 닫을 수 있습니다. 그러나 2-hop 적합이 실제 쌍 접근을 얼마나 남기는지는 검증되지 않았습니다. 우선 **bulk + HG + (b′), 현 TCF Q=0**으로 시험하고, 같은 조건의 **bulk + (b′), HG=0**을 비교하십시오. 현재 기록에서 bulk+HG는 2D Mach 10 최소 간격 \(0.062h\), Sod \(0.035h\)로 실패했고, 현 면 Q는 쌍을 막는 대신 Gresho를 감속시켰습니다. [관문 기록](/tmp/lageunha-push/Docs/design_astra_audit.md)

### 구현할 법칙

공유 면 \(f=(i,j)\)을 **한 번만** 처리하고, \(e_f=(x_j-x_i)/d_f\)로 둡니다. 기존 HG의 2-hop 가중 QR에서 얻은 기울기를 속도에 대한 선형 연산자로 저장해
\[
\delta_f=(v_i-v_j)\cdot e_f-d_f\,e_f^{T}\frac{G_i+G_j}{2}e_f
          =\sum_a\ell_{fa}\cdot v_a
\]
를 구성하십시오. 힘에는 수치적으로 구한 그 \(\ell_f\)의 **정확한 전치**를 사용합니다.

첫 시험은 양방향 잔차 감쇠로 정합니다.
\[
C_f=A_fZ_f\ge0,\qquad
F_a^{b'}=-\sum_f C_f\delta_f\ell_{fa},\qquad
H_i^{b'}=H_j^{b'}=\tfrac12 C_f\delta_f^2.
\]
그러면 \(\sum_a v_a\cdot F_a^{b'}+\sum_iH_i^{b'}=0\), 셀별 \(H_i^{b'}\ge0\)입니다. 완전 순위 QR가 아핀장을 재현하면 \(\ell_f\)는 병진·회전을 소거하므로 총 힘과 토크도 0입니다. 이 열은 *내부에너지*에 넣어야 합니다.

**제안서의 \(\delta_f<0\)은 압축 방향이 아닙니다.** 위 \(e_f\) 정의에서는 \((v_i-v_j)\cdot e_f>0\)이 접근입니다. 잔차에 단방향 조건을 시험한다면 \(\delta_f>0\)이 맞지만, 아핀 예측을 뺀 뒤의 부호가 실제 두 입자의 접근 부호와 반드시 같지는 않습니다. 따라서 초기 후보를 양방향으로 두고, 단방향형은 쌍 시험과 Gresho 결과가 나온 뒤 별도 비교하십시오.

이 항은 기존 TCF의 쌍 힘 \(F_i=-AQe,\ F_j=+AQe\)에 \(Q=Z\delta_f\)를 대입하는 방식으로 구현하면 안 됩니다. 그 방식의 운동에너지 손실은 \(AZ\delta_f d_u\)여서 음수가 될 수 있습니다. 현 Q의 힘과 \(AQd_u/2\) 열을 함께 끄고, 위 비국소 전치 힘과 열을 별도 모드로 더해야 합니다. [현 Q 경로](/tmp/lageunha-push/Exam/exam.c:7963)

### 임피던스와 HG

\(Z_f\)의 **밀도·압력·음속은 현 공유 면의 양쪽 순간 상태**를 유지하십시오. 특히 \(\rho_{if}=m_iB_{if}/(V_iW_{if})\)를 셀 밀도로 되돌리지 않습니다. 속도 인자만 새 모드에 맞춰 \(\delta_f\)로 바꾼 첫 시험식을 권합니다.
\[
\mu_f=\frac{\rho_{if}\rho_{jf}}{\rho_{if}+\rho_{jf}},\quad
c_M=\max(c_{if},c_{jf}),\quad
a_H=\frac{\gamma+1}{8}|\delta_f|,
\]
\[
Z_f=\max\!\left[
\mu_f\{c_M+\max(\delta_f,0)\},
\ \mu_f\left\{a_H+\sqrt{a_H^2+c_M^2}\right\}
\right].
\]
이는 기존 면 상태와 Hugoniot 규모를 보존한 **시험식**입니다. 현 Q의 `mm`을 그대로 쓰면 새 잔차가 커져도 강충격 임피던스가 커지지 않습니다. 면이 사라질 때 \(A_fZ_f\)가 유계 또는 0으로 가는지, 가까워지는 두 셀의 면 상태가 실제로 충분한 제동을 주는지도 계측해야 합니다. [면 상태 정의](/tmp/lageunha-push/Docs/design_astra.md:95), [현 임피던스 구현](/tmp/lageunha-push/Exam/exam.c:7980)

(b′)와 HG는 모두 아핀장을 소거하지만 같은 연산자는 아닙니다. HG는 아영역 형상 잔차 \(g_{if}R_i v\), (b′)는 면 법선 점프에서 두 기울기를 뺀 잔차입니다. 따라서 HG를 즉시 대체한다고 단정할 수 없습니다. 반대로 둘의 열이 중복되어 Gresho를 다시 감속시킬 수 있으므로, 열을 **bulk·HG·(b′)**로 분리해 위 두 조합을 비교하십시오. 아핀 충격 압축은 \(\delta_f=0\)이므로 bulk는 유지해야 합니다. [HG 연산자](/tmp/lageunha-push/Exam/exam.c:769), [5차 bulk 설계](/tmp/lageunha-push/Docs/audit_impl5_astra.md)

분리 적분은 현 HG 모드 흐름을 확장할 수 있습니다. 고정 기하에서
\[
\alpha_f=\sum_a\frac{|\ell_{fa}|^2}{m_a},\qquad
\dot\delta_f=-\alpha_f A_fZ_f\,\delta_f,
\]
\[
U_i(\delta)=U_{i0}+\frac{\delta_0^2-\delta^2}{4\alpha_f},
\quad
U_j(\delta)=U_{j0}+\frac{\delta_0^2-\delta^2}{4\alpha_f}.
\]
이 스칼라 흐름으로 \(\delta_1\)을 얻은 뒤
\[
\Delta v_a=-\frac{\delta_0-\delta_1}{\alpha_fm_a}\ell_{fa},\qquad
\Delta U_i=\Delta U_j=\frac{\delta_0^2-\delta_1^2}{4\alpha_f}.
\]
\(Z_f\)가 \(|\delta|\)와 갱신되는 \(U_i,U_j\)에 의존하므로 초기 \(Z_f\)를 고정한 단순 지수 감쇠는 정확한 모드 흐름이 아닙니다. [기존 분리 흐름](/tmp/lageunha-push/Exam/exam.c:13327)

### 시험 순서와 탈락 조건

1. **2D 기준 코드:** 불규칙 격자에서 병진·회전·전단·균일 압축의 \(\delta_f\), 힘, 열이 반올림 수준인지 확인합니다. 무작위 속도에서 에너지·총 힘·토크 장부를 검사한 다음, Mach 1·10·100·1000 쌍, 질량비 1·8·64, 체커보드와 Nyquist를 실행합니다. 특히 기존 실패 조건인 Mach 10에서 최소 간격 \(>0.1h\)와 지속적 dt 붕괴 없음이 선행 조건입니다. `audit_dyn2d.py`는 현재 bulk와 현 Q·HG 구성까지 있으므로 새 연산자와 열을 별도 경로로 추가해야 하는 *시험 명세*로 사용하십시오. [기준 스크립트](/tmp/lageunha-push/Docs/design_audit_scripts/audit_dyn2d.py)
2. **C 구현:** 공유 면당 한 모드, 양쪽 QR의 주기 영상 정합, 전역 입자 ID로 합친 \(\ell_f\), 정확한 전치, 열의 양쪽 분배를 확인합니다. 순위 결손·누락 면·소멸 면을 보고하고, 한 RHS 및 시간 간격 반분 에너지 시험을 통과시킵니다.
3. **Sod 128²:** 현 면 Q+HG와 HLLC의 동일 초기조건으로 \(t=0.1\), \(L_1\), 최소 간격 \(>0.1h\), dt 붕괴, 엔트로피 산포를 비교합니다. 쌍 또는 Sod가 실패하면 (b′)를 Q 대체안으로 기각합니다.
4. **Gresho:** 충격 관문을 통과한 조합만 \(64^2,128^2,256^2\)와 dt 반분으로 평가합니다. \(128^2\)의 \(L_1(v_\phi)<0.05\), 회전핵 평균 오차 5% 이하, 소산 열의 해상도 수렴이 필요합니다. 충격 통과만으로 전단 감속을 승인하지 않습니다. [승인 기준](/tmp/lageunha-push/Docs/design_astra.md:641)

읽기 전용으로 문서와 코드 경로를 검토했습니다. (b′)의 쌍 제동과 전단 정확도는 아직 실행 결과가 없으므로 미판정입니다.