**Astra 판정: G1 재서명 보류입니다.** 핵심 정정은 상당 부분 맞지만, 계획서에 보존을 깨뜨릴 구현 지시와 정정 전 주장이 남아 있습니다. 파일은 수정하지 않았습니다.

요청한 **반대칭 총에너지 교환식의 부호·인자는 정확합니다.** \(W=0\)에서 \(w_f=(u_i+u_j)/2+(v_i-v_j)\cdot t/d\)를 대입하면
\[
v_i\cdot F_{i,f}-P_iAw_f+\frac12AQ(u_i-u_j)
=\dot E_{i,f}^{\rm 교환}-AP_i u_i.
\]
따라서 원호 일을 포함한 셀 합계는
\[
\dot E_i
=\sum_f\dot E_{i,f}^{\rm 교환}
-P_iv_i\cdot\left(\sum_f Ae+N_i\right)
=\sum_f\dot E_{i,f}^{\rm 교환}.
\]
여기서
\[
\boxed{\dot E_{i,f}^{\rm 교환}
=-A\left[
\frac{P_j u_i+P_i u_j}{2}
+\frac{(P_jv_i-P_iv_j)\cdot t}{d}
+Q\frac{u_i+u_j}{2}
\right]}
\]
이며, 양쪽 셀을 바꾸면 정확히 반대 부호가 됩니다.

남은 정정 사항은 다음과 같습니다.

1. **계획서 §2의 원호 에너지 지시는 이중 계산을 유발합니다.**  
   [계획서 72행](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:72)은 위 교환식을 사용하면서 “부피 변화와 에너지에 \(v_i\cdot N_i\) 항을 넣는다”고 합니다. 위 유도처럼 **이 교환식에는 원호 일이 이미 흡수되어 있습니다.** 별도 \(-P_iv_i\cdot N_i\)를 총에너지 우변에 더하면 안 됩니다. 원호 항은 기하적 \(\dot V_i\)와 내부에너지·엔트로피 진단에 포함하십시오. [유도서 §2의 원호 힘·부피 정정](/tmp/lageunha-push/Docs/tcf_derivation.md:25)은 맞습니다.

2. **같은 구현 지시가 Laguerre 힘과 Voronoi 에너지식을 혼용합니다.**  
   “Laguerre면 \(\widetilde P\)”로 힘을 바꾸면서 에너지식은 \(W=0\) 그대로입니다. 공유 면이 일치하고 \(W\)가 고정이면 올바른 식은
   \[
   \dot E_{i,f}=-A\left[
   \phi P_j u_i+(1-\phi)P_i u_j
   +\frac{(P_jv_i-P_iv_j)\cdot t}{d}
   +Q\frac{u_i+u_j}{2}
   \right].
   \]
   **2단계를 \(W=0\)으로 제한하거나 이 식까지 함께 적용해야 합니다.**

3. **\(Q\)의 정정은 맞지만, 보장 범위가 본문에서 다시 확대됩니다.**  
   [유도서 51–54행](/tmp/lageunha-push/Docs/tcf_derivation.md:51)에 HLLC 전체 단조성과 “충격에서는 … HLLC 수준의 소산”이라는 단정이 남아 있습니다. 실제로는
   \[
   Q=Z_{\rm eff}\minmod(\Delta u_{\rm rec},\Delta u_{\rm cell}),
   \qquad Z_{\rm eff}\ge0
   \]
   일 때 반이산 엔트로피 부호가 보장됩니다. 부호 불일치에서는 \(Q=0\)이므로 충격 소산도 보장되지 않습니다. \(0\) 부근 연속성에는 \(Z_{\rm eff}\)의 유한하고 연속적인 극한과 영분모 처리가 필요합니다. “매끄럽게 보간하면 잡음이 점프에 들어가지 않는다”도 지나친 표현입니다. 선형장은 정확히 상쇄되지만 일반적인 매끄러운 장에는 재구성 오차가 남습니다.

4. **동적 \(W\) 조건은 정정문과 기존 문장이 충돌합니다.**  
   [유도서 65행](/tmp/lageunha-push/Docs/tcf_derivation.md:65)의 “에너지가 보존되려면 각 \(l\)에 대해 \(g_l=0\)이어야 한다”는 필요조건 주장을 삭제하거나 별도의 정류 조건으로 명시해야 합니다. 실제 조건은
   \[
   \dot E=-g\cdot\dot W,\qquad g\cdot\dot W=0.
   \]
   원판+\(W\ne0\)의 공유 면 불일치 경고는 맞습니다. 다만 [68행의 정적 Laguerre 결론](/tmp/lageunha-push/Docs/tcf_derivation.md:68)에도 **공유 면이 일치하는 유효한 분할**이라는 조건을 붙여야 합니다.

5. **계획서에 폐기된 설계와 시험 기준이 남았습니다.**  
   [서두](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:3)의 “중심력” 목표, [§0](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:21)의 격자 잡음까지 포함한 소산 소멸 주장은 정정이 필요합니다. [§1.5](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:53)의 “단일 압력 유지·합이 0인 열 재배분·셀 점프 미사용”도 현재 변분 힘과 minmod 구현에 맞지 않습니다. [G2](/tmp/lageunha-push/Docs/plan_thermo_compatible_flux.md:77)는 입자=생성점 조건에서 유체 정지와 임의 생성점 운동을 동시에 요구할 수 없으며, 일반 RK4의 엔트로피 오차를 반올림 수준으로 요구해서도 안 됩니다. 반이산 항등식과 시간간격 수렴 시험을 구분해야 합니다.

**VPH 정확성 진단 정정은 승인합니다.** 내부 셀의 균일·선형 압력 정확성과 전단 부피 잡음·Nyquist 모드 위험을 구분한 [유도서 §8](/tmp/lageunha-push/Docs/tcf_derivation.md:84)은 맞습니다. **계획서 §1.4의 두 번째 항 ‘+’, \(\Pi_{\rm ec}\), \(\sigma_f=\sum T\dot S\)라는 정정도 정확합니다.**