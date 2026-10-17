**Astra 판정: §2–3의 반이산 정식은 승인하지만, 3판 전체 서명은 보류합니다.** 엔탈피·엔트로피 장부, 정지 간격, 글래스 장부, CFL 상한, G2와 Hill/Kepler 관문은 반영됐습니다. 남은 사항은 다음과 같습니다.

1. **“정확한 Hugoniot 임피던스”라는 명칭이 틀립니다**([§6](/tmp/lageunha-push/Docs/subzonal_derivation.md:55)). 현재 δ는 좌우 전체 속도 점프입니다. 동일 상태의 대칭 두 충격에서는 한쪽 충격의 속도 변화가 δ/2이므로
   \[
   p_*-p_0=\rho\frac{\delta}{2}
   \left[\frac{\gamma+1}{8}\delta+
   \sqrt{\left(\frac{\gamma+1}{8}\delta\right)^2+c^2}\right].
   \]
   [정확 Riemann 충격식](https://www.astro.princeton.edu/~jstone/Athena/doxygen/html.with_source/exact_8c_source.html)에서 유도하면 이렇습니다. 제안 Q는 엔트로피 부호를 지키는 **Hugoniot형 인공점성**으로 사용할 수 있지만, 정확한 두 충격 해는 아닙니다. 강한 극한의 ρδ²도 Z가 아니라 Q의 척도입니다.

2. **1D 결과를 확대 서술했습니다**([§6](/tmp/lageunha-push/Docs/subzonal_derivation.md:53)). 도약 오차·울림 \(10^{-4}\)는 Noh 결과입니다. [Fable의 Sod 표](/tmp/lageunha-push/Docs/subzonal_derivation_review2_fable.md:29)는 (나)의 편향 −0.3%, 울림 0.2%입니다. 또한 원 보고의 비율은 **Φ/U_cell**이며, 본문의 Φ/U와 분모가 다릅니다. 시험별로 구분해야 합니다.

3. **고정 계수 퍼텐셜에 대한 단정이 잘못됐습니다**([§5](/tmp/lageunha-push/Docs/subzonal_derivation.md:47)). c²를 엔트로피와 무관한 상수로 영구 고정하고 전체 기하 미분을 포함하면 보존적인 다른 모형이 됩니다. 보존을 깨는 것은 상태에 따라 c²를 갱신하면서 그 미분·에너지 변화를 누락하는 경우입니다. 이는 제 2판 검토의 구별과 다릅니다.

4. **실험한 Q와 새 권고 Q를 구별해야 합니다.** [사슬 코드](/tmp/lageunha-push/Docs/subzonal_checks/chain1d.py:84)는 \(Q=Z_W\delta\)를 시험합니다. `max(Z_eff,Z_W)`와 `QLAW=2`의 결합은 별도 검증 대상입니다. 부호 정합성은 유지되지만, 기존 수치 결과가 그 조합의 충격 폭·가열·시간간격까지 검증한 것은 아닙니다.

5. **Mach 표기와 일반화도 수정해야 합니다**([§6](/tmp/lageunha-push/Docs/subzonal_derivation.md:54)). 보고된 \(P_0=10^{-10},\,\rho=1,\,|u|=1,\,\gamma=5/3\)이면 유입 Mach는 약 **\(7.7\times10^4\)**입니다. “Mach·해상도에 무관”은 “시험한 조건에서 최소 간격이 같았다”로 한정하십시오.

문서와 저장된 코드를 대조했으며, 파일 변경이나 충격 실험 재실행은 하지 않았습니다.