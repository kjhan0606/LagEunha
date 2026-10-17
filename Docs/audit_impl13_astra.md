## 결론

**(a)로 2차 적합을 먼저 복구하고, 끝내 복구되지 않는 셀에는 (b)를 적용하십시오.** 각 셀의 제한 전 2차 기울기가 유효하면 그 값을, 2차 적합이 끝내 성립하지 않으면 그 셀 자신의 제한 전 선형 LSF 기울기를 QBP에 씁니다. 선형 적합도 완전한 기울기를 식별하지 못하면 관문을 실패시킵니다. 현재의 `qin = ducell` 대체는 없애야 합니다.

면 \(i\to j\)에서 \(d=|\mathbf x_j-\mathbf x_i|\), \(\mathbf e=(\mathbf x_j-\mathbf x_i)/d\), \(u_i=\mathbf v_i\cdot\mathbf e\)라 두면:

\[
\Delta u_{\rm cell}=u_i-u_j,\qquad
\delta_{b'}=\Delta u_{\rm cell}
+\frac d2\,\mathbf e^\mathsf T(G_i+G_j)\mathbf e,
\qquad
m=\operatorname{minmod}(\delta_{b'},\Delta u_{\rm cell}).
\]

여기서 \(G_i\)는 셀 \(i\)의 제한 전 기울기입니다. 2차 적합 성공 시 \(G_i=G_i^{(2)}\), 실패 시 완전한 선형 LSF에서 \(G_i=G_i^{(1)}\)로 선택합니다. \(v(\mathbf x)=A\mathbf x+b\)이면 \(\Delta u_{\rm cell}=-d\,\mathbf e^\mathsf T A\mathbf e\)이고 양쪽 기울기가 \(A\)이므로 \(\delta_{b'}=0\)입니다. 식의 모든 항이 회전 불변 스칼라이고, 선형 fallback도 아핀 정확성을 유지합니다.

기본 법칙에서는 \(Q=Z m\), `GFS_TCF_QLAW=2`에서는 \(\Delta u_{\rm cell}\ne0\)일 때 \(Q=Z m^2/\Delta u_{\rm cell}\)를 씁니다. `GFS_TCF_FACE_Z` 또는 `GFS_TCF_ZH`가 켜졌다면 같은 \(m\)을 임피던스에도 사용합니다.

## 코드 감사

- 2차 적합은 현재 \(r_s=\min(2.5\sqrt{V_i},R_{\rm searched}\,\text{cellsize})\) 안의 후보를 모읍니다. `npt >= 8`과 normal matrix Cholesky의 `1e-8*tr` pivot 조건으로 성공 여부를 정합니다. [적합 구현](</tmp/lageunha-push/Exam/exam.c:4757>)
- QBP 캐시는 2차 적합 성공 셀의 기울기만 저장합니다. 실패하면 `qin = ducell`을 넣고, 이어지는 minmod 결과는 \(\Delta u_{\rm cell}\) 전체가 됩니다. 따라서 Q가 최대 셀 점프를 받아, 요청 문서의 실패 시 과도한 Q 진단과 일치합니다. QBP 실패 카운터도 고유 실패 셀이 아니라 해당 fallback을 만난 면 수입니다. [QBP 캐시](</tmp/lageunha-push/Exam/exam.c:4826>) [QBP 면 처리](</tmp/lageunha-push/Exam/exam.c:8114>)
- 현재 이웃 탐색 확장은 미배정 모서리나 보안 반경 부족으로 Voronoi 셀이 불완전할 때 동작합니다. **LS 적합의 rank 실패만으로는 탐색이 확장되지 않습니다.** 보안 반경은 \(2R_v\le\text{reach}\)일 때의 기하 완전성 검사이며, `GFS_PERIODIC_COMPLETE`가 켜진 경로에서만 작동합니다. [보안 반경](</tmp/lageunha-push/Exam/exam.c:4022>) [탐색 확장](</tmp/lageunha-push/Exam/exam.c:4390>)
- `shf_*`는 각 셀의 1-hop 면 이웃과 면 끝점을 기록합니다. HG 조립은 이 기록으로 이웃의 이웃까지 2-hop 표본을 만들지만, 현재 2차 적합이나 QBP에는 이 표본을 쓰지 않습니다. `SHF_K=64`를 넘는 이웃은 기록에서 잘리므로, 2-hop 적합에 재사용하려면 overflow 계측을 추가해야 합니다. 공유 기하 모드는 다중 MPI rank에서 중단됩니다. [공유 기하 기록·면 교집합](</tmp/lageunha-push/Exam/exam.c:567>) [HG 2-hop 구성](</tmp/lageunha-push/Exam/exam.c:847>)

## 순위 판정과 stencil 확장

현재 세 조건의 성격은 다릅니다.

| 코드 조건 | 판정 |
|---|---|
| \(2.5\sqrt{V_i}\) | 회전 불변인 반경 cutoff지만 경험적 국소 범위입니다. 확장의 시작 반경으로 쓰고, 실패 시 하드 정지점으로 쓰지 않습니다. |
| `npt >= 8` | 충분한 표본 수라는 휴리스틱입니다. 5계수 문제에서 필요한 최소 수는 5지만, 개수만으로 full rank를 보장하지는 않습니다. |
| `1e-8*tr` | normal matrix에 적용한 수치 pivot cutoff입니다. 기저·pivot 순서에 의존하고, normal equation은 설계행렬의 조건수를 제곱합니다. |

후보점 \(\mathbf r_k=(x_k-x_i,y_k-y_i)\)에 대해 \(h=\sqrt{V_i}\), \(X_k=r_{kx}/h\), \(Y_k=r_{ky}/h\)로 두고 다음 가중 설계행렬을 만드십시오.

\[
A_{k,:}=\sqrt{w_k}
\left(X_k,\;Y_k,\;\frac{X_k^2}{2},\;\frac{X_kY_k}{\sqrt2},\;\frac{Y_k^2}{2}\right),
\qquad
w_k=\frac{h^2}{|\mathbf r_k|^2}.
\]

대칭 Hessian에 Frobenius 정규화한 기저라 회전 시 계수 변환이 직교 변환이 됩니다. \(A=U\Sigma V^\mathsf T\)의 특이값으로 rank를 판정하고, rank 5일 때만 2차 적합을 승인하십시오. 구현 허용오차는 물리 cutoff 대신 입력 자료형과 기하 좌표의 machine epsilon으로 정합니다.

\[
\tau=\max(n,5)\,\epsilon_{\rm eff},\qquad
\operatorname{rank}(A)=\#\{\sigma_\ell>\tau\sigma_1\}.
\]

1) \(r_0=\min(2.5h,r_{\rm geom})\)에서 시작해, rank가 부족하면 2-hop 후보를 포함하고 **같은 거리의 후보를 한꺼번에** 추가하면서 반경을 키웁니다. 후보 위치는 \((i,j,k)\)의 2-hop 변위 \(\mathbf r_{ij}+\mathbf r_{jk}\)로 구하고, 입자 ID와 주기 이미지가 같은 중복만 제거합니다.  
2) 이웃 셀 검색 반경도 늘려 보안 반경 검사를 통과한 범위까지 반복합니다. 현재 코드는 이 rank 실패에 따른 확장을 하지 않으므로 별도 연결이 필요합니다.  
3) 최대 완전 탐색 범위에서 2차 rank가 여전히 부족하면 선형 기울기를 씁니다. 선형 fallback은 \((X,Y)\) 두 열의 rank 2 여부로 확인합니다. rank 2도 성립하지 않으면 관문 실패로 처리합니다.

여기서 조건수는 \(\kappa_2(A)=\sigma_1/\sigma_5\)로 기록하십시오. 2차 적합에는 rank 5, 선형 fallback에는 rank 2가 필요합니다. rank 판정은 위 machine-precision 기준을 사용하고, 별도의 `1e-8` 물리 cutoff는 두지 않습니다.

## 실패 원인 계측

각 셀의 적합 시도에 다음을 기록하십시오.

- **표본·반경:** 고유 표본 수 `npt`(1-hop/2-hop 분리), \(r_s/h\), 사용한 최대 반경, `R_searched`.
- **조건:** SVD rank, \(\sigma_1,\sigma_5,\kappa_2\), 실패 원인(`npt<5`, rank 부족, 비유한 값 등).
- **기하:** 셀 종횡비 \(AR_i=\sqrt{\lambda_{\max}(C_i)/\lambda_{\min}(C_i)}\), 여기서 \(C_i\)는 셀 다각형의 면적 중심 기준 2차 모멘트 행렬입니다. stencil 이방성도 별도로 보려면 가중 표본 공분산의 고유값 비를 함께 기록합니다.
- **완전성:** `reach`, \(2R_v\), 여유 \(reach-2R_v\), 종료한 탐색 반경. 2-hop 후보 수집에서 `SHF_K` 초과 수와 역방향 면 미발견 수도 추가합니다.
- **QBP 사용:** 셀별 기울기 출처(`Q2`, `linear`, `invalid`)와 면별 \(\Delta u_{\rm cell},\delta_{b'},m,Z,Q,\frac12 A Q\Delta u_{\rm cell}\). 기존 `qbp_nfail`만으로는 어느 셀이 어떤 이유로 실패했는지 알 수 없습니다.

실패 셀의 ID·좌표를 전부 장기 출력하기보다는 원인별 집계와 최초 실패 단계의 제한된 표본 로그를 두면, 적합 실패와 기하 탐색 실패를 구분할 수 있습니다.

읽기 전용으로 검토했으며 파일은 수정하지 않았습니다.