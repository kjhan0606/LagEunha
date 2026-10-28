# Kepler 모멘트 모형 계획 및 독립 평가 — 2026-10-18

작성: 사용자 지정 계획·평가 담당 Astra agent. 구현 및 계산 실행 담당은 root agent다. 이 문서는 승인 완료 보고가 아니라 현재 증거에 따른 작업 우선순위와 판정 기준이다.

## 범위와 현재 판정

공유 컨텍스트 전체를 읽고 저장소 `/tmp/lageunha-push`, origin `git@github.com:kjhan0606/LagEunha.git`, branch `codex/physics-validation-20261017`을 확인했다. 최초 작업 디렉토리 `/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce`와 같은 프로젝트 작업트리라고 추정하지 않았다. 연구 RHS, 적분 runner, Reynolds 설계, 누적 연구 기록과 converged angular 비교를 읽었다. 이 검토에서는 시뮬레이션을 실행하거나 기존 실행을 변경하지 않았다.

모멘트 기반 내부 회전 방향은 유지한다. 다음 작업의 중심은 전체 정상해 잔차, 퍼텐셜 폐쇄, 기존 생산 코드 통합이다. t4e-4의 coarse/fine hydro RHS에서 transport 차이가 우세하다는 결과는 항별 진단이다. 물리적 오류의 원인이나 해당 항을 억제해야 한다는 근거가 아니다. 종점 연결이 같다는 사실도 단계 중 면 교체를 배제하지 않는다.

기존 생산 A/B는 보조기억 모형이며 새 모멘트 모형과 다르다. A의 t8 이후 큰 압력·위상오차는 기존 모형의 정확도가 미확보라는 증거다. 시각이 다른 A/B 최신 snapshot을 직접 우열 비교하지 않는다. 기존 A/B 및 연구 연장은 유지한다. 짧은 연구 시간수렴은 t6–10 및 이후 생산 장시간 검증을 대체하지 않는다.

## 물리 해석의 필수 구분

고정 순질량 `dm/dt=0`은 면 유출입의 합이 0이라는 뜻이다. 면별 질량, 운동량, 내부 각운동량 수송은 남는다. 셀 질량중심의 속도도 일반적으로 `p/m`과 다르다. 생성점이 복원 유속을 따른다는 조건은 물질 경계를 보장하지 않는다. 따라서 내부 각운동량 수송 제거 또는 물질 affine 셀의 `Idot=2 tr(C)` 강제는 허용하지 않는다.

잔류 속도는 upwind 질량유속과 가중치율 풀이에도 영향을 준다. 명시적 residual 운동량 항이 작다는 결과는 residual closure의 간접 영향까지 작다는 뜻이 아니다. coarse/fine endpoint에 각기 RHS를 적용한 차이는 서로 다른 상태의 RHS 차이다. 이를 국소 truncation error 또는 불안정 원인으로 동일시하지 않는다.

중앙중력에서 원점 기준 토크는 0이지만 셀 내부 중력토크는 일반적으로 0이 아니다. `c=q/m`, `F_g=integral rho g dV`이면

`tau_internal,g = integral (x-c) cross rho g dV = -c cross F_g`.

정상 원운동의 수송, 중력, 중심 이동은 함께 실제 모멘트 변화율을 구성한다. 별도의 원심력을 추가하거나 hydro transport를 인위적으로 작게 만들지 않는다.

## 첫 구현 묶음: 전체 이동 셀 정상해 예산

정적 Eulerian Kepler 정상장의 시간에 독립인 적분밀도 `a(x)`에 대해 이동 셀의 정확한 기준은

`d/dt integral_cell a(x) dV = integral_boundary a(x) w_n dS`.

원시 RHS를 0과 비교하지 않는다. 정확한 rho, v, P와 potential에서 질량, q, Q, p, D, Pi, Ehyd, Epot를 적분하고 위 경계 이동식을 적용한다. 내부 ell 기준은 raw derivative에서 `q cross p/m`의 미분을 뺀다. dm, dq, dp를 모두 사용한다. 초기 정확 projection의 q,p,m으로 중심을 정의하고 generator를 중심으로 대신 사용하지 않는다.

하나의 full balance report에 다음을 함께 저장한다.

1. 입력 상태·기하·설정·source SHA, 셀 ID와 초기 영역, 가중치율 반복 수렴값, 실제 면속도 convention.
2. 원시 모멘트별 정확한 경계 이동 derivative, hydro 수송·압력·체적항, 중력 source, 전체 RHS 및 양자의 차이.
3. 내부 ell의 수송·압력·중력·중심 이동·전체 변화율과 정확 기준 차이. 각운동량의 원점 기준 합과 내부/궤도 분해를 함께 보관한다.
4. Ehyd와 Epot의 물리적 수송·일 및 합산 결함. runner의 유체 ledger와 별도 열로 기록한다.
5. 전체/접촉면/접촉면 밖/기존 관심 고리별 절대 및 질량가중 오차. 거의 0인 참조값으로 나눈 상대오차만 보고하지 않는다.

첫 비교에서는 수렴된 계산 면속도 w를 고정하고 정확한 정상장을 대입한다. 그 면속도에서 정확한 질량 변화율이 0이 아닐 수 있으므로 dm 항을 유지한다. 두 번째 비교에서는 정확한 정상장 질량유량으로 가중치율을 재해결하여, 고정질량 조건 자체의 오차와 모멘트 복원·수송 오차를 구별한다. exact weight solve가 불가능하거나 퇴화하면 이를 실패 증거로 남기고 임의 보정하지 않는다.

접촉면 r=.5와 기타 실제 초기장 분기점을 면 적분에서 분할한다. 체적 적분도 불연속을 반영해야 한다. 매끈한 quadrature 차수만 높인 결과를 참조로 삼지 않는다. 접촉면 분할·기하·가중치율의 오차가 관측 residual보다 충분히 작은지 확인하되, 작은 검사들을 독립 작업으로 끝없이 늘리지 않고 이 보고서의 참조 오차 예산으로 묶는다.

참해 잔차 평가는 정확한 정상장을 적분한 초기 상태에서 한다. 진화 checkpoint의 hydro+gravity 예산은 별도 진단이며 그 상태의 정확해를 알고 있다는 의미가 아니다. frozen 셀의 강제 0 RHS는 정확 이동 경계 derivative와 직접 비교하지 않는다. 활성 셀 비교와 reservoir 교환 장부를 명시한다.

root가 구현하는 첫 파일은 기존 연속 IC 적분과 수렴 initial RHS를 재사용하는 full balance evaluator다. 실제 면속도와 split convention은 `Docs/polybnd/kinetic_factor_rhs.py` 및 `QuadraticUpwindNodes`와 맞춘다. 단일 보고서를 검토하여 지배적인 정상해 결함을 결정한 뒤 후속 solver 수정을 선택한다.

## 퍼텐셜 폐쇄: 생산 승인 전 필수 관문

`Docs/run_continuous_quadratic_evolution.py`의 작은 ledger_error는 Ehyd에 중력 일을 장부로 추가한 항등식을 확인한다. 복원 밀도에서 계산한 Epot와 합친 물리적 총에너지 보존을 확인하지 않는다. 기록된 비이차 spline potential의 방향미분 결함도 미해결이다.

우선 검토할 설계는 독립 `Psi=integral rho Phi dV` 모멘트를 동일한 질량 면유량으로 수송하고 중력 일을 반대 부호로 넣으며, 밀도 복원도 Psi를 만족시키는 방법이다. 그러나 기존 추가 phi 기저 적합에서 접촉면 밀도 형상이 악화되었다. 보존량 적합 성공만으로 승인하지 않는다. 접촉면 표현, 실현가능성, 면유량 정확도까지 함께 판단한다.

대안은 동일한 이산 퍼텐셜 함수의 미분으로 힘과 일을 구성하는 방법이다. 이때 기하와 가중치 의존성까지 포함해야 한다. 셀 중심 힘만 바꾸어 에너지 수치를 맞추는 것으로 충분하지 않다. 정상해 잔차 결과에 따라 한 경로를 선택하고, 두 대안을 동시에 구현하는 확장을 피한다. 열 보정·원운동 조향·임의 clipping으로 결함을 숨기지 않는다.

## 생산 CUDA/OpenMP/MPI 연결

기존 생산 경로를 유지하고 새 모멘트 상태를 명시적인 선택 경로로 추가한다. 기존 point-state restart를 같은 물리 상태의 moment restart라고 승격하지 않는다. 비교에는 동일한 연속 초기장을 각 표현에 일관되게 투영한다.

현재 C/CUDA의 17성분 표현과 Python raw18 및 generator 포함 20성분 표현이 다르다. `Voro/Laguerre/affine_moment_flux.h`의 potential flux는 0이며 실험용 수송이다. `moment_pressure_rate.h`와 `moment_gravity_step.h`도 제한된 모형이다. 이 헤더들의 존재나 점 함수 검증을 전체 생산 solver 완료로 보고하지 않는다.

통합 순서:

1. 상태 layout, 단위, D의 인덱스 방향, Pi 대칭 표현, E/Psi 정의, 공통 면 소유권을 고정한다.
2. 승인된 Python 참조와 같은 입력을 받는 C 연산자로 면 수송, 체적 응력, 중력, factor 갱신을 구성한다. Python의 느린 조립 방식을 그대로 이식할 필요는 없다.
3. 실제 생산 geometry 및 가중치율 풀이를 포함한 전체 RHS와 수락 단계에서 CPU, OpenMP, CUDA를 대조한다.
4. MPI의 공유 면 단일 유량/작용점, ghost moment/factor/weight 상태, 전역 제약 수렴, 영역분할 수 의존성을 검증한다. 병렬 합산 순서 차이에 무조건 bitwise 일치를 요구하지 않고 보존량과 오차 예산에 근거한 허용치를 사용한다.
5. restart에 moment, factor, weight, 설정과 버전을 저장한다. nonlinear seed는 재수렴할 초기 추정으로 취급한다. 단계 거절·재시도 시 모든 확장 상태와 장부가 일관되게 복구되어야 한다.

통합 설계와 상태 계약 정리는 물리 폐쇄 평가와 함께 시작할 수 있다. 미승인 폐쇄를 장시간 생산 경로에 먼저 넣고 나중에 의미를 정하는 순서는 피한다.

## Reynolds 및 최종 장시간 승인

사용자 지정 Kepler는 Re=infinity이며 물리 전단점성·bulk viscosity는 0이다. 내부 회전, 압력 토크, 보존적 수송은 유지한다. 유한 Re 구성식 연결은 현재 Kepler 안정성 확보의 선행 작업으로 두지 않는다. 향후 연결 시 대칭 응력의 면/체적 모멘트와 에너지 일을 함께 구현하여 중복 소산을 방지한다.

기존 A/B와 두 연구 연장은 중단·재시작하지 않는다. 연구 연장의 공통 종료 시각 결과를 종합하고, 구체적인 실패가 없는 한 각 checkpoint마다 새 미세 실험을 만들지 않는다. 기존 장시간 결과는 같은 실제 시각과 같은 초기 영역으로 비교한다.

새 생산 경로는 frozen patch를 넘어 원래 전체 Kepler 영역에서 평가한다. 같은 초기 영역별 압력, 위상, 속도, 위치, 각운동량, 물리적 총에너지와 경계 수지를 보고한다. 역사적 t6–10을 통과한 뒤 더 긴 목표 시각까지 확인하며 시간·공간 해상도 의존성도 검사한다. 짧은 약 2차 시간수렴 또는 ledger 양성만으로 이 승인을 대체하지 않는다. 장시간 허용오차는 실행 전에 물리적 목적과 기존 실패 규모에 비추어 명시한다.

## 실행 우선순위

전체 이동 셀 정상해 잔차 평가 → 지배 결함에 따른 접촉면/퍼텐셜 폐쇄 결정 → 기존 생산 경로의 상태 및 승인 연산자 통합 → CUDA/OpenMP/MPI 전체 단계 정합성 → 전체 영역 장시간 Kepler 검증.

다음 평가에서는 root가 작성한 full balance evaluator의 수식, 면속도 부호, 불연속 분할, raw 인덱스, 중심 이동 및 frozen 처리부터 검토한다.

## 첫 참조 모듈 정적 검토

`Docs/polybnd/steady_kepler_ale_reference.py` 최초 구현을 읽었다. raw18+Psi 순서, `v=sqrt(gravity_coefficient) Jx`, 유체 에너지밀도 `P/(gamma-1)+rho |v|²/2`, power-bisector 면속도 및 양쪽 셀의 공유 부호는 수식상 맞는다. 실행 검증을 의미하지 않는다.

driver에 필요한 보완은 다음과 같다.

- rho의 rin/rout뿐 아니라 중력·potential 함수의 실제 분기 반경을 적분 분할에 포함한다. rc와 rout이 같은 경우 중복은 제거한다. 참조 차수 민감도, edge normal/length/bisector 잔차를 같은 보고서에 넣는다.
- 정확 initial moment와 모델 initial moment의 차이도 기록한다. 정확 ell derivative는 정확 q,p,m으로 계산하고 모델 ell derivative는 모델 상태를 사용한다.
- 모델 중력 source와 정확 rho/v의 중력 source를 분리한다. 원점 기준 중력토크0 및 내부토크 `-c cross F_g`를 확인한다.
- 모델에 Psi 진화가 없으면 해당 값을 0으로 채우지 않는다. 미구현이라고 표기하고 fitted-density potential 방향미분 결함은 별도로 다룬다.
- 독립적인 exact advective flux, pressure/volume, gravity 조립과 exact boundary derivative 사이의 항등식 검사로 참조의 부호·체적항 누락을 확인한다. 이를 모델 잔차와 구별한다.
- 외부면 wn=0은 현재 frozen-exterior 계약에 한정한다. 활성 셀 및 reservoir 장부를 구분하고 frozen 강제0를 정확 derivative0로 취급하지 않는다.

## 첫 실제 initial balance 결과 평가

`Docs/check_steady_kepler_ale_balance.py`와 `canonical_steady_ale_initial_balance_2026-10-18.json`을 검토했다. 초기 raw와 기존 split q16 projection의 bitwise 일치, 실제 gravity 입력 SHA 확인, 원시 인덱스 및 ell chain rule은 타당하다. 참고로 bitwise projection 일치는 이전 수치적 초기 투영물과의 일치이며 연속장 적분이 무오차라는 뜻은 아니다. 실제 circular split radii는 .25,.5,2.0이다.

active962에서 raw pdot 차이/m의 RMS .27751556, inner490 .00072534, .48–.55의208셀 .24304097이다. 이는 bulk 속도 가속도 오차가 아니다. 실제 w에서 exact dm가 최대2.0536e-4이며 모델의 dm=0과 다르므로, 큰 bulk momentum 수송 차이가 포함된다. 현재 geometry 잔차2.229e-11, 8→16→32 boundary 적분 차이 약5.7e-17은 이 참조 적분의 안정성을 지지하지만 독립 volume closure 검증을 대체하지 않는다.

동일 초기 raw에서 필요한 파생 변화율 차이는 다음과 같다. Delta는 모델−정확 참조다.

`Delta vdot = (Delta pdot - v Delta mdot)/m`

`Delta cdot = (Delta qdot - c Delta mdot)/m`

`Delta udot = (Delta Edot - 0.5 tr(Delta Pidot) - u Delta mdot)/m`, 여기서 `u=(E-0.5 tr Pi)/m`.

root가 저장 배열에서 계산한 `canonical_steady_ale_material_rates_2026-10-18.json`의 active bulk velocity rate 차이 RMS는 .00462720026, centroid .00257054283, specific thermal .000416386658이다. inner bulk는 .00015078978, contact bulk는 .00388146323이다. raw pdot/m .2775를 가속도로 읽으면 결함 크기를 크게 왜곡한다는 점이 실제 확인되었다. 이 파생 지표의 명칭은 moving-cell mean rates이며 물질 입자 가속도로 다시 동일시하지 않는다.

같은 실제 w에서 `Vdot=integral_boundary w_n dS`를 계산하면 정확 정상 압력의 Pdot=0과 직접 비교할 수 있다. 모델의 `Pdot=(gamma-1)(Edot-0.5 tr Pidot)/V-P Vdot/V`를 추가한다. active 최대 raw momentum 결함2.888은 .48–.55 최대.565보다 크므로 .4–.48 및 최대셀 ID/반경도 조사해야 한다. 결과를 접촉 고리 하나의 문제로 제한하지 않는다.

내부 ell 변화율 결함 RMS1.497e-8과 정확 ell 변화율 RMS1.416e-6은 비영 변화율을 참조로 삼아야 함을 지지한다. 이 비율만으로 물리 정확도 또는 장시간 안정성을 승인하지 않는다. 독립 exact volumegravity 및 exact advective/pressure/volume closure를 완료한 뒤 모델 항별 잔차를 판정한다.

다음 exact fixedmass 기준은 선형 문제다. 단일 정확 정상장의 면 양쪽 rho,v는 같으므로

`b_i = integral_face rho (v-w0) dot n dS`, `C_ij = integral_face rho ds/(2 d_ij)`

로 조립한 Laplacian `L z=b`를 active/frozen 조건으로 푼다. 이때 generator velocity는 실제 값을 고정하고 frozen z=0의 현재 계약을 유지한다. 성분별 gauge/고립 조건을 확인하고 재조립 residual과 exact boundary dm를 검사한다. 임의 residual velocity ensemble이나 nonlinear upwind solve가 필요하지 않다. 실제 w 기준과 exact-z 기준을 함께 보관한다. 모델 RHS를 exact z로 audit하면 모델의 netmass 제약이 깨질 수 있으므로 이를 실제 fixedmass 진화 RHS와 혼동하지 않는다.

## Exact-z 초기 참조 풀이 평가 및 다음 폐쇄 구현

`exact_fixed_mass_weights`, `Docs/check_steady_kepler_exact_mass.py`, `canonical_steady_ale_exact_mass_2026-10-18.json`을 검토했다. 계수와 부호, frozen Dirichlet 처리, active component anchor 검사는 이 참조 문제에 맞다. active 연결성분1, sparse 잔차7.59e-19, 독립차수48 boundary dm 최대2.34e-16, 차수16→32 z차이4.59e-17이다. 이 초기 고정질량 참조 풀이를 다음 전체 예산에 사용할 수 있다.

actual−exact z의 질량가중 RMS4.8013563e-4, 최대.00314252는 수치 closure가 선택한 면운동과 정확장의 고정질량 면운동이 다르다는 증거다. solver 개선 효과나 원인 확정을 뜻하지 않는다. z는 차원을 가진 양이므로 필요시 `Delta wn=(Delta zi-Delta zj)/(2d)`를 면 적분 가중치와 함께 보고한다. 기존 numerical RHS나 진화에 exact z를 적용하지 않은 범위를 유지한다.

다음 우선순위는 독립 exact volume+flux closure 전체 구현이다. 같은 실제 w에서 모든 원시 적분밀도 a에 대해 기본 advective rate는 `-integral_boundary a (v-w) dot n`이다. 체적 기하항은 q에 p, Q에 D+D^T, D에 Pi+P V I이다. 압력 면항은 p에 `-integral P n`, D에 `-integral x outer P n`, E에 `-integral P v dot n`이다.

중력 source는 p에 `integral rho g`, D에 `integral rho x outer g`, Pi에 `integral rho (v outer g+g outer v)`, E에 `integral rho v dot g`이다. Psi에는 advective rate와 `-integral rho v dot g`를 더한다.

현재 정확장 P는 일정하므로 Pi의 압력 gradient source는0이다. 이를 명시하여 구현할 수 있다. 압력의 보존형 면 표현을 사용할 때는 `-integral_boundary P(v outer n+n outer v)`와 `+integral_volume P(grad v+grad v^T)`를 모두 넣어야 하며 한쪽만 넣으면 틀린다.

접촉면 분할 체적 적분 두 차수로 위 source를 독립 평가하고 전체 합을 boundary derivative19와 대조한다. circular field의 E/Psi 중력 일이0인 점도 수치적으로 확인한다. 같은 identity를 exact-z 면속도에도 적용하되 동일 volume/gravity 적분은 재사용한다. 참조 identity가 닫힌 뒤 모델 hydro−exact hydro, 모델 gravity−exact gravity를 합성하여 어떤 항이 실제 잔차를 만드는지 판정한다.

## Full closure 중간 결과 및 검토 정정

root 보고에 따르면 q8 full closure의 active 최대결함은 actual/exact-mass 두 면운동 모두1.33814e-15, frozen은 약3.6e-16이다. q16은 진행 중이다. 1590셀 geometry audit에서 둘째 endpoint/edge match는0, 닫힌 법선합 최대2.3158e-11, 면적 상대차4.44e-16으로 첫 endpoint 기반 polygon 구성의 자료 계약을 확인했다. 근거는 `canonical_exact_closure_regions_2026-10-18.json`, `canonical_steady_reference_geometry_2026-10-18.json`이다.

검토 중 `g['source']`가 frozen에서0이라고 추정한 것은 잘못이었다. root가 `check_continuous_initial_gravity.py`를 확인한 결과 source 배열은 frozen에서도 실제 평가값을 저장하고, combined 진화에 더하는 것만 active에 한정한다. 따라서 이 배열은 evaluated model gravity로 기술해야 한다. active/frozen 통계를 나누는 것은 여전히 유용하지만 frozen source0와 exact nonzero를 뺀 통계라는 설명은 철회한다. q8 active evaluated model−exact gravity 최대차는8.80798e-9이다.

q16 완료 후 참조 identity와 적분 오차가 확보되면 추가 미세 참조 분기를 만들기보다 model/exact hydro 및 gravity의 bulk·thermal·ell 잔차를 분해하여 실제 수정 대상 하나를 결정한다.

## 완결 참조 이후 첫 수정 후보와 단일 교체 묶음

q16 full closure 완료: q8→16 volume 최대차1.7347e-18, active closure1.33816e-15, 초기 projection과 독립 volume moment 차6.288e-18이다. `canonical_steady_ale_model_residuals_2026-10-18.json`의 active bulkmean 잔차 RMS는 hydro .00462738, gravity3.520e-6, total .00462720이다. thermal은 hydro .00041638666, gravity5.45e-16이며, ell은 hydro1.5036e-8, gravity2.184e-10이다. r=.4–.48 bulkmean .0210104와 thermal .00329799도 거의 전부 hydro다. 현 초기 잔차의 첫 수정 후보는 중력·토크·점성이 아니라 접촉면 근처 밀도 복원과 체적/면 밀도 일관성이다. 다만 hydro 합산만으로 밀도 복원이 원인이라고 입증하지는 못했다.

최대 bulk/thermal 셀은 r=.46875로 매끄러운 저밀도 쪽이며 그곳의 체적 중력은 거의 정확하다. 기존 기록에는 인접한 접촉면 교차셀에서 지수형 밀도 복원이 크게 왜곡된다는 증거가 있다. 이로부터 밀도 trace를 첫 후보로 선택하되, 일반적인 실제 복원기를 구현하기 전에 단 한 번의 원인분리 교체 묶음으로 후보 타당성을 확인한다.

`check_density_face_pressure_effect.py`는 affine 속도와 면압력만의 과거 비교이며 massroot나 전체 RHS는 없다. 아래 quadratic 전체 RHS/수렴 massroot 비교와 중복되지 않는다. 과거에는 정확 밀도 trace만으로 면압력 RMS가 개선되지 않았으므로 이번에도 개선을 미리 가정하지 않는다.

새 audit는 baseline의 states, quadratic coefficients, factors, 압력복원, generator를 한 번 준비해 고정한다. `density_provider`로 decode/prepare를 다시 돌리면 basis/budget/velocity까지 바뀌므로 density-only가 아니다. face rhoL/rhoR만 정확 rho(point)로 교체하고 그 rho로 kinetic wl/wr, thermal specific P/[(gamma-1)rho], star pressure/work를 일관되게 재계산한다. raw moments, volume density quadrature와 velocity basis, gravity는 유지한다. 변경된 압력에 따른 F,H 및 A,a0,Pi 압력 source는 기존 RHS 식대로 반드시 재계산한다. +PV I의 P,V는 기존값이다. factor 시간진화는 하지 않는다.

한 보고서에 세 경우를 비교한다.

1. A: 원래 density와 새 공통 quadrature의 baseline. 저장 baseline rates 재현 및 quadrature 변경 효과를 기록한다.
2. B: 정확 density trace, 실제 z 고정. 비영 mass rate를 기록하고 chain rule에서 유지한다.
3. C: 정확 density trace, 기존 kinetic sample massroot 재해결. generator 고정, frozen z=0, quadratic upwind roots와 circular cuts의 합집합을 사용하며 매 outer iteration 재조립해 기존 수렴기준을 만족한다. 이는 monokinetic exact-z와 다른 문제다.

각 경우 자신의 w로 정확 boundary derivative를 계산하고 bulkmean, centroid, thermal, Pdot, ell의 signed 잔차를 비교한다. active/inner/.4–.48/contact 및 최대셀을 보고하며 압력과 advective 항도 분리한다. A의 면속도 reference를 C에 재사용하지 않는다.

이 교체는 체적/면/velocity moment가 완전히 일관된 새 물리 closure가 아니라 진단이다. 정확 원형 rho를 진화에 주입하지 않는다. B에서 개선되고 C에서 되돌아가면 가중치 제약과 결합을, 두 경우 모두 실패하면 밀도 단독 후보 기각을 보고한다. bulk 잔차만 줄고 thermal/Pdot/ell이 악화하면 생산 변경으로 승인하지 않는다. 개선되는 경우에만 원형 경계 정보를 사용하지 않는 양성·모멘트 제약 접촉면 복원 설계로 진행한다. 이 원인분리 결과 전에는 임의적인 새 밀도 기저/제한자/복원 계수를 선택하지 않는다.

## 기존 보조기억 장시간 실행의 추가 평가

root가 회수한 최신 물리 분석에서 A는 t8.8008489041, 초기 내부배경804셀의 압력 상대RMS4.21707983와 위상RMS47.1845374도다. B는 t7.6011825333, 같은 정의 영역의 압력RMS2.58147638와 위상RMS31.2080061도다. 근거는 `canonical_memory_physical_live_2026-10-18.json` 및 `canonical_corotating_memory_physical_live_2026-10-18.json`의 해당 보고 시점이다. live 파일이므로 후속 갱신 시 실제 저장시각을 함께 확인한다.

두 보조기억 후보 모두 각자의 역사적 t6–10 구간에서 큰 압력·위상오차가 지속되어 장시간 물리 정확도를 확보하지 못했다. B의 단순 회전운반 변경만으로 안정화에 성공했다는 근거도 없다. 두 최신 시각이 다르므로 이 숫자로 개선율을 계산하지 않는다. 이는 새 독립 모멘트 모형의 실패 판정이 아니다. 실행은 유지하며 이미 시작한 density A/B/C audit와 연구 fine 연장 결과를 기다린다. 통계 추가만을 위해 실행을 중단하거나 재시작하지 않는다.

## ABC 완료 판정과 일반 밀도 후보의 구현 계약

`canonical_density_trace_counterfactual_audit_2026-10-18.json`의 세 경우를 확인했다. A는 새 분할에서도 기존 rates를 최대6.919e-13 차이로 재현했다. C의 outer 변화는 .00313835,7.0192e-7,4.1554e-12이며 순질량 잔차1.0083e-17이다. active RMS A→C는 bulk .004627200259→.000125839159, centroid .002570542825→1.7498483e-5, thermal .000416386658→3.0439277e-7, ell1.4971709e-8→9.6598398e-10이다. hydro Pdot도4.1564119e-8→3.9738117e-8로 감소했다. B의 압력률 악화는 C의 질량 재제약 이후 유지되지 않았다.

따라서 밀도 trace 복원을 첫 실제 수정 대상으로 채택한다. 이것은 현재 초기 hydro 결함의 주요 원인에 대한 통제된 증거이며 새 일반 closure의 승인 또는 장시간 원인 확정은 아니다. 다음은 정확 원형정보를 사용하지 않는 한 개 후보를 구현한다. 기존 `neighbor_density_audit.py`, `moment_aligned_density_audit.py`, `discontinuous_density_audit.py`의 기하 도구는 참고하되 global phase extrema, knownphase 값, shear전용입력, angle최적화나 blend 계수 탐색은 가져오지 않는다.

### 후보: 불연속 prior에 대한 양성 exponential moment 적합

`ContactPriorDensity(poly,m,q,Q,prior_metadata,order)`가

`rho(x)=m rho_prior(x) exp(theta dot F(x))/Z`

를 제공한다. F는 기존 centered/scaled 5특징 x,y,x²,xy,y²다. 평면 cut 양쪽에서 rho_prior는 양의 상수 low/high이며 theta Newton 풀이로 현재 m,q,Q를 그대로 만족시킨다. low/high는 물리적 hard bound가 아니며 적합 후 min/max와 overshoot를 보고한다. positivity만으로 접촉면 정확도를 주장하지 않는다. clipping, moment 수정, 열 보정은 금지한다.

공통 API는 mass,area,origin,scale,points,weights,weighted_mass,value,theta,logZ,residual,iterations,cut_metadata다. polygon을 실제 cut에서 분할하고 양쪽 subpolygon에 양성 triangle quadrature를 적용한다. weighted_mass와 value는 같은 rho를 나타내야 한다. 면에서도 같은 함수와 한쪽 trace를 사용한다. 선이 면과 일치하는 퇴화에는 셀 내부에서 접근한 trace를 정의하고 임의 epsilon 이동으로 값을 선택하지 않는다.

### 지역 prior의 결정 규칙

1. 현재 mesh edge adjacency의 정확히2-hop+self를 stencil로 고정한다. frozen도 포함한 현재 mean=m/V만 사용한다. global extrema, 반경, 정확 phase값은 사용하지 않는다.
2. low/high는 이 stencil mean의 최소/최대다. rho_scale=max(localmean)>0. high−low<=64eps*rho_scale이면 constant prior(mean_i)로 둔다.
3. h_i=sqrt(V_i), geo centroid를 사용하고 d_j=(geo_cj−geo_ci)/h_i에 대해 `(d_j dot g−(mean_j−mean_i)/rho_scale)/|d_j|`의 제곱합을 최소화한다. rank2 및 |g|>64eps이면 normal=g/|g|다.
4. 위 방향이 식별되지 않으면 |mass_centroid−geo_centroid|/h_i>64eps일 때 그 방향을 normal로 쓴다. 둘 다 식별불가하면 constant prior(mean_i), reason=unidentified_direction을 기록한다.
5. fraction=(mean_i−low)/(high−low)가64eps 이내0/1이면 constant prior다. 그 외 normal 양쪽 중 높은 쪽 면적이 fraction*V가 되도록 polygon clip+bisection으로 offset을 결정한다. 각도 최적화나 phasefloor blend는 첫 후보에 없다.

constant prior는 기존 exponential family와 수학적으로 같은 명시적 퇴화 경로다. moment fit Newton 실패, feature/Hessian rank 결함, 비유한값을 기존 smooth density로 자동 우회하지 않는다. stencil LS rank부족의 위 명시적 방향 퇴화와 moment-fit 실패를 구별한다. 분할에서 임의 작은 면적 cutoff로 질량을 버리지 않는다.

### 전체 RHS 연결과 하나의 검증 묶음

모든 geometry/raw로 prior목록을 만든 뒤 density_provider로 전달한다. 기존 raw m,q,Q,p,D,Pi,E와 factors는 유지하되 새 density로 quadratic basis/예산/velocity를 다시 구성한다. generator도 새 복원 유속을 따르고, 중력은 새 밀도와 새 속도로 재적분한다. ABC처럼 velocity와 gravity를 고정하면 완결 후보 검증이 아니다. face node는 양쪽 셀 cut-line 교점과 quadratic normal-speed 영점의 합집합으로 분할하고 기존 외부 가중치율 반복을 수렴시킨다.

검증은 다음 묶음으로 제한한다.

1. uniform 및 이동/회전 planar-contact 자료에서 m,q,Q 보존, 양성, volume/face 함수 일치, 회전 공변성을 확인한다. 정확 phase/interface는 평가 자료에만 쓰고 prior inference 입력에는 넣지 않는다.
2. 실제1590 IC 전체의 inferred prior 적합과 분할 quadrature 두 차수를 비교한다. 원시 moment 불일치, min/max, cut/퇴화 빈도와 Newton 실패를 보고한다.
3. 새 density의 quadratic preparation, generator, 중력까지 일관된 초기 full RHS와 weight outer solve를 한 번 수행한다. 자신의 면속도에서 정확 정상장 reference를 구성하고 bulk,centroid,thermal,hydro/total Pdot,ell 및 기존 영역/최대셀을 baseline과 비교한다.

이 완결 후보 결과를 평가하기 전에 새 진화 실험이나 여러 prior 변형 탐색을 추가하지 않는다. 초기 잔차 개선만으로는 생산 승인하지 않으며 이후 이동기하/실현가능성, 퍼텐셜 폐쇄 및 생산 CUDA/OpenMP/MPI·장시간 관문은 유지한다.

### 실제1590셀 적합 평가

`canonical_contact_prior_actual_fit_2026-10-18.json` 확인 결과 inferred nonflat prior101셀, 차수12/24 모두1590셀 적합 성공 및 실패0이다. target normalized fit residual은 약2.997e-13, 독립 quadrature 간 normalized moment 차이는9.764e-16 이하, 독립 raw moment와 target 차이는4.202e-18 이하다. 여기서 `independent_normalized_moment_error`라는 기존 필드는 실제로 target 오차가 아니라 fit quadrature와 독립 quadrature의 차이이며, target residual과 구별하여 보고해야 한다.

표본 밀도범위 .00690626–약1.19535는 양성이지만 실제 phase 범위 밖이다. prior가 hard bound가 아니므로 clipping하지 않으며 이 범위를 full RHS 평가에 명시한다. 차수별 sampled max 차이는 표본점 변화도 포함하므로 이를 함수의 quadrature 불일치로 단정하지 않는다.

이 결과로 complete initial RHS 연결을 진행한다. 후보 최종채택은 아니다. 먼저 또는 같은 묶음에서 inferred prior의 회전/이동공변성을 지역 mesh patch의 geometry/raw/edges를 함께 변환해 확인한다. 면trace는 outward n에 대해 셀i inward=-n, 셀j inward=+n으로 평가한다. 동일분할 volume quadrature로 quadratic preparation을 다시 하고 generator와 중력까지 새밀도에 맞춘다. fullrates/own-w reference와 압력 원시 source 분해를 저장해 후속 불필요한 재실행을 줄인다.

### 공변성 실패에 따른 roundoff 판정 규칙 개정

inferred prior의 회전(.73)+이동(2.1,-1.7) 검사에서530셀이 oldflat→newnonflat으로 변했다.468셀은 거의 균일한 mean의 roundoff 대비,62셀은 phase 극값과의 아주 작은 차이 때문에64eps 경계를 넘었다. 관측 area 상대변화 최대1.4877e-14는 longdouble 누적으로도 같고 누적정밀도 차이는2.22e-16뿐이다. 주된 원인은 이미 float로 표현된 변환 좌표의 반올림이다. 진단 근거는 `canonical_contact_prior_covariance_diagnosis_2026-10-18.json`, `canonical_contact_geometry_roundoff_2026-10-18.json`이다.

고정64eps의 contrast/fraction 판정을 다음 uncertainty interval 규칙으로 대체한다. 위 초기 spec의 해당 두 판정보다 이 규칙이 우선한다. 물리 phasefloor나 threshold 튜닝이 아니며 분해할 수 없는 입력 대비를 prior에 만들지 않는 처리다.

1. optional vertex_error[nvertex,2]를 좌표성분 절대오차 상계로 받는다. 미제공 시 roundoff용 working envelope `e=gamma4*(abs(x)+h)`, `h=sqrt(A)`, `u=finfo.eps/2`, `gamma4=4u/(1-4u)`를 사용한다.4는2D affine 좌표의 두 곱·두 합 연산수다. 이 default가 전체 clipping 기하오차를 증명하지는 않는다. 알려진 geometry error가 있으면 별도 전달한다.
2. 변환 검사는 `e'=abs(R)@e+gamma4*(abs(R)@abs(x)+abs(shift))`로 입력·변환 오차를 전파한다. 관측 area 변화를 interval이 포괄하는지 보고한다.
3. 각 좌표구간을 nextafter(x−e,−inf),nextafter(x+e,+inf)로 구성한다. `r_k=I_k-I_0`, `A=0.5 sum cross(r_k,r_(k+1))`를 outward interval subtraction/product/sum으로 계산한다. product는 네 끝점 조합의 min/max이며 각 연산을 nextafter로 바깥 반올림한다. 명목A는 기존 계산을 유지한다. area 구간이0을 포함하면 명시적으로 실패한다.
4. mass 구간은 m의 양쪽 nextafter와 실제 알려진 projection 오차만 사용한다. rho 구간은 `[mlo/Ahi,mhi/Alo]`다. 임의의 projection 오차를 만들어 추가하지 않는다.
5. stencil rho_j=[lj,uj]가 공통교집합 `max lj<=min uj`를 가지면 constant prior, reason=unresolved_local_contrast다.
6. low 구간=[min lj,min uj], high 구간=[max lj,max uj]다. 공통교집합이 없으면 high_lo−low_hi>0이다. fraction 구간을 `(rho_i-low)/(high-low)`의 interval 연산으로 계산한다. 이 구간이0 또는1을 포함하면 constant prior, reason=unresolved_phase_fraction이다. 그 외 nominal low/high/fraction과 기존 normal/offset 방식을 사용한다.
7. metadata에 area/rho/fraction 구간과 오차모델을 기록한다. 기존 mixed101셀 유지 여부는 검증 결과이지 ID나 정확 원형마스크로 강제하는 조건이 아니다.

이번 실패는 contrast/fraction 판정에서 발생했다. true mixed 영역의 LS 방향 판정 문제는 아직 보고되지 않았으므로 rank/64eps 방향fallback까지 임의 확장하지 않는다. 개정 후 inferred 공변성과 실제 fit를 확인하고 complete RHS 연결을 재개한다. 종전 fit 결과와 실패 결과는 보존한다.

구간공변성 재검증 `canonical_contact_prior_interval_covariance_2026-10-18.json`은1590셀 전체 classification mismatch0, nonflat101 유지, 모든 관측 area 변화의 구간포함을 확인했다. 지정 회전·이동에서 최대 normal 오차1.901e-14, offset4.574e-14이며 선택3개 비영theta 적합의 density transform 상대차는2.411e-10 이하다. 정확히 cut선상인 양쪽 trace1/4도 확인했다. 이는 지정변환과 선택density점 검증 범위이며 모든 상태의 공변성 증명은 아니다. 이 결과로 complete RHS 연결을 재개한다. outward contrast 분모가0이하인 극한 roundoff에서는 unresolved_local_contrast로 처리하는 보호조건을 추가하며, 기존 실행을 재시작하지 않는다.

## Complete RHS 결과: 동적 연구 후보 채택

`canonical_contact_prior_complete_rhs_2026-10-18.json` 확인: outer 마지막 z 변화1.7347e-17, mass 잔차1.1933e-17이다. active baseline→candidate RMS는 bulk .00462720→.000685368, centroid .00257054→.000275380, thermal .000416387→1.05747e-5, ell1.49717e-8→2.38546e-9, Pdot4.15641e-8→1.73120e-8이다. r=.4–.48의 큰 결함도 줄고 inner r<.4는 대부분 유지되었다. 정확한 원형 정보를 prior에 넣지 않고 전체 density/velocity/generator/gravity를 일관되게 바꾼 결과다.

`canonical_contact_prior_rhs_budget_2026-10-18.json`의 active 정규화 최대결함은 p/(mv)9.49e-15, D/(mhv)1.347e-13, Pi/(mv²)1.774e-14, implied thermal fraction2.310e-10이다. pressure+remaining raw hydro 합의 잔차는2.118e-22다. 이 근거로 ContactPriorDensity를 다음 동적 검증에 사용할 연구 후보로 채택한다. 기존 production/default 교체 또는 장시간 성공 승인은 아니다. prior 변형과 계수탐색을 멈추고 다음 묶음으로 진행한다.

1. 같은 fitted density theta/cut 및 quadratic coefficients를 고정한 독립 중력 적분24/48을 수행한다. 차수를 바꾸면서 재fit하지 않는다. inferred cut subpolygon과 spline rc/2,rc 분기를 함께 분할하고 p,D,Pi,E source와 bulk/ell/thermal 기여를 비교한다. 이전12차 적분 오차가 현재 개선량에 비해 작은지, thermal source 상쇄와 원점 토크를 확인한다.
2. 기존 midpoint+Sylvester runner의 새 선택 경로에서 모든 accepted/half stage의 움직이는 geometry,2-hop prior,density,quadratic reconstruction,generator,face cut/speed split,z outer solve,중력 및 covariance source를 재구성한다. t0 geometry에 prior를 캐시하지 않는다. frozen mask는 원래IC ID 기준이며 frozen raw는 유지하되 density/prior는 실제 stage geometry에서 복원한다. 실패는 숨기지 않고 모든 상태·seed·장부를 일관되게 처리한다.
3. 동일한 연속 t0 raw/factors와 r<.60 활성영역에서 dt1e-4 및5e-5, 종료t1e-3의 두 분기만 실행한다. oldclosure의 중간 checkpoint에서 새 후보를 시작하지 않는다. 기존 완료 oldclosure 결과를 재사용하고 공통 종료시각 비교에 집중한다. 두dt로 차수를 주장하지 않는다. stage 양성·moment 적합·massroot·장부, prior 분류 변화/theta 성장, factor rebuild 열분배와 실제 오차를 함께 평가한다.

이 묶음은 동적 실현가능성 검증이며 t6–10을 대체하지 않는다. 계산 대기 동안 기존 생산 C/CUDA의 상태layout,공유면,ghost,restart 계약 구현을 진행할 수 있다. 정상해 초기 잔차 감소만으로 장시간 통과를 추정하지 않는다.

hydro ledger 통과를 Ehyd+복원Phi 보존으로 표현하지 않는다. 가능한 경우 같은 face massflux와 midpoint 가중치를 재사용해 진단용 Psi_transport 변화율 `-integral_face Phi*massflux-gravitywork`를 누적하고 복원 integral rhoPhi와 차이를 기록한다. 이 누적량은 열 보정이나 물리상태 수정을 하지 않는다. 다음 물리 설계 관문은 새 prior에서의 potential-moment closure이며 장시간 생산 승인 전에 해결해야 한다.

### 중력 정련 완료와 동적 두 분기 시작 판단

독립 중력 정련은 q12→24 source 최대차2.3126e-16, bulk RMS2.855e-13, ell RMS4.493e-19, thermal RMS7.652e-16이다. q24→48 source 최대차1.735e-18, bulk RMS5.143e-16으로 수렴했다. thermal source 최대9.04e-19, 원점 토크2.168e-19이고 기존 source bitwise 재현도 확인했다. 근거 `canonical_contact_prior_gravity_refinement_2026-10-18.json`. 이 동적 검증에서 q12 중력 사용에 동의한다.

`run_continuous_quadratic_evolution.py`의 선택적 contact 경로와 `contact_upwind_face_nodes.py`를 정적으로 검토했다. stage별 prior/density/quad 재작성, stage 내 outer 동안의 density 재사용, 새 density gravity/covariance source, 원래 ID 기반 frozen mask와 closure restart guard가 계약에 맞다. 지정한 두 분기를 새 root에서 시작할 수 있다. 초기 raw/factors 동일성, immutable module 경로/해시 및 적분 설정을 launch 기록에 남긴다.

Psi_transport 진단의 미구현은 이번 admissibility 시작의 차단사유로 삼지 않지만 물리적 총에너지 보존 통과도 주장하지 않는다. stage 분류 전환은 직전 evaluate stage 대비이며 accepted-step 간 변화와 구별한다. frozen geometry 변화 등으로 moment fit가 불가능해지면 실패를 보존하고 smooth fallback이나 열 보정으로 숨기지 않는다.

## 생산 통합의 첫 구현: 상태 layout adapter와 통신 계약

연산 canonical layout은 Python과 동일한 raw18로 정한다: m0,q1:3,Q3:7(row-major),p7:9,D9:13(row-major),Pi13:17(row-major),Ehyd17. Psi는 별도 scalar와 semantic mode(absent / diagnostic_reconstruction / transported_moment)다. generator x[2], additive power W[길이²], factor L[4], stable global ID,frozen,epoch도 raw와 분리한다. 기존 wm17 bytes와 생산 particle struct는 변경하지 않는다.

첫 구현은 `Voro/Laguerre/moment_state_layout.h`의 host/device index enum과 checked raw18↔wm17 adapter다. wm17의 Psi를 absent일 때0으로 채우지 않는다. 모든 입력의 finite/layout/필수 Psi mode를 검사하고 Qxy=Qyx,Pixy=Piyx의 정확 equality가 성립할 때만 full→compact 변환한다. 평균·clipping은 하지 않는다. Dxy,Dyx는 독립이며 절대 대칭화하지 않는다. 검사가 실패하면 output이 변하지 않는 transactional API로 한다. rate도 같은 layout을 사용하므로 adapter에는 m>0/열양성 같은 state 조건을 넣지 않고 별도 물리state 검증과 구분한다. Psi semantic mode는 호출자가 명시하고 유지하며 진단량을 자동으로 수송모멘트로 승격하지 않는다.

geometry의 canonical weight는 additive W다. `laguerre.h` 설명의 weight²와 연구 RHS의 W를 혼동하지 않는다. 실제 생산 식을 확인한 뒤 변환이 필요한지 정하고, 이름만 보고 sqrt/square 또는 음수W 거부를 하지 않는다. 첫 state adapter는 geometry 변환을 수행하지 않는다.

공유면 key는 stable ID 두 개의 정렬쌍 및 periodic image/subface 식별자다. minID cell을 소유한 rank가 공통 면좌표·normal·wn·quadrature·raw18 flux 및 필요한 Psi flux를 한 번 조립하고 양쪽에 반대 부호로 전달한다. 양쪽 rank가 각각 반올림한 독립 유량을 적용하지 않는다. ghost payload는 raw18,L,x,W,ID,frozen,stage/time/geometry epoch를 포함한다. halo는 density2-hop, quadratic velocity의 radial support3sqrt(area), geometry 필요범위의 합집합이다. MPI 경계를 frozen Dirichlet로 취급하지 않고 massroot는 전역 동일 문제로 구성한다. partition과 무관하게 같은 stable ID stencil에서 prior를 결정한다.

restart는 C struct memory dump를 쓰지 않는다. schema version,field order,dtype/endianness,ID,raw18,L,x,W,frozen,Psi mode/value,closure/version,integrator/time/dt,weight seed 및 boundary/gravity/energy 장부를 명시한다. seed는 재수렴할 초기추정이며 stale geometry/prior cache는 재생성한다. point-state checkpoint를 암묵 변환하지 않는다.

첫 검증은 하나의 C/CUDA adapter harness에 정확 roundtrip,비대칭D 보존,Q/Pi 비대칭 거부,absent Psi/NaN 거부,실패 시 출력불변을 묶는다. 실제 캡처 raw18과 wm_atom fixture를 사용해 CPU/OpenMP/CUDA 입력·출력을 대조한다. 아직 MPI 통신이나 물리kernel까지 한꺼번에 확장하지 않고 adapter와 계약 산출물을 먼저 평가한다. 계산은 허용된 lageunha에서 수행하고 기존 실행은 유지한다.

### Adapter 완료와 동적 초기 분류 변화의 해석

`canonical_moment_state_adapter_comparison_2026-10-18.json` 및 harness를 검토했다. 실제1590 raw fixture와 legacy atom,비대칭/absent/invalidmode/NaN/signed-zero/transaction 검사가 CPU1/OpenMP4/CUDA에서 통과하고 출력 SHA가 모두 `29721734110750b88c4f1da78c11e5f53063f7a8fd0d34d2b6dfc8625ece923b`로 같다. 해당 adapter 계층의 통과로 수용한다. full hydro 및 MPI 통신 통과는 아니다.

동적 후보의 첫 midpoint에서 nonflat이101에서 약1369–1382로 증가했다. 분류 수만으로 실패나 실제 접촉면의 급증을 판단하지 않는다. 매끄러운 영역의 작은 mean 차이도 roundoff 구간보다 커지면 local fraction/cut을 정의할 수 있다. contrast→0에서 prior가 constant로 수렴하면 큰 cut 개수가 큰 밀도 불연속을 뜻하지 않는다.

현재 실행을 유지하고 최종 공통시각 후처리에서 새 nonflat 셀의 `(high-low)/mean` 분포, 기존 nonflat과의 구분, theta/face trace/압력 변화 관계를 한 번 조사한다. 작은contrast에서 유한한 trace 변화가 남거나 dt 감소에도 진폭이 줄지 않는다면 시간매끄러움 문제로 검토한다. count를 줄이려는 phasefloor나 새 계수를 추가하지 않는다.

### 다음 생산 연산자: raw18 kinetic face-node flux

adapter 다음은 `Voro/Laguerre/moment_kinetic_face_flux.h`의 host/device 순수 node flux 함수다. 현재 지배적인 hydro 연산을 실제 공용 경로로 옮기는 단계이며 새 물리closure나 Psi 설계를 동시에 바꾸지 않는다. geometry/stencil/massroot 계산은 이번 범위 밖이다.

입력은 point x[2], outward normal n[2], scale=face length*quadrature weight, left/right velocity samples 각4개, signed sample massrates wl>=0/wr<=0, thermal specific eL/eR, Pstar,ustar_n이다. wl/wr에는scale이 아직 없다. fm=sum wl+sum wr, fp=sum wl vl+sum wr vr+Pstar n으로 두고 raw18 flux를

`(fm, x fm, xx^T fm, fp, x outer fp, sum wl vlvl^T+sum wr vrvr^T, sum wl(.5|vl|²+eL)+sum wr(.5|vr|²+eR)+Pstar ustar_n)*scale`

로 계산한다. Pi의 압력/체적 source는 별도 연산이며 이 node flux에 추가하지 않는다. D는 독립성분을 유지하고 Psi항은 이 함수에 없다. finite,scale>=0,wl/wr 부호를 검사하고 실패 시 output을 바꾸지 않는다. 단위normal/물리압력 admissibility는 상위caller의 계약이며 임의 허용치를 helper에 추가하지 않는다.

완료 contact RHS의 최종 수렴z와 동일 initialraw/factors/prior를 사용해 audit_weight_rates 고정 단일RHS로 face records를 캡처할 수 있다. full outer solve는 반복하지 않으며 mass residual/gen/fullrates가 기존 저장값과 일치하는지 확인한다. input SHA와 closure/state/z를 저장하고 기존 Python 산식으로 독립 expected flux를 조립한다.

한 harness로 CPU/OpenMP/CUDA 및 Python을 스케일 기반 오차로 대조한다. backend간 bitwise 일치는 달성 시 별도 보고하며 필수로 강제하지 않는다. 공유면 ±singleflux의 m,p,D,E 보존, Q/Pi 대칭과 D의 독립성, 균일속도 및 zero-relative transport의 압력 일을 함께 검사한다. 실제 node 순서의 셀별 면flux 합도 기존 records 합과 비교한다. 이 단계 통과를 전체 solver/MPI 통과로 확대하지 않으며 이후 massroot/factorcovterms/전체RHS 연결을 별도로 평가한다.

## Kinetic face-node 결과 평가 및 다음 생산 묶음

`canonical_kinetic_face_backend_comparison_2026-10-18.json`, 비교기와 harness를 읽었다. 실제213064 node/1590셀에서 CPU1/OpenMP4/CUDA synthetic 및 actual 검사가 통과했다. Python 대비 node 최대차3.46945e-18, 셀 합 최대차2.86229e-17이다. node absolute elementary-term bound의64eps 기준 대비 최대비율 .04349, 셀 누적 roundoff bound 대비 .009442다. CPU/OpenMP는 bitwise 같고 CUDA는 다르지만 선언한 오차범위에 있다. capture의 전체 RHS 재현차6.94e-17 및 mass 잔차1.09e-17도 입력 재현성을 지지한다.

해당 node algebra를 통과로 평가한다. 다만 비교기의 셀별 합산은 각 backend에서 출력한 node flux를 **Python에서 같은 순서로 합산**한 것이다. GPU의 실제 cell reduction이나 MPI 공유면 통신이 검증된 것은 아니다. synthetic의 `f+(-f)=0` 역시 동일 owner flux의 부호 적용을 확인하며 독립 MPI 보존을 증명하지 않는다. node bound와 cell bound는 현재 유한값·정상범위 입력에 대한 수치 비교 기준이며 모든 overflow/underflow 상태에 대한 보편적 오차 정리는 아니다.

### 다음 구현 묶음: backend 셀 조립과 pressure/volume raw18 source

새 중력/Psi 모형을 추가하기 전에, 통과한 node 함수를 **실제 CPU/OpenMP/CUDA 셀별 raw hydro RHS**로 연결한다. 이미 저장한 node input/owner/raw 및 complete RHS 배열을 재사용하며 Python RHS나 시간진화를 다시 실행할 필요가 없다.

1. owner 목록에서 cell별 CSR incidence를 만들고 각 행 안에서는 기존 전역 node 순서를 유지한다. node flux는 한 번 계산하고 CPU/OpenMP는 셀별 루프, CUDA는 우선 한 thread가 한 cell의 incidence를 순서대로 처리한다. atomics나 새 병렬 reduction 최적화는 이번 범위에 넣지 않는다. i쪽은−flux,j쪽은+flux다.
2. 같은 incidence에서 pressure traction을 사용해 `F=Σ sign scale Pstar n`, `Hface=Σ sign (point-c) outer (scale Pstar n)`를 조립한다. c=q/m이며 generator를 대신 쓰지 않는다. pressure work도 같은 node 입력에서 누적해 pressure 원시 분해를 반환한다.
3. 새 host/device `moment_hydro_cell_source.h`에 raw18,F,Hface,PV를 받는 순수 source 함수를 구현한다. `S=Q-qq^T/m`, `H=Hface+PV I`, `A=solve(S,H)^T`, `a0=F/m-Ac`를 구한다. source의 q=p,Q=D+D^T,D=Pi+PV I,Pi=a0 p^T+p a0^T+A D+D^T A^T이며 m,p,E source는0이다. 면 유량에 이미 pressure momentum/work가 있으므로 source에 F나 pressure work를 다시 더하지 않는다.
4. 2x2 solve는 양의 S를 확인한 scaled Cholesky 등 명시적인 안정 풀이로 구성한다. S를 임의 clipping/regularization하거나 실패 시 다른 행렬로 바꾸지 않는다. finite/mass/S admissibility 실패 시 output을 쓰지 않는 transactional 방식이며 A,a0 및 solve residual을 진단으로 반환한다. raw18 source와 legacy wm17/Psi adapter를 혼합하지 않는다.
5. pressure 원시 분해는 p=F,D=H+c outer F,Pi=위 pressure Pi source,E=면 pressure work,q/Q/m=0이다. 전체 hydro raw−pressure raw를 transport+kinematic-volume로 저장한다. 모든 셀의 unmasked RHS를 먼저 얻은 뒤 frozen exchange 장부를 수집하고 frozen raw rate를0으로 만든다. MPI 경계를 frozen으로 대체하지 않는 기존 계약은 유지한다.

### 한 번의 통합 검증과 통과 범위

같은 captured node 자료, 초기 raw/factors 및 PV를 사용해 CPU/OpenMP/CUDA에서 실제 CSR cell assembly+source를 수행하고 Python의 전체 hydro raw rates, pressure raw, transport remainder, F/H 및 frozen/exterior ledger와 비교한다. PV는 실제 `P*area` 값과 입력 gamma/열에너지 정의를 명시한다. 기존 저장물에서 필요한 값을 독립 산식으로 재구성할 수 있으면 그렇게 하고, 단순 capture 편의를 위한 새 full RHS 반복은 피한다.

검증에는 닫힌 polygon의 일정 압력 정지 상태에서 면압력과 D의 PV 체적항이 상쇄되는 fixture를 포함한다. source solve 잔차, Q/Pi 대칭과 D 독립성, pressure+remaining=hydro, frozen 처리 전후 장부를 같은 묶음에서 확인한다. tolerance는 node bound·셀 합산 오차·S solve의 조건수/잔차를 구분하여 제시하고 단일 큰 absolute tolerance로 숨기지 않는다. CPU/CUDA의 bitwise 일치를 강제하지 않는다.

이 단계가 통과하면 **주어진 geometry/reconstruction/face split/수렴 weight-rate 입력에서 raw hydro RHS를 세 backend가 재현한다**고 판정한다. factor의 K/incoming Qfactor, 중력/covariance source, nonlinear massroot, geometry/reconstruction, MPI와 시간진화는 여전히 별도다. 다음 단계는 그 결과를 본 뒤 factor evolution에 필요한 K/Qfactor와 중력 source를 연결하여 완전한 stage operator로 나아가는 것이다. 현재 동적 연구 실행들은 그대로 유지한다. 이번 문서는 다음 구현 지시이며 새 계산은 root가 수행한다.

### CSR/source 검토 및 비교 경계 수정

계획·평가 담당 Astra가 `moment_hydro_cell_source.h`, `moment_hydro_cell_harness.cu`, capture 준비 코드와 비교기를 정적으로 검토했다. scaled Cholesky의 전치 방향, pressure Pi source, p/E 중복 방지, CSR 부호 및 frozen/exterior 장부는 일관된다. 첫 보고의 GPU masked hydro 차3.361e-18, pressure 차6.615e-22, source 차1.231e-19는 작지만 최종 판정은 Pi 대칭 거부를 추가한 v2와 수정 비교기의 결과를 기준으로 한다. v1 산출물은 보존한다.

비교기의 `g/(1-eta)` 성분별 경계를 `solve(I-B,g)`로 수정했다. B=abs(invS)@dS, eta=norm_inf(B)<1이며 g에는 H와 S 섭동 및 열별 norm으로 배분한 solve 반올림 경계를 넣는다. 성분별 abs(A)에 비례하는 상대오차 가정은 작은 성분으로 유입되는 오차를 놓치므로 사용하지 않는다. H에는 centroid/point 계산 오차와 traction 곱의 절대합을 추가한다. 노드 elementary-term 오차, 순서 있는 합산, solve 섭동 및 backward residual 관문을 계속 분리한다. 이는 현재 정상범위 입력의 수치 비교이며 모든 부동소수점 입력에 대한 엄밀한 interval 인증은 아니다.

## 다음 생산 묶음: covariance transport와 중력을 연결한 고정 입력 stage RHS

v2 최종 artifact `canonical_hydro_cell_backend_comparison_2026-10-18.json`를 확인해 CSR/source 묶음을 통과로 판정한다. 실제1590셀/213064노드에서 CUDA masked hydro 최대차3.36103e-18, source1.23124e-19, pressure6.61499e-22, ledger1.18829e-16이며 backward residual1.71306e-16이다. 수정 경계 대비 hydro 최대비율 .009181, source .0001084다. CPU/OpenMP는 bitwise 일치하고 CUDA는 허용범위에서 다르다. binary SHA는 `55c89b21dd2149457ed80f3674721e3aafaaafb2b2c372f65e776249c1a683d3`다. 이 판정은 실제 backend 셀 조립을 포함하지만 아래에 남긴 전체 stage 및 진화 검증을 대신하지 않는다.

v2 CSR/source 통과 후 아래를 한 묶음으로 구현·대조한다. 목표는 기존 raw hydro 조립에 factor RHS 및 중력을 연결하는 것이다. 시간적분, geometry 재생성, massroot 및 MPI 구현은 이번 묶음에 포함하지 않는다. 기존 장시간 및 contact 두 진화는 그대로 유지한다.

### 입력과 구현 계약

1. 기존 node 순서/owner/CSR/raw18/PV를 재사용한다. 추가 입력은 각 면의 좌우 sample coordinate xi 각4×2, 각 receiver의 해당 point에서 affine velocity, cell factor L, 그리고 중력 체적 quadrature의 point/positive massweight/mean velocity/acceleration이다. xi는 quadratic reconstruction의 실제 좌표이며 상수 ALPHA로 대체하지 않는다. affine velocity도 quadratic sample 평균으로 대체하지 않는다. 기존 33-double node capture에는 이 정보가 없으므로 immutable 초기 상태와 동일 reconstruction에서 재구성하고 기존 velocity samples와 대조한다. 저장 자료만으로 복구되지 않는 항에 한해 기존 수렴 z를 고정한 단일 capture를 허용한다. outer 재수렴이나 새 진화는 하지 않는다.
2. 별도 host/device covariance transport helper에서 owner i의 K에 scale/m_i * sum(wl xi_l xi_l^T), owner j의 K에 scale/m_j * sum((-wr) xi_r xi_r^T)를 더한다. incoming Qfactor 열은 i에 sqrt(-scale wr_k)*(vr_k-vaff_i(point)), j에 sqrt(scale wl_k)*(vl_k-vaff_j(point))다. receiver별 원래 node/sample 순서를 유지하며 exterior j<0에는 수신 cell을 만들지 않는다. 0 rate는 0 열로 허용한다. signed rate를 임의 abs 처리하지 않는다.
3. K와 Qfactor 자체를 반환하고 Qfactor@Qfactor.T도 독립 합산한 양의 incoming covariance와 비교한다. covariance만 저장하고 factor를 임의로 재분해하여 기존 factor gauge를 바꾸지 않는다. hydro covariance generator는 기존 B=solve(S,C)^T, C=D-q p^T/m에서 -B다. pressure source에서 구한 A와 이 B를 혼동하지 않는다.
4. 중력 node source는 dm을 포함하여 p=dm*a, D=dm*x outer a, Pi=dm*(v outer a+a outer v), E=dm*v dot a, 나머지0이다. 같은 cell에서 정해진 순서로 합산한다. 가속도 값은 검증된 기존 q12 quadrature 입력으로 공급하며 이번에 spline force 함수를 별도로 이식하지 않는다. frozen도 평가값은 보존하되 적용률과 gravity 장부에는 active만 포함한다.
5. 중력 covariance source는 dC=gD-q outer gp/m, dT=gPi-(gp outer p+p outer gp)/m, G=dT-B dC-dC^T B^T다. 물리적으로 대칭인 연산은 대칭 결과를 구성하고 대칭화 전 잔차도 진단한다. R=L L^T에 대해 대칭 Ag R+R Ag=G를 푸는 host/device helper를 연결한다. 엄격한 SPD를 요구하고 clipping/regularization은 하지 않는다. scaled 3×3 대칭 Sylvester 풀이 등 명시적인 작은 선형계를 사용할 수 있으며 L의 SVD gauge와 같은 Ag를 낼 필요는 없고 고유한 대칭 해와 residual을 비교한다.
6. stage 출력은 raw hydro+active gravity, cov_A=-B+Ag, K, Qfactor, 분리된 exterior/frozen/gravity 장부다. frozen의 raw/cov_A/K/Qfactor는 적용 단계에서0으로 만들고 unmasked 진단과 구별한다. generator 위치율과 additive weight율은 기존 고정 입력 metadata로 보존한다. Psi 상태나 총 물리 에너지 보존을 암묵적으로 추가하지 않는다.

### 단일 대조 묶음과 판정

CPU1/OpenMP4/CUDA가 실제 CSR 및 volume reduction을 수행한 결과를 Python 독립 산식과 저장 stage 자료에 비교한다. K/Qfactor/incoming covariance, raw gravity, G/Ag 및 최종 raw/cov_A/장부를 함께 저장한다. primitive 절대합과 합산 길이, S와 R의 조건수 및 source/Sylvester residual을 분리해 허용치를 정한다. 작은 성분에 단일 상대오차를 적용하거나 큰 절대 floor를 넣지 않는다.

추가 물리 항등식은 acceleration-only에서 E−trace(Pi)/2의 source가0, 중심력 quadrature의 원점 토크가0, affine velocity의 gravity residual covariance source가0임을 확인한다. 내부 중심 토크 자체는0을 요구하지 않는다. frozen 및 exterior 처리도 현재 장부 정의로 검사한다. zero incoming/zero acceleration과 일반 비대각 SPD fixture를 같은 harness에 포함하고 rank-deficient R은 명시적 실패로 다룬다.

이 묶음 통과는 **고정 reconstruction/면 운동/중력 quadrature에서 raw 및 covariance stage RHS를 세 backend가 재현함**을 뜻한다. factor exponential/QR/Sylvester 시간 update와 Pi rebuild, 전체 moving geometry·massroot, MPI 및 장시간 안정성은 여전히 미검증이다. 다음은 그 고정 stage 출력을 입력으로 기존 midpoint factor update+Pi rebuild를 공용화한 뒤, geometry/reconstruction/massroot 연결 및 실제 분할 MPI 일관성으로 전진한다. 짧은 research 진화의 결과는 후보 admissibility에 사용하고 t6–10 생산 목표를 대체하지 않는다.

### Stage 정적 검토에서 발견한 장부 부호 문제

`moment_stage_harness.cu`의 초기 gravity ledger와 `prepare_moment_stage_reference.py`가 함께 active gravity source의 음수를 기록했다. 기존 runner는 `dy[-4:]+=sum(active g)`이고 진단은 `totals-initial-ledger`이므로 연결 계약상 중력 장부는 양수다. 이는 reference/backend가 같은 실수를 공유하면 수치 대조만으로 발견되지 않는 사례다. 초기 결과는 보존하고 양쪽을 수정한 뒤 재평가한다. exterior/frozen 음수 부호는 그대로이며 중력 부호와 혼동하지 않는다.

최종 관문에 `sum(active raw p,E,Dxy-Dyx)=exterior+frozen+gravity ledger`를 추가한다. active masked K/incoming 출력도 unmasked 결과와 별도로 대조한다. transport xi/receiver 선택, G 및 대칭기저 Sylvester 연산은 정적 검토상 일관된다. affine G0 fixture는 source를 대수적으로 구성한 항등식 검사이며 독립 체적적분 검증과 구별한다.

### 기존 smooth 연구모델의 t=.001 연장 완료 평가

`canonical_quadratic_extended_comparison_2026-10-18.json`의 ready=true를 확인했다. dt1e-4/5e-5 두 분기 모두 동일 t=.001에 도달했고 pressure와 factor covariance가 양수다. fine minP2.77560838e-5, minR7.51852053e-13, hydro ledger 오차1.2172e-16이며 frozen 상태는 두 분기에서 동일하다. 두 해의 active mass-weighted velocity 차8.52427e-9, centroid 차4.81088e-9, specific thermal 차9.55015e-12다.

이는 기존 smooth 모델의 짧은 연장 admissibility 및 두 시간간격의 근접성을 보여준다. pressure RMS는 coarse1.17198e-6, fine1.14100e-6으로 비슷하므로 두 값만으로 시간차수 또는 공간오차 지배를 확정하지 않는다. 두 분기는 t1e-4 restart에서 연장되었으며 steps9/18은 연장 구간의 step 수다. 재구성 potential 변화 약−5.7e-9에는 열린 reservoir 교환이 섞이며 총에너지 보존오차가 아니다. 진행 중인 contact prior 모델의 동적 결과와 섞지 않고 별도 기준 자료로 보존한다. 이 결과를 근거로 추가 smooth 분기를 시작하지 않는다.

## Stage v2 통과와 다음 factor update·Pi rebuild 묶음

`canonical_moment_stage_backend_comparison_2026-10-18.json`의 v2 결과를 직접 확인해 고정 입력 stage RHS 구현을 통과로 판정한다. 1590셀/213064면노드/1431792중력노드에서 raw total 최대차1.12757e-17, K5.68434e-14, Qfactor0, G2.61355e-17, Ag2.00156e-7이며 각 선언 경계 안에 있다. Sylvester backward residual 최대2.17131e-16, 독립 전체 장부 항등식 차3.73535e-16이다. masked K/incoming도 직접 대조했다. CPU/OpenMP는 bitwise 같고 CUDA는 다르지만 수치 관문을 통과한다. binary SHA는 `604468ca5a6bf3cc18bd0f2a9e487931f4bfbd433ef6cb4c8615e42ec448c822`다.

Ag 차2e-7은 숨기지 않는다. R 조건수 최대5.11, Sylvester 조건수6.60이므로 큰 조건수에 의한 불안정으로 단정할 수 없다. 작은 R의 절대 척도와 G 원시모멘트 상쇄오차의 증폭을 구분해야 한다. 다음에는 dt*deltaAg와 최종 R/Pi/열분할 차이로 영향을 평가한다. 현 stage 장부가 총 물리 에너지 보존을 증명한다는 뜻은 아니다.

### 구현 대상과 고정 입력 계약

다음 한 묶음은 host/device factor update와 Pi rebuild다. 기존 captured stage와 backend 결과를 재사용하며 새 RHS·massroot·진화를 실행하지 않는다. 입력은 Lbase, Llinear, cov_A, K, 순서 있는 Qfactor 열, h>=0, 그리고 provisional raw18이다. Llinear는 stage의 K를 정의한 factor이고 Lbase는 갱신 출발 factor다. 실제 midpoint에서 half update는 둘 다 L_n, full update는 Lbase=L_n이고 Llinear=L_half다. 두 인자를 합쳐 버리면 기존 integrator와 달라진다.

1. Rlinear=Llinear Llinear^T, J=Llinear K Llinear^T에서 대칭 loss T가 T Rlinear+Rlinear T=J를 만족하도록 계산한다. 기존 Python SVD 식의 loss와 동일한 해다. 현재 SPD 입력에서는 검증한 대칭기저 Sylvester 풀이를 재사용할 수 있다. Ac=cov_A−T다. rank-deficient Llinear는 별도 구현 전까지 명시적 unsupported로 거부하고 Python의 rank 처리까지 이식했다고 주장하지 않는다. frozen은 계산 전에 bypass하여 state/factor를 정확히 보존한다.
2. Ehalf=exp(h*Ac/2), Efull=Ehalf*Ehalf, U=[Efull*Lbase, sqrt(h)*Ehalf*Qfactor]를 구성한다. incoming은 Ehalf로 곱하고 loss를 이미 Ac에 넣었으므로 별도 right exp(−hK/2)를 다시 곱하지 않는다. h=0은 정확 복사한다. 행렬지수는 비대칭 2×2를 지원해야 한다. 구현 선택으로 norm_inf를 .5 이하로 줄이는 scaling, 18차 Taylor 다항식, 반복 squaring을 사용할 수 있다. 항 누적과 squaring 오차·finite를 진단하고 과도한 scale/overflow를 조용히 허용하지 않는다. scipy expm이 독립 기준이다.
3. U^T의 두 열에 안정적인 thin QR를 적용하고 Lnew=Rqr^T로 저장한다. CSR incoming 순서를 유지하며 CPU/OpenMP/CUDA에서 실제 QR를 수행한다. QR 대신 U U^T를 먼저 만든 뒤 일반적인 임의 제곱근을 선택하면 factor gauge가 달라지므로 이번 이식에서는 그렇게 하지 않는다. 최종 factor를 하삼각·양의 대각으로 정규화하는 것은 열 부호 변경만 허용한다. Python QR도 같은 열 부호로 정규화하여 비교한다. 현재 네 sample은 ±축 쌍이므로 열 부호 변경은 sample 순열에 해당하지만 일반 직교회전은 같은 이산 sample ensemble을 보장하지 않는다. 아직 진행 중인 Python 진화 factor를 변환하지 않는다.
4. provisional raw의 S=Q−qq^T/m, C=D−q p^T/m, B=solve(S,C)^T를 새로 구하고 Pi=pp^T/m+B S B^T+Lnew Lnew^T로 재구성한다. 나머지 raw 성분과 Ehyd는 유지한다. 이전 stage의 B를 재사용하지 않는다. 기존 rebuild와 같은 의미이며 raw Pi의 provisional 값과 차이, 그에 따른 −trace(deltaPi)/2의 열분할을 저장한다. m/S/finite/열 admissibility 실패는 보존하고 clipping·열 재주입을 하지 않는다. 출력은 transactional로 한다.

### 한 번의 검증과 다음 연결 기준

실제 저장 stage를 사용한 h=5e-5 및1e-4의 고정 계수 map을 하나의 harness에서 CPU1/OpenMP4/CUDA와 Python에 대조한다. 이 계산은 새 midpoint stage를 평가하지 않으므로 실제 한 timestep 또는 시간차수 시험이라고 부르지 않는다. 고정 입력에서 loss/Ac/Ehalf/Efull/Lnew/Rnew/Pi/thermal partition을 저장하고 backend stage 오차를 제거한 공통 입력 비교와 각 backend stage 입력을 이어 붙인 비교를 분리한다. 후자의 Ag 오차 전파를 최종 R/Pi/열 단위로 보고한다.

provisional은 고정 mass 계약을 따른다. 실제 runner는 별도 mass 배열과 dm 없는 flat17 갱신을 사용하므로 raw18 deck도 m을 초기 입력 그대로 유지하고 다른17성분에만 h*rate를 적용한다. 수치 massroot 잔차는 stage 진단에 남기며 h*dm을 상태에 적용하지 않는다. 처음 모든18성분을 갱신한 deck/결과는 보존하고 고정 m deck을 v2로 대조한다. 또한 CPU↔CUDA 연결 차이만으로 독립 Python stage 대비 Ag 약2e-7의 영향을 평가할 수 없으므로, 저장된 원래 Python stage에 동일 Python map을 적용해 CPU stage 기준 R/Pi/열 변화도 별도로 후처리한다.

fixture는 h0/frozen 정확 보존, 비대각 SPD와 비대칭 Ac, zero incoming, Lbase≠Llinear, 회전 아닌 열 부호 정규화, 실패 transactional을 같은 묶음에 둔다. loss Sylvester residual, expm 오차, QR 재구성 및 orthogonality, Rnew=U U^T, Pi 대칭과 보존 성분 불변을 검사한다. 하나의 큰 absolute tolerance 대신 입력 섭동·조건수·행렬지수 및 QR 오차를 구분한다. 현재 고정 계수에서 열양성을 통과해도 evolving geometry에서 열양성이 보장된다고 일반화하지 않는다.

이 묶음 뒤에는 reconstruction/면 운동과 massroot를 Python에서 제공하더라도 compiled stage 및 update를 실제 두-stage 호출에 연결할 수 있다. 그 연결은 별도 승인된 결과 기반 계획으로 정하고, 현재 contact 두 진화를 중단하거나 동일한 미세 진화 분기를 추가하지 않는다. 생산 moving geometry·전역 weight solve·MPI 공유면 및 장시간 Kepler 관문은 계속 남는다.

## Factor v2 최종 판정과 실제 midpoint 한 step 연결

`canonical_factor_update_backend_comparison_2026-10-18.json`과 `canonical_factor_original_stage_propagation_2026-10-18.json`을 직접 확인했다. 고정 mass v2의 공통 입력 및 backend 고유 입력 map이 통과했고 CPU/OpenMP는 bitwise 일치한다. loss backward 최대2.19085e-16, Ehalf 차4.44089e-16, Efull 차8.88178e-16, Lnew1.49078e-19, Rnew4.80799e-24, Pi5.42101e-20, 열분할1.08420e-19, QR orthogonality residual2.39693e-14이며 최종 열 최소3.86742e-8이다. binary SHA는 `faa57b9fefd236d879f98751cacfd2019bae2e2fb0352a4bca62dca673037c23`다. 이 범위에서 factor update와 Pi rebuild를 통과로 평가한다.

원래 Python-stage Ag 차2.00156e-7의 별도 전파는 h5e-5/1e-4에서 R 최대차1.30653e-21/2.61251e-21, Pi2.71051e-20/2.16840e-19, 상대열4.16252e-14/2.67835e-12다. CPU↔CUDA 차이만 비교한 결과와 구별했으며 두 고정 계수 map에서 실제 상태 영향이 작다는 근거다. 누적오차·시간차수·장시간 Kepler 안정성을 보장하지 않는다.

### 다음 한 묶음의 목적과 실행 범위

다음은 contact closure의 동일 초기 상태, active r<.60, dt1e-4에서 **실제 explicit midpoint 한 step**을 CUDA stage와 factor update가 구동하도록 연결하는 단일 검증이다. 기존 연구 진화와 A/B 실행을 유지한다. Python provider는 geometry/prior/density/quadratic reconstruction/face split/weight massroot/중력 quadrature를 맡고 CUDA가 raw/covariance stage 및 factor/Pi update를 실제 수행한다. 이를 전체 생산 CUDA solver나 완전한 geometry 이식으로 부르지 않는다.

각 stage의 공통 입력을 CPU/OpenMP/Python에도 제공하는 shadow 비교를 같은 작업에 넣는다. 서로 다른 네 backend trajectory를 새로 진화시키지 않는다. shadow는 CUDA가 선택한 실제 stage 상태에서의 연산 비교이며 backend별 독립 midpoint 궤적 검증은 아니다. 필요 입력을 모두 한 번에 capture하여 후속 수식 진단 때문에 full RHS를 반복하지 않는다.

### 두 stage와 상태 commit 계약

1. 시작 상태 raw18/L/x/W/mass/IDs/frozen/Psi mode와 code/input SHA를 기존 contact t0와 대조한다. stage0는 동일 초기 입력의 기존 capture를 재사용할 수 있으나 실제 입력 SHA와 generator/weight seed 일치를 확인한다. cache key에는 state/factor/geometry epoch를 포함하고 이후 midpoint에는 초기 capture를 재사용하지 않는다.
2. CUDA stage0에서 raw hydro+gravity, cov_A/K/Qfactor, exterior/frozen/gravity 장부를 조립한다. generator와 Wdot는 provider의 무조향 결과다. half provisional은 고정 m, 기타 raw17 및 x/W를 시작 상태+h/2*stage0로 갱신한다. factor는 Lbase=Llinear=L_n으로 h/2 갱신하고 새 Pi를 rebuild한다. frozen raw/L/x/W는 원래 값 그대로 복사한다.
3. CUDA half-state를 Python provider로 전달하여 t+h/2의 moving geometry→local prior/density→quadratic samples→generator/face cuts→weight outer solve를 실제로 새로 수행한다. 수렴 seed만 재사용하며 초기 z나 face nodes를 고정하지 않는다. density는 해당 stage 내에서만 재사용하고 outer 변화<1e-11, mass 잔차 및 fit/admissibility를 기록한다. gravity quadrature도 새 상태에서 구성한다. 이 시점이 고정 계수 map과 실제 midpoint의 핵심 차이다.
4. CUDA midpoint stage 출력을 이용해 full provisional을 **초기 상태**+h*midpoint rate로 계산한다. factor는 Lbase=L_n, Llinear=L_half, midpoint cov_A/K/Qfactor로 h 갱신한다. half 상태에 h*rate를 더하거나 stage0 generator/Wdot를 full step에 쓰지 않는다. 최종 raw Pi는 full provisional에서 새 B와 Lnew로 재구성한다. 장부도 초기 장부+h*midpoint 장부로 갱신한다.
5. 최종 finite/SPD/thermal/mass/frozen/장부 관문을 통과한 뒤 한 번 commit한다. 실패 시 입력·stage·provisional·factor·상태를 저장하고 clipping, timestep 변경 재시도, smooth fallback 없이 중단한다. 정확한 평가 시각과 accepted 시각을 구분한다. 최종 진단에 geometry/density가 필요하면 계산하되 새 최종 RHS나 두 번째 step을 시작하지 않는다.

### 연결 검증의 판정과 산출물

두 실제 stage에서 기존 stage 항목의 backend shadow 차이·오차 경계, massroot convergence, 재구성 p/D/Pi budget, prior fit/theta/nonflat, generator 및 additive Wdot, pressure와 residual covariance를 저장한다. x/W는 provider 입력이므로 그 transport가 CUDA kernel에서 계산됐다고 주장하지 않는다. frozen은 bitwise 동일, m은 원래 값과 bitwise 동일이어야 한다.

각 half/full update는 **동일 Lbase/Llinear/provisional/raw stage 입력**으로 독립 Python SVD-loss/expm/QR/rebuild와 비교한다. 여기서 half-state가 달라지는 독립 Python trajectory를 추가로 계산하지 않는다. 두 update에서 Pi rebuild의 열분할과 provisional 최저열, factor positivity, 잔차 및 backend 오류 전파를 보고한다. QR 부호 정규화는 기존 계약을 유지한다.

stage별 `sum(active pdot,Edot,Ldot)=exterior+frozen+gravity` 및 accepted 상태의 `total-initial-ledger`를 검사한다. 이 장부는 hydro E와 외력 일에 관한 것이며 Psi transport가 없는 현재 총 물리 에너지 보존 판정은 아니다. 이미 진행 중인 Python contact coarse의 t1e-4 checkpoint가 있으면 최종 물리량을 후처리로 비교할 수 있지만 provenance와 factor sign 정규화 및 독립 trajectory 차이를 명시하고 단순 bitwise 관문으로 사용하지 않는다.

통과하면 연구 Python provider와 compiled stage/update를 연결한 실제 moving midpoint 한 step이 확인된다. 이후 우선순위는 이 인터페이스를 유지한 생산 geometry/weight-rate 경로 연결과 MPI shared-face/ghost/global massroot 분할 일관성이다. 물리 closure 개선 판단은 진행 중 contact 두 분기의 동일 종료시각 결과를 기다리며, 이 연결 검증을 새 장시간 안정성 증거로 대신 쓰지 않는다.

## 기존 auxiliary-memory A의 t10 최종 판정

원래 실행 A는 중단·재시작 없이 t10.000069146956758, runner exit0으로 완료됐다. `canonical_memory_physical_final_2026-10-18.json`의 최종 snapshot SHA는 `b22a2a00cb2b84760769447debcbbbcf2bb6bf3f054b120be71a37ef20ab4639`다. 51개 완전 frame이 분석됐고 최종 mass는 불변, weights는0이다.

정확도 목표는 실패로 판정한다. inner_background804셀에서 상대 pressure mass RMS6.31317, wrapped phase RMS64.076도, radial displacement RMS .0309785, 정확 원궤도 위치 RMS .349667이다. inner_shear_band521셀은 pressure7.30911, phase74.0666도, 위치 .418067로 더 큰 지역 오차를 보인다. 반면 disk_core10304셀의 pressure RMS .00472991와 phase .00386631도는 작다. 큰 core 영역의 좋은 값으로 내부 저밀도/전단 영역 실패를 평균 내어 숨기지 않는다.

phase는 modulo2pi로 감긴 값이므로 궤도 위치·속도 오차와 함께 읽는다. 실행 완료와 finite/양성 유지가 전체 장시간 Kepler 정확성 통과를 뜻하지 않는다. A는 기존 auxiliary-memory 모델이며 새 독립 moment/contact 모델의 물리적 검증 결과로 전용하지 않는다. 이 결과를 보존하고 A를 새 파라미터로 재시작하지 않는다. B와 이미 진행 중인 contact 및 compiled midpoint 실행은 기존 계획대로 유지한다.

## 실제 midpoint 연결 최종 판정

`canonical_compiled_midpoint_run_2026-10-18.log`와 `canonical_compiled_midpoint_endpoint_2026-10-18.json`을 확인했다. CUDA가 실제 stage/update를 구동한 .60 active1194셀의 dt1e-4 한 midpoint step이 완료됐고 stage0/midpoint 및 half/full shadow 관문을 통과했다. 최종 minR7.51772837e-13, minthermal3.86742499e-8, minP2.77762092e-5, pressure relative RMS6.02093e-8이다. accepted 장부 최대잔차2.12047e-16이다.

이미 진행 중인 독립 Python contact coarse의 동일 t1e-4 checkpoint와 저장 상태만 비교했다. active 최대차는 raw6.50521e-19, velocity4.44089e-16, centroid5.55112e-17, ell1.08420e-19, R3.27840e-21, thermal3.25261e-19, pressure2.21940e-16이다. generator는 동일하고 W 차1.21232e-20이다. factor 원소 대신 LL^T를 비교했다. 해당 한 step의 연결 구현을 통과로 판정한다. geometry/reconstruction/massroot를 Python이 공급한 범위이며 전체 GPU/MPI나 장시간 정확성 통과는 아니다.

## 다음 생산 통합 우선순위: 기하 → weight 제약 → MPI 공유면

계획·평가 담당 Astra가 다음 순서를 결정한다. 먼저 생산 기하에서 정확한 moment 입력을 내보내고, 그 위에서 kinetic massroot 연산을 공용화한 다음, 전역 weight 제약과 shared-face를 MPI 분할에 연결한다. 현재 contact 두 진화와 B는 계속 유지하며 추가 미세 진화 분기는 만들지 않는다. Re=infinity, 무조향 generator, 고정 mass, 상대 W의 동일 기하 의미를 유지한다.

### 첫 구현 묶음: 생산 polygon export와 엄격한 moment geometry adapter

대상은 실제 생산 `Exam/exam_gpu.cu`의 `voronoi_tessellate_kernel`, `Exam/exam_gpu.h`의 PolygonOut/SoA 및 해당 생산 host 호출 경로다. CPU의 관련 경로는 `Voro/voro.c:Voro2D_FindVC`다. `Voro/Laguerre/laguerre_cuda.cu`의 example용 결과만 검증하고 생산 연결로 주장하지 않는다. 첫 작업에서 실행 대상 Makefile/전처리/실제 호출 symbol을 기록해 이 경로를 재확인한다.

현재 생산 GPU 경로를 읽어 확인한 차이는 pw2 입력 float, 3×3 neighbor 검색과 최대256개, 최대48 corner, 상대 작은 edge 병합 및 volume floor다. 병합 주석과 실제 상수도 다르며 실제 기준은 maxedge의1e-3이다. 이 처리를 거친 polygon/volume을 정확 moment geometry로 그대로 쓰지 않는다. 새 moment 선택 경로는 double x/W/Wdot sidecar와 명시적 status를 가지며 legacy particle bytes를 암묵 변경하지 않는다. W 필드 이름이 w2여도 실제 power plane은 `2*(xj-xi)·(point-xi)=|xj-xi|²+Wi-Wj`이므로 additive length² 자체를 전달한다. sqrt(W), 음수 W 거부, clamp, radius derivative로의 자동 변환은 금지한다.

첫 구현은 다음 계약으로 한다.

1. 입력: stable global IDs, double generator x[2], signed additive W, 해당 snapshot의 무조향 generator velocity와 Wdot, physical frozen/active 역할, geometry epoch, global exterior generator registry 또는 완전성이 확인된 candidate CSR. local index와 global ID mapping을 저장한다. patch 밖 실제 이웃의 global ID를 버리고 단순 -1만 남기지 않는다. physical reservoir와 MPI remote neighbor는 서로 다른 종류다.
2. 생산 clipping 산식을 추출해 host/device에서 재사용할 수 있는 helper 또는 생산 kernel의 명시적 moment mode로 구현한다. 검사 가능한 CPU/OpenMP/CUDA polygon export API를 만든다. 결과는 CCW vertices, 각 edge의 neighbor global ID, outward normal, length, 정확 polygon area·centroid·geometric second moment, status와 epoch다. double 산술과 additive W를 유지하고 밀도/압력 갱신은 이 geometry API에 넣지 않는다.
3. moment 경로에서는 상대 tiny-edge 병합 및 volume floor를 적용하지 않는다. 정확 중복점 제거처럼 같은 polygon을 표현만 바꾸는 처리는 근거와 함께 가능하다. ambiguous/empty/coincident/hidden/capacity 실패는 status로 반환한다. neighbor/corner cap에 도달했는데 일부 polygon을 성공으로 반환하는 경로는 금지한다. 필요한 경우 cap을 검사하고 입력 크기에 맞는 재할당·재실행으로 해결하며 물리 상태를 바꾸지 않는다. 실패 후 caller 출력은 commit하지 않는다.
4. neighbor 완전성은 기존 3×3라는 이유만으로 가정하지 않는다. 첫 두 저장 snapshot에서는 전체 generator를 대상으로 모든 출력 vertex의 power halfspace 위반을 tiled 방식으로 독립 확인한다. candidate 누락이 발견되면 누락 ID를 추가해 동일 snapshot의 기하만 재구성한다. 기존 global regular triangulation neighbor 목록을 참조 입력으로 사용하는 비교는 허용하지만 그 경우 native neighbor search 통과라고 주장하지 않는다. 향후 생산 검색은 공간 bound와 W 범위를 이용한 완전성 조건을 별도로 충족해야 한다.
5. 공유면 key는 정렬된 global ID 쌍+periodic image/subface다. 양쪽 polygon의 plane/segment 일치와 opposite normal을 확인하고 owner가 한 geometry를 선택해 양쪽 사용자가 공유한다. 불일치한 양쪽 segment를 평균해 숨기지 않는다. floating-point 표현 차이는 geometry uncertainty와 plane residual로 평가하고 실제 topology 차이는 따로 기록한다. zero-length contact를 유한 면으로 바꾸지 않는다.
6. face 정상속도는 actual point에서 `wn=((xj-xi)·vj-(vj-vi)·(point-xi))/d+(Wdot_i-Wdot_j)/(2d)`로 계산한다. Python provider와 기존 생산 `dev_get2dUpqradRk4`의 GFS_W2_LINEAR 분기를 수식으로 비교하되 float 입력과 tangential convention 차이를 분리한다. legacy sqrt/clamp 분기를 새 moment 경로에서 호출하지 않는다. Wdot는 실제 물리 weight derivative이며 restart seed와 구분한다.

### 기하 검증 자료와 중단 관문

완료한 compiled midpoint의 stage0와 실제 half-state 두 저장 입력을 사용한다. half-state에는 비영 signed W가 있으므로 zero-weight 정규격자만 확인하는 검사를 대신한다. global exterior generator는 동일 IC registry를 사용한다. 새 RHS·outer solve·진화는 필요 없다. CPU/OpenMP/CUDA polygon과 연구 reference polygon을 비교하고 area/centroid/second moment, power plane/공유면, normalsum 및 전역 halfspace 완전성을 한 보고서에 묶는다. 상대 W를 공통 상수만큼 바꾼 입력의 plane/geometry 불변성은 같은 harness fixture로 검사한다. geometry tolerance는 입력좌표 roundoff 및 clipping 연산 규모로 정하고 기존1e-3 merge 규모를 허용치로 사용하지 않는다.

cap 초과, 누락 neighbor, 비양의 면적, unresolved edge ownership, 일치하지 않는 유한 공유면은 provider 연결 전 중단 관문이다. Python clip으로 조용히 fallback하지 않는다. 실제 production callsite와 선택 mode를 구현하고 snapshot harness가 그 동일 helper를 호출해야 한다. 검증을 위해 별도 clipping 알고리즘을 새로 만들어 생산 연결처럼 제시하지 않는다. 해당 묶음이 통과한 뒤에만 Python provider의 geometry callback을 native export로 바꾸는 연결을 검토한다.

### 두 번째 묶음: kinetic weight-rate 연산과 전역 제약

기하 계약 후 `vector_kinetic_root.py`/`frozen_weight_root.py`의 현재 정의를 공용 host/device 연산으로 옮긴다. 같은 reconstructed face samples/rho/scale/basewn/d를 입력으로 wn(z), wl/wr, outward mass residual 및 Laplacian coefficient를 한 owner 면에서 계산한다. `d residual/dz = -L`이므로 Newton은 L*step=residual, z+=step이다. kink의 반미분(.5), physical frozen z0와 j<0 외부면의 z 비의존성도 그대로 유지한다. 면 quadrature의 velocity-zero split이 z에 의존하므로 기존 outer 수렴을 제거하지 않는다.

첫 native 단계에서는 sparse linear solve를 기존 검증 Python solver에 남길 수 있으나 residual/J·업데이트가 실제 native 경로인지 명확히 구분한다. 이후 sparse backend를 선택할 때 anchored component와 nullspace를 검사하고 미지의 gauge를 임의 pinning하지 않는다. 저장된 두 stage에서 입력/J/수렴 z/wn/mass 및 전체 RHS를 비교해 결과 기반으로 연결한다. arbitrary fixed-z replay를 실제 mass 제약 solve 통과로 부르지 않는다.

### 세 번째 묶음: MPI partition 일관성

native geometry/face/weight contract를 기준으로 stable ID 기반 owner/ghost와 epoch 검사를 구현한다. 같은 face flux 및 Jacobian 항을 owner가 한 번 생성하고 두 cell에 ±로 전달한다. halo는 geometry 완전성, density2-hop, quadratic radial stencil을 함께 충족한다. MPI 경계를 physical frozen으로 대체하지 않으며 전역 active massroot를 같은 문제로 푼다. rank별 독립 local massroot를 금지한다.

첫 MPI 관문은 동일 저장 stage의1/2/4 rank 분할(접촉면과 active/frozen 경계를 가로지르는 분할 포함)에서 ghost/state/면 key·원시유량·전역 residual/J·stage 장부를 대조하는 한 묶음이다. fixture를 위해 frozen/exterior를 rank 경계로 재정의하지 않는다. 단일 rank 결과와 tolerance 내 일치, duplicate/missing face 없음, collective 실패 전달을 확인한 뒤 실제 두-stage 연결로 나아간다. 이때도 geometry/reconstruction을 global Python reference에서 제공한 범위와 native 분산 구현 범위를 구별한다.

Psi transport closure와 contact 장시간 물리 검증은 별도로 미완료다. 생산 통합 작업으로 이 물리적 미폐쇄를 해결했다고 주장하지 않는다. 현재 contact .001 두 분기의 완료 평가가 나오는 즉시 다음 장시간 후보 판단을 갱신하되, 위 순서의 첫 구현은 지금 지정한 geometry adapter다.

### 외부면 속도 계약 보충

현재 Python `kinetic_factor_rhs.py`의 j<0 hydro 면은 boundary_velocity=0의 mirror wall이며 z에 의존하지 않는다. 이는 외부 generator와의 power-plane 기하 속도 정의와 별개다. export는 geometry_wn과 hydro_bc_wn, 양쪽 global IDs 및 역할을 구분해 보존한다. 내부 공유면은 위 plane 속도식으로 직접 비교한다.

이번 .60 patch에서는 외부면의 patch쪽 cell도 모두 frozen, 외부 registry의 v/Wdot도0인지 확인하고 그 결과 geometry_wn=hydro_bc_wn=0을 검사한다. active에 닿는 외부면 또는 비영 geometry_wn이 발견되면 기존 고정 wall BC와의 등가 연결을 중단해 평가한다. 외부면을 비교에서 제외한 채 moving-domain/GCL 불일치를 통과시키지 않는다.

## 기하 실패의 확정 원인과 API v2 결정

`canonical_production_geometry_exact_all_registry_2026-10-18.json`을 확인했다. 원래 binary64 x/W를 정확한 dyadic 값으로 해석한 유리수 clipping에서 owner3은 후보 누락이 없고 owner5는 global ID2의 halfspace가 누락됐다. ID2를 추가한 뒤 두 cell의 모든 vertex가 전체36864 generator의 정확 halfspace를 만족한다. 3–5 면은 정확히 양의 길이이며 근사 길이4.860850672431795e-16, global 유리수 endpoint가 양쪽에서 역순으로 정확히 일치한다. 따라서 이상적인 원형 IC의 공점이라고 가정해 이 면을 없앨 수 없다.

v1에는 서로 다른 두 문제가 확인됐다: Qhull sparse candidate의 정확 topology 누락, double SH 반복 교차점의 역방향 미세 segment 생성이다. 기존64eps*좌표규모*plane수 tolerance로 위반을 허용한 all-registry 검사는 정확 topology 인증이 아니었다. v1 실패·입력·라이브러리는 보존한다. 기존 연구의 수치 RHS 대조와 한-step 결과는 해당 연구 geometry에 대한 결과로 유지하지만 정확 power topology를 보장했다고 해석하지 않는다.

### 결정한 구현: support-plane 표현과 filtered exact resolution

1. 별도 API/schema v2를 만들고 vertex마다 두 support-plane IDs를 저장한다. cell i의 j plane은 원래 x/W에서 정의한 정확 `n=xj-xi`, `b=(|xj-xi|²+Wi-Wj)/2`, `n·(point-xi)<=b`다. 부동소수점으로 한 번 계산한 dx/d²를 exact 원본으로 재해석하지 않는다. box plane도 명시적인 별도 ID를 갖는다. outgoing edge의 neighbor ID와 vertex support를 함께 유지하고 후보 순서는 stable ID 순서다.
2. topology 판정은 이전 단계에서 반올림한 vertex 좌표를 다음 plane에 넣어 결정하지 않는다. 두 support plane의 intersection을 세 번째 plane에 대입한 determinant 부호로 결정한다. D=n1x*n2y−n1y*n2x, X=b1*n2y−n1y*b2, Y=n1x*b2−b1*n2x이면 검사량은 `(n3x*X+n3y*Y−b3*D)/D`다. D와 numerator 부호를 모두 다룬다. 원래 binary64 입력에서 시작한 outward interval arithmetic으로 부호가 인증될 때만 fast path를 통과시키고, interval이0을 포함하면 unresolved다. 임의eps cutoff나 sign flip은 하지 않는다. interval 연산은 fast-math 없이 각 기본연산의 outward rounding을 포함한다.
3. GPU unresolved cell/plane은 목록과 상태를 반환한다. 명시적인 native host exact-resolution 단계에서 같은 원래 입력을 유리수로 해석해 topology를 확정한다. C++ arbitrary-precision rational 구현을 별도 host translation unit에 둘 수 있다. Python Fraction은 독립 oracle이며 runtime 연구 clip fallback으로 쓰지 않는다. exact-resolution 여부·cell/plane 수·시간을 report/manifest에 남긴다. 전체 API는 모든 cell이 확정된 뒤에만 commit한다.
4. 정확한 길이0 edge/중복 vertex만 symbolic equality로 제거한다. 양의 미세면은 유지한다. export에는 support IDs, local endpoint 근사, local midpoint, 별도 directed edge vector, length, normal, midpoint wn을 둔다. midpoint와 edge vector는 exact 합/차를 먼저 만든 뒤 각각 binary64로 반올림한다. 이미 반올림된 endpoint의 차로 미세면 길이·방향을 다시 만들지 않는다. 길이는 exact squared length에서 충분한 정밀도로 계산하고 필요한 경우 결과 구간을 저장한다. 정확 판정된 polygon의 area/centroid/second moment도 일관된 exact 또는 검증된 고정밀 적분에서 export할 수 있다.
5. double endpoint가 exact 양의 면을 표현할 수 없는 경우와 정상적인 endpoint 표현을 구분한다. downstream density/face quadrature가 midpoint+edge-vector 표현을 사용하기 전에는 그런 cell을 일반 polygon provider에 연결하지 않는다. density API를 우회해 미세면을 없앤 polygon을 전달하지 않는다. 현재3–5 사례가 이 출력표현 관문을 통과하는지 수치로 확인한다.

### 후보 완전성과 공유면 확정

Qhull 목록은 initial candidate로만 사용한다. 저장 stage 두 개의 각 출력 vertex에 대해 전체 registry plane을 filtered predicate로 검사하고, exact 양의 위반이면 해당 ID를 추가해 재구성한다. unresolved 부호는 native exact resolver로 보낸다. candidate completion이 끝나 모든 vertex가 정확히 feasible이면 bounded convex polygon의 모든 점이 전체 halfspace를 만족한다. 이후 missing face를 길이 cutoff로 허용하는 규칙은 사용하지 않는다.

양쪽 cell의 shared segment는 global exact endpoint의 역순 equality 또는 동등한 exact support/intersection 관계로 확인한다. 유한·양의 면의 상대 cell 누락은 크기와 무관하게 실패다. exact zero contact는 finite-face ledger에 넣지 않는다. canonical owner가 midpoint/edge vector/normal/wn을 한 번 export해 두 사용자가 반대 방향으로 사용한다. 양쪽 반올림 polygon을 평균하여 새 면을 만들지 않는다.

### 다음 구현·검증 묶음의 범위

먼저 동일 native v2 구현 안에서 cell3/5 witness와 ID2 completion을 재현하고, 그 구현을 저장 stage0/midpoint1590셀에 적용하는 한 geometry 검증 묶음으로 진행한다. 기존 status6을 통과로 바꾸거나 v1 tolerance를 키우는 패치는 하지 않는다. independent Python Fraction 결과는 두 witness cell과 새로 발견되는 unresolved 유형의 근거로 재사용하고 모든 plane의 비싼 Fraction 반복을 무조건 추가하지 않는다.

보고서는 filtered 판정 수, exact-resolution 수, candidate 추가 IDs, 정확 shared topology, output representation 실패, 면적/모멘트 차이, negative common-W gauge 및 transactional 실패를 구분한다. 각 snapshot의 native/참조 수치차이는 보고하되 Qhull 참조의 누락된 미세 topology를 native 결과에 강제로 맞추지 않는다. 정확 topology가 완료되기 전 native provider/RHS/weight solve/MPI 연결은 보류한다. 이번에도 새 진화나 기존 실행 변경은 없다.

gauge fixture는 encoded binary64의 상대 W가 정확히 유지됐는지 먼저 유리수로 확인한다. 일반적인 W−.125는 덧셈 반올림으로 실제 Wi−Wj를 바꿀 수 있으므로 그 경우 exact topology 불변을 요구하지 않는다. 정확한 공통상수 이동과 반올림된 입력 perturbation을 별도로 분류한다. W=0 전체에 같은 표현가능 상수를 적용하는 fixture 등 상대 차이가 정확히 유지되는 사례는 exact gauge 검사로 사용할 수 있다.

## Native exact witness 통과와 전체 snapshot/filter 관문

`canonical_native_exact_geometry_witness_2026-10-18.json`을 확인했다. native host resolver는3/5의 Python Fraction endpoint·length²·적분 RNE 기준을 재현하고 ID2 completion을 찾아냈다. CPU1/OpenMP4, 정확한 negative common-W gauge와 registry 순서 반전에서 bitwise 일치한다. coincident 실패는 caller output을 보존한다. library SHA는 `c05ba86e6e3dbf383b8e5d1d56160327784180b3788589d7406e608aecd306d7`다. 이 두 cell의 exact-resolution 관문을 통과로 판정한다.

endpoint collapse/reversal은0이지만 반올림 endpoint 차분의 최대 상대 edge-vector 오차는 .007155/.004715다. 정확 edge 차분에서 각각 반올림한 별도 vector는 보존된다. 따라서 no-collapse만으로 기존 endpoint-only quadrature를 승인하지 않는다. 작은 absolute face measure와 topology/표현 정확도를 구분한다.

### 다음 단일 묶음: 전체1590 두 snapshot의 host completion과 GPU sign filter

1. 완료 midpoint의 stage0/half 두 저장 registry/candidate 입력을 그대로 사용한다. 각 snapshot 전체1590셀을 CPU1/OpenMP4 native exact resolver로 처리하고 exact all-registry completion, 공유면 global endpoint equality, 양의 area와 양의 finite-face measure를 확인한다. missing candidate ID 목록 및 witness, box/hidden/coincident/capacity/representation failure는 구분해 저장한다. failure를 작은 면 제거로 해결하지 않는다. native search가 Qhull 없이 완료됐다는 의미는 아니며 initial candidate는 여전히 Qhull이다.
2. 원래 x/W와 support-plane IDs로 정의되는 실제 clipping predicate 목록을 기록해 동일 `moment_geometry_interval.h` filter를 CPU/GPU에서 평가한다. 전체 registry feasibility filter는 exact vertex의 outward binary64 enclosure와 동일 plane 입력을 사용한다. 같은 scalar 입력을 두 backend에 전달하고 정확 부호 기준을 유지한다. nonzero filter 부호가 정확 기준과 다르면1건이라도 실패다. CPU/GPU가 둘 다 resolved이면 같은 부호여야 하며 어느 한쪽만 unresolved인 것은 허용한다. GPU가 resolved이고 host filter가 unresolved인 항은 native exact arithmetic으로 반드시 확인한다. 양쪽 unresolved는 exact resolver가 결정한다. 전체 부호 query의 count/hash를 기록하고 불일치·unresolved witness를 보존한다.
3. 이 단계는 GPU에서 interval predicate를 실행하는 관문이다. host에서 이미 확정한 topology를 GPU가 생성한 것처럼 부르지 않는다. GPU가 unresolved queue를 생성하고 host exact resolver가 확정하는 역할과 처리 수/시간을 구분한다. GPU topology builder 연결은 이 결과 뒤에 정한다.
4. 원래3–5/ID2 witness, exact equality, nearly parallel support, determinant zero, subnormal/overflow 구간이 whole-line 또는 unresolved로 가는 사례를 같은 filter harness에 포함한다. fast-math를 끄고 binary64 outward arithmetic을 지키며, GPU/CPU resolved 비율을 높이기 위한 임의 tolerance는 넣지 않는다. FMA/연산 순서 차이가 구간 폭에 영향을 줄 수 있으나 잘못된 resolved sign은 허용하지 않는다.
5. export 보고서는 실제양의 face 수, 짧은면 길이 분포, 후보 completion 수, endpoint 반올림차이/collapse/reversal, edge-vector 상대오차, 정확 midpoint/length bracket/적분의 재현성, 외부면 geometry_wn=wall0를 분리한다. shared owner API가 읽을 canonical 면 데이터도 저장한다. legacy 연구 polygon과의 area/centroid/secondmoment 차이는 참고로 보고하고 exact topology를 Qhull reference에 맞춰 삭제하지 않는다.
6. 전체 resolver 성공과 filter sign 인증이 끝나도 density/velocity/face quadrature provider 연결은 별도 representation 관문 전까지 보류한다. 다음 interface는 midpoint+별도 edge vector로 quadrature point와 measure를 만들고, 내부 shared face를 owner가 한 번 생성해야 한다. 그 연결은 전체 snapshot에서 실제 발생한 표현 문제를 보고 결정한다.

이 묶음에는 새로운 RHS, massroot 또는 진화가 없다. 기존 contact coarse t.001 완료 결과는 보존하고 fine의 동일 종료시각을 기다린다. corotation B도 원래 계획대로 유지한다.

### GPU filter의 all-unresolved 결과는 통과가 아님

전체 bundle stage0의 초기 binary에서 GPU support42829/42829 및 feasibility351719424/351719424가 전부 unresolved였다. CPU는 각각1/28623만 unresolved다. wrong_nonzero_signs=0은 이 경우 공허한 부호 안전성 결과이며 GPU filter의 유효한 이식 통과로 판정하지 않는다. host exact completion과 shared topology 결과는 별도로 보존한다.

nvcc의 `nextafter(double,INFINITY(float))` host constexpr overload 경고와 실제 연산의 관계는 중간구간과 binary64 nextafter 결과로 진단한다. 명시적 double 방향 인수의 `::nextafter(...,±HUGE_VAL)` 판본은 새 binary hash로 관리하며 기존 저장 predicate와 Cartesian 입력만 replay한다. 전체 exact completion이나 RHS를 반복하지 않는다. 수정 뒤 nonzero 부호의 정확성뿐 아니라 명확히 부호가 떨어진 query에서 GPU가 실제 resolved 판정을 내리는지도 필수 관문으로 둔다.

## 사용자 목표 수정에 따른 우선순위 재설정 — 이후 지시가 앞선 계획보다 우선

사용자가 목표를 **long-term Kepler 진화의 안정성 확보**로 수정했다. Voronoi의 Laguerre 확장과 시간에 따라 변하는 weight는 필요하면 허용된다. CUDA+OpenMP+MPI 정합성은 더 이상 독립적인 필수 완료 조건이 아니다. 따라서 위 문서의 GPU topology builder→native massroot→MPI 확장을 자동으로 진행하라는 순서는 여기서 해제한다. 이미 확보한 backend/기하 검증 결과는 재사용하되 성공 기준은 장시간 물리 진화다.

현재 sign replay와 compact-input hash 확인은 실제 발견한 오류의 해결 근거를 완결하는 작업으로 마무리한다. 그 뒤 새로운 GPU/MPI 구현은 시작하지 않는다. B 및 진행 중 contact fine은 유지하고 완료한 A/coarse 자료는 보존한다. exact 미세면을 임의 삭제하거나 기존 연구 geometry의 오류를 없던 것으로 취급하지 않는다.

### 다음 우선순위1: 후보의 동일시각 물리 평가와 기하 영향 확인

contact coarse/fine의 공통 t.001 완료 결과를 먼저 평가한다. pressure, 열에너지, residual covariance 양성, m/frozen 보존, radial/azimuthal velocity, internal angular momentum, 영역별 오차 및 두 dt 차이를 같은 시각에서 보고한다. inner<.4, contact 인근, outer active 경계와 전체 active를 분리한다. nonflat 수만으로 contact 증가나 실패를 판정하지 않고 contrast/theta/tracejump 분포와 dt 의존성을 함께 확인한다. 두 dt로 시간차수를 확정하지 않는다.

새 moment 모델의 moving fixed-mass cell 중심은 material particle이 아니므로 초기 centroid의 원궤도를 정확해로 강제하지 않는다. 정확 Kepler Eulerian field를 **현재 polygon에 적분한 질량·운동량·위치/속도 모멘트**를 물리 기준으로 쓴다. 정확 밀도 integral과 고정 numerical m의 차이 자체도 진단이다. 압력의 일정값, 평균속도/내부 ell의 정확 cell projection, 지구중심이 아닌 문제 원점에 대한 총 angular momentum을 구분한다. auxiliary-memory A/B의 particle phase 분석을 이 모델에 그대로 이식하지 않는다.

기하 작업의 다음 실제 산출물은 GPU topology 이식이 아니라 **기존 저장 stage0/half의 old/exact adjacency가 local density prior에 미치는 영향**이다. 원래/current raw를 고정한 채2-hop ID 집합·prior low/high/normal/offset/flat 판정과 face trace를 비교한다. exact provider의 midpoint+별도 edge vector를 소비할 수 있는 표현 검사를 이 평가에 필요한 만큼만 진행한다. geometry 모멘트 차이가 작아도 adjacency 기반 prior가 불연속적으로 바뀔 수 있으므로 그 영향을 숨기지 않는다. prior/reconstruction의 차이가 유의미하면 해당 저장 상태의 완결 RHS 한 묶음으로 귀속하고 후보 closure를 수정한 뒤 장시간 분기를 정한다. 차이가 미미하면 geometry precision 추가 이식을 목표 자체로 확대하지 않는다.

### 다음 우선순위2: 실제 장시간 후보의 물리 계약과 비용

현재1590 patch의 frozen reservoir는 장시간 전체 Kepler 검증을 대신할 수 없다. 이웃 운동으로 frozen polygon이 변할 수 있고 mass/thermal fit·경계 수송이 물리해와 달라질 수 있다. 최종 장시간 후보는 필요한 유체영역 전체를 evolve시키거나 물리적으로 정의하고 검증한 경계 조건을 가져야 한다. reservoir와 active 경계에 생긴 오차를 bulk 물리오차에서 분리하며, 경계를 바깥으로 옮겼을 때 결론이 유지되는 근거를 최종 검증에 포함한다.

Psi transport closure는 여전히 물리적 미폐쇄다. hydro ledger 통과를 전체 Ehyd+potential energy 보존으로 부르지 않는다. 새 장시간 실행 전에 실제 face mass flux와 potential work를 사용하는 일관된 에너지 budget을 정의하고, 저장 stage에서 continuity/potential work 및 frozen/exterior 교환을 확인한다. 진단적 potential 재적분과 수송 상태의 차이가 축적되는 경우 열을 보정해 감추지 않는다. 원인과 closure 수정 방향을 평가한다.

dt1e-4에서 t10은10만 step이며 현재 Python provider의 stage당 수분 비용을 그대로 적용하면 현실적인 검증 계획이 되지 않는다. stage별 측정시간과 예상 full-domain 크기로 t1/t10 비용을 먼저 산정한다. 필요한 성능 작업은 이 비용을 낮추는 지배 항만 대상으로 한다. 이미 검증된 compiled stage/update는 유지하고, density fit/reconstruction/geometry/massroot 중 측정된 병목에만 native 또는 병렬 구현을 투입한다. MPI 또는 GPU 이식 자체를 다음 완료 목표로 삼지 않는다. dt를 물리 근거 없이 키워 비용 문제를 숨기지 않는다.

### 다음 우선순위3: 한 장시간 후보를 연속 검증

위 물리 평가와 비용 관문 뒤 한 후보 설정을 고정해 새 immutable 실행을 계획한다. 새 미세 분기를 반복하는 대신 같은 실행의 t.01/.1/1 및 원래 장시간 t6–10에서 적절한 진단 시점을 둔다. 각 시점은 성공 기준이 아니라 진행 중 물리오차/성장률·실패 원인을 확인하는 checkpoint다. 정확 horizon·dt·물리영역·출력 비용과 pressure/궤도/열/각운동량 허용 기준은 실행 전에 문서화하고 결과를 본 뒤 유리하게 바꾸지 않는다. 시간·공간 정련 또는 기준해 비교는 실제 확인하려는 오차를 분리하는 최소 묶음으로 정한다.

최종 성공은 단순 exit0/finite/양성이 아니다. 원하는 장시간에서 영역별 압력과 궤도/속도·모멘트 오차가 허용 기준 안에 있고, 고정질량·각운동량·에너지 교환이 정의된 물리 budget과 일치하며, timestep/경계/기하 선택에 대한 결론의 신뢰성이 확보돼야 한다. 기존 A의 t10 정확도 실패는 기준 자료로 유지한다. 현재는 새 장시간 분기를 승인한 것이 아니며, 진행 중 결과와 물리 계약을 평가한 뒤 계획·평가 담당 Astra가 다음 실행 묶음을 지정한다.

## 수정 GPU interval filter 저장 replay 최종 평가 및 adjacency 귀속 범위

계획·평가 담당 Astra가 `canonical_geometry_v2_filter_replay_v2_2026-10-18.json`을 확인했다. 수정 라이브러리 SHA는 `9c7042941014d85796949d5dcef5f6706715cbcb15fdb11a6ca4984bb8df6d9b`다. stage0/midpoint clipping 42829/42724개에서 CPU/GPU sign 배열 SHA가 같고, 각각 42828/42724개를 실제 판정했으며 exact oracle 대비 잘못된 nonzero sign은 0이다. 각 snapshot의 351719424개 전체 registry 질의에서도 sign SHA가 같고, 351690801개 판정·28623개 unresolved·9541개 self exact zero를 확인했다. 실제 로드 입력의 hash 검증과 전체 tile coverage가 기록됐다. 따라서 앞선 all-unresolved 상태를 해소했고, 저장된 두 snapshot에서 GPU interval filter의 유효성을 통과로 판정한다.

clean probe는 double nextafter bit, interval add/mul, plane 및 D/N을 정상 재현했다. nearly-parallel/singular/subnormal/overflow fixture는 판정 가능한 경우의 정확 부호와 판정 불가능한 경우의 명시적 unresolved 계약을 통과했다. subnormal/overflow에서 0 반환은 exact zero의 증명이 아니라 fallback 요청이다. 앞선 failed probe에 legacy 호출이 남아 수정 필터까지 잘못 진단했던 문제와 실제 필터의 explicit double nextafter 수정은 별개로 기록한다. 이 결과는 저장 predicate replay의 검증이며 GPU topology builder, 실제 provider 연결 또는 장시간 정확도 검증으로 확대하지 않는다. 이 묶음을 종료하고 사용자 수정 목표에 따른 물리 우선순위를 유지한다.

다음 묶음은 저장 stage0/midpoint의 adjacency-only 영향 평가다. old polygon/raw 및 동일 working vertex error와 mass error 입력을 고정하고 stable global ID의 old/exact 무방향 graph만 교체한다. 양의 exact shared face는 길이 cutoff 없이 포함하고 reservoir 이웃 규칙은 동일하게 유지한다. patch 밖에 필요한 이웃 정보가 없으면 명시적으로 실패한다. 직접 이웃 및 정확 2-hop+self 집합 변화, low/high 구간, flat 사유, 방향·offset·cut fraction을 기록하고 phase label/normal 동치 표현을 정규화한다. 저장 fit 객체가 없으므로 baseline은 셀당 한 번 복원하되, prior가 같은 분기는 객체를 재사용하고 실제 prior 변경 셀만 추가 fit한다. 동일 평가점에서 A/B density 및 양쪽 face trace를 비교하고 old/exact geometry별 길이 가중 통계는 분리한다. 반올림으로 중복되는 미세면 점은 ID와 개수를 기록하며 exact 면적분 검증으로 부르지 않는다. mass_error=0은 저장 raw를 이번 비교의 공통 입력으로 고정한다는 뜻이다. 새 RHS·provider·진화 변경은 없으며 결과 뒤 추가 물리 연산의 필요성을 결정한다.

## Adjacency 귀속 종료와 장시간 비용 관문

`canonical_adjacency_prior_provenance_v2_2026-10-18.json`에서 stage0/midpoint 각각 1590셀의 volume quadrature points/order/offset과 dm이 원래 저장값과 bitwise 일치했다(1431360/2169936 nodes). 면 inventory 대응이 모호하지 않고 모든 segment가 사용됐다. stage0 내부 pair [3,4]의 두 record는 끝점을 공유한 연속 segment이고 같은/역방향의 동일 encoded segment 중복은 없었다. 원래 vertex bits 자체는 저장되지 않아 직접 재현을 주장하지 않는다. 이 범위에서 adjacency-only 결과를 채택한다: stage0의 13셀 stencil/metadata 변화에도 실제 prior 변화는 0이며 midpoint에서는 graph/prior 변화가 없다. 동일 객체 재사용으로 발생한 trace 차이 0을 독립 flux 검증으로 세지 않는다. 추가 adjacency RHS 실험 없이 이 묶음을 종료한다.

`canonical_contact_stage_cost_2026-10-18.json`의 stage wall median은 coarse 310.6816초, fine 301.3211초다. 같은 patch와 dt 및 2 stages/step의 단순 t10 외삽은 약 719일/1395일이다. 동시 부하를 포함한 관측이며 full-domain 비용이나 개별 병목의 증거는 아니지만, 현재 경로 그대로 장시간 실행하는 계획은 실용적이지 않다.

다음 두 작업을 우선한다. 첫째, 진행 중 fine를 유지하고 완료 후 t=.001의 coarse/fine 물리 비교를 수행한다. raw, centroid/mean velocity, ell, R, thermal/P, generator/W 및 prior contrast/trace jump를 동일 영역으로 비교하며 두 dt만으로 수렴차수를 확정하지 않는다. 정확 정상장 기준은 현재 움직이는 polygon에 적분한 Kepler장이다. 기존 smooth 결과와 비교할 때도 시각·영역·정의 일치를 요구한다.

둘째, 단 한 번의 저장 midpoint provider 성능 재생을 준비한다. 원래 입력·초기 z seed·outer 종료조건을 재현하고 geometry, prior/fit, velocity reconstruction, outer별 face nodes·mass residual/Jacobian·linear solve, flux accumulation, gravity, diagnostic의 exclusive wall/process CPU time 및 호출 수를 계측한다. seed를 복원하지 못하면 full-stage 재현으로 부르지 않고 범위를 재평가한다. fine 종료 뒤 실행해 동시 부하를 줄이는 것을 선호한다. 새 timestep 진화나 튜닝 반복은 없으며, 측정 후 지배 연산 하나의 최적화만 선택한다. 이 1회 성능 재생은 종료된 adjacency 정확도 실험과 구분한다. 새 장시간 실행은 물리·경계·에너지·비용 관문 평가 뒤 결정한다.

## A/B 종료 및 coarse 참조·계측 실행 순서의 묶음 평가

Astra가 `canonical_memory_terminal_bundle_2026-10-18.json`, `canonical_contact_coarse_polygon_reference_2026-10-18.json`, `canonical_contact_coarse_reference_components_2026-10-18.json`, `canonical_provider_profile_preflight_2026-10-18.json`을 확인했다. A/B 모두 원래 실행을 보존해 t10까지 exit0로 완료했지만 장시간 물리 정확도 목표는 실패다. B의 inner background 상대 P RMS6.2008/위상68.2367도, shear P7.0885/위상78.1574도는 큰 오차다. A t10.00006914696와 B t10.00041349724는 서로 다른 시각이며 우열 비교를 하지 않는다. 추가 연장·재시작은 지정하지 않는다.

coarse t.001의 active 상대 P RMS7.4218e-7, 속도 RMS1.3100e-6, 상대 질량 RMS4.6319e-5는 선택 q32 참조 차이보다 충분히 크다. 관측된 오차를 적분 차수 오차로 설명할 근거가 없다. 짧은 시간의 압력 오차가 작다는 것만으로 장시간 성공을 판정하지 않는다. active-edge/boundary의 밀도·각운동량 오차와 inner-transition의 속도·열 오차를 fine 결과와 함께 평가한다. 새 파생량 최대점의 q32 미선택 여부는 기존 기록대로 유지하며 필요성은 통합 평가에서 판단한다.

B 종료로 자원이 확보되었으므로 이전의 fine 종료 후 계측 선호를 변경한다. OMP/BLAS 각1thread로 승인된 저장 provider 계측 1회를 fine와 병행할 수 있다. fine 종료 guard가 있는 원본은 보존하고 변경 실행본 SHA와 동시 실행 fine ID·thread 환경을 기록한다. 무부하 walltime이나 독점 CPU 성능이라고 부르지 않는다. fine의 설정·우선순위를 바꾸지 않는다. 입력/outer/수치 결과 재현 관문을 유지하고 단순 walltime 변동 때문에 반복 계측하지 않는다.

후속 순서는 fine 원래 종료 → 정확장 후처리 및 공동 비교 → profile/공통시각 물리/참조 신뢰성/지역별 오차의 통합 평가다. 구현 담당 root는 승인 범위의 사소한 코드·provenance 보완을 자율 처리한다. 사용자 요청에 따라 건별 Astra 호출을 줄이고, 수치 관문 실패·새 물리 가정·진행 불가 외에는 완결 묶음으로 평가를 요청한다. 새 장시간·geometry·MPI·튜닝 분기는 지정하지 않았다. 다음 지배 병목 하나의 최적화와 경계/Psi 작업 순서는 통합 결과로 결정한다.

## 공통 t=.001 물리·성능 통합 평가와 다음 실행 묶음

Astra가 `canonical_contact_common_time_physical_comparison_refined_2026-10-18.json` 및 `canonical_provider_profile_summary_2026-10-18.json`을 확인했다. 두 dt의 동일 초기자료·물리상수·현재 polygon 참조와 상태 보존/양성 관문을 통과했다. active의 coarse→fine 상대 질량 RMS4.63185e-5→4.59752e-5, 상대 P7.42182e-7→5.79862e-7, 속도1.30998e-6→1.27693e-6, ell9.57547e-12→9.52452e-12, factor R9.89481e-14→9.85851e-14다. 참조 차수 차이는 이 오차보다 훨씬 작다. 전체 질량/속도/ell/R은 dt에 약하게 반응하지만 inner-transition은 더 뚜렷이 개선된다. 압력 개선은 active-edge에서 크고 inner/contact는 거의 같다. 따라서 시간오차 또는 경계오차 하나로 전체 원인을 확정하지 않는다. 세 번째 dt나 장시간 연장은 현재 지정하지 않는다.

계측 재생은 outer3회 및 z/hydro/gravity/raw source를 bitwise 재현했다. 총 wall363.04243/process CPU361.06222초 중 density trace25.34%, velocity samples18.42%, scalar flux17.98%, exact-star10.82%, face record 나머지10.80%다. mass linear solve나 geometry가 지배 병목이라는 근거는 없다. 호출 계측 overhead와 동시 fine 부하 때문에 이 시간을 무계측 성능 기준으로 쓰지 않는다.

### 한 성능 연산자: face 상태 표본 일괄 평가

동일 point/owner/inward normal에서 density trace와 quadratic mean/sample/xi를 반환하는 batch 경로를 구현한다. 이 연산자에 해당하는 두 비용의 합은 약43.76%이며 이것만으로 장시간 비용이 해결된다고 주장하지 않는다. 변경 위치는 contact_prior_density의 batch trace, quadratic_velocity_reconstruction의 batch samples+mean, kinetic_factor_rhs의 node state 준비 연결이다. scalar 경로를 oracle로 보존하고 geometry/prior fit/면 분할점과 가중치/record 순서/exact-star/massroot/flux accumulation/gravity는 이번에 유지한다. 원래 순서로 scatter하며 cache는 해당 outer 안에서만 유지한다. cut/tangent/flat 분기, mirror BC, shared mean 정의를 보존하고 작은 부호차를 임의 tolerance로 덮지 않는다.

저장 stage0/half 표본 및 cut/tangent/flat/mirror fixtures를 묶어 대조한다. 이어 동일 실제 midpoint 입력에서 무계측 baseline1회와 batch1회를 같은 thread 환경으로 실행한다. 이는 구현 교체 검증 두 호출이며 새 진화가 아니다. 큰 구간 timer만 사용하고 호출 수/batch 크기/peak memory/총 wall와 process CPU를 기록한다. 최종 z·mass residual·outer 종료·hydro·gravity·raw source에 기존 성분별 수치 관문을 적용한다. 분기 불일치나 관문 초과는 결과를 보존하고 보고한다. 추가 자동 튜닝 재생은 지정하지 않는다.

### 물리 우선 작업: 저장 midpoint step의 Psi와 경계 장부

기존 stage0/half face inputs에서 outward J=scale*(sum wl+sum wr)를 구성한다. 같은 point의 potential 및 중력 E source gE를 사용해 dPsi_transport=−sum(Phi*J)−gE를 계산한다(a=−grad Phi). 각 shared face는 한 번 계산해 반대 부호로 배분하고 active-active 상쇄, active-frozen, patch exterior를 분리한다. frozen evolution source0과 실제 flux/geometry 교환을 구분한다. hydro E와 합칠 때 gE 상쇄를 확인한다. 실제 midpoint step의 full budget은 half rate에 h를 곱한다.

필요한 endpoint의 수치 density를 한 번 재구성해 Psi0/Psi_end를 계산한다. 동일 density parts의 차수24/48 비교로 적분 불확실성을 분리한다. 재구성 Psi 변화와 수송 장부의 차이 및 Δ(Ehyd+Psi_reconstructed)−h*(hydro+Phi boundary exchange)를 저장한다. 해석적 Kepler rho를 수치 Psi 대신 쓰거나 열을 보정하지 않는다. 이는 진화된 Psi 보존 증명이 아니라 재구성 일관성 결함 진단이다. 새 RHS/weight solve/진화는 없다.

같은 준비를 재사용해 t=.001 두 endpoint의 prior contrast/theta/양쪽 density jump를 active-boundary/contact/inner별로 수집한다. 저장 checkpoint에 없어 필요한 fit 재구성은 각 상태당 한 번 허용한다. 새 해석적 reservoir 주입이나 경계 확장 진화는 아직 지정하지 않는다. root는 일반 구현·검증을 자율 처리하고 수치 관문 실패/새 물리 가정 이외에는 건별 호출을 하지 않는다. 다음 평가는 성능 실측 개선, 한 step Psi 장부, endpoint prior/trace 및 지역별 물리오차를 모은 완결 묶음으로 요청한다.

## Face-state batch 채택 및 Psi 결함 평가 — 다음 완결 묶음

Astra가 `canonical_face_state_provider_comparison_2026-10-18.json`, `canonical_saved_midpoint_Phi_flux_2026-10-18.json`, `canonical_numerical_Psi_endpoint_bundle_2026-10-18.json`을 직접 확인했다. face-state batch는 실제 midpoint의 outer3회와 최종 z/hydro/gravity/raw를 bitwise 재현했다. wall325.38449→168.75066초, CPU323.41319→167.73817초, RSS1826.87→1533.67MiB를 관측했으므로 이 경로를 채택한다. 한 번씩의 관측이며 이전 per-call profile 비율에서 계산한1.778은 엄밀한 Amdahl 상계가 아니다.

수송 Psi를 정의한 hydro+Psi−boundary 잔차 약−1.35e-18은 대수·구현 장부 검증이다. 수치 density로 재구성한 Psi와 수송의 차이+6.46016e-11은 q24/48 차이2.77332e-16보다 크지만, 한 step으로 시간오차·재구성 폐쇄·branch 전환을 분리하지 못한다. active-boundary의 직접 기여가 작아도 경계의 간접 영향을 배제하지 않는다. frozen raw가 고정돼도 움직이는 polygon에서 Psi와 trace가 변하므로 현재 frozen reservoir를 장시간 정상 경계로 채택하지 않는다. 별도 Psi 장부만 추가해 재구성 일관성 문제가 해결됐다고 부르거나 E/열을 보정하지 않는다.

nonflat 수 자체는 실패 기준이 아니다. t.001 inner contrast 중앙값은 약1+1.9e-6, theta RMS약3.3e-5다. 반면 frozen-boundary theta RMS약.164, 상대 face jump RMS약.0283은 coarse/fine에서 거의 같다. contact jump RMS약.276도 dt 반응이 작지만 실제 밀도 접촉에는 물리적 점프가 있으므로 이를 줄이기 위한 smoothing/점성/cutoff를 도입하지 않는다. Re∞, 무조향, 고정mass 계약을 유지한다.

### 다음 성능 대상: moment-flux accumulation

기존 `Docs/polybnd/vector_moment_flux.py`의 accumulate를 현재 quadratic/contact 계약에 맞춰 검토·연결한다. 반환 계약은 dr18, exterior ledger, F/H, K, 순서 보존 incoming Q, area rates 전체다. frozen mask와 pressure-volume source의 중복 적용을 금지한다. j<0의 right state를 음수 index로 읽고 폐기하는 관행에 의존하지 않고 valid mask로 처리한다. incoming은 기존 face/sample 순서를 유지한다. 기존 kinetic/CSR stage의 저장 oracle로 stage0/half 반환항·장부·K/incoming을 검증한다. face-state batch를 유지한 scalar-accumulation baseline1회와 새 accumulation1회를 동일 midpoint 입력/seed/thread 환경에서 비교한다. 기존 성분별 roundoff 관문, 총 wall/process CPU 및 RSS를 사용한다. 새 큰 absolute tolerance, exact-star 변경, GPU/MPI 확장 또는 진화는 지정하지 않는다.

### 다음 물리 대상: Psi 재구성의 반이산 일관성

저장 자료에서 P2=a+b·x+1/2 xᵀHx(H대칭)에 대해 Psi2=a m+b·q+1/2 H:Q 및 그 raw-rate 미분을 계산한다. face −sum(P2 J)와 volume b·p+H:D를 비교한다. constant/linear/quadratic 여섯 기저 및 local Taylor P2를 검사한다. local Taylor 계수는 midpoint cell의 고정 중심에서 실제 Phi·중력계수의 해석적 도함수로 구성하고 initial/half/end에 동일 계수를 쓴다. m은 실제 update대로 고정하고 dm roundoff residual은 분리한다. 각 셀의 Phi=P2+remainder에서 endpoint 재구성·flux·volume work를 직접 분해한다. 셀마다 다른 P2의 face term만으로 전역 shared cancellation을 요구하지 않으며 실제 Phi 합에서만 상쇄를 요구한다. 중력장을 교체하는 실행은 아니다.

추가로 저장 midpoint의 x,W,q,Q를 저장 actual stage rate 방향으로 ±epsilon, ±epsilon/2 이동한 고정4상태를 재구성한다(epsilon=기존 h/8). m과 frozen 상태는 원래 규칙대로 고정한다. factor/velocity/Pi 업데이트, 새 RHS/weight solve/시간진화는 없다. geometry→prior→m/q/Q density fit 및 Psi 적분만 수행한다. 중심차분 두 값과 한쪽차분, fit 잔차, graph/prior branch 전환을 저장한다. 동일 branch에서 차분이 적분/재구성 오차보다 안정적으로 일치하면 C=directional dPsi_recon−(−sum(Phi J)−gE)를 반이산 closure defect로 읽는다. branch가 바뀌거나 scale 간 불안정하면 미분 미확정으로 보류하고 epsilon sweep을 추가하지 않는다. C와 실제 step defect/h를 비교하되 나머지를 자동으로 순수 시간오차라 부르지 않는다. remainder 직접 적분으로 큰 다항식 항의 상쇄를 줄이고 q24/48, fit 제약 잔차 및 합산오차를 함께 보고한다.

이 묶음의 provider 실행은 성능 비교2호출, 물리 재구성은 고정4상태 범위다. 새 장시간·경계 확장 진화는 없다. root는 일반 구현/검증을 자율 처리하고 관문 실패 또는 새로운 물리 가정 외에는 중간 호출하지 않는다. 성능과 polynomial/방향미분 진단을 묶어 평가한 뒤 실제 density/Psi closure 변경 및 장시간 경계 계약을 결정한다.

## 2026-10-19: moment-flux 채택과 potential 제약 density 후보

Astra가 `canonical_potential_consistency_bundle_2026-10-18.json`의 성능·다항식·고정4상태 결과를 직접 확인했다. moment-flux batch는 actual provider outer3회와 z/hydro/gravity/raw를 bitwise 재현했으므로 채택한다. 이번 측정에서 wall167.77452→109.03531초, CPU166.68135→108.41303초, RSS1530.78→1721.35MiB다. 성능 개선과 메모리 증가를 함께 기록한다.

P2 여섯 기저 및 local Taylor 항등식은 통과했고 actual full-step P2 결함 signed1.41026e-18, remainder 결함6.46016e-11이다. 따라서 이번 저장 step에서 보존된 m/q/Q의 다항식 부분보다 실제 potential의 비다항식 remainder 재구성과 수송의 불일치가 확인됐다. candidate directional C≈6.45406e-7와 step defect/h≈6.46016e-7의 근접은 지지 증거다. 하지만 prior branch 전환 및 두 scale 차이4.26132e-9가 수치 오차 추정보다 크므로 매끄러운 전역 반이산 미분의 증명으로 부르지 않는다. 추가 epsilon sweep이나 classification threshold 변경은 하지 않는다.

### 물리 후보: m/q/Q/Psi를 함께 맞추는 양의 상대 entropy reconstruction

새 독립 ContactPotentialDensity API를 준비한다. 기존 raw18은 유지하고 Psi는 별도 sidecar와 명시적 mode로 취급한다. 첫 실행은 저장 accepted polygon/raw m/q/Q/기존 local prior를 고정한 active1194의 offline 실현 가능성 평가다. target은 저장 수치 density의 Psi_initial+h*dotPsi_transport_half이며 해석적 Kepler density를 주입하지 않는다. 초기 self-target control은 원래 reconstructed Psi를 쓴다. E/Pi/thermal/m/q/Q를 수정하지 않는다.

기존 다섯 quadratic feature F에 potential 독립 feature를 추가한다. 동일 positive quadrature의 baseline density 확률을 가중치로 {1,F}에 Phi를 weighted QR 투영해 P2proj를 만들고 chi=Phi−P2proj, sigma²=<chi²>를 계산한다. 추가 feature는 chi/sigma, target은 (Psi_target/m−<P2proj>target)/sigma다. P2 target은 m/q/Q에서 계산하고 이 제약이 원래 Psi 제약과 동등함을 검증한다. 차감과 remainder는 longdouble/compensated 계산 및 직접 적분으로 처리한다. sigma가 평가·투영 오차와 구별되지 않으면 unresolved-rank로 보고하며 임의 floor/regularization으로 통과시키지 않는다. 정확 degree≤2 potential fixture는 종속 제약임을 확인하고 target consistency까지 요구한다. 상수 potential gauge 변화는 Psi_target에도 C*m을 반영하며 물리 density가 바뀌지 않아야 한다.

fit용 q24 및 독립 q48을 사용하고 기존 prior cut과 potential spline 분할을 유지한다. accepted active1194에 constrained fit을 한 번씩 수행해 m/q/Q/Psi, positivity, rank/condition, dual residual, entropy 변화와 기존 density/face trace 대비 변화를 기록한다. constant/linear/quadratic 종속 potential, 비다항식 potential, 부적합 target의 명시적 실패 fixture를 묶는다. Newton 실패를 물리적 infeasibility로 단정하지 말고 quadrature moment hull/dual witness가 있는 경우와 unresolved를 구분한다. 실패 셀을 원래 closure로 자동 대체하지 않는다. 이 단계는 backward target의 재구성 가능성이며 self-consistent 진화 개선 증명이 아니다. 후보 full RHS·진화 연결은 아직 지정하지 않는다.

### 성능 후보: face-record 생성과 exact-star batch

승인된 두 batch 경로를 유지하고 face-record 생성 연산자를 배열 중심으로 구성한다. node geometry/pressure/BC를 배열로 준비하고 exact Euler star를 CPU native batch API로 평가해 원래 node 순서에 기록한다. scalar oracle를 보존한다. 기존 vacuum reject, pressure scale, bracket 확대, normalized root tolerance8eps 및 maxiter 계약을 유지한다. HLLC/점성 fallback, 다른 massroot/fit/geometry 변경은 없다.

저장 stage0/half의 실제 Riemann 입력과 shock/rarefaction/equal-state/near-vacuum/invalid fixtures를 묶어 검증한다. bracket/residual에서 pstar·speed 오차를 전파하고 기존 provider 수치 관문을 적용한다. actual midpoint 같은 input/seed로 baseline1회/candidate1회만 실행하며 큰 구간 timer로 record 준비/exact-star/누적/capture를 분리한다. capture 포함·제외 시간은 같은 호출에서 보고한다. GPU/MPI 및 새 진화는 지정하지 않는다.

### 장시간 경계·비용·dt 계약 준비

장시간 후보는 검사 영역과 동적으로 진화하는 buffer를 포함해야 하며 현재 r=.60의 frozen raw reservoir를 채택하지 않는다. 원래 전체 IC registry의 외곽/주기성/기하 BC와 rho/P/v를 소스·저장자료로 조사해 정상 Kepler장과 외곽 조건의 양립 여부를 기록한다. 전체 유체를 evolve할 때 셀/node/메모리·시간 비용을 현 측정으로 범위 추정한다. 현재109초/stage와 dt1e-4를 고정하면 patch만 t10 약252일이므로 아직 장시간 실행 가능 상태가 아니다.

다만 dt1e-4가 장시간에 반드시 필요한 간격이라고 입증된 것도 아니다. 저장 stage0/half의 outgoing kinetic mass rate/m, 상대 face 속도와 sound speed, Omega, strain/factor generator norm에 h를 곱한 무차원 수 및 최소 cell geometric scale을 집계한다. 이는 timestep 설계 자료이며 positivity/realizability 보증이나 dt 증가 승인이 아니다. 실제 outer BC는 정상장의 유입·유출과 맞는 명시적 물리조건을 가져야 한다. 해석적 profile로 내부셀/생성자를 재설정하는 조향은 허용하지 않는다. closure 실현 가능성·성능·BC 계약을 다음 완결 묶음으로 평가해 실제 경계 선택과 dt/장시간 horizon을 결정한다. 이번에는 저장자료 분석 외 새 경계 확대/진화 계산을 하지 않는다.

root는 일반 구현·검증을 자율 처리한다. 이번 provider 호출은 성능 비교2회이며 closure 후보는 저장 상태의 fit/적분만 수행한다. 수치 관문 실패나 새 물리 가정 이외에는 건별 Astra 호출 없이 완료 묶음으로 평가를 요청한다.

## 2026-10-19: Ψ-coupled midpoint 및 원형 slip-wall 장시간 경계 후보

Astra가 `canonical_closure_cost_boundary_bundle_2026-10-19.json`의 closure fixtures/accepted fit, native saved-node/provider 및 boundary/timestep 자료를 직접 확인했다. native face-record는 실제500040 node 및 actual provider의 bitwise 재현을 통과했으므로 채택한다. wall108.5778→68.5439초와 RSS1723.77→1949.20MiB를 함께 기록한다. face export만 제외한51.6063초는 다른 audit를 포함한 이전 closure의 관측이며 새 solver production 비용으로 부르지 않는다.

ContactPotentialDensity는 active1194의 offline positive fit과 독립q48 검증을 통과했다. maxPsi오차1.4093e-17, absolute합6.2334e-16은 이전 reconstruction−transport absolute합8.7407e-11보다 작다. 이는 추가 제약의 실현 가능성을 입증하며 self-consistent 진화 정확도는 아직 입증하지 않는다. self-control scaled remainder3.65e-9는 기존5feature fit의 유한 raw 제약 잔차를 작은 sigma로 정규화한 결과다. 초기 coupled 상태에서도 raw m/q/Q와 기존 numericalPsi0를 고정한 새 constrained fit을 실제 수행한다.

### 주요 실행: 실제 Ψ-coupled midpoint 한 step

기존 .60 patch, 같은 t0 raw18/factors/generator/W/m/IDs, h1e-4를 사용한다. Psi0는 기존 수치 density의 저장q48 적분이다. active의 새 density를 volume moments, velocity basis/mean, generator, face trace, massroot, gravity 전체에 일관되게 사용한다. stage0에서 동일 face J와 gravity E source로 dotPsi=−sum(Phi J)−gE를 계산한다. halfPsi=Psi0+h/2 dotPsi0와 실제 half raw/geometry/L로 constrained fit 및 full provider를 다시 계산한다. finalPsi=Psi0+h dotPsi_half이며 raw/L은 기존 midpoint/factor map이다. 이전 offline accepted target이나 이전 half flux를 새 coupled step에 재사용하지 않는다. Psi를 재적분값으로 덮어쓰거나 E/thermal을 보정하지 않는다.

이번 patch의 frozen 조건은 이전 short reference와 동일하게 유지하되 psi_constrained_mask=active를 명시한다. frozen sidecar 고정과 실제 재구성Psi 변화는 분리하며 이 patch를 장시간 경계 후보로 취급하지 않는다. provider 호출은 stage0/actualhalf 총2회다. final density 재구성 및 독립q48·저장 endpoint 물리 후처리는 허용한다. 6feature density가 모든 batch/face API에서 실제 평가되는지 확인하고 기존5feature shortcut을 우회 재사용하지 않는다.

기존 m/frozen 보존, Q/Pi 대칭, R/thermal 양성, fit/rank/dual 상태, outer massroot 및 stage/accepted active boundary 장부 관문을 적용한다. q48 m/q/Q/Psi와 scaled remainder를 별도로 검증하고 반올림·적분·제약 오차를 energy 판정에 반영한다. initial/half/final의 sidecar/target/mode/SHA/rates/원시 장부/fit metadata를 저장한다. 같은t1e-4 기존 unconstrained 저장 endpoint와 비교하되 각 현재 polygon의 exact-reference P/v/m/ell/R를 함께 보고한다. Psi 일치가 좋아져도 궤도·열 정확도는 악화될 수 있으므로 closure 전체 채택은 결과 뒤 결정한다. 실패 시 원자료 보존·명시적 중단이며 자동 fallback은 없다.

### 장시간 BC 주 후보: R_B=2.5 정지 원형 slip wall

현재 물리상수 rin=.5, rout=2, floor=.01, 일정P0의 정상 Kepler장은 원점 중심 원형벽에 접한다. 주 후보를 R_B=2.5의 정지 slip wall, 내부 모든 fluid 및 rout2→2.5 buffer를 동적으로 evolve하는 문제로 정한다. 벽 mass/E/Phi flux는0, radial wall traction의 원점 torque는0이며 내부 analytic reset/생성자 조향을 요구하지 않는다. 이 문제는 원래 A/B의 외곽BC와 다른 명시적 benchmark다. 이후 buffer/반경 민감도는 확인해야 한다.

원래 rotating floor는 opposite square edge의 tangential flux가 다르고 square wall normal velocity가0이 아니므로 naive periodic/정지 square wall을 사용하지 않는다. 고정영역의 open square boundary는 일반적으로 총 boundary mass flux가0이 아니어서 모든 cell fixedmass massroot와 양립하지 않을 수 있으므로 이번 주 후보로 선택하지 않는다.

이번 BC 작업은 reference geometry/flux 계약과 정적 fixture까지만이다. power cell∩disk를 straight shared face와 실제 circular arc wall로 표현한다. arc를 chord로 몰래 대체하지 않는다. arc node에서 x/R_B가 outward normal, wn=0, ds=R_B dtheta이며 mirror는 normal 성분만 반전한다. closed-system massroot의 constant-W gauge nullspace를 명시적으로 제거하고 기존 frozen anchor 계약과 분리한다. generator는 mean velocity대로 움직이며 wall projection/clamp를 하지 않는다. empty/unresolved cell은 명시적 실패다.

정적 fixture는 partition areaπR², q0, Qxx=Qyy=πR⁴/4, shared cancellation, 각 cell의 divergence identities(∮n=0, ∮x⊗n=VI), zero wall mass/E/Phi flux와 origin torque를 검사한다. 향후 domain IC는 실제 clipped cell에 정상장을 적분해 m/raw/Psi를 초기화해야 하며 이전 boundary mass를 그대로 재사용하지 않는다. 이번에 전체 약2만cell domain geometry/projection/fit/RHS/진화를 시작하지 않는다.

### 비용 및 dt 순서

이번에는 추가 최적화 연산자를 시작하지 않는다. 새 constrained density가 실제 stage 비용을 바꾸므로 coupled2provider에서 density/prior/potential fit, velocity prepare, face construction, massroot, accumulation, gravity, capture 및 독립q48 비용을 큰 구간으로 측정한다. h1e-4의 작은 Courant/strain/Omega 지표는 후속dt 설계 근거지만 새로운Psi moment realizability를 보장하지 않는다. 이번 step은 h1e-4 고정이다. coupled 물리/에너지·새 비용·원형벽 정적 계약을 모아 다음 domain 구성과dt 검증 및 장시간 경로를 결정한다. 일반 구현/fixture 보완은 root가 자율 처리하고, 수치 관문 실패나 새 물리 가정 외에는 건별 호출하지 않는다.

## 2026-10-19: Ψ-coupled 연결 채택과 전체 동적 원형 영역의 첫 dt 비교

Astra가 `canonical_psi_coupled_boundary_bundle_2026-10-19.json`, 원형벽 계약 및 `run_psi_coupled_midpoint_once.py`, `circular_power_reference.py`, positive quadrature 구현을 확인했다. Ψ-coupled actual midpoint는 기존 상태 보존/양성/수송 장부와 독립 재구성 검증을 통과했으므로 다음 full-domain 물리 후보로 채택한다. transported energy residual−1.63537e-19, reconstructed residual−4.44647e-16, measured constraint+empirical quadrature 합1.06985e-15 및 보수적 summation bound5.0612e-12를 구분한다. 기존 defect6.46016e-11보다 작지만 새 trajectory이며 density closure와 volume quadrature가 함께 바뀌었다. 같은 시각 P/v/m/ell/R 변화가 작고 일부 지역은 악화됐으므로 물리 정확도 개선 또는 장시간 성공을 주장하지 않는다.

현재 circular reference의 all-constraint plane-pair loop는 작은 fixture용이다. 이 구현을 전체 약2만cell에 그대로 적용하지 않는다. 다음 완결 묶음은 확장 가능한 full geometry/새 IC/closed nonlinear root 연결과 동일 새 IC의 두 dt 비교다.

### 전체 물리영역 및 기하 API

원래 IC의 |x|<2.5 generator를 stable global ID로 선택한다(예상19945, 실제 count/hash를 확인). 외부 generator/periodic image/frozen reservoir를 사용하지 않는다. 내부 모든 fluid와 2<r<2.5 buffer가 active이며 W0=0이다. power cell∩disk를 실제 straight shared face와 circular arc로 표현한다. area/centered moments/halfplane clip/positive volume quadrature/boundary quadrature API를 polygon과 curved cell에 명시적으로 제공하며 prior cut과 potential spline 분할도 동일 curved domain을 사용한다.

candidate search는 sparse seed와 누락 halfplane completion/인증을 사용한다. Qhull 목록 자체를 완전성 증거로 삼지 않는다. 기존 native filtered/exact plane 경로 재사용을 우선한다. 전체 n²/n³ Python loop, 미세면 cutoff, arc-to-chord 치환으로 비용이나 topology 문제를 숨기지 않는다. uncertain topology/empty cell은 명시적 실패다. 내부 shared face는 global pair key로 한 번 구성한다.

line/arc 공통 face 출력은 node별 point/normal/ds weight/wn/global owner다. arc에 기존 a+g(b−a), constant normal을 적용하지 않는다. prior-cut/arc 교점을 분할하고 normal-speed/traction 적분의 부호 분할 또는 독립 차수 확인을 포함한다. x/R_B, ds=R_B dtheta, wn=0 및 mirror BC를 실제 kernel 입력으로 연결한다. 기존 native face/stage/factor kernels는 계속 재사용하며 별도 GPU/MPI 이식을 필수 목표로 만들지 않는다.

### 새로운 보수적 IC

새 disk-clipped cell에 정상 rho/v/P0의 raw18/Psi를 투영해 m을 정의한다. analytic field는 IC와 독립 reference에만 사용하고 evolve 중 내부 reset은 하지 않는다. rho interface .5/2, potential spline .25/.5 및 원호를 분할한 positive quadrature와 독립 차수를 사용한다. 기존 square boundary mass/raw에 새 area만 대입하지 않는다.

R는 centered velocity residual의 직접 positive integral/QR로 안정적으로 초기화한다. 부정/불확정 rank를 floor/noise로 감추지 않는다. 전체 area/q/Q identity, 질량 pi*[.01*2.5²+.99*(2²−.5²)], 대칭 p/q 및 독립 radial energy/Psi/L 적분을 검사하고 projection/factor consistency 오차를 기록한다. 새 IC SHA는 다른 명시적 BC benchmark의 초기조건이며 기존 A/B나 patch와 bitwise 일치를 주장하지 않는다.

### Closed nonlinear massroot

internal face residual/Jacobian을 쌍별 누적하고 wall massflux0을 확인한다. connected graph의 constant-W nullspace 하나를 제거한다. sparse reduced gauge solve 후 mean z=0 환원 또는 동등 sparse augmented solve를 사용할 수 있으나 이는 frozen cell 생성이 아니다. ridge, 비물리 residual 보정, boundary mass 유출은 없다. sum(r) 호환성과 J*1=0을 실제 절대 flux/J항의 roundoff로 검사한다. 첫 full-stage의 전체 area-rate0, mass0, wall E/Phi0 및 central gravity와 radial wall의 origin torque0도 확인한다. 정적 fixture 통과를 실제 closed root/dynamic boundary 검증으로 확대하지 않는다.

### 다음 실행 예산과 dt

먼저 full IC/geometry/closed stage0/양성 관문을 통과한다. 통과하면 동일 immutable IC와 stage0를 공유해 h=1e-3 한 step, h=5e-4 두 step의 두 분기를 공통 t=.001까지 수행한다. unique provider는 최대5회다: 공통 stage0, coarse half, fine 첫 half, fine 두 번째 step의 stage0와 half. endpoint density 및 독립 reference 후처리는 별도다. h1e-3은 이 검증의 후보이며 장시간 dt 승인이나 positivity 보장이 아니다. 실패 상태를 보존하고 자동 h축소/추가 분기는 시작하지 않는다.

각 현재 cell의 정상장 reference에 대해 P/mean velocity/radial·azimuthal velocity/mass/ell/origin L/R/thermal/Psi를 비교하고 global/wall mass·energy·torque budget을 보고한다. 기존 inner 구획과 함께 disk core/outer density contact/dynamic buffer/wall-touching을 분리한다. 실제 wall 및 density-interface 교차mask도 저장한다. 두 dt로 시간차수나 장시간 성공을 확정하지 않는다. 이번 검증은 장시간을 대체하지 않고 다음 dt와 실행 비용을 정할 전체영역 근거다.

### 자원·비용·출력

full geometry와 positive quadrature의 실제 셀/node 수로 peak memory·저장량을 먼저 산정하고 원격 가용 자원과 대조한다. 각 stage의 독립q48/full backend shadow를 production timing에 섞지 않는다. cell별 density/quadrature/velocity basis/gravity 준비는 bounded workers로 독립 처리하고 stable ID 순서로 모은다. 우선4 workers, worker별 BLAS/OMP1을 사용하며 실제 메모리·자원이 허용하면8까지 가능하다. immutable 입력 공유/chunk 전달 및 필요 없는 quadrature/trace 배열의 수명 제한으로 전셀·전stage 중복 보관을 피한다. 동일 stage의 동일 입력 prepare는 재사용할 수 있으나 stage간 cache는 명시적으로 무효화한다. 계산법과 물리 제약을 바꾸지는 않는다.

출력은 새IC, accepted3states, cell rates/ledger/fit상태, node hash/count·요약 중심이다. full node dump는 실패 witness 또는 사전에 크기를 산정한 자료만 보존한다. 기존 결과는 삭제하지 않는다. 독립 검증 비용과 실제 provider exclusive 비용을 분리하고 full step 실측 및 두 dt 정확성으로 t1/t10 비용·출력량을 보고한다. 이번 t1/t10 실행은 지정하지 않는다. root는 일반 구현/자원 보완을 자율 처리하며 관문 실패·새 물리 가정·명시 범위 초과 때만 중간 요청한다.

## 2026-10-19: 전체 stage0 원호 K 적분 실패와 제한된 동일 상태 재시도

Astra가 원격 `Codex_full_disk_initial_stage0_20261019/report.json`과 `Codex_saved_full_wall_quadrature_diagnostic_20261019/report.json` 및 현재 `curved_face_nodes.py`, `arc_velocity_roots.py`를 확인했다. 최초 full stage0는 closed massroot outer 3회에서 z 변화4.49414e-13, mass 잔차1.09267e-17로 수렴했으나 wall Gauss8/16 비교의 K가437셀에서 실패했다. 최대차2.219680777e-6(ID20110)는 기존 관문을 크게 넘는다. dr/F/H/incoming/area는 통과했다. wall sample normal speed의 실제 1–3회 부호교차와 기존 prior-only arc 분할은 K의 upwind 부호 꺾임을 적분하지 못한 구체적 원인 후보다. 질량·에너지의 mirror 상쇄가 K 적분 정확성을 보장하지 않는다. stage0는 미승인 상태이며 실패 산출물과 시도 횟수를 보존한다.

### 수정 범위와 수학적 관문

실제 네 kinetic sample velocity 모두의 arc normal-speed 영점을 prior-cut과 함께 분할하는 변경을 승인한다. 고정 stage의 Cartesian quadratic velocity를 n(theta)와 내적하면 삼각다항식 차수는 최대3이다. arc 중심으로 회전한 tangent-half-angle chart에서는 양의 분모(1+s²t²)^3와 최대6차 분자가 된다. 실제 affine/표준 quadratic basis/projection/residual map/L/m을 사용하며 mean velocity만의 영점으로 대체하지 않는다. 정지벽 wn=0이고 정확 mirror의 오른쪽 속도 영점은 왼쪽과 같다는 계약을 확인한다. density prior 교점 분할은 계속 유지한다.

chart는 pole을 포함하지 않도록 짧게 나누고 angle wrap 및 endpoint를 명시적으로 처리한다. 현재 최대 pi/2 chart 선택은 이에 맞는다. 계수는 실제 reconstruction에서 대수적으로 구성하고 샘플 polynomial fitting으로 대체하지 않는다. Cartesian 성분 자체의 cancellation까지 포함해 실제 evaluator와 대조한다. 현재 17점 대조와 2048eps primitive envelope는 구현 점검용 보수적 검사이며 전체 구간의 엄밀한 root 인증으로 부르지 않는다. 허용오차를 K gate로 전파해 완화하지 않는다.

현재 derivative recursion + critical point에서 fa==0 검사는 일반 접선·중근의 완전 검출을 입증하지 않는다. 근사 critical 위치의 작은 잔차로 even root를 놓칠 수 있고 거의 접하는 두 crossing도 부호 불확실성에 민감하다. derivative critical points를 추가 분할점으로 유지하고, critical 값이 polynomial 평가의 반올림 범위 안이면 unresolved로 분리해 구간 또는 더 높은 정밀도의 판정으로 해결한다. 단순 np.roots의 imaginary cutoff, sign-grid, 작은 계수/근 간격의 임의 제거는 사용하지 않는다. identically zero, 차수 저하, endpoint root, 접선/중근, 가까운 crossing, chart wrap 및 cancellation 사례를 fixtures에 포함한다. 동일근의 중복 제거와 실제 서로 다른 가까운 근의 삭제를 구별한다. 불확실 구간을 남긴 채 '모든 영점 인증'으로 보고하지 않는다.

### 저장 wall 검증과 stage0 재사용 계약

저장된 동일 physical state에서 새 분할의 Gauss8/16을 적용해 dr/F/H/K/ordered incoming/area를 모두 비교한다. 기존 성분별 gate `1e-14 + 1e-10*abs(order16)`를 유지한다. roots 및 critical/unresolved 수, chart별 분할점, 최대오차 셀을 저장한다. K의 smooth pieces 확보는 Gauss8 정확성을 자동 보장하지 않으므로 실제 관문 통과가 필요하다. 분할 변경으로 밀도·traction 적분도 달라질 수 있어 K만 바꾸고 다른 wall rates를 재사용하지 않는다. mirror mass/E/Phi 및 원점 torque, symmetry/finite/positive node weights도 확인한다.

valid initial geometry/density fit/velocity/gravity와 converged z는 입력 SHA를 고정해 재사용할 수 있다. 새 wall outputs 전체와 internal outputs를 합쳐 stage를 다시 조립하고 전체 mass residual, sum residual 호환성, J*1/gauge, area 및 energy/angular ledger의 기존 관문을 다시 적용한다. wall mass가 이론적으로0이라는 이유만으로 검사를 생략하지 않는다. 이 재검사 통과 시 outer solve를 다시 수행할 필요가 없다. 실패 시 old wall flux와 혼합하거나 residual을 보정하지 않고 중단한다.

### 이번 추가 예산 및 다음 실행

허용 범위는 저장 wall 수정/fixtures 묶음과 **동일 physical stage0 재조립·관문 재시도 1회**다. v1/v2 evaluator-envelope 실패 및 v3 저장 진단은 각각 보존하며 physical provider state 수와 implementation attempt 수를 별도로 기록한다. 현재 root 처리의 일반 보완은 root가 자율 완료한다. 재시도에서 관문이 다시 실패하면 무한 차수 증가, 반복 provider, dt 축소를 하지 않고 실패 묶음으로 평가한다.

통과한 stage0만 앞서 지정한 공통 stage0로 사용한다. 후속 h1e-3 한 step 및 h5e-4 두 step, unique physical provider 최대5개의 기존 예산과 종료 t=.001은 유지한다. 이번 동일 상태 수정은 추가 시간 stage가 아니다. 새로운 초기 fit/중력/전체 root의 불필요한 반복은 하지 않는다. 완료된 짧은 결과와 실제 비용으로 이후 dt/장시간 경로를 평가하며, 이번 repair나 stage0 통과 자체를 장시간 안정성 확보로 부르지 않는다.

## 2026-10-19: 전체 두 dt 묶음 평가 및 t=1까지의 다음 실행 경로

Astra는 `canonical_full_disk_two_dt_bundle_2026-10-19.json`, `canonical_full_disk_equal_time_comparison_2026-10-19.json`, fine-final physical comparison, saved pressure budget 및 원격 `Codex_full_disk_midpoint_prior_metadata_20261019/report.json`을 확인했다. 실제 19945셀/5개 unique provider/3개 accepted checkpoint가 완료되었고 질량 bitwise 보존, 양성, Ψ와 각운동량 장부가 통과했다. fine의 transported E+Ψ 변화−1.953e-18, reconstructed q48(initial q48 대비)−2.641e-16, origin L 변화−4.446e-18은 보존 검증이다. 각 셀의 정상장 정확도와 별개로 보고한다.

압력 RMS fine/coarse 비는 전체.32725, disk core.23307, inner.96180, outer contact.99332다. 전체 오차 감소를 공간 전체의 시간 수렴으로 확대하지 않는다. inner 압력 상대 RMS6.75e-7→6.49e-7에는 dt를 줄여도 남는 성분이 있다. 현재 오차 크기는 다음 물리 진화를 막을 근거가 아니다. 반대로 t=.001 오차를 t10까지 선형 외삽해 안정성을 승인할 수도 없다. endpoint pressure budget의 volume/raw thermal/factor RMS는 서로 상쇄한다. 직접 gravity thermal trace가 작다는 사실도 중력의 궤도·geometry·factor 경유 영향을 배제하지 않는다.

### 1. 저장 stage의 순간 압력 미분: 추가 RHS 없이 한 번의 분석

저장 initial geometry/gen/z와 raw stage/L/A/K/incoming을 사용한다. 필요하면 이미 저장된 4개 후속 stage에도 같은 분석을 일괄 적용한다. fit/root/provider/epsilon sweep은 추가하지 않는다. shared face의 실제 midpoint/normal/length로 Vdot=Σ(ds*wn), 정지 원호 wn=0을 구한다. gen/z를 실제 Laguerre 면속도 식에 대입한다. global ΣVdot=0과 local orientation을 확인하며 finite-step area change/h를 순간 Vdot의 정의로 쓰지 않는다.

고정 m에서 S=Q−qqᵀ/m, C=D−qpᵀ/m, R=LLᵀ를 쓴다. saved loss T는 TR+RT=LKLᵀ, Ac=A−T이고, Uin은 ordered incoming columns의 Σuuᵀ다. 실제 factor map의 h→0 tangent는 Rdot=AcR+RAcᵀ+Uin이다. Pi 재구성은 ppᵀ/m+CᵀS^-1C+R이며, 그 미분은 (pdot pᵀ+p pdotᵀ)/m+CdotᵀS^-1C+CᵀS^-1Cdot−CᵀS^-1 Sdot S^-1C+Rdot다. Sdot/Cdot에는 저장 qdot/Qdot/pdot/Ddot를 그대로 넣고 dm=0이다. inverse를 새로 명시 구성하기보다 solve를 사용한다.

Ethermal=E−tr(Pi)/2, P=(gamma−1)Ethermal/V에서 Pdot=(gamma−1)[Edot−tr(Pidot)/2]/V−P Vdot/V를 구한다. raw Pi RHS를 사용한 Pdot, factor tangent Pdot, 두 Pi 미분의 signed 차이를 각각 저장한다. hydro/직접 gravity/pressure source 분해는 기존 배열 범위에서 하며 gravity subset을 중복 합산하지 않는다. 정상장 Pdot_exact=0와 비교하고, m/q/Q/p/D의 정확 moving-cell derivative가 필요하면 기존 exact reference의 경계 적분만 이용한다. raw RHS=0 또는 물질점 궤도를 참조로 삼지 않는다.

지역별 signed 평균/RMS/max ID와 성분 간 signed 내적·공분산을 함께 기록한다. initial h*Pdot와 실제 endpoint pressure pattern의 상관을 보고하되 초기 tangent를 유한시간 진화의 증명으로 취급하지 않는다. Pi_rawdot−Pi_factordot가 반올림보다 크면 먼저 모델상 의도된 차이인지 문서/식으로 구분한다. 대수 구현 오류, 누락 항 또는 부호 오류가 확인될 경우 진화 전 수정·검증하고 보고한다. 단지 정상장 잔차가 nonzero라는 이유로 이번 진화를 무기한 보류하지 않는다. 이 분석만으로 analytic source 상쇄나 열 보정을 도입하지 않는다.

### 2. 실행 비용을 낮추는 하나의 구현 묶음

우선순위는 native prior bisection, outer 내부 불변 데이터 재사용, 감사와 실행 경로 분리다. curved_contact_priors의 local stencil/interval/LS normal/64회 이분법/phase 규칙 및 q24/fit tolerance를 유지한다. tiny contrast를 flat으로 바꾸지 않는다. coarse/fine predictor의 contrast median3.16e-8/8.01e-9와 비≈4는 forward predictor와 geometry의 영향이며 제거할 잡음이라고 확정할 수 없다.

* native 이분법은 우선 arc 없는 셀의 exact rational halfplane/area 경로를 공용 GMP 기반으로 옮긴다. 64개 binary64 midpoint와 area→binary64 비교의 기존 의미를 보존하고 exact branch witness를 남긴다. arc 셀은 기존 true-arc 경로를 유지해도 된다. 이것은 명시적 shape dispatch이며 물리 실패 fallback이 아니다. 새 오차 범위로 비교 방향을 바꾸거나 cut을 생략하지 않는다.
* 같은 stage의 priors/line endpoints/normals와 velocity sample의 t=0,.5,1 기본값은 z와 무관하므로 cache할 수 있다. 각 outer의 wn(z), relative speed, roots, 분할점과 새 node의 density/sample/flux는 다시 계산한다. 이동 quadrature node에서 이전 face-state 배열을 재사용하지 않는다. 원호 roots는 같은 stage의 고정 velocity/geometry일 때만 재사용한다. stage가 바뀌면 모두 무효화한다.
* Python provider와 검증된 native stage/factor 사이의 gravity-node 배열 등은 같은 process의 bounded lifetime buffer로 연결한다. 5.14GB/stage를 매번 파일로 왕복하지 않는다. geometry/fit/velocity/gravity/q24 입력은 변경하지 않는다. workers·BLAS thread·peak RSS를 실제로 확인하고 stage당 전체 복사를 중첩하지 않는다.
* 매 stage q48/모든 backend shadow/정상장 전체 reference를 생산 시간에 넣지 않는다. fit/rank/residual/양성/geometry completeness/root/gauge/질량/장부/원호8–16 gate는 그대로 매 stage 유지한다. 독립 q48와 정상장 reference는 지정 checkpoint에서 수행한다. CUDA 경로는 유지하고 OpenMP/MPI 개발을 새 병목 과제로 추가하지 않는다.

검증은 저장 초기 및 대표 moving stage 입력의 component oracle와 대표 moving provider의 baseline/candidate 각1회까지다. 가능하면 기존 immutable baseline 자료를 재사용해 baseline 재생을 생략한다. prior branch·offset·quadrature inventory, massroot 종료·z, full raw/Ψ/K/incoming/factor 결과를 기존 bitwise 또는 이미 정한 성분별 roundoff 관문으로 비교한다. root/fit 수렴 tolerance를 바꿔 속도를 얻지 않는다. 같은 stage·동일 seed·동일 thread 조건에서 exclusive 준비/provider/factor/I/O/audit 시간과 peak memory를 분리한다. native 이식 일부가 미완성이어도 검증된 나머지 개선을 이용할 수 있다. 성능 작업을 끝없는 사전 과제로 만들지 않는다.

### 3. 다음 물리 진화: t=.05 두 dt, 조건부 t=1 연속 진화

기존 immutable full-disk IC에서 h=.005와 h=.0025 두 경로를 t=.05까지 각각10/20 step 진행한다. 기존 t=.001 endpoint에 서로 다른 dt를 이어 붙여 시작하지 않는다. 검증된 동일 initial stage0는 재사용할 수 있다. 이것은 dt 확대의 검증 실행이며 사전 안정성 보장이 아니다. 같은 midpoint/Sylvester/exp/QR, Re∞, 무조향, fixed mass, Ψ constraint, true circular wall, q24 및 기존 모든 numerical gates를 유지한다. 자동 dt 축소·retry evolution·heat repair는 없다.

두 경로의 공통 checkpoint는 t=.01,.025,.05다. 모든 현재 cell 정상장 reference와 기존 region mask를 사용하고 mass/P/velocity(radial/azimuthal 포함)/ell/R/Psi를 평가한다. phase는 finite-volume reference mean velocity 방향과 비교하며 속도가0에 가까운 중심에서 무의미한 angle을 만들지 않는다. trajectory ID 차이와 current-cell exact-reference 오차를 분리한다. q48 constraint와 전역 E+Ψ/L 장부, wall torque도 검사한다. 두 dt 비교를 시간차수 증명으로 부르지 않는다.

진화의 임시 물리 감시 기준을 다음과 같이 명시한다. 이는 장시간 성공 허용오차를 사용자 대신 확정한 것이 아니라, 명백한 악화를 조기에 포착하는 운영 관문이다. 각 기존 비어 있지 않은 region에서 pressure relative RMS<.01, maximum<.1; velocity error RMS/reference velocity RMS<.01; 위의 유효 방향 cell에서 phase RMS<1 degree를 요구한다. 중심의 reference velocity RMS가 매우 작아 정규화가 불안정한 별도 subset은 absolute radial/velocity 값을 병기하고 전체 inner region 정규화와 혼동하지 않는다. positivity/fit/geometry/root/ledger 관문은 이보다 먼저 적용한다.

coarse 경로를 t=1까지 연장하는 조건은 세 공통 checkpoint 모두 위 감시 관문을 통과하고, coarse−fine current-cell pressure relative error 차이의 region RMS<.002, velocity error-vector 차이 RMS/reference speed<.002, 유효 phase error 차이 RMS<.2 degree인 것이다. fine도 전체 감시 관문을 만족해야 한다. 이는 초기의 거의0인 오차를 분모로 한 비율 대신 장시간 감시 기준의20%를 시간 오차 예산으로 둔 임시 기준이다. 불연속적인 성장·새 fit/root 문제는 원자료와 함께 보고한다. 이 조건이 실패하면 두 경로를 t=.05에서 끝내고 임의로 fine만 장시간 연장하지 않는다.

조건을 통과하면 root는 추가 건별 평가 요청 없이 coarse accepted t=.05를 h=.005로 t=1까지 이어간다(추가190 step). t=.1,.25,.5,1에서 같은 독립 진단을 수행하고 각 accepted step의 저비용 양성/pressure/global 장부를 기록한다. 다음 독립 checkpoint에서 감시 기준을 넘으면 상태를 보존하고 더 진행하지 않는다. 매 step reference/fit을 추가로 반복하지 않는다. 이 조건부 승인 전체는 최대220 accepted step, 공통 initial stage를 포함한 unique physical provider 최대439회(기존 initial stage 재사용 시 신규438회)다. 저장 자료를 읽는 진단과 성능 replay는 별도 횟수로 기록한다. timestep 실패 시 자동 축소나 추가 분기는 없다.

출력은 IC 원본 링크/hash, 두 경로 공통 checkpoint 및 coarse .1/.25/.5/1의 compact state/factor/Psi/seed, 매 step cell 요약·장부·수렴 기록과 실패 witness다. 모든 step의 수GB node 배열을 남기지 않는다. 시작 전에 실제 checkpoint bytes/free space 및 전체 예상 저장량을 계산한다. production 모드 실측 stage 시간을 바탕으로 439 provider와 factor 비용의 예상 walltime, 같은 h로 t6/t10에 필요한 비용을 보고한 뒤 지정 실행을 시작한다. 현실적인 비용이 아직 높다는 사실은 숨기지 않으며 이미 검증한 workers/cache/메모리 경로 안에서 자원 배치를 정한다.

### 4. 장시간 목표의 유지와 다음 평가 지점

이번에는 micro-test 추가만으로 종료하지 않고 조건부 t=1까지 실제 전체 동적 영역을 진화시킨다. t1 결과·장시간 비용·두 dt 비교를 모아 다음 평가를 한 번 요청한다. 후속 경로는 같은 후보의 t2→4→6→8→10이며, 오래된 A/B가 실패했던 t6–10 구간의 물리 정확도 확인이 목표다. t1 통과조차 장시간 성공으로 부르지 않는다. t6–10에는 pressure/위상/방사 운동/각운동량/Ψ와 경계 buffer 민감도 및 필요한 dt 검증이 남는다. 현재 stage 비용이나 짧은 endpoint 개선을 근거로 이를 생략하지 않는다. root는 위 범위의 일반 구현·검증·실행을 자율 처리하고 numerical gate 실패, 새로운 물리 가정 또는 명시 예산 초과 때만 중간 평가를 요청한다.

#### 같은 평가 중 도착한 tangent 완료 자료 반영

원격 `Codex_full_disk_saved_pressure_tangent_20261019/report.json`을 직접 확인했다. initial inner relative Pdot RMS는 raw6.485097991424e-4, factor tangent6.485097991387e-4이고 factor−raw thermal contribution RMS1.2102e-12/s다. initial ΣVdot3.3312e-20는 저장 면속도 roundoff 관문을 통과했다. initial Pdot*t=.001의 inner RMS6.48510e-7은 실제 fine6.49270e-7과 가깝다. 따라서 기존 instantaneous factor map이 raw thermal RHS와 서로 모순된다는 근거는 없고, 반이산 공간/closure 정상장 잔차가 초기 압력 drift를 설명한다는 근거가 강해졌다. 이 tangent 계산은 SVD/solve/LD 후처리이며 별도 certified bound가 아니다. root가 추가 보고한 per-cell alignment cosine.9992496/잔여RMS2.515e-8도 같은 해석을 뒷받침하며 해당 원자료를 최종 묶음에 포함한다. 이미 완료된 tangent 분석은 재실행하지 않는다. 새로운 RHS/epsilon 탐색이나 정상장 forcing은 추가하지 않는다. 다음 t=.05/t1 진화의 inner signed drift와 성장률을 추적하는 기준 자료로 사용하며 위 실행 계획은 유지한다.

## 2026-10-19: 사용자 지시 — 기존 병렬 본체에 같은 full-moment 모델을 연결해 병행 실행

사용자가 ‘기존 병렬 본체로 함께 실행’을 명시했다. 따라서 앞선 MPI 작업의 우선순위 제한을 이번 범위에서 해제한다. 실행 중 `Codex_full_disk_continuous_two_dt_to_t1_20261019`의 두 경로는 중단·재시작·모듈 교체하지 않는다. 기존 경로도 CPU workers8+CUDA를 사용하는 병렬 계산이며, 이번 작업은 기존 실행 파일의 MPI 제어 아래 같은 모델을 실행하는 별도 연결이다. 구형 보조메모리/RK4 모델의 독립 실행만으로 요청 완료를 선언하지 않는다.

Astra는 eunha2.c 진입점, Exam/Makefile, reversible_memory_runtime.h 및 Docs/full_disk_stage_memory.py, full_disk_factor_memory.py, run_full_disk_continuous_bundle.py, moment_backend_bridge.cu를 확인했다. 기존 본체는 현재 full-moment stage/factor를 호출하지 않는다. auxiliary memory 경로에는 MPI1/stress-sized/square-periodic 가드가 있다. 이를 제거해 full-moment 지원이라고 주장할 수 없다. 새 모델을 legacy particle 구조에 억지로 넣거나 RK4/dual-energy/periodic 후처리를 적용하지 않는다.

### 본체 진입점·빌드 정체성

`eunha2.c`의 MPI 초기화/checkarg 뒤, 기존 Hydro3D 분기와 같은 수준에 명시적 `FullMomentKepler` 모드를 추가한다. 구현은 `Exam/full_moment_runtime.{h,cpp}` 등 별도 translation unit으로 둔다. 기존 top-level 실행 파일 및 Exam 라이브러리에 연결하고 같은 MPI communicator를 전달한다. 기존 physics 진입점은 그대로 보존한다. 단순히 본체에서 전체 Python continuous runner를 exec하고 MPI rank들을 대기시키는 launcher 구현은 완료로 인정하지 않는다.

새 mode의 runtime이 authoritative state/ID 소유권, stage 순서, MPI 통신, accepted commit, checkpoint와 실패 처리를 제어한다. 물리 stage는 현재 검증된 provider와 native stage/factor ABI를 사용한다. midpoint half는 initial base/initial linear, full은 initial base/actual half linear라는 구분을 그대로 유지한다. production RK4는 legacy mode에 보존하며 새 모드의 적분기는 midpoint+Sylvester/exp/QR라고 명시한다.

본체/Exam/새 helper/bridge/Python 서비스의 commit 및 SHA, 링크 명령·라이브러리 경로, 실제 실행 binary의 SHA와 MPI/CUDA 버전을 manifest에 기록한다. 기존 Exam CUDA Makefile에 --use_fast_math가 있으므로 새 moment/geometry helper는 별도 명시적 검증 flags로 빌드한다. 미검증 fast-math/FMA contraction 설정을 기존 통과 helper에 상속하지 않는다. 기존 실행 바이너리를 덮어쓰지 않고 새 build/run root를 사용한다.

### 실용 provider 서비스: 전체 C++ 이식을 선행조건으로 삼지 않는다

typed C ABI의 persistent provider service를 둔다. CPython embedding/shim 또는 사전 기동된 rank-local service와 공유 메모리 IPC 중 환경에서 안정적인 방식을 root가 선택한다. 매 stage interpreter 재시작이나 GPFS에 수GB node 배열을 기록하는 RPC는 피한다. 서비스가 native runtime에 요청받는 것은 geometry/owned-cell prepare/face-root 자료 처리이며, 전체 timestep을 몰래 자체 수행하지 않는다. MPI는 native main thread에서만 호출한다. MPI/CUDA 초기화 이후 fork로 workers를 복제하지 않는다. 필요한 worker는 thread 또는 명시적 spawn/exec lifecycle로 구성한다.

계약 필드는 ABI version, run/step/stage/outer epoch, global stable IDs, x/W/raw18/m/L/Psi, z seed, 물리 설정 SHA, shape/stride/dtype/byte count, owner map, request status다. Psi의 현재 longdouble 정확도를 double로 낮추지 않는다. 실제 플랫폼의 itemsize/byte order를 확인하고 동종 노드에서 동일 표현을 사용한다. C caller가 버퍼를 쓰는 동안 Python 배열을 해제하지 않는 명시적 release/token 계약을 둔다. stale epoch/hash, 불완전 응답, rank 실패는 전체 stage 실패로 처리한다.

### 첫 MPI 분할: 전역 root는 하나, 셀 작업은 분산

최초 구현은 한 노드 1 rank 및 2 ranks다. stable-ID 정렬 순서의 deterministic contiguous ownership을 사용하고 소유권을 step 중 변경하지 않는다. 각 rank는 local authoritative cells를 소유하며 rank0는 stage transaction용 전역 compact state를 모은다. 전역 registry/geometry/neighbor prior stencil은 모든 rank가 읽을 수 있는 immutable stage 자료다. 전체 원형벽은 하나이며 rank 경계는 물리 경계가 아니다.

1. rank0에서 global moving Laguerre/true-arc geometry와 completion을 수행한다. 결과를 stage epoch와 함께 공유한다. 각 rank의 서로 떨어진 독립 원형 mesh를 만들지 않는다.
2. global area/mass/centroid/neighbor 정보를 동일하게 사용해 owned cell의 prior, constrained Psi density, quadratic velocity와 gravity quadrature를 rank별로 준비한다. global IDs로 compact fit/reconstruction/gen을 모으고 필요한 ghost 자료를 제공한다. ghost는 현재 stage의 읽기 전용 상태다. rank별 density phase/평균의 min/max를 쓰지 않는다. 처음부터 가장 큰 CPU 비용을 분산하므로 GPU 계산 몇 초만 분산하는 구성에 머물지 않는다.
3. 초기 연결에서는 rank0가 기존 ordered shared-face construction과 closed nonlinear massroot를 유지한다. internal shared face는 global key로 정확히 한 번 계산하고 전역 z/gauge와 모든 cells의 residual을 함께 푼다. outer마다 새 z를 공유하며 prior/geometry/velocity cache의 불변 범위는 기존 규칙과 같다. rank별 root, frozen anchor, reservoir 또는 boundary -1 치환은 금지한다. 중앙 root의 비용/직렬 비율을 숨기지 않는다.
4. 수렴 face의 native primitive flux는 canonical node당 한 번 계산한다. 최소 구현은 rank0의 기존 node kernel로 계산해 각 receiver owner rank에 같은 node result를 전달하는 방식이다. 각 rank에서 두 번 독립 flux를 구해 상쇄를 기대하지 않는다. 글로벌 양측 상태/xi/receiver affine 및 K/incoming에 필요한 원래 입력도 필요한 rank에 전달한다.
5. native stage ABI에 ‘precomputed node flux + owned cell CSR’ 경로를 별도로 추가하여 기존 stage_cell 연산을 재사용한다. 기존 lagmoment_stage는 전체 incidence를 요구하므로 partial 배열을 그대로 넣지 않는다. ghost raw를 업데이트하거나 내부 cross-rank face를 wall로 바꾸어 ABI 검사를 통과시키지 않는다. local CSR는 canonical global node 순서와 signed incidence를 보존한다. 각 rank의 GPU에서 owned cell hydro/gravity/covariance, 이어 factor 업데이트를 실제 수행한다. incoming columns의 순서와 Lbase/Llinear를 보존한다.
6. 모든 rank의 candidate/gates를 모은 뒤 collective 성공에만 accepted commit한다. mass/IDs와 psi mask, 전역 closed graph/gauge/장부를 확인한다. partial accepted checkpoint는 만들지 않는다. native runtime이 gen/W/raw/L/Psi/time을 업데이트하고 다음 실제 half stage를 요청한다.

전역 장부는 면마다 한 번의 canonical flux와 owned gravity로 계산한다. MPI reduction 순서 변경은 별도 기록하고 기존 roundoff 상계로 검사한다. 가능한 첫 구현은 stable-ID 순서 gather 후 동일 합산을 사용해 진단 비교를 간단히 한다. 각 rank의 ledger를 중복 전역 합산하지 않는다. global output은 stable-ID sorted checkpoint로 써서 1rank 저장값을 2rank로 읽을 수 있게 한다. legacy restart 바이트와 섞지 않고 format/model/version을 명시한다.

### 검증·실행을 하나의 묶음으로 수행

기존 저장 initial/moving 입력에서 1rank/2rank 분할·reassembly, CSR와 cross-rank sharedface, ordered incoming 및 owned gravity 합을 검사한다. ABI 실패/epoch mismatch/누락 ghost가 global no-commit으로 끝나는지도 검사한다. 필요 데이터는 현재 오래된 capture를 재사용하며 기존 진행 trajectory를 중단해 새 dump를 요구하지 않는다.

다음 실제 본체에서 동일 full-disk IC, h=.001 한 midpoint step을 1rank/2rank로 각각 수행한다. 각 실행은 실제 stage0와 실제 half provider를 호출하며 최대4개의 validation provider다. 저장된 이전 coarse h=.001 endpoint와 현재-cell reference도 비교한다. rank 수마다 동일 provider model/IC/z seed를 사용한다. stage/root/factor/보존/양성/원호8–16 관문 및 기존 성분별 roundoff 상계를 유지한다. MPI 순서 때문에 bitwise 일치가 없으면 그 차이를 정량화하고 허용오차를 결과에 맞춰 확대하지 않는다. 2rank의 실제 owned cell 수, 실제 GPU kernel 호출 수, face 단일계산 수 및 MPI 통신량으로 두 rank가 일을 했음을 확인한다.

위 관문 통과 후 2rank 본체에서 같은 IC의 h=.005 coarse를 t=.05까지 10step 실제 병행 실행한다. 현재 CPU8+CUDA 경로의 .01/.025/.05 checkpoints와 비교한다. 아직 기존 run이 그 시각에 도달하지 않았다면 해당 비교를 pending으로 남긴다. 새 MPI run만 더 먼저 통과했다고 두 dt 물리 관문을 건너뛰지 않는다. 현재 계획의 fine 경로와 같은 물리 감시·시간 오차 관문이 모두 완료되면 병렬본체 coarse도 t1까지 연장할 수 있다. t=.05 후 checkpoint에서 rank 개수를 바꾸는 자동 재시작은 하지 않는다. 신규 병렬 경로는 coarse 최대200 accepted steps/400 physical providers이며 validation4개는 별도다. 별도 CPU8+CUDA 실행의 기존 예산을 이 숫자로 숨기거나 바꾸지 않는다. 실패 시 원자료 보존, 자동 timestep 축소/legacy fallback은 없다.

### Slurm H200 실행과 진도·효율 보고

H200 자원에서 별도 Slurm job으로 실행한다. 실제 요구량은 1rank 검증에 `--nodes=1 --ntasks=1 --gres=gpu:H200:1`, 2rank 검증/진화에 `--nodes=1 --ntasks=2 --gres=gpu:H200:2`다. cpus-per-task는 우선8, 각 rank 서비스 workers는 CPU affinity 안에서 main/control 자원을 포함해 배치한다. BLAS1 및 phase별 OpenMP thread 수를 명시하고 두 thread pool을 중첩해 과할당하지 않는다. 메모리는 actual per-rank peak와 root gather peak로 산정하고 노드 가용량과 비교한다. 부족하면 kernel/prepare 배열 chunking 또는 lifetime을 고치며 old run 자원을 회수하지 않는다.

rank당 GPU 한 개를 scheduler binding으로 매핑한다. CUDA_VISIBLE_DEVICES의 실제 가시 UUID를 기록하고 local rank 숫자를 무조건 물리 device 번호로 해석하지 않는다. 기존 Ada 실행과 H200 job을 혼동하지 않는다. Slurm pending을 실행 시작으로 보고하지 않으며 jobid/allocation/hostname/rank/visible GPU/PID/binary SHA를 실제 시작 로그에 저장한다.

효율은 같은 H200 노드의 1rank와 2rank, 같은 initial/half 입력에서 geometry/fit/velocity/gravity/face-root/GPU stage/factor/통신/감사 exclusive wall과 rank별 max/평균, 총 CPU/GPU hours를 보고한다. speedup=T1/T2, parallel efficiency=speedup/2는 같은 resource-per-rank 설정일 때만 계산한다. Ada 대비 차이는 hardware까지 달라져 MPI speedup으로 부르지 않는다. 새 본체의 rank0 face/root 병목을 별도 보고하고 speedup을 사전 보장하지 않는다. 2rank가 느려도 동등성 결과는 보존하며 장시간에 무조건 더 많은 GPU를 쓰지 않는다.

매 step accepted time/dt/cells/root iterations/physical gates 및 최근 step wall을 기록하고 두 진행 경로의 실제 accepted 진도를 별도로 보고한다. routine 구현/빌드/Slurm 제출/검증에는 재승인 질문이 필요 없다. 본체 연결의 물리 관문 실패·예산 초과·새 물리 가정만 묶어서 평가 요청한다. 목표는 같은 full-moment 모델의 병행 진화와 장시간 물리 정확도이며, binary 실행 또는 MPI_Init 성공 자체가 완료 조건은 아니다.


### 2026-10-19 묶음 정적 평가 — native runtime RELEASE 경쟁
기존 Astra agent의 검토(추가코딩/실행/polling없음)는 full_moment_evolution/native_factor/native_stage/state_io/provider/face/factor/global_audit를읽고 P1 한건을확인했다. v7 RELEASE에서각rank가meshSHA를검사하는동안rank0가자기RPC에서sharedgeometry를삭제가능;RPC후collective로는경쟁을막지못한다. ranklocalACK 전체collective후 별도rank0삭제단계필요. 수식/소유권에서 Lbase/Llinear,고정mass, native x/W/raw/L/Psi와LD,canonicalflux/CSR/incoming,전관문후collectivecommit은중대오류발견없음. 이는실행검증완료주장이아님. restart새token import/format검사유지 및durablefailurecandidate보존,fsync/atomicacceptedpublish는t1전필수. native단독으로시간오차gate없으므로현재t.05제한및외부coarse/fine비교후연장유지.
Root는유효v71rank한step완료를보존하고후속2rank초기준비에서oldbatch를취소했다. v8는공유releasehelper allrankACK→Allreduce→rank0DISPOSE_GEOMETRY→Allreduce로수정했고rank1RELEASE1.5sec지연시험427897통과. fsync/atomicpublish/포맷roundtrip및초기보존도동시에통과. 새로운2rank검증427898이afterok로실행중. 상세근거/범위/남은restart계약은Docs/reversible_residual_memory_2026-10-18.md 마지막절참조. 완료1rank를다시돌리지않고v8결과와정확비교하며controlbinary차이때문에MPI성능벤치라고부르지않는다.


## 2026-10-19 묶음 평가: cross-host seed와 restart v9/v10
기존Astra정적평가(새계산/테스트/폴링없음)에서새P1없음. samehost originalglobalstage+Pythonmidpoint가nativeowned/timecontrol을우회하는bridge exactvalue검증으로타당. 단crosshostseedfailure해결/호스트산술원인확정/공유물리연산자독립정확성증명으로부르지말것. np.array_equal은signedzero/LDpadding byteidentity가아니며 exactvalue로표현. Restartpayloadbyteidentity는별도실제hash시험으로입증한것과구분.
최소추가증거: IC/physics/seed/provider/operator 실제SHA,host/GPU,제어binaryv7/v8차이,original두RHS+half/fullfactor관문,oldcrosshostfailure/newsamehostpass별도기록. Root는새RHS없이427901allocation내명시된originalPID695089/직계자식713702만읽는runtimeidentity캡처. 실제mappedv4operatorlibSHA51014e.../Pythonexe/provider4모듈SHA/ICsha/v7v8v10snapshotsha/hostsyn104/실제부모GPUcontextUUIDc5feb54b72313bc24cffb61c4cf3bac3기록. Docs/canonical_full_disk_original_live_runtime_identity_2026-10-19.json. Nativecoarse실행뒤같은readonly도구로nativepid/실제mappedoperator/provider를추가연결할것. 고정OPERATOR파일SHA만으로실제로loadedlibrary/provider전체동일성증명아님(P2).
Astra는provenance→nativeimport→실제x/Wgeometry→q48/reference/admissibility-beforeRHS 및v10추가200step/400RHS/globalindex구분/최종checkpoint설계적절판단. 실제restartphysicalstep은필수. 다음순서유지: samehostgate→v8coarse.05;native.01에서v10.015한step;uninterruptedv8.015와state/factors/Psi/mass/IDs/seed전필드exactvalue,importpayloadhash/oldsource불변/newRHS2/finalq48물리관문검증;Ada두dt.01/.025/.05와nativecoarse비교후만.05→1연장. 비교실패시tol증가/seed덮어쓰기금지. t>1budget기능을t6–10안정성증거로확대금지. 추가Astra평가는다음의미있는수치경계/물리실패시묶음으로만요청.

## 2026-10-19: 밀도 drift와 접촉 압력 성장 — t1 자동 연장 조건 보완

Astra는 `compare_full_moment_checkpoint_errors.py`, `diagnose_saved_full_moment_density_volume.py`, `diagnose_full_moment_subcell_energy.py` 및 원격 checkpoint-error 진단의 saved_density_volume_drift/coarse_t0p025_subcell_energy/fine_t0p025_subcell_energy report를 읽었다. 계산·테스트·job polling은 하지 않았다. 새 자료에 따라 이전 자동 .05→1 승인을 다음 조건으로 보완한다. **진행 중 Ada fine/native .05 비교와 지정된 restart 동등성 검증은 그대로 완료한다. .05→1 신규 child는 아래 완료 평가 전에 자동 시작하지 않는다.** 기존 accepted 상태를 취소하거나 진행 중 궤적을 재시작하지 않는다.

### 물리 평가

초기 uniform cavity754셀의 밀도 상대 RMS1.77e-15와 현재 밀도−고정질량 체적 항등식 오차~4e-15는 초기 밀도 투영 오류가 현재 .025의 RMS.00549/.00441을 설명하지 못함을 보인다. 초기와 현재 모두 cavity에 속하는 셀에서는 m/(rho_floor V)−1이 밀도 오차와 같다. 이는 fixedmass의 운동학적 결과이며 체적 적분기, 생성자 속도, 수치 face density, massroot 또는 stress closure 중 하나를 원인으로 특정하지 않는다. initial Vdot 선형 예측의 낮은 alignment는 초기 압력 tangent의 성공을 밀도에 그대로 적용할 수 없음을 보여 준다. fixedmass cell은 material cell이 아니므로 실제 rho/velocity divergence를 함께 확인해야 한다.

subcell 에너지 식 E=Kcm+Kaffine+tr(R)/2+thermal은 맞으며 ell²/(2trS)는 Kaffine의 하한이다. 이를 별도 회전에너지로 더하지 않는다. 큰 R 상대오차와 작은 unresolved-energy/thermal 오차는 동시에 가능하다. 작은 thermal 척도 기여가 R의 양성/방향/각운동량 중요성을 없애지는 않는다. affine-energy 오차가 상대적으로 더 크지만 signed E/운동/체적/열 상쇄를 보지 않고 압력 성장의 단독 원인으로 정하지 않는다.

coarse inner-contact 압력 RMS는 .025의.0003348에서 .05의.0075707로 약22.6배 증가했다. 기존 감시 기준을 통과했다는 사실과 이 급격한 성장은 함께 보고해야 한다. 현재로서는 장시간 실패가 확정되지 않았지만, 남은190step을 기존 절대 기준만으로 자동 승인할 근거는 부족하다. fine .05 결과를 기다리며 두 경로의 지역별 signed pattern/max ID/분포/시간 변화와 정상장 오차를 종합한다. 한 endpoint나 두 dt만으로 차수·성장률의 장시간 지속성을 확정하지 않는다.

### 밀도와 연장 관문

기존 numerical/보존/양성/P-v-phase 기준은 유지하고, 같은 actual-cell reference의 e_rho=(m−m_ref)/m_ref=(rho−rho_ref)/rho_ref를 추가한다. 모든 기존 비어 있지 않은 region에서 numerical-mass-weighted RMS<.01, max<.1을 임시 물리 감시 기준으로 삼는다. 기존 P/v 감시와 동일한 규모를 사용한 운영 기준이며 최종 과학적 정확도 허용치를 확정한 것이 아니다. 초기 고정 cavity mask와 현재 cavity mask의 교집합/이탈 cell도 별도로 보고하고 contact와 bulk를 섞어 밀도 오차를 숨기지 않는다.

동일시간 coarse/fine e_rho 차이 RMS<.002를 각 region에 추가한다. .01/.025/.05 모두에 적용하고 기존 pressure.002/velocity.002/phase.2degree 관문을 유지한다. 기존 완료 결과는 새 진단으로 표시하며 이전 accepted를 다시 쓰지 않는다. .05 density/temporal 관문 통과만으로도 이번에는 자동 t1 연장하지 않는다. 접촉 압력 성장과 cavity 밀도 drift를 포함한 완료 묶음을 한 번 평가해 연장 길이·dt 또는 수정 후보를 정한다. 기준 실패 시 tolerance 확대나 fine 경로만의 자동 장시간 대체는 없다.

### 실패 시 최소 분리 진단과 새 접근의 순서

1. 우선 저장 .01/.025/.05의 같은 ID와 공통 mask에서 density/pressure/radial velocity/ell/R/affine energy의 signed 변화 및 최대 오차 cell을 한 번에 모은다. reference quadrature 차이와 mask 변화도 분리한다. 초기해 재주입이나 새 density fit을 진단 때문에 반복하지 않는다.
2. 저장 stage 자료가 남아 있는 구간에서 실제 geometry Vdot와 accepted V 변화의 budget을 작성한다: 각 step δV_map=V_next−V_base−h*Vdot_half. 같은 구간의 mass/thermal/Pi factor budget과 함께 본다. 실제 stage의 face wn으로 Vdot를 계산하며 endpoint difference/h를 Vdot라고 정의하지 않는다. 전체 ΣδV와 cell별 누적 signed 값을 모두 확인한다. δV_map은 유한시간 geometry update와 midpoint 적분의 불일치이며 그 자체만으로 기하 구현 버그는 아니다. smooth 구간의 기대 차수와 topology/prior branch 변화 구간을 구분한다.
3. 다음으로 현재 reconstructed density/mean velocity/sample flux의 질량 제약과 geometric volume rate를 연결한다. uniform-reference cavity에서도 numerical trace는 상수가 아닐 수 있다. zero net mass residual만으로 Vdot=0이 보장되지 않는다. 같은 node의 mean advective·residual kinetic·mesh-motion 항과 pressure/thermal flux를 signed 합으로 분리한다. gravity의 직접 thermal 상쇄와 중력을 통한 궤도·geometry feedback도 구별한다.
4. 필요한 stage가 성공 후 이미 회수되어 저장 budget으로 원인을 분리할 수 없다면 ‘자료 없음’으로 명시한다. 이 경우 temporal gate 실패 뒤 **저장 coarse t=.025 상태의 일관된 provider 재평가 최대1회**를 승인한다. 현재 geometry/raw/L/Psi/seed/동일 closure를 고정해 volume·flux·pressure·factor tangent를 동시에 capture한다. 원래 stage와 동일 seed가 없으면 동일 재생이라고 주장하지 않는다. trajectory/epsilon sweep/복수 closure counterfactual은 추가하지 않는다. 한 instantaneous RHS만으로 .025→.05 누적 오류를 전부 귀속할 수 없다는 한계를 유지한다.

새 접근은 위 분리 결과에 따라 하나를 선택한다. finite-step geometry-volume defect와 시간 오차가 지배적이면 시간에 일관된 면 sweep/면속도와 전역 massroot를 묶는 geometric conservation 경로를 우선 설계한다. smooth oscillatory coupling의 시간 증폭이 지배적이면 전 상태/geometry/W/Psi와 factor를 함께 평가하는 implicit midpoint 등 대칭 시간결합 후보를 검토한다. explicit midpoint의 선형 test equation iω에 대한 |1+iωh−(ωh)²/2|²=1+(ωh)^4/4는 장기 증폭 가능성의 수학적 경고이지 현재 nonlinear solver의 원인 증명은 아니다. factor exponential/QR 부분만으로 전체 방법이 symplectic이라고 부르지 않는다.

반대로 instantaneous mass/pressure equilibrium residual이 지배적이면 reconstruction/face mass·stress·thermal 구성식의 호환성을 먼저 수정한다. 단순 integrator 교체로 공간 잔차가 사라진다고 가정하지 않는다. 어떤 경우에도 원형 정상장을 강제로 유지하는 source subtraction, analytic velocity steering, heat repair, tiny-contrast flat화, W clamp를 넣지 않는다. 고정mass·Re∞·무조향·mean-free 가중치율과 공통 상대W의 Voronoi 환원 조건을 유지한다. 지금 승인하는 것은 분리 진단까지이며 새 integrator/closure 진화는 완료 증거로 한 후보를 결정한 뒤 시행한다. 생산 RK4 보존은 계속하되 검증 없이 새 모델에 대입하지 않는다.

### Native/Ada 동일 h 비교의 기준

새 checkpoint error 비교기의 actual-cell error difference와 mask union은 물리 재현성 진단으로 적절하다. union뿐 아니라 양쪽 mask cell 수/교집합/한쪽에만 속한 ID 및 phase-defined 제외 cell 수를 기록한다. 현재 load가 ‘physical gate passed’만 허용하므로 실패 결과의 원인 진단에는 실패 상태도 읽되 passed와 분리하는 경로가 필요하다. 이를 acceptance 우회로 사용하지 않는다.

동일 호스트 original/native exact-value 및 native1/2rank 관문, v10 restart .01→.015 exact-value 관문은 그대로 유지한다. Ada와 H200 사이의 seed 실패를 완료 처리하거나 모든 checkpoint에 bitwise 일치를 강요하지 않는다. 새 코드의 같은-h 비교는 독립적 물리 재현성 관문이며 roundoff 인증이 아니다.

연장 평가에는 두 플랫폼 각각의 numerical/P-v-phase/density 관문과 .01/.025/.05의 같은-h error difference가 필요하다. 플랫폼 간 pressure/density error 차이 RMS<.0002, velocity error 차이/reference-speed RMS<.0002, phase error 차이 RMS<.02degree를 각 region의 임시 **플랫폼 오차 예산**으로 둔다. 이는 기존 temporal 예산의10%를 할당한 물리 기준이며 기존 seed/backend 수치 허용오차의 완화가 아니다. 실제 coarse/fine 차이와 비교한 상대 크기도 병기한다. 플랫폼 차이가 시간 오차만큼 커지거나 특정 cell에서 급성장하면 고정 예산 안이라도 연장 전에 원인을 평가한다. raw/state/seed/geometry 및 R/ell/열 정규화 energy 차이도 자료에 포함하며 새 물리 예산만으로 ‘동일 구현 인증’이라고 표현하지 않는다.

목표는 여전히 t6–10에서의 물리 정확도와 안정성이다. 이번 .05 비교를 조기 종료하지 않고 문제를 드러낸 상태로 끝까지 보존한다. 병렬 bridge의 성공이나 전역 에너지 보존이 밀도·접촉 압력의 정확도를 대신하지 않는다.

## 2026-10-19: 보존된 native .015→.02 구간의 bounded volume-map 진단 승인

Astra가 canonical saved_t0p025_signed_budget/tangent_growth/volume_rate_crosscheck 및 coarse_t0p025_to_t0p05_finite_error_budget JSON을 확인했다. .025 instantaneous face-plane Vdot가 capture와 fractional rate2.68e-14 이내로 일치하는 것은 순간 기하식의 교차검사다. .025 압력 tangent가 후속 .05 변화와 반대 방향인 사실, 작은 순간 factor−raw 차이, 유한구간 energy/kinetic 항의 큰 상쇄는 추가 시간정보가 필요함을 보여 준다. unresolved-energy 오차의 높은 alignment도 독립적인 원인 증명은 아니다.

보존 native accepted3(.015), actual half e8(.0175), accepted4(.02)의 **3 geometry 재구성/0 RHS/0 fit/0 root/0 진화** 한 구간 진단을 승인한다. 이 결과를 다른 Ada 궤적의 .025→.05 구간 또는 장시간 전체의 attribution으로 확대하지 않는다. 기존 fine/native .05 진행 및 t1 보류는 그대로다.

### 입력 및 역산 계약

세 상태의 원본 파일 hash/IDs/mass/time/epoch/run lineage와 실제 h=.005의 binary64 값을 고정한다. h를 .02−.015의 부동소수점 차로 새로 정의하지 않는다. 실제 full_moment_build_candidate 및 실행 binary의 연산 옵션을 확인하여 x_next=RN(x_base+h*g_half), W_next=RN(W_base+h*z_half) 구현임을 확인한다. half가 단순 endpoint 평균이 아니라 실제 provider 입력이었다는 lineage를 보존한다. accepted next의 weight_rate_seed가 후속 root seed나 gauge 재중심화된 값으로 덮어써지지 않았으며 실제 full update의 z_half라는 계약을 확인한다.

각 좌표의 g_half는 하나의 역산값이 아니라 enclosure로 취급한다. binary64 b의 RN preimage I(b)는 prev(b)/b/next(b) 사이의 midpoint 경계로 구성한다(외향 포함은 ties 처리에 보수적이어도 됨). FMA이면 g∈(I(x_next)−x_base)/h. multiply/add가 각각 반올림되면 p=RN(h*g), x_next=RN(x_base+p)이므로 먼저 p∈I(x_next)−x_base 안의 가능한 binary64 p 범위를 구하고, 그 p들의 RN preimage 합집합을 감싸 h로 나눈다. 정확 rational 또는 올바른 outward 연산으로 구현한다. signed zero/subnormal을 포함하고 overflow/nonfinite는 거부한다. 실제 contraction을 확정하지 못하면 두 연산 계약의 enclosure 합집합을 사용한다. 단순 relative eps*|g| 또는 (next−base)/h 한 값만으로 엄격 bound를 주장하지 않는다.

중심 추정치는 저장 좌표의 정확 차/h에서 얻되 추정임을 표시한다. 별도 저장된 z_half로 W update의 rounding-cell 포함 여부를 검사한다. 이 검사가 실패하면 seed lineage 또는 update 가정이 맞지 않으므로 진단을 중단하고 원자료를 보존한다. x displacement가0인 경우에도 g interval을0으로 강제하지 않는다. 생성자 조향/후처리/reset이 있었던 cell에 이 역산을 적용하지 않는다.

### 독립 geometry와 오차 전파

세 실제 x/W로 같은 true-disk/complete Laguerre geometry를 구성한다. half의 각 internal face에서 d=xj−xi, d_dot=gj−gi, b=(|xj|²−|xi|²+Wi−Wj)/2, b_dot=xj·gj−xi·gi+(zi−zj)/2를 사용한다. face outward i→j인 경우 wn=(b_dot−d_dot·x_face)/|d|이다. 직선 face의 affine speed 적분은 length*wn(midpoint)로 정확하며 정지 circular wall은0이다. global canonical face를 한 번 계산하여 두 owner에 반대 signed 기여로 넣는다. 저장 half geometry의 x/W는 고정 입력이고 g만 역산 구간이며 z는 저장 실제 값이다.

cell Vdot는 g에 선형이므로 가능한 경우 global generator별 계수를 먼저 모아 enclosure를 전파한다. 같은 g가 여러 면에 등장하는 의존성을 무시해 큰 상쇄를 잃지 않도록 한다. global ΣVdot의0은 shared flux 반대부호 구조의 검사로 유지한다. exact rational support/edgevector를 재사용하여 미세면을 삭제하거나 rounded endpoint 차로 길이를 대체하지 않는다. geometry coefficient 및 계산 합산 오차도 포함한다.

δV_map=V_next−V_base−h*Vdot_half의 interval을 계산한다. base/next area 오차를 제외한 채 generator enclosure만으로 전체를 ‘엄격 bound’라고 부르지 않는다. non-arc cell의 rational area는 정확하게 계산할 수 있다. arc area/길이/중점 등 transcendental 값은 검증된 outward bound가 있으면 사용하고, 없으면 해당 cell은 `geometry_bound_unavailable`로 명시해 전체 certified 통계에서 제외한다. 경험적 차수 차이를 엄격 enclosure로 대체하지 않는다. 3geometry 전체를 재구성하되 엄격 결과와 경험적 결과의 적용 cell 수를 구분한다.

보고는 cell별 δ 추정/하한/상한, δ/V 및 밀도 척도 기여, zero 포함 여부, bound 폭의 주요 원천, max ID, 기존 regional signed합/RMS를 포함한다. 실제 ΔV와 hVdot를 별도로 보존한다. 여러 topology 변화가 있어도 endpoint와 실제 midpoint의 map residual은 정의되지만 smooth O(h³) 법칙이나 중간 경로의 적분오차 분해를 자동 가정하지 않는다. 한 step에서는 시간차수를 판정하지 않는다. shared cancellation과 disk 전체 area 불변도 교차검사하되 그 항등식 자체를 local δ의 정확성 증거로 삼지 않는다.

### 판정과 충분하지 않을 때의 종료 조건

interval이0과 명확히 분리되는 cell에서는 반올림 enclosure를 넘어선 finite-step map residual을 검출했다고 할 수 있다. 이것이 곧 버그 또는 현재 압력·밀도 drift의 주원인이라는 뜻은 아니다. 관련 region에서 bound가 residual과 비슷하거나 geometry bound가 없으면 그 범위는 unresolved로 보고하고 이 한 구간 진단을 종료한다. 추가 RHS/epsilon/dt sweep이나 누락 generator를 얻기 위한 과거 step 재실행은 승인하지 않는다.

불충분한 경우 다음 새로 승인될 진화 구간에 직접 journal을 설계한다: base/half/next state hash와 실제 h, global IDs, actual half gen/z 배열, cell Vdot 및 계산 error estimate, base/half/next area, canonical geometry hash/face count/topology 변화, operator/build SHA를 저장한다. 이는 cell-size 자료이며 대형 face-node dump를 요구하지 않는다. 기록은 update 전에 실제 solver 값에서 수행하고 endpoint 역산값으로 덮어쓰지 않는다. 지금 journal을 위해 진행 중 immutable job을 교체하거나 새 trajectory를 시작하지 않는다.

## 2026-10-19 완료 .05 묶음 평가: 자기일관 midpoint 후보와 의미 보존 성능 작업

계획·평가 담당 Astra의 통합 결정이다. completed_density_temporal_budget, fine_t0p025_to_t0p05_finite_error_budget, coarse/fine_t0p05_subcell_energy, native_restart_equivalence, native_ada_t0p05_platform_budget, runtime_cost_baseline, saved_density_cost 및 앞선 tangent/volume-map 자료를 근거로 한다. 이번 검토에서는 계산·폴링·과학 코드 변경을 하지 않았다. 아래 구현과 실행은 운전자 root가 한 묶음으로 진행하며 routine 변경마다 추가 평가를 요청하지 않는다.

### 1. 완료 결과와 연장 판정

full disk fine h=.0025는 t=.05까지 20step/40RHS를 완료했다. 기존 양성·보존·q48·P/v/phase 관문 통과는 유효하다. 그러나 새 밀도/시간오차 관문까지 통과한 것은 아니다. cavity density RMS는 coarse .010219127, fine .007606562이며 coarse는 절대 RMS .01 기준도 실패한다. 두 경로의 cavity density 차이 .00410605111, inner-contact density .00260498440, pressure .00756272938은 각각 .002 시간오차 예산을 초과한다. .025 cavity 시간오차 .00200579도 이미 실패했다. **기존 coarse/fine 또는 native accepted .05에서 t1로 자동 연장하지 않는다.** fine만 골라 장시간 진행하는 것도 승인하지 않는다. 완료 궤적과 coordinator 보류 상태는 보존한다.

contact P RMS가 .00757073589에서 .000137585742로 약55배 작아진 것은 시간 처리에 민감한 오류가 있음을 보여 준다. 두 dt로 차수 또는 특정 선형 불안정을 확정하지 않는다. fine에서도 cavity density drift가 남는다. initial spatial Pdot bias, 한 구간의 유한 geometry-map 잔차, 뒤 구간에서 tangent와 반대 방향인 압력 변화, 큰 signed energy 상쇄가 함께 관찰됐다. 따라서 ‘시간적분만 바꾸면 해결’, ‘geometry만 원인’, ‘four-stream만 원인’이라는 결론은 모두 이르다. 초기 conditional residual energy의 작은 규모만으로 kinetic transport 영향을 배제하지도 않는다.

native/Ada 같은 h의 .05 물리 차이는 기존 플랫폼 예산보다 매우 작고, native 실제 restart .01→.015는 seed 포함 전필드 exact value 및 여섯 payload byte-preservation을 통과했다. 이 범위에서는 병렬 소유권 또는 restart 오차를 현재 큰 물리 오차의 우선 원인으로 삼을 근거가 없다. 이전 cross-host seed 수치상계 실패는 별도 실패 기록으로 유지한다. bridge/restart 통과는 물리 정확도 또는 t6–10 안정성의 대체 증거가 아니다.

### 2. 이번에 선택하는 물리 후보 하나

**기존 공간 연산자를 유지한 자기일관 midpoint 결합**을 선택한다. 현재 explicit predictor의 half 상태를 고정하고 끝내는 대신, endpoint와 그 midpoint의 geometry/W/raw moments/Psi를 함께 수렴시킨다. density closure, 네 stream, pressure 구성식, gravity, true circular wall, massroot 및 quadrature 규칙은 바꾸지 않는다. finite-step geometry와 pressure/kinetic coupling을 함께 시험할 수 있다는 점이 선택 이유다. geometry-map 잔차의 정확 소거 또는 공간 평형 잔차 제거를 보장하는 방법으로 소개하지 않는다.

구체적 정의는 다음과 같다. endpoint trial Y*는 x,W,raw18,Psi를 포함한다. mass는 원본 고정값이며 W gauge는 기존 mean-free 계약을 따른다. provider 입력은 Y_mid=(Y_base+Y*)/2의 **물리 raw moments와 x/W/Psi**다. L 자체를 평균하지 않는다. midpoint raw에서 S,C 및 R=Pi−ppᵀ/m−CᵀS^-1C를 계산해 기존 convention의 Cholesky factor를 복원한다. 실제 midpoint geometry 위에서 여섯 feature density와 quadratic velocity를 다시 준비하고 global massroot를 수렴시킨다. midpoint admissibility/fit/rank 실패는 실패이며 floor나 analytic target 보정으로 통과시키지 않는다.

그 midpoint provider로 기존 full-step map M_h(Y_base,F_mid)를 만든다. x/W와 비-Pi raw, Psi 및 장부는 동일 midpoint rate로 갱신한다. covariance는 기존 exponential/ordered QR map을 사용하되 **Lbase는 원래 base, Llinear는 trial midpoint**다. endpoint Pi를 새 q/Q/p/D와 갱신 R로 재구성한다. 풀어야 할 방정식은 Y*=M_h(Y_base,F((Y_base+Y*)/2))이다. endpoint Pi/R/thermal도 잔차에 포함한다. 한 번의 추가 corrector만 적용하고 이를 implicit convergence라고 부르지 않는다.

이는 positive covariance map을 결합한 자기일관 midpoint이며 표준 implicit RK midpoint 또는 symplectic integrator와 동일하다고 주장하지 않는다. covariance injection의 유한-step 오차는 남을 수 있다. 구현 문서에는 frozen-stage covariance map의 작은 h 전개와 raw tangent 일치, h=0, 회전/zero-incoming/비대각 covariance, Lbase≠Llinear 사례를 포함한다. 공간적 steady source subtraction, 정상 궤도 조향, 열 보정은 없다. production RK4는 보존한다.

비선형 해법은 결정적 Picard를 기본으로 하고 필요하면 제한된 Anderson 또는 damping을 같은 잔차 방정식에 사용할 수 있다. 변경 이력을 기록하며 해법이 물리 state를 clip하지 않는다. 큰 step 실패 때 dt 자동 축소, explicit fallback 또는 다른 closure로의 전환은 없다. density/prior branch 변화 수를 기록한다. massroot의 내부 수렴과 이 외부 시간결합 수렴은 별개로 검사한다.

### 3. 수렴 및 장부 계약

수렴 판정은 연속 iterate 차이만으로 하지 않는다. 마지막 trial midpoint에서 새 provider를 평가해 실제 map residual을 확인하고, accepted 값은 검증된 M_h 출력으로 commit한다. 그 출력에서 다시 정의되는 산술 midpoint와 실제 마지막 provider midpoint의 불일치를 저장한다. 따라서 마지막 확인 provider도 예산에 포함한다.

기존 성분별 수치 roundoff 관문을 유지하고 외부 시간결합에는 추가 운영 기준을 둔다. 각 cell의 residual을 raw 원점 모멘트만으로 정규화하지 말고 centroid, mean velocity, centered S/C/R, thermal, geometry 및 Psi remainder로도 검사한다. 길이는 base sqrt(area), 속도는 base mean-speed와 sound-speed 중 큰 값, covariance는 각각 m*l², m*l*v, m*v²로 정규화한다. x/centroid와 이 centered moment/velocity 정규화 최대 residual은 1e-10 이하를 요구한다. W는 공통 gauge 제거 후 l²로 정규화한다. thermal residual은 base thermal로, Psi는 density fit에 쓰는 실제 potential-remainder scale과 그 기존 rank 분기로 정규화하여 각각 최대 1e-8 이하를 요구한다. dependent potential 분기는 기존 dependent constraint를 직접 검사한다. 양수 물리 scale을 임의 floor로 키워 residual을 숨기지 않는다.

실제 provider midpoint와 accepted 산술 midpoint 사이의 cell density/pressure 차이도 현재 geometry로 확인해 각 region RMS 1e-8 및 최대 1e-7 이하를 요구한다. 이 숫자들은 nonlinear solve 오염을 기존 .002 물리 예산보다 충분히 작게 제한하는 추가 운영 기준이며 엄격한 해 오차 상계나 기존 tolerance 완화가 아니다. 작은 residual만으로 nonlinear conditioning을 인증하지 않는다. 기존 수치 관문을 통과했더라도 residual 정체/branch 순환이면 실패로 종료한다.

canonical face/global gravity의 동일 midpoint ledger로 E+Psi, p, origin angular momentum을 검사한다. midpoint 산술 평균을 맞추기 위해 endpoint E/Psi/L을 사후 수정하지 않는다. 반올림과 nonlinear residual의 보고된 기여를 구별한다. 기존 fixed m, symmetry, R/thermal positivity, wall8/16, root gauge, q48 및 current-cell reference 관문은 유지한다.

### 4. 유한 실행 범위와 의사결정

새 후보 두 경로는 **동일한 보존 fine t=.025 checkpoint**에서 출발한다. h=.005로 .05까지5step, h=.0025로 .05까지10step이다. 이는 모델의 integrator를 명시적으로 바꾼 fork이며 과거 궤적 재현 또는 같은 integrator restart 검증으로 부르지 않는다. 양쪽 initial state/factors/Psi/mass/seed와 IDs를 동일하게 import한다. 기존 explicit fine .025→.05가 같은 초기값의 대조 경로이므로 이를 그대로 사용한다. old coarse .025는 다른 초기값이므로 새 coarse와의 직접 시간오차 비교 기준으로 삼지 않는다.

먼저 새 coarse의 .025→.03 한 step과 fine의 .025→.03 두 step을 같은 승인 안에서 수행한다. 이 세 step이 모든 수치/물리/비선형 관문과 기존 시간오차 예산을 통과하면 나머지 .05까지 진행한다. 첫 구간 실패 시 remainder를 시작하지 않는다. 각 step은 base predictor가 필요하면 포함하여 **실제 provider 최대10회**, 두 경로 합계 **15 accepted step/150 provider 이하**다. 실패 trial/line-search/final residual 확인도 모두 센다. 재사용할 수 없는 과거 RHS를 재사용했다고 기록하지 않는다. 기존 .025 단일 진단 RHS 승인은 이미 소진됐으며 이번 후보 계산과 섞지 않는다.

두 후보를 통틀어 계산 wall time 합계 24시간을 상한으로 둔다(동시에 실행하면 job별 경과시간을 합산; queue 대기는 제외). 처음 세 step의 실제 평균 비용으로 잔여 예산을 예측하고, 상한을 넘을 것으로 보이면 새 step 진입 전에 종료·자료 보존한다. provider 한도나 wall 한도를 맞추려고 수렴을 생략하지 않는다. step 도중 hard failure 외에는 검증된 checkpoint까지 보존하고 미완료를 완료로 표시하지 않는다.

모든 새 step에서 actual provider gen/z, cell Vdot, base/provider/accepted area, h, geometry hash/topology/prior 변화, raw/factor/Psi signed budgets를 작은 cell journal로 직접 저장한다. δV_map=Vaccepted−Vbase−h Vdot_provider와 그 누적값, 고정질량 density 변화, instantaneous density-trace mass/mesh balance, thermal/affine/R 기여를 같은 ID로 비교한다. endpoint difference/h를 Vdot라고 정의하지 않는다. 이 자료는 후보 진화의 부수 출력이며 별도의 diagnostic RHS를 요구하지 않는다.

.05 후보 채택에는 기존 모든 지역의 절대·시간오차 관문이 필요하다. 추가로 공통 initial-cavity ID에서 e_rho(.05)−e_rho(.025)의 질량가중 RMS가 기존 explicit fine의 해당 증분 RMS보다 최소2배 작아야 이번 후보를 density 문제의 유효 개선으로 판단한다. 이는 결과를 보고 고르는 tolerance가 아니라 현재 가설의 사전 개선 기준이다. P/velocity/phase/ell/R 및 접촉 signed pattern의 다른 지역 악화를 함께 보고하며 절대오차가 초기부터 상속된 부분을 구별한다. 밀도 개선 없이 contact P만 좋아지면 현재 장시간 blocker를 해결한 것으로 채택하지 않는다. 기준 실패는 후보 기각 또는 증거 부족으로 보고하고 자동 추가 dt sweep은 하지 않는다.

이번 성공도 t1/t6–10 연장을 자동 승인하지 않는다. 다음 한 번의 완료 평가에서 비용당 오차 개선, cavity drift와 δV 누적의 관계, 남은 spatial residual을 함께 판정한다. 그 뒤에만 .05 이후 유한 길이의 장시간 경로를 정한다. 목표를 .05 검증으로 축소하지 않으며 정확하지 않은 18일 계산을 먼저 시작하지 않는다.

### 5. 동일 묶음의 성능 작업: face geometry/transport-root native 경로

현재 stage261–273초 중 compiled 약5.2초만 보고 나머지를 Python 인터프리터 시간이라고 부르면 안 된다. outer4회의 face-quadrature 포함 geometry 약96초, face-state 약40초, massroot56–67초가 이번 구현 우선순위다. q24 약9160만 node와 cell당 median4608, nonflat96.4%도 실제 비용이며 GPU 개수만 늘려 해결될 근거가 없다. worker별 density wall 합은 job wall과 구별한다.

`Docs/curved_face_nodes_fast.py`의 straight-face speed polynomial/cut union/ordered quadrature fill과 `Docs/closed_root_transport_arrays.py`의 transport array 및 root residual/Jacobian 구성에 대응하는 native batch 경로를 구현한다. 이미 native인 star/flux를 다시 이식하는 작업으로 대체하지 않는다. 기존 within-stage line/prior cache를 그대로 확장하고 각 outer에서 z에 의존하는 항만 갱신한다. cache key는 전체 physical-stage input/geometry/density/velocity hash 및 epoch를 포함하며 다른 nonlinear midpoint 사이에 값이 같다고 가정하지 않는다. canonical face/node/incidence 순서와 accumulation 순서를 유지한다. sparse solve, global gauge, stopping tolerances는 우선 기존 것을 그대로 둔다.

true arc의 인증된 speed root isolation과 미세면/support 표현은 유지한다. straight-face 이식은 기존 polynomial validation과 exact/prior cut 계약을 보존해야 한다. q24/48 차수, density feature/positivity/rank, tiny contrast 판정, wall8/16, exactstar의 root 정책을 바꾸지 않는다. production와 audit timer를 분리하고 반복 q48/전체 face dump를 불필요하게 production inner loop에 넣지 않되 기존 checkpoint의 정확한 q48 감시는 유지한다. 대형 GN/face 배열은 기존 공유메모리 및 chunked native 소비로 다루고 입력/출력 lifetime을 명시한다. 물리 cache와 shape-only workspace를 구분한다.

검증은 저장 initial 및 moving의 두 대표 입력에서 component oracle를 수행하고, 하나의 저장 moving physical state에서 기존/새 provider를 각각 최대1회만 재생한다(합계2회; 위 후보150회와 별도). 기존 재생 가능한 baseline 전체 출력과 실행 identity가 완전하면 baseline 새 실행을 생략할 수 있다. 원래 순서/연산을 보존한 경로는 exact-value 비교, 불가피한 산술 차이는 기존 성분별 roundoff bound와 기존 seed 관문을 그대로 사용한다. 같은 physical output인데 bounds를 넓혀 속도를 채택하지 않는다. root outer 수와 root 입력, z/hydro/gravity/Psi/K/incoming 및 conservation을 함께 검사한다. 최초 실패는 immutable 보존하며 oracle 또는 방법 의미를 바꿔 통과시키지 않는다.

성능 재생의 목표는 같은 입력·같은 host에서 face geometry+transport-root 합산 wall을 최소2배 줄이는 것이다. 이것은 물리 acceptance tolerance가 아니라 작업 목표다. whole provider와 전체 step의 실제 개선, peak RSS, CPU/GPU 시간 및 데이터 이동도 함께 보고한다. 2배에 미달해도 동등성 결과는 보존하되 성공으로 과장하지 않는다. 물리 후보는 위 24시간/150provider 예산에 들어갈 때만 진행하며 native 이식의 예상 속도를 실제 측정 대신 사용하지 않는다. 이번 묶음에서는 density quadrature 축소, 새로운 spatial closure, angular ensemble 변경, 추가 MPI 확장 또는 장시간 dt 확대를 동시에 하지 않는다.

### 6. 한 번의 후속 완료 보고

운전자는 성능 동등성/실측 비용, 후보 nonlinear 잔차 및 실패 시 witness, 공통 .03/.05 물리·시간오차, direct volume journal과 signed pressure/energy budget을 하나의 완료 자료로 제출한다. 수치 관문 실패·예산 초과 시 해당 분기를 종료하고 다른 수정 후보를 자율 추가하지 않는다. 일반 구현/빌드/검증 세부는 root가 소유하며 별도 반복 평가를 요구하지 않는다. Re∞, 고정mass, 무조향 및 공통 상대W의 Voronoi 환원은 전 과정에서 유지한다.
