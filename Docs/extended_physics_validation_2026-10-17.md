# 확장 물리 검증 — 2026-10-17 진행 기록

사용자 요청: 충격파, 장시간, 3차원, 케플러 검증.
실제 repo `/tmp/lageunha-push`, 모든 계산은 lageunha foreground SSH에서 수행.
무조향, 원래 셀 질량 고정, 가변 Laguerre 가중치 유지. 기존 승인 경로 보존.
현재 네 가지 전체 유체 검증이 모두 통과한 상태는 아니다.

## 현재 판정 — 2026-10-17 14:34 갱신

| 범위 | 완료에 필요한 범위 | 현재 판정 |
|---|---|---|
| 충격파 | T=.15, 공간 비교, 물리 L1<.05 및 수치 경로 검사 | 미완료. 새 피팅16² T=.02 통과; 전체 구간 계산 중 |
| 장시간 와류 | 10중심주기(T=10.4685979156) | 미완료. 기존 중간 궤적의 밀도 기준 초과 확인; 전체기간 예산 대조 진행 |
| 3차원 | 비선형 물리장, 독립 전체 경로, 시간·공간 수렴 | 미완료. 6³ coarse 압력 기준 미달; 8³/10³ 및 시간 정밀화 진행 |
| 케플러 | 공식192² 유체 원반100궤도 | 미완료. 초기 적분 완료; 셀별·전체 이웃 재구성 감사 중; 공식 유체 진화 미시작 |

아래 내용은 순차 연구 기록이다. 초기 시점의 구현/실행 상태를 현재 완료 상태로
해석하지 않는다. 최신 체크포인트 캡처와 제한은 JSON의
`current_completion_audit`에 있으며 이전 판정도 보존한다.

## 충격파

현재 2D 모멘트 해법에 Sod adapter를 추가했다. 이전 `stagSod.py`는 다른
Delaunay 영역 모델이므로 대신 실행해 현재 모델의 통과로 처리하지 않았다.
현 모멘트 코드와 같은 gamma=5/3, rhoL=1/rhoR=.125, PL=1/PR=.1, u=0,
단위 정사각형 x=.5에 정렬된 초기 셀이다. 초기 셀 생성 위치를 정했으며 진화 중
조향하지 않는다. 해석해의 각 파동 전선에서 다각형을 분할해 참조 적분한다.
충격 Rankine–Hugoniot, 엔트로피 증가, 전역 질량/에너지/경계압력 운동량 검사를
통과했다. 추가 후 전체 Python407개 검사가 통과했다(KT). 정확한 접촉 압력의
공유 면 보존 장부 및 공분산 모멘트 항등식도 포함한다.

사전 판정: 첫 T=.02, 이어 T=.15 및 격자 수렴. rho/u/P 장 및 셀평균 L1<.05.
이 수치는 본 Sod 시험의 rhoL=PL=1 및 속도 단위로 정규화한다.

| T=.02 | density field | velocity field | pressure field | 판정 |
|---|---:|---:|---:|---|
| 4차 압력8² | .0231404 | .0916667 | .0404139 | 속도 미달 |
| 2차 압력8² | .0185781 | .0594079 | .0291546 | 속도 미달 |
| 4차 압력16² | .0153369 | .0326060 | .0237808 | 첫 T=.02 통과 |
| 2차 압력16² | .0135297 | .0296741 | .0220005 | 첫 T=.02 통과 |

두8² 모두 수치 보존/양성 gate는 통과했다. 수치 종료 성공을 물리 통과로 세지 않는다.
4차8² 연장은 t=.0541331756974829,402승인 단계에서 HLLC 음의 중간 압력으로 실패했다.
실패는 승인 상태의 첫 RHS에서 발생했다. 오른쪽 벽 재구성 rho=.0257460,
P=.00115133, 반사 전후 법선속도 ±.1754137인 비진공 상태이며, 정확한 Euler
접촉 압력은 .000345003>0이다. KN에 입력과 checkpoint hash를 보존했다.
정확한 접촉 압력 옵션은 명시적 모델 변경으로 구현하고405개 회귀 검사를 통과했다.
실패 상태의 새 RHS는 유한하며 질량 잔차6.96e-17, 조향0(KP)이다.
KQ는 동일 승인 상태부터 T=.06까지의 별도 진단 분기다. HLLC 실패나
기존 물리 오차를 없던 것으로 처리하지 않으며, 단일 해법 수렴 시험도 아니다.
4차16²는 T=.02의 수치·물리 관문을 통과해 동일 상태에서 T=.15로 연장 중이다.
2차 모델은 동일 원래 IC 배열로 별도 비교한 것이며4차 진화를 재시작한 것이 아니다.
초기 압력 L1은 N8/16/32에서 P2=.0231211/.0115898/.00580217,
P4=.0425142/.0146901/.00732844. 불연속에서도 고차가 항상 우수하지 않다.

## 장시간

완료된4차 vortex32 JB의 T=.05에서 JU로 승인 배열/질량/누적량/카운터를
그대로 이어 첫 연장 T=.2를 실행 중이다. T=.2는 아직 한 중심 궤도보다 짧으므로
장시간 통과로 부르지 않는다. 최종 의도는10중심주기이며, 짧은 구간의 정확도
관문을 먼저 확인한다. 초기 실행 환경의 mass-root backend 누락으로 첫 시도는
계산 전에 중단됐고 로그를 보존했다. 올바른 환경 설정 후 같은 체크포인트에서
이어갔다. 승인된 시간 구간을 반복하지 않았다.

## 3차원

이번 확장 작업을 시작할 당시 Python 모멘트 유체 진화는2D 전용이었다. 기존 C3D 기하부터 검사해 두 결함을 발견했다.
1. 단위 정육면체 부피4, 길이2 정육면체 부피16: 면적을 부피처럼 반환하고 있었다.
2. 여섯 축방향 이웃으로 셀을 만드는 절단 코드가 무한 순회했다(15초 timeout).

`Voro/Laguerre/laguerre.c`의 CPU 부피 및 절단을 수정했다. 부피는 각 면에서
내부점까지의 사면체 적분, 절단은 모든 교차변을 모은 뒤 연결관계를 재구성한다.
기존 `Voro/voro.c`의 EPS guard나 CUDA 코드는 수정하지 않았다.

- 크기변환/이동/비직교변형 정육면체, 사면체, 원점–면 pyramid 통과.
- ASan/UBSan 기본 기하와 절단 검사 통과.
- 무작위 가중 셀64건을 SciPy 반공간 교집합 및 convex hull과 독립 비교 통과.
- 가중치 gauge 이동과 이웃 순서 바꾸기 통과.
- 27/64셀 전체 상자 분할 및 생성점 이동 후 고정 부피 역문제 통과.
  최대 상대 셀질량 오차1.25e-12, 전체 부피 합 오차3.34e-16,
  독립 셀 부피 차이1.1e-16 이내. 공유 면적 상호 일치4.2e-16 이내.
  가중치 Jacobian은 독립 중앙차분과6.1e-10 이내로 일치했다.

추가로 생성점이 셀 밖에 있는 경우, 빈 셀, 저장 용량 부족, 중복 생성점 계약을
ASan/UBSan으로 검사했다(KK). HLLC fan의 차원 일반화 후3D 회전·갈릴레이
변환·이동 접촉·균일 상태 검사를 통과했다(KM). 기존2D1000개 입력은 수정 전과
비트 단위로 일치했다.

이는 기하/용량 및 Riemann 구성요소 검증이다.35개 원시 모멘트·양의 밀도 복원·
시간적분의 고차 결합은 아직 미완료다. 아래 LB–LD에서 초기3D 결합을 구현했으며,
전체3D 유체 정확도 통과로 판정하지 않았다.
현재 C vertex API는 꼭짓점에3개 면만 저장한다. 기존 꼭짓점을 지나는 비단순
절단은-1로 명시 거절하며, 일반 동적3D 셀의 모든 퇴화를 처리한다고 주장하지 않는다.

## 케플러

공식192² 비진공 원반 정의는 변경하지 않았다. 기존 C`Exam/exam.c`에서 실제
중력 함수를 추출해100궤도 시험입자 대조를 수행했다. 이 결과는 유체 원반과 구분한다.
Spline rc=.5, 중심0, r=1의 원궤도, 각 궤도당512/1024/2048단계다.

2048단계/궤도,100궤도:
- midpoint RK2: 최대 상대E오차2.27e-6, L오차1.13e-6, 위상오차-.004518rad.
- RK4: E4.74e-12, L2.37e-12, 위상3.53e-9rad.
- KDK: E2.22e-11, L2.58e-14, 위상-.001971rad.

원궤도에서의 관측값이며 일반 비원궤도 차수 또는 유체 성능을 보장하지 않는다.
확장 시작 시2D 모멘트 진화에는 중력이 연결되지 않았다. 아래 LH/LJ에서
중력 kick과 유체 결합을 추가했으며 공식 동결 원반은 아직 미구현이다. 평균 운동량에 중력만 더하면
D/Pi/E와 양의 공분산·열에너지의 일관성이 깨질 수 있다. 이를 위한
`moment_gravity_kick.py`는 양의 질량 quadrature의 속도 표본을 갱신하고 QR로
잔여 공분산을 구성하는 구성요소다.2D/3D 균일 가속도 모멘트, 중심력의
각운동량, 내부에너지 불변, 비선형 힘의 양의 잔여 공분산, 역방향 kick을 검사한다.
유체 driver에는 아직 연결하지 않았다. drift/면 플럭스/중력 퍼텐셜 에너지 장부를
결합해야 하며,100궤도 공식 유체 원반은 아직 검증되지 않았다.

정확한 실행 경로와 진행 상태는 `extended_physics_validation_2026-10-17.json` 참조.

## 독립 충격파 재적분

8² 두 모델의 저장 끝점에서 참조 quadrature를8차 높여 재적분했다(KJ).
장 L1 차이는4차 압력에서 최대1.47e-6,2차에서5.54e-8이었다.
두 모델 모두 전체 에너지 상대 변화1.35e-16이지만, 정확한 경계 충격량
.9t 대비 x운동량 오차는 각각 .00324845와 .0000614197이었다.
보존 장부의 일치만으로 실제 경계 힘의 정확도를 보장할 수 없다.

## 08:35 후속 판정

- KQ 정확한 접촉 압력 전환 분기는 T=.06 도달, 수치 관문 통과이나 물리 관문 미달.
  밀도/속도/압력 장 L1=.0453388/.0742135/.0780827. Riemann 근사 중단을
  해소해도 거친 격자의 물리 오차가 해소되지는 않는다. 이 진단 분기는 완료했다.
- P2 N16은 T=.02의 모든 관문을 통과해 원래 승인 상태부터 T=.15로 연장(KU).
  P2 N32는 별도 해상도 수렴 계산(KV). 두 작업은 서로 독립적이다.
- N16 P4 독립 고차 참조 적분(KR): 장 L1 차이 최대1.09e-7,
  전체 에너지 상대 변화4.04e-16, 정확한 경계 충격량 대비 오차2.24e-6.
- JU의 저장 T=.0633816338,335단계에서 별도 중간 평가(KS): 정규화 밀도/속도/
  압력 장 L1=.000475276/.000344012/.000320889, 원래 질량 배열 완전 일치.
  참조 적분 오차1.72e-11. 이는 중간 평가이며 장시간 완료 판정이 아니다.

현재 우선순위는 충격파의 연장 및 해상도 수렴, 소용돌이 장시간 오차의 확인이다.
3D 완전 유체 검증에는3D 양의 밀도 복원·모멘트 수송·시간적분 결합이 필요하다.
케플러 유체 검증에는 중력 모멘트·퍼텐셜 에너지·면 플럭스 결합 및 공식 원반
경계 구현이 필요하다. 해당 구현이 없는 상태를 검증 통과로 기록하지 않는다.

## 최종 회귀 검사와 보존

KT407개 검사는 정순 및 역순 모두 통과했고 건너뛴 검사는 없다.
P2 N16 독립 적분(KW)의 장 L1 변화는 최대2.14e-7, 전체 에너지 상대 변화
4.04e-16, 정확한 경계 충격량 대비 오차2.53e-8이다.
전체 물리 검증은 미완료이며 JU/KL/KU/KV가 계산 중이다. 각 실행의2시간 제한에
도달하면 승인 checkpoint와 pending cache를 보존한 상태에서 종료한다.
이 제한은 물리 목표 달성으로 간주하지 않는다. 기존 승인 구간은 재계산하지 않는다.

소스·결과 검토 보관본: `/gpfs/kjhan/LagForce/Codex_KX_extended_review_20261017_083534`.

## 후속 개발: 3D 밀도와 케플러 퍼텐셜

- `moment_density_3d.py`: 기존2D와 같은 양의 지수형 이차 밀도를3D의
  질량·3개1차·6개 독립2차 모멘트에 맞춘다. ConvexHull 면–내부점 사면체의
  양의 Duffy quadrature를 사용한다. 불가능한 모멘트는 clipping 없이 거절한다.
  균일/비균일 밀도, 아핀 변형 셀, 독립 고차 적분, 불가능한 중심/공분산을
  포함해 전체411개 검사가 통과했다(KY). 아직3D 유체 driver에 연결하지 않았다.
- `kepler_spline_gravity.py`: 공식 C spline 질량식과 일관된 퍼텐셜을 도입했다.
  중심은 유한하며 softening 밖에서는 정확한 -GM/r이다. 힘의 독립 적분 및
  퍼텐셜 중앙차분을 통과했다. 실제 C 함수와153개 위치에서 최대4.44e-15로
  일치했다(KZ). 전체415개 회귀 검사 통과.

남은 결합 순서:3D 밀도→35개 원시 모멘트 해석/수송→양의 공분산 적분→
가중치 경로·고정질량 제약→3D Euler 수렴 시험. 케플러는 중력 kick과
유체 수송의 결합, 면을 통한 퍼텐셜 에너지 수송 장부, 공식 원반 IC·경계,
시간/공간 수렴 및100궤도 검증까지 필요하다. 구성요소 통과로 대체하지 않는다.

- `moment_state_3d.py`:35개 원시 모멘트(m,q,Q,p,D,Pi,E)의 배열 계약,
  양의 밀도를 사용한 아핀 속도·잔여 공분산·열압력 해석, 동일 면에서 양쪽에
  반대 부호로 쓰는 운동학적 모멘트 플럭스를 구현했다. LA에서 아핀 회복,
  균일 Euler 플럭스, 면 반전 시35개 성분의 상쇄를 통과했다.
  전체418개 검사 통과, 건너뜀0(LA). 기하·시간적분 결합은 아직 미완료다.

## 3D 실제 기하–모멘트–시간적분 결합

LB는 실제 C Laguerre 셀과3D 양의 밀도,35개 원시 모멘트, 공유 운동학적
플럭스, 순간 고정질량 가중치 방정식을 연결했다. 균일 상태의 정확한
Reynolds RHS와 최대7.37e-15로 일치한다. 비균일 상태에서 보존 장부 오차
6.34e-16, 공분산 방정식 항등식 오차2.30e-17, 순간 질량 잔차9.05e-19,
생성점 조향0을 확인했다.

LC는 같은 연산자에3D 양의 Sylvester midpoint를 연결했다. 기존 factor
함수의 두 군데2D 단위/영행렬을 실제 인자 차원에 맞게 일반화했다.
3D rank-deficient loss와 covariance injection 검사를 포함해419개 검사 통과.
균일 병진 유체를 실제 움직이는 셀에서 T=.1까지 dt=.05/.025/.0125로
진화시켜 원시 모멘트 오차5.28e-16 이하, 원래 질량 배열 완전 일치를 확인했다.
이는 균일해 보존 검사이며 비균일 물리 정확도나 연속 경로 질량 제약의 증명이 아니다.

LD는3D 셀평균 보존 선형 압력 복원을 추가했다. 아핀 압력의 정확한 Euler
모멘트 가속도를 먼저 확인하고, gamma5/3·진폭1e-5·파수2pi의 반사 경계
정상 음파를27셀에서 T=.02까지 계산한다. 참조는 선형 음향해로 비선형
오차는 O(A²)임을 명시한다. rho/P/u L1을 각각 A/gamma*A/c*A로 정규화해
기준 .05와 비교한다. 현재2D의 고차 압력·속도 복원과는 별도 명시적 모델이며,
3D 고차 확장·공간/시간 수렴·조밀 경로 감사는 여전히 미완료다.

### LD 음파 첫 결과와 후속 수렴

LD27셀(dt=.002,T=.02)은 끝까지 진화했으나 물리 오차 기준에 미달했다.
정규화 rho/u/P L1=.0164089/.0253946/.1419404이며 초기 압력 오차도 .1373332다.
단계 질량 잔차 최대1.75e-13, 원래 질량 배열 일치, 에너지 장부 오차2.23e-16.
기하·밀도·양의 모멘트·선형 압력·실제 시간적분이 결합된 첫 비균일3D 시험이며,
성공 종료를 물리 통과로 처리하지 않는다. 사전 기준 .05는 유지했다.
아핀 압력의35개 모멘트 RHS 정확도는6.15e-16로 통과했다.

LE는27셀 dt=.001, LF는64셀 dt=.002, LG는216셀 dt=.002의 독립 수렴 시험이다.
세 계산 모두 T=.02, 동일 진폭/방정식/오차 기준을 사용한다. 원래 LD 결과는 보존했다.

## 09:12 판정 갱신

- 3D 선형 압력 수렴: N3(dt=.001)/N4/N6의 최종 압력 L1은
  .141924/.099575/.061160으로 모두 .05 기준 미달. N3 dt반감 변화가 작고
  초기 압력 오차가 지배적이다. 에너지 장부 잔차는4.45e-16 이하.
- 이차 압력의 전역 다항식 회복, 셀평균 보존, Bernstein 양성 제한 및 정확한
  Euler 모멘트 가속도 검사는 통과(LI). 그러나 N3 실제 음파는 압력 .177680,
  속도 .061524로 실패했다. 이차 복원을 무조건 개선으로 채택하지 않는다.
- 동일 초기 배열의 정적 비교(LK): N6 선형/이차 압력오차 .060191/.059693,
  N8 .032529/.025163. 이에 따라 LL(N6,.002), LM(N8,.002), LN(N8,.001)의
  이차 압력 진화를 별도 수렴 시험으로 시작했다. 기준은 그대로 .05다.
- LH는 현재2D 아핀/이차 속도 모멘트 상태에 중력 kick을 연결해422개 검사 통과.
  LJ는 유체 midpoint와 Strang 결합 후 균일 자유낙하의 움직이는 셀·위치/운동량·
  Egas+Epotential 검사까지 포함해423개 검사 통과. 공식 케플러 원반과 구분한다.
- 공식192² 고리 IC는 실제 binary 36864점, [0,6]²임을 확인하고 hash를 저장했다.
  동결 경계/초기 불연속 적분/에너지 장부의 남은 구현은
  `kepler_moment_coupling_plan_2026-10-17.md`에 명시했다.

## 동결 경계와 공식 원반 초기조건

LO427개 검사 통과. `frozen_weight_root.py`는 동결 weight-rate=0을 두고
활성 셀 질량 방정식만 sparse solve한다. 동결 셀의 생성점·모멘트·factor는
유체/중력 전체 단계에서 정확히 유지되며, 버린 반대쪽 P/E/L 변화는 외부
저장소 교환 장부로 기록한다. 활성 생성점 조향0, masked 보존 장부,
전체 동결, 잘못된 audit-mode 조합 거절을 검사했다.

LP431개 검사 통과. `kepler_disk_reference.py`는 원형 밀도 경계 및 spline
분기(.25,.5,2)에서 양의 polar quadrature를 분할한다. 전체 상자의 정확한
원반 질량, 대칭, 원운동–중력 평형, 경계 셀16/24차 독립 적분을 검사했다.
고리 생성기는 공식192² binary와36864개 위치가 비트 단위로 같다.
작은 명목 n12/16/24/32에서는 원본 규칙상 점수가 n²보다 적다. 공식 배치는
그대로 두며, 보조 시험은 실제 미패딩 점수를 표기한다(n16은246셀).

LQ의246셀 초기화는 시간 진화 전에 zero polygon edge로 실패했다. 중복된
동일 꼭짓점만 제거하도록 polar_rule을 수정했다. 모양이나 면적을 변경하는
허용오차 병합은 하지 않는다. 해당 회귀 검사 후 LR에서 첫 단계를 실행하도록
연결했으며, 현재 회귀 검사 중이다. LQ 실패 로그는 보존했다.

## 09:32 케플러 첫 RHS 실패 진단

LR432개 검사가 통과했고246셀 IC 생성·밀도 복원도 끝났다. 독립16/24차
초기 원시 모멘트 상대 차이2.57e-15, 초기 압력 상대 오차1.67e-9,
밀도 맞춤 잔차2.92e-13, 동결 셀114개다. 그러나 첫 유체 RHS는 진공 조건으로
실패해 승인된 유체 시간단계는 아직 없다. 생성된 초기 상태는 보존했다.

LS는 동일 저장 IC를 사용해 중력 전/반 kick 후를 따로 확인했다. 활성 셀0–7,
r=.56633인 면에서 중력 전 법선속도 차이 .03559209 > 진공 문턱 .03543482.
반 kick 후도 같은 면이다. 좌측 rho=1.8463, 우측 rho=1이며 실제 기준장의
속도는 양쪽에서 연속이다. 따라서 시간단계 감소만으로 제거할 수 있는 문제가
아니다. 거친 공간 재구성 오차가 만드는 국소 Riemann 진공이며, 현재 정확한
압력 solver의 진공 처리도 미구현이다. 실제 vacuum 지원이 없는 한 일반적인
진공 Euler 검증을 통과했다고 주장할 수 없다.

해상도 영향을 확인하기 위해 동일 물리·고정질량·무조향·동결 경계의 명목
n32/n64 미패딩 고리 제어(LT/LU)를 독립적으로 시작했다. 공식192² 초기 파일과
100궤도 목표는 바꾸지 않는다.246셀 실패나 진공 조건을 floor로 숨기지 않았다.

## 09:47 보조 첫 단계 통과 및 새로운 실패 진단

- KU P2 Sod16은 정상 저장 시각 t=.04841389062313528에서 HLLC 별 압력으로 중단됐다. LZ 실패 재현에서 정확한 별 압력은 .07231228>0이고 진공이 아니다. 마지막 상태의 rho/v/P L1은 .0229334/.0460362/.0210308로 기존 평균 기준 안이지만, 최대 셀 평균 횡속도는 .345784이다. MB에서 T=.02에도 P2 .06197, P4 .09262의 횡속도가 확인됐다. 평면 Sod의 대칭 오차를 별도로 추적해야 하며 장시간 충격파 검증 완료가 아니다.
- LV의 n32 압력 오차 .374112는 극좌표 적분이 셀 면적 .0395691을 .0247658로 잘못 계산한 것이었다. 원점의 반올림 크기 inward vertex loop가 가짜 halfspace를 만들었다. convex power cell의 hull 면적과 원래 면적이 반올림 오차 범위에서 일치할 때만 정리하고, 별도 polygon 면적 검사를 추가했다. 넓은 좌표에서 shoelace 상쇄를 피하도록 면적은 셀 내부 기준점 좌표로 계산한다. 진짜 비볼록 셀은 거절한다.
- n64는 동결 셀–벽 면에서만 strict Euler vacuum이 발생했다. 두 쪽 모두 진화하지 않는 면의 압력 traction을 계산하지 않도록 수정했다. 활성–동결 면의 물리 solver는 유지한다. LW는 기존 n64 초기 상태를 그대로 사용하여 dt=1e-5 첫 단계를 통과했다. 질량 유량 잔차6.52e-14, 조향0, 동결 상태와 factor 정확 일치.
- LZ n32는 잘못된 기존 IC를 보존한 별도 수정 초기화 후 첫 단계를 통과했다.1009셀, 압력 상대 오차3.04e-9, mass residual5.23e-17, 조향0. Egas+재구성 potential 변화−6.54e-9는 진단값이며 완전한 potential 경계 장부 보존 판정이 아니다.
- LZ 회귀436개 정순 통과. 역순 검사 진행 중. LW/LX 첫 회귀 실행은 snapshot에 p2_geometry_omp.so를 누락하여 각각1개 오류가 났다. 해당 로그를 보존하고, 바이너리를 포함한 LY/LZ에서 전체436개 정순 통과를 확인했다.
- 3D acoustic 재개 기능은 원래 질량/설정 일치와 유한값을 검사한다. LY에서 중단 없는 n3 두 단계와 저장 후 재개 결과의 배열이 비트 단위로 같았다. 과거 체크포인트에 없던 stage mass peak는 실제 측정치 대신2e-13 상한으로 명시한다. 실행 중인 LM/LN 소스는 수정하지 않았다.
- MA는 저장된 n32/n64 초기/첫 단계 상태의 실제 Eulerian 오차를 독립16/24차 적분으로 조사한다. n32 evolved polygon의 면적 검사에서 멈춰 MC에서 적분 해상도 진단 중이다. 오류 허용치를 완화하지 않았다.

전체 요청(충격파 T=.15, 장시간 와류, 3D 비균일 물리 정확도와 경로 질량, 공식192² 케플러100궤도)은 아직 미완료다. 초기/첫 단계 통과를 전체 물리 검증으로 해석하지 않는다.

### 09:48 확정 결과와 바로 다음 작업

LZ436개 역순 검사도 통과했다. MA의 n32/n64 field audit은 모두 evolved
near-origin polygon의 독립 면적 검사로 실패했다. MC에서 n32의6개 중심 셀은
16,24,32,48,64,128,256차 모두 실패한다. 원점에서 약1e-8 떨어진 짧은
변이 생긴 상태이며, 초기 동일 꼭짓점 문제와 구별해야 한다. 감사 적분기를
개선하고 기존 저장 상태를 재평가해야 한다. 시뮬레이션을 처음부터 반복할
이유는 없고, 단순 차수 증가나 허용치 완화로 통과 처리하지 않는다.

다음 작업: (1) near-origin polar 적분의 독립 면적/모멘트 검사를 해결,
(2) Sod 벽 근처 횡속도 성장의 경계/재구성 기원 진단, (3) LM/LN 종료 결과
평가 및 시간 제한 때 명시적 resume, (4) JU/KL/KV의 원래 정상 prefix 유지.
JU timeout은 약09:56, LM/LN은10:09이다. 실행 중 프로세스에 소스를 덮어쓰지
않는다. 수정된3D resume driver는 LZ snapshot에 있고, LY resume_check에서
비트 일치가 검증됐다. 장시간 실행용 source 변경은 프로세스 종료 확인 후
별도 재개 snapshot에서 원래 checkpoint/IC/ledger를 보존해야 한다.

### 09:48 3D 비균일 음향 첫 물리 게이트 통과

LM n8(512셀), quadratic pressure, dt=.002가 T=.02에 도달했다.
최종 정규화 L1 rho=.00630597, v=.01099084, P=.02667823으로 모두 .05
이하다. 원래 셀 질량 정확 유지, stage mass flux1.63e-13, 에너지 장부
오차4.44e-16. 이 결과는 기존 n3/n4/n6 실패를 대체하거나 지우지 않는다.
LN의 dt=.001 시간간격 비교는 현재 t=.010으로 계속 실행 중이다.
3D 전체 검증에는 시간 수렴 비교, 연속 경로 질량 검사 및 비선형 시험이 남았다.

## 09:54 기준해 적분 및 반사벽 압력 개선

MD에서 원점 근처 면의 polar r=b/(n·direction) 함수가 급변하는 각도 구간을
기하급수적으로 분할했다. 면적 허용치를 변경하지 않았고, 기존 실패 중심6셀의
16/24차 적분 모두 독립 Cartesian 면적·1차·2차 모멘트와5.25e-15 이내로
일치한다. MF437개 정순/역순 회귀 통과. MF n32/n64 field audit은 원래 저장된
IC와 첫 단계 상태에 대해 진행 중이며 시뮬레이션 prefix는 반복하지 않았다.

ME는 Sod 원래 초기 상태부터 경계 셀의 가짜 횡방향 힘을 재현했다. P2/P4의
최대 y가속도6.1288/9.1975, 내부 셀 y가속도는 반올림 수준이다. P2/P4 압력
재구성에 명시적 y반사벽 ghost-cell 평균/기하 확장 옵션을 추가했다. 셀 평균
압력·에너지·질량·생성점은 변경하지 않는다. MG에서 가속도는3.11e-15/1.60e-14로
줄었고, 동일 x열의 재구성이 모든 y행에서 일치하는 회귀 검사도 통과했다.

MH439개 정순 회귀 통과 후 n16/n32 수정 경계 모델 비교를 독립 실행했다.
원래 t0 초기 배열을 SHA로 확인하여 재사용하며, 잘못된 기존 경계 모델의
진화 상태를 되감거나 횡속도를 강제로 지우지 않는다. 기존 실패와 진행 중인
원래 모델 결과는 별도로 보존한다. 두 새 설정은 pressure_reflection_y=
"0,1", T=.02, 기존 오차예산 그대로다. 수정 모델의 장시간 성공은 아직
미검증이다. 설정/환경 및 저장 RHS의 reflection 설정 일치도 검사한다.

활성 실행: JU 와류(약09:56 timeout), KL P4 Sod16(10:19), KV 원래 P2 Sod32
(10:33), LN 3D N8 dt=.001(10:09), MF 두 field audit(약10:10 timeout),
MH n16/n32(약11:53 timeout). JU/LN timeout은 실패한 물리 판정과 구별하고
실제 종료 확인 후 승인된 checkpoint 및 pending cache에서 재개할 것.

### 09:55 반사벽 수정의 한계 확인

MH439개 역순 회귀도 통과했다. 하지만 MI에서 초기 셀 평균 횡방향 가속도가
반올림 수준이어도 셀 내부 y압력 곡률은 P2 .007563/P4 .010295, Dyy 변화율은
3.0774e-5/4.1076e-5로 남는다. 반사벽만으로 완전한1D 보존이 해결됐다고
주장하면 안 된다. MH 계산은 이 명시적 부분 수정의 비교 증거로 유지한다.

다음 재구성 진단 후보는 tensor-product 공간 가중치와 자유 절편의 FV fit 후
원래 셀 평균으로 정확히 재중심화하는 방식이다. 정렬된 Cartesian stencil의
1D 데이터가 y기저에 잘못 투영되는지를 먼저 정적/첫 RHS 시험으로 확인해야
한다. 현재 ridge/셀평균 고정/radial stencil 조합이 가짜 y² 항을 만들 수 있다.
이는 아직 구현하거나 검증한 결과가 아니라 다음 실험 가설이다.

## 10:02 tensor pressure 초기 게이트와 와류 재개

MJ는 tensor-product compact weights, target cell 포함, 자유 절편 FV fit 후
원래 열에너지 cell mean으로 재중심화하는 명시적 대안을 검사했다. 반사벽과
결합하면 P2/P4 초기 y압력 곡률<=2.22e-16, Dyy 변화율<=3.73e-17, 셀 평균
y가속도<=3.25e-14다. 원래 radial/fixed-mean 구현은 기본으로 남는다. 설정과
저장 RHS는 fit mode를 확인한다. MK443개 정순/역순 검사 통과(1D 횡방향
불변성, quartic 재현, mean/양성 보존, 설정 거절 포함). MK/sod16,sod32는
각각 원래 t0배열을 사용하는 독립 모델 비교로10:00 시작했다. 장시간 성능은
아직 미검증이며 MH/원래 모델의 결과를 보존한다.

JU는09:56:37 timeout exit124로 종료됐다. PID1538678 부재를 확인하고
t=.09939543512386753,582steps checkpoint, pending.json/source_rhs.json의
SHA를 resume_0957_receipt.json에 기록했다. 같은 폴더·같은 코드·같은 설정으로
새 execution_resume_0957.txt/run_resume_0957.log를 사용해 재개했다.
583–585번째 승인 단계를 확인했으며 완료 prefix를 반복하지 않았다.

MF n32 field audit 완료:24차 상대L1 rho=.0319413290, 질량가중 v=.0004219745,
P=1.0624971e-7 (t=1e-5).16차 density와의 차이는1.06e-5다. n64 감사는 계속
진행 중. MM은 n32 기존 첫-step checkpoint에서 dt1e-5로 T=.001까지 이어
계산한다. 매 승인 단계 저장, frozen 상태 정확 일치와 stage mass/무조향을
검사한다. dense mass path와 완전한 potential 경계 장부는 미검증으로 표기한다.

ML은 convex polygon의 vertex angle로 나뉜 빈 각도 구간을 한 번 판별하여
Gauss node 처리를 생략했다. 기존 실패6개+무작위24개 polygon에서 적분점과
가중치 비트 일치;1.088초→.608초. 실행 중 MF 소스는 바꾸지 않았다.

## 10:10 3D 실제 경로 검사와 반사벽 endpoint

MN의 3D RHS에 명시적 audit_geometry_rates를 넣었다. 진화 단계는 원래대로
질량 root를 풀고, 진단은 제공된 생성점/가중치 속도를 그대로 사용한다.
교란된 속도 입력은 mass residual=.005878을 남겨 재해결하지 않는 것을
확인했다. 실제 C 기하8셀의 비균일 단계 dt=.01/.005/.0025에서 적분 상대
질량 오차8.479e-7/1.050e-7/1.304e-8, Simpson 중첩 차이9.12e-11/
1.63e-11/2.81e-12. 단계 mass residual은1e-18 수준이어도 경로에는 별도
오차가 있다. 한 단계 오차의 약3차 감소이며 전체 음향 궤적 인증은 아니다.

MP 새 단위 검사2개는 fixture의 Pi=.01과 factor .1의 제곱이1.7e-18 차이,
너무 작은 차분 간격의 상쇄 때문에 실패했다. Pi를 factor에서 구성하고
정확한2차식의 중앙차분 간격을1e-5→1e-3으로 바꿨다. 물리 코드나 오차
허용치를 완화하지 않았다. MQ445개 정순 검사 통과, 역순 확인 예정.

MH/n16 반사벽 control은 T=.02,131steps/12reject,exit0로 완료됐다.
rho/v/P L1=.01324895/.02952490/.02198305로 기존 기준 통과. MR에서
max mean vy는 원래 .0619734→.00077303, transverse kinetic6.303e-6→
7.347e-7로 감소했다. 내부 transverse RMS .010547이 남으므로 tensor-fit
모델(MK)의 결과와 비교해야 한다. MH/n32는 계속 진행 중(누락된 PID를
잘못 연결하지 말 것: 반사벽n16이 먼저 완료됐다).

MO에는 LN3D acoustic 재개용 소스만 준비했다. 핵심8개 물리 파일과
libgate.so의 SHA가 LN원본과 동일함을 확인했다. 기존 프로세스 종료 후
initial_state/initial_score/final checkpoint만 복사하고 --resume으로
진행할 것. 경로 검사 추가 전 물리 코드와 동일한 LY의 검증된 resume
driver를 사용한다.

케플러 감사기는32셀마다 원자적으로 진행량을 저장하고 입력/config/source
manifest hash를 확인하도록 개선했다. 반경2.5 원과 만나지 않는 bounding
box는 적분에서 건너뛴다. MS는 기존 N16 IC에서 중간 저장 후 강제 중단과
재개 결과가 원래 감사값과 같은지 검사 중이다. 실행 중 MF는 변경하지 않았다.

### 10:11 시간 제한 작업의 재개와 감사 분할

LN은10:09:42 exit124, 최종 승인 t=.015000000000000006으로 종료됐고 PID
1631545 부재를 확인했다. MO로 원래 initial_state/initial_score/checkpoint를
복사하고 SHA receipt를 저장한 후 --resume을 실행했다. 원래5개 남은 단계
진행이며 기존 .015 구간을 반복하지 않는다. physics8파일/라이브러리는 원본과
비트 동일하다. 최대 stage residual은 구 checkpoint에 없으므로2e-13 상한으로
기록된다.

MF n64 감사는 tool session99530에서 exit124를 확인했다. 기존 감사기가
중간 결과를 저장하지 않아 완료된 부분 적분값을 회수할 수 없었다. MS에서
새 감사기의 강제 중단/재개 결과가 원래246셀 scorer와 정확히 같음을 확인했다.
MT는 저장된 n64 IC/첫 단계 checkpoint의 initial/evolved×16/24차 네 항목을
독립 병렬 실행한다. 각 항목60분 제한,32셀마다 진행값과 입력/설정/소스 SHA
저장. 이것은 read-only 감사 적분 재실행이며 시뮬레이션 prefix는 그대로다.

MQ445개 역순도 통과했다. 준비/진행/완료를 혼동하지 말 것: MO와 MT는 현재
진행 중이며 전체 목표는 미완료다.

## 10:16 전체 기간 예산 회계 수정과 새로운 시간제어 시험

MK/n16은 T=.02에 완료(exit0):rho/v/P L1=.01318826/.02942935/.02225530.
max mean vy=.0001739633, transverse kinetic1.31568e-10으로 원래 P2보다
각각 약356배/4.79만배 작지만 정확히0은 아니다. n32 수렴 비교는 진행 중.

중요한 장시간 검증 회계 문제를 발견했다. assess_values의 rate는
duration*max|mass_flux/mass|이다. end를 .02→.15로 바꿀 때 원래 prefix의
peak.rate도7.5배 해야 한다. MU에서 원래P2/반사벽/tensor의 재정규화값은
모두약1.9e-4로1e-4 목표를 넘는다. 기존 duration-extension들의 peak를
그대로 쓴 전체 구간 주장에는 추가 감사를 해야 한다. 또한 final numerical
flag가 누적 mass/quadrature만 검사하던 부분을 수정하여 path rate, generator,
ledger, stage mass, positive P/factor도 명시적으로 포함했다. 기존 허용치를
변경하지 않았다. MU448개 정순 회귀 통과(역순 확인 예정).

따라서 MV는 기존 tensor/reflected 공간 모델과 정확히 같은 원래 t0 배열을
쓰되 처음부터 T=.15 전체 기간의 시간제어 예산을 적용하는 독립 시험이다.
기존 T=.02 계산은 그대로 보존한다. 해당 prefix의 재정규화 rate가 목표를
넘기 때문에 그대로 이어서 전체 기간 gate를 통과했다고 주장할 수 없다.
기존 정상 진화 상태를 덮어쓰거나 되감지 않는다. MV시작10:16,2h 제한 후
필요하면 같은 설정에서 checkpoint/pending cache를 이어갈 수 있다.

### 10:19 케플러 n64 field audit 완료 및 dt 비교

MT initial/evolved×16/24 네 감사가 모두 완료됐다.24차 첫-step 상대L1은
rho=.0225747955, 질량가중v=4.9706287e-5, P=1.8997242e-8이다. n32보다
줄었지만 두 해상도의 첫 극단적으로 짧은 단계만으로 장시간 정확도나 일반적인
수렴차수를 주장하지 않는다. density16/24차 차이는1.92e-5다.

MW는 n32/n64의 기존 첫 단계에서 dt1e-4, T=.001로 진행하는 독립 시간간격
비교다. n32는 MM dt1e-5와 동일 endpoint에서 비교할 예정. 사전점검용
초소형 dt를 장시간 계산에 무조건 적용하지 않고 안정성과 오차를 측정하여
선택한다. 각 단계checkpoint, stage mass/조향/frozen 상태 검사 유지.

MU448개 역순 검사도 통과했다. 기존 KL 장시간 baseline은 곧2h 제한이며,
물리적으로 큰 횡속도와 재정규화 prefix gate 문제를 이미 확인했다. 종료 후
자동으로 성공 판정하거나 같은 baseline을 무조건 연장하지 말고, 결과를
보존하여 MV 전체-period 수정 모델과 비교하는 진단 자료로 취급한다.
KV 원래 n32 T=.02 비교, MH n32 반사벽 비교, MK n32 tensor 비교는 각각
독립 분해 실험이다. 전체 요청의 충격파 장시간 정식 gate는 MV가 담당한다.

### GitHub 저장 시점 추가 상태

- MX 회귀검사: 정순 451개(38.571초), 역순 451개(39.872초) 모두 통과.
- KL 기존 P4 Sod16 연장 시험은 10:19:26에 exit 1로 종료했다. 마지막 승인 상태는 t=.06261794533199819, 1573단계이며 half-stage HLLC의 nonpositive star pressure가 원인이다. timeout 또는 성공으로 분류하지 않는다. 기존 checkpoint와 결과를 보존한다.
- 새 3D 비선형 와류 실행기와 실제 dense record의 무손실 NPZ 저장/복구를 추가했다. MX n6/n8 진화 시험은 이 기록 시점에 아직 시작하지 않았다. 저장 경로 전체의 독립 mass audit도 남아 있다.
- 충격파 장시간, 3D 비선형, 장시간 케플러 물리 검증은 미완료다. 회귀검사 통과는 이들 물리 검증의 완료를 의미하지 않는다.

### 10:35 후속 검증

MO 3D 음향파 n8, dt=.001이 T=.02에 완료(exit0). 진폭 정규화 L1은
rho=.00629938395, v=.01099444769, P=.02667693760으로 모두 .05 미만이다.
원래 mass 배열은 정확히 유지, energy ledger error=2.22e-16이다.
이 실행의 stage mass 값2e-13은 과거 prefix를 포함한 상한이며 측정 peak가 아니다.
전체 dense path는 아직 인증하지 않았다.

MM/MW 케플러 n32 dt1e-5/1e-4와 n64 dt1e-4가 모두 T=.001 완료.
MY는 해당 endpoint를 별도로 복사하여16/24차 field audit를 진행한다.
n32 두 dt의 최대 상태 차이9.969e-11, 고정질량 배열은 동일하다.
24차 상대L1(n32 fine/coarse): density .03198107952/.03198107919,
v .000422040992/.000422040982, P1.06080049e-5/1.06072491e-5.
짧은 구간의 시간간격 민감도는 작지만 장시간 궤도 검증은 미완료다.
NA에서 기존 T=.001 상태를 이어서 T=.01까지 n32 dt1e-4/1e-3와
n64 dt1e-3 독립 비교를 시작했다. 각 단계 checkpoint와 frozen/mass/steering gate 유지.

MX n6/n8 비선형 3D 와류는10:31:37에 시작했다. n6 첫단계T=.002와
실제 dense record를 저장했다. MZ는 해당 기록의 prescribed geometry
derivative로17node mass audit를 시작했다. mass-root를 재해결하지 않는다.
감사 코드는 node별 checkpoint와 입력/소스/library hash를 기록한다.
새 Simpson 적분 검사와 기존 dense path/직렬화 검사5개 통과.
전체 물리 검증의 완료 판정은 보류한다.

### 10:36 실패 진단 및 공간 비교

NB는 KL의 승인 checkpoint, source_rhs 및 당시 dt를 복사하고 hash를 확인하여
실패한 half-stage를 복원했다. t=.06262489138, face232–248,
위치(.5605422,.9385684)에서 HLLC의 nonpositive star pressure를 재현했다.
속도차1.15146은 진공 임계5.75392 미만이고 exact star pressure=.02159340135>0.
좌/우 횡속도 .47949/1.65841로 이미 심각한 비물리 횡방향 운동이 있다.
따라서 solver 실패와 공간 재구성/진화의 횡방향 오차를 구분하며,
exact solver 하나로 해당 원래 모델이 검증됐다고 주장하지 않는다.
진단원본: /gpfs/kjhan/LagForce/Codex_NB_halfstage_sod_diagnosis_20261017.

MY n64 T=.001 field audit도 완료.24차 relative L1은
rho=.02259628025,v=4.97073218e-5,P=1.90158006e-6.
n32 값 .03198107919/.000422040982/1.06072491e-5보다 작다.
짧은 시간의 두 해상도 비교이며 공식192²100궤도의 증거가 아니다.
NA n32 dt=.001은 T=.01까지 종료(exit0), 작은 dt 및 n64는 진행 중이다.

KV 기존 n32 Sod는10:33:48 timeout124, t=.01746306122에서 중단됐다.
PID1588567 종료 확인 후 같은 디렉토리/설정의 checkpoint와 pending을 유지하여
재개했다. resume_1034_receipt.json에 원본 hash, 별도 실행 로그를 기록했다.

MX n6/n8은 각각t=.006/.002 진행. MZ 감사는 face-order4의 독립 고차
면 적분을 사용하며, 실제 진화는 order3이다. NC는 face-order를 명시적으로
설정하고 identity에도 넣어 동일 기록의 order3 검사를 시작했다. 두 결과를
비교하여 면 적분 오차와 실제 경로의 질량 오차를 구분한다. 아직 첫 기록의
감사도 진행 중이므로 전체 경로 통과를 주장하지 않는다.

### 케플러 분할 적분 경로 기록 추가

ND는 gravity_coupled_2d의 선택적 hydro_record 기록을 추가했다. 기본 수치
결과는 그대로이며 기록 유무의 배열 동일성을 검사했다. frozen mask가 양
RK 단계에서 같은 경우에만 dense extension을 허용한다. mask 전환은 명시적
실패로 남긴다. 기록된 geometry derivative를 그대로 사용하는 mass-only
검사는 frozen generator/weight rates=0을 요구하며 root를 다시 풀지 않는다.
기록은 첫 중력 kick 이후부터 둘째 kick 직전의 실제 유체 subflow이다.
중력 kick은 density moments/geometry를 바꾸지 않는다. 이 검사는 전체
퍼텐셜 에너지/경계 교환 장부를 대신하지 않는다.

ND 정순456개 회귀 통과(45.963초), 역순 진행. NE는 보존된 n32 T=.01
상태에서 dt=.001/.0005/.00025 독립 trial을 실행하며 실제 경로17node
적분과9node 비교를 저장한다. trial endpoint는 자동 승인 checkpoint로
승격하지 않는다. 각 node를 저장하여 중단 후 동일 계산을 재사용한다.

이 시점 비교 진행: JU 장시간 와류t=.1173277, MV whole-interval Sod
t=.0278209, KV 기존 n32 t=.0177286, MH 반사벽 n32 t=.0097679,
MK tensor sod32 t=.0084588. MK 하위폴더 이름은 n32가 아니라 sod32다.

### 10:41 케플러 실제 유체 경로 결과

NE 세 trial 모두 계산 완료. dt=.001/.0005/.00025에서 활성 셀의
최대 적분 absolute 상대 mass defect는 각각9.6044e-7/4.6876e-8/5.8496e-9,
최대 상대 mass rate는 .00379573/.000192494/4.79843e-5이다.
최대 generator interpolation defect는6.8545e-6/1.7096e-6/4.2689e-7.
17node와9node 적분 차이는1.0325e-9/5.1609e-10/7.9523e-16.
작은 두 dt의 누적 오차는 약8배 감소하지만 큰 dt까지 단일 수렴차수를
단정하지 않는다. quadrature/topology 민감도 추가 검사가 필요하다.
endpoint time 비교만으로 내부 mass path를 통과시킬 수 없음이 확인됐다.
이 결과는 한 hydro subflow의 진단이며 전체 과거 구간 인증/공식100궤도
통과 또는 trial endpoint의 자동 승인을 의미하지 않는다.
ND456개 역순 검사도 통과(45.103초).

### 3D 비선형 첫 완료와 적분 정밀화

MX n6 T=.02 완료(exit0), pressure normalized L1=.0555262>.05로 물리 gate
실패. density=.00453980, velocity/amplitude=.0426413이나 velocity의 참조
적분 대비 상대L1은 .259607이다. 질량 배열 정확 유지와 energy ledger0은
물리 정확도 통과를 대신하지 않는다. n8은 계속 진행한다.

첫 n6 실제 경로의 mass audit에서 면order3(NC)/order4(MZ) 적분 absolute
relative defect는8.59181e-7/1.37456e-5로 일치하지 않았다. 실제 진화의
면order3 수렴이 미검증이며 face quadrature 오차를 분리해야 한다.
NG에서 pressure traction/공분산 미분을 제외한 독립 mass-only evaluator를
구현하고 실제 기록 s=0,.5,1의 기존 full RHS mass rate와 비교했다.
최대 절대 차이4.0e-18, 질량 정규화 차이8.6e-16 이하, 각 평가2.3–2.6초.
이것은 해당 세 표본의 동등성 증거이며 모든 상태 동등성의 증명은 아니다.
NH에서 동일기록의 면order3/4/6/8/12 수렴을 검사 중이다.

NF는 NE의17개 값을 그대로 재사용하여33node 적분 완료. 세dt absolute
relative defect9.603745e-7/4.687559e-8/5.849642e-9로 주요 결론은 유지됐다.
중간dt의 worst quadrature-difference cell29는 끝점에서 상대 유량이
직전1.257e-5에서 -3.658e-5로 급변했다. 끝점 근접 추가 표본을 검사 중이며
원인을 단정하지 않는다. 모든 정밀화는 기존 계산/표본을 보존했다.

### 10:49 3D 면 적분 원인 확인

NH의 동일 초기 상태/기존 weight-rate에서 unsplit 면order3/4/6/8/12의
최대 상대 mass rate는 약0/.0070062/.0062481/.0053061/.0049680이었다.
NI는 아핀 속도 표본의 face-relative speed가0인 직선으로 실제 face
polygon을 분할하여 각 양/음 부분을 독립 양의 적분했다. 분할order3/4/6
결과는 .004875259097336571/.004875259097336238/.004875259097336047로
일치했다. 압력항은 mass equation에 없으며, 참조장을 주입하거나
weight-rate를 다시 풀지 않은 비교다.

따라서 기존 order3 mass root의0 잔차는 이 상태에서 면 적분 오차에
의존한다. 이전 stage gate만으로 연속 면 유량 보존을 주장할 수 없다.
다음은 split 적분과 그 Jacobian으로 weight root를 비교하고, 같은 적분을
다른 모멘트 flux에도 일관되게 연결하는 것이다. 현재 진화 코드는 아직
변경하지 않았고 기존 n8 실행은 baseline 증거로 보존한다.
NI의 선형 양/음 부분 해석적 적분과 경계/영속도2개 검사는 통과했다.
NG mass-only evaluator의 실제3개 표본 일치와 함께 독립 진단에 사용했다.

케플러 끝점 근접 probe에서 cell29 유량은 s=.99에서2.184e-6,
s=.9999에서-3.619e-5, s=1-1e-8에서-3.657669e-5, s=1에서
-3.657672e-5다. 정확한 끝점만의 불연속이 아니라 짧은 구간의 빠른
변화가 관측됐다. 이것의 공간 재구성/geometry 원인은 아직 미확정이다.

NA n64도 T=.01 완료, NJ에서 n32/n64 실제 field error 감사 중이다.
NK는 최신 소스 전체의 정순·역순 회귀를 실행 중이다.

### 10:53 分割 면 적분 weight root 검증

NL은 동일 n6 초기 상태, 동일 물질 생성점 속도, 동일 원래 질량/모멘트에서
부호 분할 face integral과 그 Jacobian으로 weight-rate만 풀었다.
기존 rate에서 시작한 최대 절대 flux는2.22037e-5→1.25516e-6→2.70530e-10
→1.86038e-15로 Newton3회에 감소했다. 독립order6 재평가도1.86029e-15,
질량 정규화 잔차4.03283e-13이다. gauge를 맞춘 weight-rate 최대 변화는
8.07583e-5. 물질 생성점 속도와 고정질량은 변경하지 않았다.

이는 기존 적분 오차를 weight-rate의 올바른 제약식 풀이로 줄일 수 있다는
한 상태의 증거다. 아직 split kinetic moment flux를 진화 RHS에 연결하지
않았다. 기존 mode3진화가 이 결과로 소급 검증되는 것도 아니다.
Jacobian 유한차분/상수 gauge/쌍별 질량 보존/root 수렴과 해석적 face
적분 검사5개 통과. 다음 필수 작업은 같은 분할 적분의 운동량·D·Pi·E와
양의 covariance factor 소스 연결이다.

NJ field audit는 n32/n64 모두 독립 polygon area gate에서 실패했다.
NM에서 문제 polygon을 보존: angular scale1.53e-4 및1.58e-4의 면이
이전 graded-angle 적용 범위(<1e-4) 밖에 있었다. 기존 적분16/24/32차
실패,64차 통과를 확인했다. 적용범위를1e-2로 넓혀 좁은 각도구간을
분할했고, 실제 n32 polygon의16/24/32차 area error<2e-10 및 양의
weights를 검사했다. 면적 gate를 완화하지 않았다. NN은 같은 저장된
T=.01 상태의 field audit 및 전체회귀를 진행 중이다. 진화 재시작이 아니다.
NK460개 정순/역순 회귀는 모두 통과했다.

### 10:58 공유 split moment RHS와 케플러 field 결과

NN464개 정순/역순 검사는 통과했다. T=.01 field audit도 두 해상도 완료.
24차 relative L1(n32/n64)은 rho .03232853558/.02278131067,
질량가중v .000422846085/4.97404848e-5, P .000104635112/2.13475436e-5.
이 결과는 짧은 시간의 공간 비교이며 공식100궤도 검증은 아니다.
NA n32 작은dt=.0001도T=.01 완료하여 동일 endpoint 시간 비교가 가능하다.

kinetic_moment_rhs_3d에 선택적 face_integration='split'을 추가했다.
각 내부 공유 면을 양쪽12개 velocity-sample 상대속도의0선으로 나누고,
모든 모멘트/압력 traction/공분산 source를 같은 양의 quadrature로 계산한다.
fixed reflecting wall도 표본 부호로 분할한다. 움직이는 normal wall은
아직 split mode에서 명시적으로 거부한다. 기본point 모드는 유지했다.

NO의 실제8셀 정적 RHS 검사에서 coincident cuts가 만든 정확히 zero-area
fan triangle 때문에 예외가 발생했다. NP에서는 split piece의 정확한
zero-area triangle만 제외하며 작은 양의 면적은 버리지 않는다.
중복 절단선과 zero-area piece에 대한 검사 추가 후 동일8셀 보존식/공분산
검사를 다시 실행 중이다. 시간진화에는 아직 연결하지 않았다.

NP의 실제8셀 검사 통과: order3/4에서 mass flux≈8e-18, 독립order6
평가도8.6e-18 이하. 운동량/에너지/각운동량 ledger error<2.34e-15,
raw covariance 미분과 primary factor 소스식의 차이<5.14e-17.
7개 targeted 검사도 통과했다. NQ는 NL의 기존 split root를 재사용하여
실제 n6 nonlinear vortex의 quadratic pressure RHS를 검사 중이다.

NA n32 T=.01 시간간격 비교: dt=.0001/.001 최대 상태 차이9.366e-8,
생성점 위치 차이9.368e-8, 셀 평균 속도 차이3.0243e-6. 원래 질량 배열
정확히 동일. 이 endpoint 민감도가 전체 dense path 인증을 대체하지 않는다.

### 11:03 split 3D 시간적분 시작

NQ 실제 n6 nonlinear vortex quadratic-pressure RHS 검사 통과:
mass=1.86125e-15, ledger=3.46885e-15, covariance identity=5.96918e-16,
generator defect=0. NL root를 그대로 재사용한 독립 공유모멘트 적분 검사다.

moment_evolution_3d의 두 RHS 단계에 선택적 split mode를 연결했다.
dense record는 face integration/order/pressure reconstruction을 저장하고
무손실 NPZ 복구와 후속 audit에 전달한다. 기본point+order3 config는
그대로 유지하여 과거 checkpoint 식별을 바꾸지 않았다. NR 정순467개
검사 통과, 역순 검사 진행 중. NS 실제8셀 split time step은 fixed mass
정확, stage residual1.38636e-16, energy ledger0, dense endpoint 일치 통과.

NS/n6은 기존 MX/n6의 원래t0 sites/weights/raw/factors/mass를 정확히
동일하게 복사하고 provenance hash를 남긴 독립 face-integration 비교다.
재초기화 적분을 하지 않고 기존 IC/initial score를 재사용한다.
face integration만 split/order4로 명시하고 dt=.002,T=.02를 유지했다.
11:02:47경 시작,2h 실행 제한 및 단계별checkpoint/dense NPZ 저장.
기존 point-mode 승인 prefix/전체 결과는 보존한다. 이 실행 완료와
전체 dense mass path 검사가 끝나기 전에는 새3D 방식 통과로 분류하지 않는다.

### 전체3D 경로 누적 감사 준비

NT의 audit_3d_history는 시작시점 checkpoint를 캡처하고 그 승인 단계의
모든 dense NPZ를 순서대로 읽는다. 시간/상태 연결, 원래 mass 배열의
정확한 일치, 최종 checkpoint 일치를 확인한다. root를 다시 풀지 않고
분할 면 적분으로 각 node의 실제 유량을 평가한다. 셀별 signed/absolute
오차를 먼저 시간누적하며, duration-scaled rate는 항상 설정의 전체
목표기간으로 정규화한다. prefix만 존재하면 목표기간을 완료했다고
표시하지 않는다. node cache는 source/library/record/config hash로 고정한다.

누적회계/연결/고정질량 검사와 split tests7개 통과. 실제1셀 정지 유체의
저장된 time step을 감사하고 두번째 호출에서 evaluator를 강제로 실패하게
해도 모든node cache가 재사용되어 동일report를 얻었다. 이 stationary
검사는 경로 감사 구현의 검사이며 nonlinear3D 물리 검증을 대신하지 않는다.
NT 전체 회귀 진행 중, 직전NR467개 정순/역순은 모두 통과했다.

MX 기존point 방식 n8은T=.02 완료. 12차 field L1은 density=.00325169,
velocity/amplitude=.0241335, pressure/amplitude²=.0286741로 해당
정규화 기준을 통과했다. 속도의 참조 적분 대비 L1은 .146928이며,
면 적분 mass 제약의 수렴은 별도다. NU는 n8 첫실제경로의 시작/끝
유량을 split order4/6으로 검사 중이다. 기존point 결과를 새split 방식
검증으로 혼동하지 않는다. NS의split n6 진화는 계속 진행 중이다.

### 11:11 baseline n8 독립 mass audit와 potential 교환 항

NU 완료: 기존point n8 첫 경로의 s=0 split order4/6 상대 mass rate는
.007778775468310845/.007778775468310694, s=1은
.007981507286457596/.007981507286456833이다. 절대 flux는 각각
1.50359e-5/1.54277e-5. 물리장 오차가 n6보다 감소했어도 면 적분에
의존한 mass root 문제가 해소된 것은 아니다. 따라서 MX n8의 field
pass만으로 전체3D 통과를 선언하지 않는다.

NT469개 정순/역순 검사 모두 통과. split n6 NS는 실제 계산 중이며
아직 첫 승인 checkpoint를 기다리는 상태다. 승인되지 않은 경로를
전체 history 감사 대상으로 사용하지 않는다.

smooth_mass_path에 선택적 potential callback을 추가하여 active–frozen
면 및 외벽에서 같은 mass flux의 Φ-weighted incoming 교환량을 구한다.
내부 active–active 면은 두 셀에서 정확히 상쇄한다. 기존potential 없는
계산 경로는 유지한다. NV에서 상수Φ의 교환량=Φ×active mass flux,
내부면상쇄 및 기존gravity dense tests5개를 통과했다. 아직 이 항을
전체 potential-energy ledger 검증에 연결한 것은 아니다. 다음에는
동일 saved trial의 Δ(Egas+Uactive)에서 gas boundary ledger와 시간적분된
potential mass exchange를 빼고, 독립16/24차 U적분 차이를 함께 검사한다.

### 11:18 첫split 경로 결과와 전체 energy ledger 검사

NS 첫 승인단계t=.002의 독립split order6,17node 감사 완료.
absolute 누적 상대mass defect2.44316e-7, signed2.25130e-7,
nested temporal quadrature difference3.72886e-9, 전체T=.02로
정규화한 peak relative rate9.15632e-6이다. generator interpolation
defect1.72489e-5이며, 이를 stage별 실제 조향과 구분한다(stage 조향0).
전체10단계 중1단계 결과이므로 covers_target=false를 유지했다.
NX는 동일split 물리소스와 같은 원래초기 배열, 같은T=.02에서
dt=.001로 시간세분화 비교를 시작했다. NS dt=.002는 계속 보존/진행한다.

NW는 저장된NE trial 세개에서 실제 hydro subflow의 potential mass
exchange를17node로 적분하고 active Egas+U 변화를 평가한다.
U는 polar16/24차 비교, node마다 원래mass audit flux와의 일치 검사,
frozen mask 전환 금지, frozen raw 배열 보존 검사를 포함한다.
모든 에너지 보정 없이 Δ(Egas+U)−gas_boundary_exchange−potential_mass_exchange
잔차를 기록한다. cold disk의 initial thermal energy로 나눈 잔차도 별도 출력한다.

공통initial U를 먼저 한 번 계산한 뒤 hash가 같은 경우만 세 trial에서
재사용한다. 노드/cell 적분은 checkpoint 가능하다. 대상은t=.01의
단일trial 세개이며 전체과거구간 에너지 보존을 소급 인증하지 않는다.
NW targeted4개 검사 통과 후 독립세trial 감사를 진행 중이다.

### 11:26 Kepler energy residual and density normalization diagnosis

NW completed all three saved trials. Energy residuals for dt=.001/.0005/.00025
are -6.21308847e-7/-3.10855550e-7/-1.55954859e-7. They are approximately
linear in dt; polar16/24 differences are <=1.78e-15. Gas ledger/gravity work
identity errors are <=8.12e-16. No complete energy-consistency pass is claimed.
NW regression: 473 tests normal and 473 reverse passed.

NY independently integrates every active reconstructed density at saved t=.01.
Worst cell410: native fit residual4.16e-17 but continuum relative mass error
-1.37610083e-5. Polar16/24/32 agree. Native fitted mass is exact to roundoff;
this distinction invalidates using native fit residual alone as physical proof.
Largest centered-moment discrepancy is1.97547e-5 (cell387).

NZ refits exactly the same cell410 polygon and original raw moments using
triangle Gauss orders6/8/12/16/24, without changing any accepted history.
Independent polar32 mass errors are -1.3761e-5/-9.5325e-8/-8.1224e-13/
-1.3656e-14/-1.3989e-14. This identifies underresolved density-fit quadrature
in this cell. It does not yet establish the cause of the entire energy residual.
Next: explicit quadrature-controlled density closure, followed by same-state
coupled trial comparison and independent energy/mass audits; no mass rescaling.

NS accepted t=.004 (2 steps). NT independent order6 history audit resumed
from its existing cached first step, now auditing step2. NS/NX processes were
confirmed live; no trajectory restarted. Full shock/long-duration/3D/official
192-squared100-orbit Kepler requirements remain incomplete.

### 11:28 Explicit density quadrature controls and completed short Sod32

Added density_quadrature_order configuration (6..64); default6 retains the
original quadrature function exactly. Higher orders affect density fit points,
normalization, and all existing consumers of density.weighted_mass, including
gravity and quadratic-velocity budgets. This is an explicit numerical-model
comparison, not a silent change to accepted trajectories. Same raw masses and
moments are retained. Default identity, configuration validation, saved-cell
independent polar mass convergence, and energy/exchange targeted tests:7 passed.

OA immutable snapshot runs same saved t=.01 Kepler state for q16 dt=.001,
q16 dt=.00025, and q24 dt=.001. Each trial must finish its17node hydro audit
before its energy audit runs. Independent cases run concurrently on lageunha.
Full476-test regression normal/reverse launched. No complete pass claimed yet.

NT NS independent history now covers t=.004 (two accepted steps): accumulated
absolute relative mass4.71138181e-7, signed4.51020696e-7, quadrature estimate
3.72929869e-9. NX first accepted t=.001 now has its own independent order6
history audit running; compare at equal physical times once both are available.

KV original P2 Sod32 completed T=.02,475steps,exit0. Driver numerical and
physical gates passed. Field L1 rho=.01071721,v=.01795574,P=.01256974;
absolute mass-path defect5.23020e-6, energy ledger8.84709e-17. This is the
short baseline, not the required long shock interval. MV long shock t=.048375,
JU vortex t=.135310, MH Sod32 t=.016988, MK Sod32 t=.015555 remain active.

### 11:31 Higher density quadrature fixes normalization, not energy consistency

OA476 tests passed normal and reverse. q16 dt=.001/.00025 full saved-step
energy audits completed. Independent relative mass errors <=1.93e-14.
Energy residuals remain -6.22333492e-7/-1.56019596e-7, comparable to q6.
Polar16/24 residual differences <=7.11e-15; relative to initial thermal energy
these are7.42231e-4/1.86078e-4. Thus density normalization underintegration is
real but is not the main cause of the energy residual. q24 comparison remains
running. Original fixed masses, material sites, and accepted history preserved.

OB now checks the actual saved hydro-stage directional derivative of U using
symmetric perturbations h=1e-5,5e-6. It compares this with independent native
and polar gravity work rates and the frozen-reservoir potential mass flux.
This discriminates semi-discrete closure consistency from time integration.
It does not apply any energy correction or advance an accepted trajectory.

NX first-step t=.001 independent order6 path audit completed: absolute mass
5.27184804e-8, signed4.38562873e-8, nested temporal difference1.15600e-9,
whole-target-scaled peak rate3.87203e-6, generator interpolation4.31441e-6.
NS first step was dt=.002: differing endpoints preclude calling this an
equal-duration convergence result. Need NX second step for that comparison.

### 11:33 Semi-discrete potential compatibility diagnosis

OA q24 dt=.001 energy audit complete. Residual -6.22333491813e-7 agrees
with q16 -6.22333491811e-7. Independent density mass errors <=2.07e-14.

OB on the actual saved q16 dt=.00025 initial hydro stage:
physical gravity work rate .00162809416064705 (polar), native differs4.3e-15.
Frozen-reservoir potential rate -1.97700365712e-7. Symmetric U derivatives
h=1e-5/5e-6 are -.00225292273726/-.00225292478007, giving residuals
-.00062463087624/-.00062463291905. This matches the approximately linear
single-step energy defect. This is semi-discrete, not solely timestep error.

OC checks fixed potentials1,x,y,x²/2,xy,y²/2 directly from stored q,Q,p,D
and the same face flux. All six conservation identities close <=9.66e-16.
Thus the basic moment transport identities pass on this actual state; the
nonpolynomial potential compatibility remains unfulfilled.

OD independently differentiates the exponential quadratic density closure.
For fixed local polynomial basis b and A=int rho b b^T, solve
A lambda=int rho Phi b. Reynolds transport yields
Udot=lambda.Mdot+integral_boundary rho(Phi-lambda.b)wn.
Compare against OB finite differences before assigning the residual entirely
to closure projection. No force modification or energy correction is applied.
A future compatible closure must address potential outside the retained
quadratic moment space; simply increasing quadrature or reducing dt is not
supported as a remedy by these controls.

NS reached t=.006 (3 accepted stages); NT audit resumed with previous cached
steps preserved. NX remains active after first accepted t=.001 step.

### 11:36 Analytic derivative and potential-enrichment feasibility

OD analytic Udot=-.00225292411258, agreeing with OB symmetric differences
within8.1e-10 at h=5e-6. Face quadrature16/24 agrees6.4e-15. Components:
lambda.Mdot=-.00167755149117 and moving-boundary remainder=-.000575372621407.
This confirms the actual closure derivative and the semi-discrete defect.

OE adds diagnostic cellwise Phi-weighted mass flux. Constant-Phi cellwise
identity and global internal-face cancellation pass targeted tests (7 total).
For each cell, orthogonalize r=Phi-lambda.b against the retained six moments.
An extra density tangent rho_dot=rho*r*missing/(m*E[r²]) can change Udot
while retaining m,q,Q to first order. No dynamics is changed by this audit.
Missing rate summed over active cells=.000624632251572. Maximum RMS relative
density-rate requirement=2.9264481 (cell13), corresponding to7.31612e-4 for
dt=.00025. This is local tangent evidence, not finite-step realizability.

OF tests positive exponential tilts in the six independent normalized
features (five spatial + potential remainder), holding the original cell
geometry,m,q,Q fixed and varying only the extra U moment by the missing rate.
Trials dt=.00025/.000125/.0000625 use independent polar24/32 checks.
No energy repair is applied to a simulation; no accepted history is replaced.
Full regression and actual worst-cell feasibility test are running.

OF actual worst-cell13 control completed: all three finite perturbations
converged in2Newton iterations. Independent polar24/32 checks retain mass
within2.23e-16 and five scaled spatial moments within7e-15; target potential
mean error<=1.03e-14. Minimum sampled density>.00959, and the exponential
representation is positive wherever finite. This supports local representability
of the missing potential direction, not global stability or coupled evolution.
Next: check all active cells, then propagate an actual extra potential moment
with the same mass-flux quadrature and gravity work. Recompute all density-
dependent velocity/gravity closures consistently; a diagnostic scalar alone
would not establish physical conservation.

### 11:40 All-active-cell enrichment feasibility

OG completed all572 active cells at fixed geometry with dt=.00025 missing-U
perturbations. All fits converged. Independent polar24/32 maxima: relative
mass9.10383e-15, retained scaled spatial moments1.99606e-12, potential mean
7.98916e-13. Cellwise progress and original hashes are preserved.
This is not a coupled timestep: geometry and m,q,Q were held fixed deliberately.

Implemented opt-in PotentialDensity interface exposing positive density values,
quadrature points/weights/weighted masses, and moment-fit diagnostics. It fits
supplied original m,q,Q and extra U; it rejects numerically dependent potential
bases and never changes targets on fit failure. Current evolution still uses
its original density until explicit stage-provider plumbing is implemented.
OH tests actual cell13 against independent polar32, weighted-density identity,
and rejection of redundant constant potential, plus full regression.
OF prior source passed478 tests normal and reverse.

NS third-step independent history covers t=.006: integrated absolute mass
6.91521e-7, signed6.87042e-7; fullT=.02 still pending.

### 11:43 Explicit density providers and first coupled potential-moment trial

Added opt-in volume density providers to kinetic RHS, gravity source and dense
mass-path evaluator, plus explicit density argument to decode. Existing original
providers reproduce RHS, gravity state/factors, and path flux arrays exactly.
OI targeted8 and483 normal passed; reverse pending. OH480 normal/reverse passed.

Potential-moment coupling now evolves per-cell U with Phi-weighted shared mass
flux minus physical gravity work. Every fluid stage and gravity kick reconstructs
its own density from original m,q,Q and its stage U. Frozen states and U remain
fixed. Physical acceleration is unchanged; gas energy is not reset from U.
A Strang kick/RK2 drift/kick trial uses independent same-state t=.01,dt=.00025.
OJ checks original masses/frozen data exactly and independently integrates the
final enriched density with polar24. No old trajectory or checkpoint is replaced.
Full history, dense enriched-path, temporal convergence and canonical longrun
remain required even if this first trial succeeds.

NS accepted t=.008; NX accepted t=.002. NT resumed both independent audits
from prior caches after confirming terminal receipts. NX will now permit the
first equal-duration dt=.001 vs.002 mass-path comparison at t=.002.
MV shock t=.054544; JU vortex t=.140387; MH/MK short Sod32 t=.018264/.017470.

### 11:46 First enriched coupled trial substantially improves energy balance

OJ dt=.00025 completed. Egas+U residual after gas/potential boundary exchange
is6.09726573e-11, compared with baseline q16 -1.56019596e-7 (about2559-fold
smaller absolute residual). Stage mass5.20417e-17, generator defect0. Independent
polar24 relative mass9.32587e-15, maximum per-cell U error2.46955e-14.
Gravity kick work2.45324062233e-8 and hydro potential-work integral
2.44714335797e-8 differ by the residual to roundoff, consistent with a remaining
split-time error; its scaling is not yet established. No gas-energy correction.
This is one saved-state step at t=.01, not a retroactive correction of earlier
accepted history or canonical longrun proof.

OI483 normal/reverse tests passed. Added optional actual hydro record including
U,du0,dum and its quadratic dense extension. Frozen U remains fixed and
extrapolation is rejected; OK targeted7 tests passed. OK dt=.000125/.0000625
trials run independently, each followed only on success by a17-node enriched
dense mass/potential boundary audit. Their common initial U was reused from
OJ after AST equality of initialization function and byte equality of relevant
density/geometry/potential dependencies; provenance/hash recorded per case.
No already completed fluid trial is replayed. Old snapshots stay immutable.

### 11:50 Temporal controls reveal remaining source-integration issue

OK485 normal/reverse tests passed. Enriched single-step residuals at
h=.00025/.000125/.0000625:6.09726573e-11/1.52394583e-11/3.80994439e-12.
They decrease by4, not8, per halving. This is a dt² local energy defect;
do not call the complete coupled method second-order based on its Strang/RK2
labels. Independent density mass errors remain <1e-14 and cell U error<2.8e-14.
Dense audits of the smaller trials are still running; all trial outputs kept.

ON isolates the gravity source with fixed density/geometry over duration.001
using1/2/4/8/16 substeps. Same enriched densities are cached by exact geometry
and original density moments; velocity reconstruction is recomputed at every
kick. Purpose: determine whether the nonlinear velocity closure makes the
current finite kick only first-order for its actual moment source flow.

OL field audit now detects stored potential moments and reconstructs enriched
density for reference L1 comparisons. Original checkpoint handling is retained.
OM continuation driver saves each trial before its17-node path audit, caches
nodes, updates accepted checkpoint only after cumulative mass/rate/quadrature,
energy/thermal and generator gates. Interrupted-audit unit test verifies no
fluid trial replay and no replay of a completed prefix. Targeted6 passed.
Physics continuation is prepared but awaits these path/source-flow controls.

At equal3D physical time t=.002, NS dt=.002 absolute path defect2.44316e-7
versus NX two dt=.001 steps6.05961e-8: reduction4.032. This supports temporal
convergence of this early path diagnostic, not complete nonlinear3D accuracy.

### 11:53 Source-flow first-order convergence confirmed in partial subdivision data

ON fixed-duration.001 source results1/2/4/8 subdivisions: successive summed
energy differences9.75451742e-10,4.87725832e-10,2.43862983e-10. Maximum state
differences1.30515046e-10,6.52575229e-11,3.26287616e-11. Exact factors of2
identify first-order convergence of the finite kick for the current nonlinear
velocity closure. n16 remains running; source files and cached endpoints kept.

OO experimental implicit midpoint gravity source reconstructs midpoint velocity
and updates p,D from fixed force moments, Pi from midpoint cross(v,g), E from
physical kinetic work. Cholesky requires positive-definite residual covariance;
rank-deficient/indefinite trials are rejected, never clipped. This prototype is
not yet a general replacement. Affine-force comparison to exact moment kick and
original provider tests4 passed. Same full-state source subdivision launched.
OO is not enabled in coupled production; original ON remains the comparison.

OL enriched physical-field audit completed at t=.01025. Relative L1 density
.0323340251, velocity .000422868464, pressure .000107238644, close to source
values. Sampled minimumdensity8.98695e-5 is below background.01, and maximum
relative density error47.4092 near the disk contrast remains large. Bulk L1
agreement alone does not certify pointwise accuracy or long-run disk stability.
OM486 normal/reverse tests passed. Dense enriched paths remain in progress.

11:53 MH n32 timeout124 confirmed at11:53:09, originalPID1685764 absent. Entire source manifest verified. Resume at t=.0188102991538,348steps with original checkpoint/pending caches; execution_resume_1153.txt and resume_1153_receipt.json preserve evidence. No initial restart.

### 11:55 Midpoint source improvement and completed enriched paths

ON16 subdivisions completed; energy differences halve at every refinement,
confirming finite-kick first-order source convergence. OO midpoint1/2/4 source
results show max state differences2.61228e-13 and6.57755e-14 (ratio3.97).
Summed energy differences7.17e-15/2.53e-15 approach roundoff and are not a
reliable energy order estimate. Larger subdivisions continue; this is a short
fixed-geometry source test, not a full time-integrator proof.

OP explicitly selects gravity_scheme=midpoint in the coupled extra-U method.
The default finite scheme and immutable old trials remain available. Targeted7
passed. Same-state dt=.00025/.000125 coupled trials started; common initialU
reused with unchanged initializer AST and verified dependency bytes. Results
will determine whether the complete split step recovers expected order.

OK finite-kick path audits completed: absolute relative mass1.16037149e-9
at h=.000125 and1.45035257e-10 at h=.0000625 (ratio8.00). Temporal nested
quadrature estimates2.51e-15/1.07e-15; generator defects1.08582e-7/2.71378e-8.
Integrated potential boundary exchange differs from midpoint-stage ledger by
1.62349e-13/2.02876e-14, consistent with dt³ scaling of that quadrature term.
It does not explain the dominant finite-kick dt² energy defect. No mass root
is re-solved at path audit nodes. Both remain single-step controls only.

### 11:57 Midpoint continuation preparation

OQ continuation identity now includes explicit gravity integration scheme;
default finite and selected midpoint cannot silently share a resume identity.
Path audit checks dense endpoint geometry and density moments against the
saved post-kick endpoint, and checks U as well. Targeted5 passed including
interrupted-audit reuse; full regression running. No physical continuation
launched before midpoint coupled/path evidence. OP trial PIDs1837441/1837442
verified active at11:56:59, both ~99% CPU; OO source PID1834097 also active.

### 11:58 Midpoint coupled energy residual reaches roundoff against stage ledger

OP completed dt=.00025/.000125 saved-state coupled trials. Energy residuals
against the midpoint-stage potential exchange are -1.89318e-15/-3.14676e-15.
Stage mass<=6.99e-17, generator defect0. Independent relative mass<=9.11e-15,
maximum per-cell U error<=2.80e-14. These residuals are at roundoff and do not
provide an energy convergence slope. Independent continuous-path boundary
exchange still needs evaluation: the prior finite-kick tests showed a dt³
stage-versus-path boundary difference, so no whole-path roundoff claim is made.
OQ now audits both actual saved midpoint paths with endpoint consistency checks.

JU timeout124 confirmed at11:57:39 and original PID1691812 absent. Complete
source manifest verified before resume. t=.145917239502,934steps retained,
as were pending stage and quadrature caches. New execution_resume_1158.txt,
run_resume_1158.log and resume_1158_receipt.json record the continuation.

### 12:00 Higher-resolution split3D control and midpoint source completion

OO completed all1/2/4/8/16 source subdivisions. Maximum state differences
2.61228e-13,6.57755e-14,1.64703e-14,4.12084e-15 decrease approximately4x
per halving until floating-point scale. Global summed energy differences are
near roundoff; use maximum state differences for this order observation.
OQ488 normal/reverse tests passed; midpoint dense path audits remain live.

NS t=.012, NX t=.004; independent audits were resumed from their existing
cached prefixes after terminal receipts. OR starts the independent n8 split
face-integration control at dt=.002,target.02, using byte-verified NS physics
snapshot and the original MX n8 physical initial arrays/score. Only experiment
configuration changes; prior n8 point-mode accepted trajectory is preserved.
OR timeout4h permits this heavier reference calculation to complete while
saving every accepted step and exact dense record. It is not a restart of MX.

12:00 MK sod32 timeout124 at12:00:02, PID1694969 absent. Verified complete immutable source manifest and resumed t=.018595780775,344steps. Checkpoint and pending quadrature retained; execution_resume_1200.txt/run_resume_1200.log/resume_1200_receipt.json record the action.

### 12:03 Completion requirements and equal-time3D field comparison

Refreshed the four top-level scope statuses and an explicit completion audit
in the JSON. ShockT=.15, ten vortex periods, complete nonlinear3D validation,
and canonical192²100-orbit fluid disk all remain unfulfilled. Unit counts or
single-step energy conservation do not substitute for those requirements.
Current shock t=.0606774; long vortex t=.1471585. Midpoint path audit handle
44705 verified running and both controls reached node12/16 at12:02:22.

OS reads preserved NSstep2/NXstep4 endpoints at identical t=.004, verifies
original mass arrays match, and computes physical reference scores at volume
orders8/12 for both. It also records raw/site differences and source hashes.
This comparison reuses saved endpoints without evolving either trajectory.

### 12:09 Independent midpoint path results and completed Sod32 reflection comparison

OQ both saved midpoint paths completed (exit0). At dt=.00025/.000125,
independently integrated potential-boundary energy residuals are
-1.3010934523e-12/-1.6549740554e-13. Absolute relative mass-path defects are
9.284028626e-9/1.160371925e-9. This approximately eightfold local reduction
supports the time refinement result; it is not a long-time convergence proof.
OT actual n32 continuation is live (PID1850161 observed12:07:44), seeded from
the already completed .01025 trial and its cellwise path accumulators, with
original source .01 and full target .02. No first-step replay or budget reset.

OS completed equal-time t=.004 n6 comparison. Volume-order12 density/v/P
amplitude-normalized L1 values for dt=.002 are .00094074635/.0386212833/
.0469039501; dt=.001 gives .00094075952/.0386214675/.0469034895.
Velocity reference-relative L1 remains .235132/.235133. Time refinement does
not materially improve this spatial reconstruction error; OR n8 remains
needed. NS independent path audit covers .012 (6steps), absolute relative
mass1.2974242e-6; NX covers .004 (4steps),1.1681250e-7. Different-time values
are not a temporal convergence ratio. Neither covers target .02 yet.

MH n32 completed T=.02,388steps,37rejections,exit0 at12:07:46. Existing
numerical/physical L1 gates passed; density/v/P field errors are
.01073227185/.01761985600/.01243192199. OU read-only saved-endpoint comparison
shows original P2 max mean |vy|=.21259948 versus reflected .0014733383;
transverse kinetic1.3304255e-5 versus3.0477577e-7. Reflected internal transverse
RMS .01231119 remains. Original n16/n32 max mean |vy|=.06197/.21260 worsens
with resolution, while reflected .000773/.001473 also fails a monotone
pointwise decrease: integrated L1 passage alone does not establish planar
invariance. Tensor-pressure n32 comparison remains necessary. No acceptance
tolerances were changed. Long shock .15, ten vortex periods, full nonlinear3D
and canonical Kepler192²/100orbits all remain incomplete.


### 12:11 Initial 3D representation and new independent path coverage

Initial order12 reference-relative velocity errors are .22955734 (n6) and
.13344684 (n8), already present before any evolution. Amplitude-normalized
velocity errors .03770561/.02191912 imply a two-resolution slope about1.89.
This diagnoses a substantial initial affine-reconstruction contribution, not
a proven evolved convergence order. Initial pressure errors .04535658/.02044214
also decrease. The evolved t=.004 n6 velocity error .23513 includes this
initial error; attributing all of it to time integration would be incorrect.

NS/NX independent order6 dense-path audits resumed after previous exit0,
reusing completed node caches (sessions94176/67255, execution_resume_1210.txt).
NS fluid checkpoint .014,7steps and NX .006,6steps were observed beforehand.
The MV long-shock calculation reached .06324772,1397steps and thus passed the
prior P4 failure time .06261795; this is progress, not completion of T=.15.


### 12:13 First additional enriched midpoint step accepted

OT reached t=.0105 with two accepted steps; original first-step state and
accumulators were retained. Cumulative absolute relative mass1.85288243e-8,
nested quadrature difference2.57746261e-10, full-target-scaled rate7.35994e-7,
latest generator defect4.36049e-7, energy/initialthermal3.10126e-9 all satisfy
the unchanged auxiliary continuation gates. These do not prove field accuracy.
OV captures this checkpoint byte-for-byte and evaluates density errors in
fixed radial bands [0,.45,.55,1.95,2.05,2.5], at independent orders16/24.
All bands contribute to the full-domain norm. This localizes density errors
near both jumps versus the bulk disk without masking them from acceptance.
MV shock t=.0643550,1429steps; JU vortex t=.1510875,976steps. MK tensor Sod32
t=.0193792,373steps. Full four-part validation remains incomplete.

### 12:16 Radial density comparison completed

OV evaluated the unchanged captured t=.01/.0105 enriched states at orders16/24.
All five radial bands sum to the analytic area pi*6.25 and reference mass
pi*3.775 within1.1e-14; full-domain L1 equals the sum of band absolute errors
normalized by total reference mass. This checks coverage and normalization,
not the numerical solution's pointwise accuracy.

Order24 full-domain density L1 .03232730684 -> .03233824175. The two jump bands
(.45,.55) and(1.95,2.05) contain83.27% of the evolved absolute density error.
Band-local relative L1 errors at .0105: inner background3.084%, inner jump48.04%,
disk interior .4886%, outer jump38.20%, outer background15.89%.
The outer minimum density8.94156e-5 is far below physical background.01.
The starting state already has minimum8.92696e-5 and outer-background L1
15.8865%; the new two-step midpoint evolution does not originate that error.
It also does not repair it. Order16/24 full-domain L1 differ1.0052e-5; sampled
extrema depend on quadrature points and are not rigorous continuum bounds.
No band was omitted from acceptance or the full-domain norm.

12:17 MV timeout124 confirmed at12:16:31 and original PID1715438 absent.
All386 immutable source manifest files matched. Resumed t=.0653657690,
1469steps with checkpoint and pending caches preserved, session65744,
execution_resume_1217.txt/run_resume_1217.log/resume_1217_receipt.json.
This resumes the original T=.15 budgeted run, without restarting from t0.


### Tensor Sod32 completed; remaining wall-localized transverse growth

MK/sod32 completed T=.02,390steps,exit0 at12:18:20. Existing numerical and
physical gates passed: density/v/P field L1=.010708526/.017752911/.012432669,
absolute relative mass-path defect2.94221666e-6. OX read-only endpoint audit
shows max mean |vy| increases n16->n32 from1.73963e-4 to8.73399e-4 and
transverse kinetic1.31568e-10 to1.65691e-9. These remain far below uncorrected
P2 but do not establish pointwise symmetry convergence.

OX sampled accepted archives1,10,50,100,200,300,390 with original masses,
accepted commit states and hashes verified. First-step max |vy|2.64e-19;
step10 1.10e-9; step50 6.78e-6; final8.73e-4. Fromstep50 onward the worst
cell17 is at y~.015625 (bottom wall), x~.547. Final maximum transverse
acceleration~.3084. This localizes growth but does not yet distinguish
reconstruction/boundary/flux/time error. No evolved state was repaired.

OW launches identical radial density diagnostics on original n32/n64 saved
states at equal physical t=.01, orders16/24, with unchanged original closure
configurations. No new trajectories or enriched/original mixing. Full physical
validation remains incomplete; new comparisons do not change acceptance goals.

### Velocity curvature diagnosis and explicit tensor fit alternative

OY saved-step1 has longitudinal velocity's artificial transverse curvature
4.93944e-8 while pressure curvature4.44e-16 and column pressure spread1.78e-15.
Step10 velocity curvature2.06808e-5 precedes substantial later pressure/mean-vy
error. This identifies a reconstruction defect, without claiming that it is
the only source of subsequent transverse growth.

An explicit velocity_fit=tensor_free_mean option replaces the candidate fit
with separable weights and a free intercept. The intercept is discarded; the
same density-orthogonal basis and existing covariance budget retain stored
m,p,D,Pi. Default radial_fixed_mean formula is bit-identical. No steering,
state projection or modified mass is introduced. OZ targeted20 tests and all17
existing velocity tests with the new mode passed. Full default regression
is running. OZ same saved step1 transverse curvature falls to3.92395e-15;
step10 to3.62245e-9; step50 to3.81827e-5. Existing evolved pressure errors are
unchanged in this read-only comparison. PA n16/n32 independent controls are
prepared from original byte-verified t0 arrays; launch awaits full regression.

OW original n32/n64 same-time radial audit completed at orders16/24. Order24
full density L1 .03232727 -> .02277358. Outer-background band relative L1
.158868 -> .00147766 improves strongly, but inner-background .0298436 ->
.100127 worsens; min inner density at n64 is .000115958 versus physical .01.
Disk-interior error .00487661 -> .00427214 decreases only slightly. All bands
still contribute to the full norm. Added polar splits change quadrature of
absolute error compared with earlier unsplit radial scorer; these are sampled
numerical estimates, not uniform convergence or continuum error bounds.

OZ491 normal regressions passed. PA independent n16/n32 controls launched
(sessions81887/37229) with original t0 physical arrays and explicit new
velocity_fit setting. Reverse-order suite is running (48629). Earlier MK/MV
trajectories continue unchanged; PA is a distinct model comparison, not their
restart. Long shock/long vortex/full3D/canonical Kepler remain incomplete.


### 12:27 3D spatial comparison and refreshed dense-path audits

OZ491 reverse-order regressions also passed. PA n16/n32 accepted their first
steps with the explicitly selected tensor velocity fit; evolution comparisons
remain pending. OT enriched midpoint Kepler reached .01075,3accepted steps:
cumulative absolute relative mass2.77263e-8, energy/initialthermal4.66256e-9.

NS fluid t=.018,9steps; NX t=.007,7steps; OR n8 firststep t=.002 available.
Previous NS/NX audit exit0 receipts confirmed before resuming. Independent
order6 audits now capture these prefixes, reuse all completed node caches,
and add the first OR n8 path (sessions60502/52849/48201). PB compares physical
fields of original NS/OR firststep endpoints at equal t=.002, verifying each
mesh's original masses separately. Cross-resolution arrays are not subtracted.
Full four-part physics validation is still incomplete.

PB completed. At equal t=.002, order12 n6->n8 density/v/P amplitude-normalized
L1 .000473037/.038149375/.046082011 -> .000343816/.022065751/.020825988.
Velocity reference-relative L1 .232259 -> .134340 decreases. Sampled maximum
vertical velocity .000307511 -> .000405582 increases, so this comparison
supports early integral-error reduction but not uniform transverse symmetry
convergence. First-step dense-path n8 audit remains separate and running.


### 12:29 Equal-time evolved symmetry does not yet improve

PC initial comparison failed safely because its binary selector treated
accepted report t as a start time; report t is the endpoint while record t
is the start. PD corrects this selection and adds a passing regression for
unequal step widths and exact endpoints. No simulation states were changed.

At identical t=.004 using preserved accepted longdouble dense extensions:
MK n16 max mean |vy|6.47690318e-7, transverse kinetic4.13372e-15; PA tensor
velocity n16 max mean |vy|6.48134085e-7, transverse kinetic7.88454e-15.
There is no demonstrated evolved-symmetry benefit at this time despite the
large reduction of reconstructed transverse curvature on earlier saved states.
The new mode remains explicit and experimental; default physics is unchanged.
Both full short comparisons continue. Additional sources or amplification of
transverse error remain unresolved. This result is evidence against treating
the candidate fit alone as a completed symmetry repair.


### 12:32 Reused 3D path convergence evidence and long-vortex field audit

PE verifies identical original n6 physical arrays, audit identities and exact
saved-record hashes, then accumulates cached per-cell absolute/signed fluxes
at identical times. Coarse dt=.002 versus fine dt=.001 absolute mass ratios
at t=.002/.004/.006 are4.03188/4.03329/4.01951. At .006 defects are
6.91521e-7/1.72041e-7, nested differences3.72963e-9/1.15840e-9. This supports
early second-order accumulated path-error reduction, not full-interval or
physical-field completion. No previous quadrature or fluid step was replayed.

PF independently scores captured JU t=.1558767064,1026steps against the
stationary reference using its immutable source configuration. Density/v/P
normalized field L1=.001301275/.000651476/.000359208; minimum density.942644,
minimum pressure.685133. Original mass array matches exactly, reference
quadrature difference1.94e-11 and reference total-mass error6.66e-16.
Ten central vortex periods and historical full-duration path accounting remain
unverified; this is an intermediate field audit only.


### 12:34 Nonlinear split-face n6 evolution completed; pressure gate failed

NS reached T=.02,10accepted steps. Order12 density/v/P normalized L1:
.004549497/.042642571/.055557249. Pressure exceeds the unchanged .05 goal,
so physical_passed=false. Original mass array exact, max stage mass flux
1.77323e-13, energy ledger error0. Reference-relative velocity error.259614;
maximum sampled vertical velocity.00257459. Completion of evolution is not
completion or passage of the full 3D validation.

Independent NS order6 audit completed through .018 with absolute relative
mass1.90362e-6; after its exit0, finalstep audit started (72933) reusing nodes
for all9earlier steps. OR n8 first-path audit completed at t=.002: absolute
relative mass3.154845e-7, quadrature difference8.86216e-10, normalized peak
rate9.69050e-6, generator path defect1.48475e-5. This is a first-step prefix,
not coverage of .02. Ongoing NX time refinement and OR spatial refinement
remain required. All failed/narrow results are retained with their scope.


### 12:36 Long-shock refinement and tensor-velocity n16 endpoint

PG captured MV t=.07091958,1666steps: density/v/P field L1
.02747232/.04784140/.02471261. Velocity approaches .05 target; original mass
array exact. Mean transverse maximum.00718606, transverse kinetic1.54175e-7.
This is an intermediate score, not completion of .15. PH starts the independent
n32 resolution control with whole T=.15 budgets from t0 and established
pressure-tensor/radial-velocity physics from MV, original n32 physical IC.
All prior short trajectories remain preserved. Session29011, timeout4h.

PA n16 completed .02 at12:35:06,127steps,exit0, numerical and L1 gates passed.
Density/v/P errors .01320747/.02938712/.02228145; absolute relative path mass
2.80088e-6. PI compares actual endpoints: max mean |vy| .0001739633 ->
.0001018828 (41.4% lower), transverse kinetic1.31568e-10 ->7.37754e-11
(43.9% lower). This later-time modest benefit contrasts with the nearly equal
.004 diagnostic, but does not establish exact planar invariance or improved
spatial convergence. PA n32 remains running; default model stays unchanged.

NS final independent order6 audit completed T=.02,10steps,covers_target=true.
Absolute relative mass2.07609926e-6, signed2.06423125e-6, nested quadrature
difference3.72985510e-9, duration-scaled peak9.15631982e-6, generator path
defect1.72489193e-5. All10actual saved paths are covered; this does not
reverse the pressure L1 failure or establish continuum enclosures.


### 12:42 Canonical192 initial integration launched with resumable chunks

Added initialize_kepler_cached.py: immutable source/site identity, cached
geometry and physical moment integrations in32-cell chunks, independent
order16/24 initial quadrature check, and no overwrite of different initial
states. It performs no fluid evolution. Reference potential integrals are
stored under reference_potential_moments, distinct from evolved enriched U.

PJ cache unit test passed. Actual n16(246cells) stop after32/resume preserves
the first chunk and yields bit-identical mass/state/factors to a direct
unchunked integration of all cells.493 normal regressions passed; reverse
running. PK starts official n192 initial integration with the unchanged ring
layout and physical fields, source manifest verified, one CPU thread on
lageunha, timeout4h, session35559. Completed chunks survive timeouts. Density
reconstruction gates, time evolution, and100orbits are not certified by this
initial-only calculation. All original completion requirements remain active.


### 12:46 Early transverse resolution comparison remains inconclusive

PJ493 reverse regressions passed. PK official192 process1896632 was verified
live during initial geometry construction; no initial integration completion
claim is made. PL uses preserved dense trajectories at exactly t=.003:
old/new velocity-fit n16 max mean |vy|2.61110e-7/2.43596e-7; n32
3.05883e-6/2.96911e-6. Transverse kinetic old/new n16
6.63463e-16/1.10851e-15; n32 1.26789e-14/2.31630e-14.
The small maximum-mean improvement coexists with larger transverse kinetic
energy, and the finer mesh has larger symmetry error. The optional fit is not
a verified symmetry repair. Same physical time, archived-record hashes, and
original mass arrays are checked; no fluid evolution was replayed.


### 12:48 Saved-state mass-root tolerance sensitivity

PM evaluates the same accepted MK/sod32 step10 endpoint t=.00038406982 with
mass-root tolerances2e-13 and2e-15 in an isolated read-only process. Both
converge in12iterations to3.96631e-19 residual; derivative and weight rates
are identical. Maximum transverse acceleration remains1.08121757e-5,
transverse kinetic rate7.95031e-17, residual-covariance y rate4.47962e-21.
No physical state or production solver tolerance was changed. The root's
stopping tolerance alone does not explain this instantaneous transverse force;
this is not an exclusion of cumulative error at every earlier state.

PK canonical192 geometry completed and began persisting32-cell integration
chunks; at least768/36864 cells observed. Existing initial and control work
is retained, and completed chunks will be reused on any required resume.


### 12:57 Free-affine velocity fit: counterexample and independent controls

PN/PO force decompositions show pressure traction dominates the early
transverse acceleration. At the same first accepted state, the original fit
has column pressure-rate spread 6.178851e-6; the tensor/free-constant fit
still has 6.79093e-6. A truncated tensor stencil with a planar odd residual
produces false xy curvature (0.2928317) under the free-constant fit.
The optional tensor_free_affine mode fits nuisance constant and linear terms,
then discards them: the stored momentum and mixed moments still determine
the physical affine velocity. No mass rescaling or generator steering occurs.

PP targeted tests (4), existing velocity invariants under the new option (17),
and normal regression (494 passed, one Sod-reference log warning) completed.
On identical archived MK states, PP reduces first-step pressure-rate spread
to 5.339773e-11; at step10 it is 4.837546e-5 versus the original 0.00245476.
The pressure traction already encoded in later archived states remains;
this read-only experiment does not prove evolution or spatial convergence.

PQ independent n16/n32 T=.02 controls were launched on lageunha using the
verified PP source manifest and the original initial-array/configuration hashes.
Only the velocity fit setting changes relative to PA; all error thresholds
are retained. Two-hour bounded CPU runs use sessions33324/55313 and preserve
existing trajectories in their original directories. Root:
`/gpfs/kjhan/LagForce/Codex_PQ_velocity_affine_sod_controls_20261017`.
The default remains radial_fixed_mean. Full shock, long-time, nonlinear3D,
and canonical192 100-orbit validation remain incomplete.


### 12:59 Early evolved planar-symmetry comparison

PR compares preserved accepted dense states at exactly t=.0004 for n16.
Original radial, tensor/free-constant, and tensor/free-affine velocity fits
respectively give max mean |vy| 6.376876e-10, 6.670424e-10, 2.908398e-17;
transverse kinetic energies are 1.156560e-21, 1.103417e-21, 1.013085e-34.
Trial hashes and unmodified initial mass arrays were verified. This is
actual evolution evidence for early symmetry improvement, not a long-time
or mesh-convergence certificate. PQ n32 is still performing its initial
adaptive step selection; no failure or restart is inferred from rejections.

The independent n8 3D audit completed through t=.004 (2 steps): accumulated
absolute relative mass 6.175911e-7, nested quadrature difference8.870148e-10,
duration-scaled rate1.053320e-5, generator path defect1.484746e-5.
Full target .02 is not covered. PS separately evaluates n6/n8 physical fields
at t=.004 from their actual saved endpoints; no evolution is replayed.


### 13:00 Extended 3D temporal mass-path evidence

NX independent order6 audit completed through t=.01,10steps, exit0.
PT compares identical initial arrays and independently cached mass-path nodes
at matching times: coarse/fine absolute mass ratios4.017124 at t=.008 and
4.016554 at t=.01. At t=.01 absolute relative mass is1.099813e-6 versus
2.738201e-7; nested quadrature differences3.729846e-9 versus1.158399e-9.
NX generator defect4.314414e-6 and duration-scaled rate3.872029e-6 are
reported separately. This extends evidence of second-order path convergence;
it does not resolve the known n6 pressure-field failure at T=.02.

PS first invocation exited2 before computation because the comparison
script was absent from the source snapshot. The repository script was copied
into the standalone audit directory and a local launcher imports immutable
PP dependencies. The corrected bounded read-only audit is session3754;
no simulation trajectory was restarted.

PS corrected audit completed exit0. At t=.004, n6→n8 order12 density L1
.000940746→.000682323, velocity/amplitude L1 .0386213→.0222433,
pressure/amplitude² L1 .0469040→.0213843. Sampled maximum vertical velocity
increases .000602517→.000789471. Integral field errors improve with refinement,
but vertical symmetry does not improve uniformly. No full spatial convergence
certificate is inferred from two resolutions and this early time.


### 13:02 Vertical pressure reconstruction and evolved planar symmetry

PU read-only audit and two derivative unit tests completed. The analytic
vertical derivative of each quadratic pressure polynomial is affine; its
absolute maximum is evaluated at all convex-cell vertices. Original files
and masses are checked. At t=0, n6/n8 maximum reconstructed |dP/dz| is
0.09274172/0.14774216, despite exactly zero initial vertical mean velocity
and a vertically invariant physical reference. At t=.004 the gradients are
0.09369875/0.14495708 and maximum mean vertical velocities are
0.000245971/0.000373882. No pressure positivity limiter was active.
This identifies a pre-evolution spatial reconstruction defect; it is not a
complete momentum-force decomposition or proof that every later error has
one cause. Sites are irregular, so tensor-product planar invariance cannot
be assumed from the 2D Cartesian control. No state, mass or generator changed.

PV same-time actual dense-path comparisons extend the free-affine evidence.
At n16,t=.004 original/free-constant/free-affine max |vy| is
6.476903e-7/6.481341e-7/2.122667e-15; transverse kinetic energy is
4.133719e-15/7.884541e-15/9.916858e-32.
At n32,t=.0001 max |vy| is2.095339e-11/2.170079e-11/4.561392e-17,
and transverse kinetic energy3.609127e-25/3.840665e-25/9.767878e-35.
These two resolutions are compared at different times, so they support
within-resolution improvement, not a cross-resolution convergence order.
Full-duration shock and physical-field criteria remain required.

OT midpoint Kepler has accepted9steps to t=.01225: absolute mass8.393648e-8,
quadrature mass1.273482e-9, rate1.130149e-5, generator4.464017e-7,
energy/initial thermal1.422086e-8. Official192 initialization and all longer
validation targets remain incomplete; no canonical completion is inferred.


### 13:05 3D candidate fit tradeoff and long-vortex gate review

PW evaluated only preserved initial states, with prototype modules saved beside
the audit and the original production implementation unchanged. Free-constant
radial fitting reduces max vertical pressure gradient n6 .0927417→.0428911,
n8 .147742→.0608217; tensor/free-constant gives .0496972/.0752265.
However order6 pressure L1 slightly worsens: n6 .0453460→.0459774 (radial)
or .0469269 (tensor); n8 .0204311→.0204953 or .0212563. All cell pressure
means are preserved to1e-15 and no limiter activates. This is a tradeoff,
not a validated model improvement; no new evolution with these models launched.

PX checks PF's saved vortex score against actual per-field configuration.
At t=.1558767064 density_field=.00130127493 exceeds its .001 limit;
the other five checked field/cell measures pass their respective limits.
The reference quadrature is valid. This is an intermediate physical-gate
failure, not a prediction of the future endpoint. JU remains preserved/running
for diagnostic evolution; it cannot certify all-time physical accuracy.
The full10-period target is10.4685979156. Even scaling the current stored
rate peak by full_target/.2 gives a lower bound exceeding the1e-4 rate
budget; earlier shorter prefix durations can make the exact bound stricter.
Therefore merely extending the target cannot certify the existing trajectory
under the full-duration rate gate. PX captures the actual checkpoint and
calculation without rewriting any accepted history or relaxing requirements.


### 13:08 Full-horizon allocation and completed short free-affine shock

PY adds optional validation_end, an absolute target no earlier than execution
end. Every step mass/quadrature allocation and rate gate uses the duration
from original source time to validation_end. Default behavior is unchanged;
shorter inspection intervals cannot enlarge their per-step budgets. Invalid
or nonfinite horizons reject. Report includes covers_validation_target;
configuration hashes still guard resumes, no existing history is migrated.
Six targeted and497 normal regressions passed (existing Sod log warning).

PZ starts an independent n32 vortex temporal control with original JU initial
arrays and model, execution end=.01, validation_end=10.468597915605091,
max_dt=.00025, unchanged thresholds. This isolates stricter full-horizon
time allocation; it does not promise to repair spatial density error.
Session86530, lageunha CPU threads1/path workers4, timeout3h. Prior JU runs
and caches remain intact. Full10period completion is explicitly false.

PQ n16 free-affine shock reached T=.02 in117accepted steps, exit0 at13:05:50.
Numerical and physical gates pass for that short interval: rho/v/P field L1
.013303149/.029607562/.022291339; absolute relative mass3.417785e-6.
QA compares actual accepted endpoint paths: original/free-constant/free-affine
max mean |vy|1.739633e-4/1.018828e-4/8.498615e-12; transverse kinetic
1.315684e-10/7.377544e-11/2.455946e-25. Global field errors are not uniformly
smaller, although symmetry is strongly improved. PQ n32 remains running.

QB starts the independent n16 free-affine whole T=.15 test, using original
initial arrays and full-duration allocation fromt0. Source PY manifest and
original initial/configuration hashes verified. Session45359, timeout4h.
This is a new candidate control; completed and active older trajectories are
preserved. All full-scope validation targets remain incomplete.


### 13:10 Full-horizon initial-step regularity and continued 3D audits

PZ rejected its initial .00025 step and adaptively reached an accepted
first endpoint t=3.90625e-6. QC evaluates actual archived rejected samples,
not newly solved geometry. Full-duration rate decreases .00170398→.000425979
→.000106494 as dt halves initially (worst cell336), then5.892825e-5→
4.232966e-5→3.018376e-5 (worst cell877). Initial-node relative rate remains
2.593558e-11 and stage mass residual4.061993e-15. Later rejection ratios
are dominated by quadrature (9.59457,6.79213,4.80419), with additional
mass-allocation violations. No solver tolerance or budget was relaxed.

Cell877 starts with zero residual covariance-factor trace. Across these
six saved trials its endpoint trace increment/dt tends to1.23033e-11;
cell336 starts at6.39339e-9 and has increment/dt about1.33e-8.
This is evidence of O(dt) covariance generation at an initially rank-zero
cell, consistent with square-root velocity-scale regularity; a causal
attribution of the full path error is not yet established. Values are
preserved factor norms, with no clipping, seed covariance or mass correction.

At13:08 NX was live at t=.012 and OR live at t=.006. Earlier independent
order6 audits had terminal exit0; resumed from cached nodes with new receipts
execution_resume_1308.txt (NX session65329, OR session53040), timeout1h.
No old path nodes or simulation steps are replayed. Official192 initial
integration was live with10592/36864 cells persisted. QB accepted7steps to
.000177492 under full .15 allocation; all full targets remain incomplete.


### 13:14 Third 3D spatial resolution and pipelined canonical IC audits

QD uses the verified OR source manifest and unchanged split4/affine-velocity/
quadratic-pressure nonlinear vortex model at n10, dt=.002,T=.02. Initial-only
run exited0 at13:12:37, preserving the exact normal initializer output.
Order12 initial rho/v/P amplitude-normalized L1 is9.639214e-12/
.013800178/.012762244; velocity relative to reference speed is.084017529.
Initial sampled vz=0 and positive pressure/density. The same saved initial
state now proceeds into evolution, session73509,timeout4h; no initialization
replay or altered original masses. lageunha observed32 CPUs, load19.63,
937GiB available before adding the single-thread control.

QE audits already completed official initial-integration chunks while PK
continues. Source/site/geometry/chunk/script identities protect cache reuse.
It reconstructs the actual factor-completed cell state with density order16,
and independently integrates mass with polar orders16/24. It neither edits
physical arrays nor solves alternative geometry. A32-cell n16 real-data
preflight gives max independent relative mass8.379257e-14, quadrature
difference3.349748e-15, minimum pressure2.777777773e-5. A second run replacing
decode with a function that fails proves completed cache reuse without
reconstruction; cache and report bytes remain identical.

Canonical QE session49974,timeout2h, processes available chunks and stops at
the first missing chunk. Report explicitly distinguishes audited prefix from
all36864 cells. Density-order16 validity does not validate the native order6
model, enriched potential closure, intercell pressure fit, or100orbital periods.
All full-scope completion claims remain false.


### 13:16 Third-resolution vertical-gradient evidence

QF n10 initial-state audit completed: max reconstructed |dP/dz|=.13807905,
worst cell549 near the upperz wall, initial meanvz=0, no pressure limiting.
The corresponding n6/n8 values were .09274172/.14774216. Thus the third
resolution lowers the n8 maximum slightly but does not establish monotonic
maximum-gradient convergence, despite improved integrated field errors.

At13:15:43 NX remained live at t=.013, so no premature restart was made.
Its independent order6 audit completed through captured t=.012 with exit0
at13:14:43. QE canonical initial reconstruction reached928cells with no
reported fit failure; this is only the audited prefix, not completion.


### 13:18 Verified NX timeout and checkpoint continuation

NX timeout exited124 at13:16:44; oldPID1790182 absent at13:17:16. Verified
every immutable source-manifest entry before resuming the actual checkpoint
t=.013000000000000005,13steps. No pending next-step dense record existed;
all accepted records and original mass/state arrays retained. Receipt
resume_1317_receipt.json includes checkpoint hash and terminal evidence.
Session31155 resumes for2h with unchanged parameters; no accepted-step replay.

QG independent cached-node comparison at t=.012 gives coarse/fine absolute
relative mass1.297424e-6/3.230168e-7, ratio4.016584. Nested quadrature estimates
are3.729851e-9/1.158400e-9. Temporal path evidence remains consistent with
second order; it does not erase the n6 physical pressure failure atT=.02.


### 13:20 Equal-time two-resolution symmetry and canonical reconstruction prefix

QH evaluates all three velocity fits at the same physical time t=.004 for
both n16/n32. Max mean |vy| original/free-constant/free-affine:
n16 6.476903e-7/6.481341e-7/2.122667e-15;
n32 8.061891e-6/7.434164e-6/5.280132e-14.
Transverse kinetic energy under free-affine is9.916858e-32/1.185391e-29.
Both meshes show strong improvement. These near-roundoff values do not
establish a conventional spatial convergence order. Actual archived paths,
original mass arrays and record hashes were checked without replay.

QI captures the first2624 canonical density-order16 audited cells and checks
each cache against its original IC chunk hash and common audit identity.
Worst independent mass relative error1.257942e-12 (cell2423), maximum polar
order16/24 difference3.320217e-15, minimum pressure2.777777710e-5.
This prefix includes only2624/36864 cells; no complete initial/evolution
validation is claimed. PK had15488 integrated cells at13:18:23 and remains
live. PQ32 t=.00616421, QB16 t=.01599194, PZ vortex t=.000616781 were
observed advancing. QJ separately captures the original long-shock current
state for an exact-reference physical score, preserving source files.


QJ completed the exact-reference original MV shock score at t=.0845738673,
2244steps. Density/velocity/pressure field L1 .02949335/.04199748/.02711791
and all three cell errors are below .05 at this captured time. Reference
quadrature difference5.392618e-16 and original masses unchanged. Max mean
transverse velocity.00642806, transverse kinetic1.350978e-7 remain nonzero.
This is an intermediate field assessment, not fullT=.15 or every-time proof.


### 13:23 Missing n8 temporal control added

QK starts the n8 dt=.001,T=.02 split4 nonlinear control from exact original
OR sites/weights/raw/factors/mass arrays. Only configuration dt changes;
initial integration and score are reused, with original IC hash recorded.
All source manifest entries verified. Session58672, lageunha one thread,
timeout4h. This fills the missing temporal comparison at the finer physical
resolution; no result is claimed before actual paths/fields are audited.
Refreshed the explicit completion matrix from current checkpoint artifacts;
all four full-scope requirements remain unfulfilled.


### 13:26 Runtime velocity-fit identity guard

QL verifies the actual environment velocity-fit mode matches the recorded
configuration, and rejects nondefault curvature fits with affine-only velocity.
These checks run before initial/checkpoint creation. Two driver-level negative
cases plus four stencil tests pass; normal regression499passed, with the
existing Sod-initializer logarithm warning. Existing running snapshots were
not edited; this guard changes configuration validation, not the equations.
At13:24:46 PQ32 had86steps,t=.00733461; QB16 had234steps,t=.02392065;
PK had18592 integrated initial cells and QE5472 reconstructed/audited cells.
All counts are intermediate observations, not full validation claims.


### 13:31 Time accuracy versus spatial error in3D and vortex

QM compares exact same-time n6 saved endpoints at t=.012, dt=.002/.001,
with identical original IC arrays. Order12 rho/v/P errors are coarse
.002769770/.040599882/.050888944 and fine
.002769755/.040600279/.050887405. Both pressure errors exceed .05;
maximum raw endpoint difference5.084309e-9 and site difference3.960232e-8.
Thus halving time step improves the separately measured mass path but does
not remove this physical pressure error. No evolution was replayed.
OR reached t=.008; its prior order6 audit exited0 at13:16:25 through .006.
Resumed existing caches at13:27,session71849,timeout1h, for the new path.

QN evaluates identical original n32 vortex states under three velocity fits.
Initial density equals1 to geometric precision. Max spurious mean-density
rate radial/free-constant/free-affine=.03162034/.03485959/.03604841;
volume-integrated absolute rates=.001718594/.001849956/.001771351.
Global area rates are around1e-17 and mass-root residuals4–6e-15; original
state/factor/mass bytes remain unchanged. The new planar fit is not a measured
improvement for this irregular vortex and was not substituted in active runs.
This instantaneous consistency diagnostic does not prove a later-time error.

QO prepares n16/n64 initial-only spatial consistency controls using the same
vortex definition and baseline radial velocity fit. Exact initializer arrays
and initial physical scores are saved for later reuse. Sessions29930/94011,
one thread each,timeout1h. No evolution has been launched; purpose is to
measure how the initial spurious density rate scales with spatial resolution.
All original long-time and canonical targets remain required.


### 13:34 Spatial improvement of initial vortex consistency; tight path events

QO initial-only n16/n64 completed with original generated arrays saved.
Combined with QN's unchanged n32 baseline, volume-integrated absolute initial
density rates are .01087597/.001718594/.000289014; max rates
.10458178/.03162034/.006430708. This is clear initial consistency improvement
with spatial refinement, not a long-time density-error certificate.
n64 initial pressure/velocity field errors6.441056e-6/1.556644e-5 and density
field3.345010e-11. Reference moment/score data and initial states are retained.

PZ reaches t=.00109775476 but rejects some very small steps. QP reads only
actual rejection archives: rejections17/18 at the same source time have dt
3.304679e-6/1.652339e-6 and nearly identical quadrature ratio15.833464.
Worst flux is at endpoint s=1,cell719; initial-node relative rate1.20655e-11
and stage residual3.68802e-15 are small. Halving again lowers the ratio;
there is no evidence of a terminal failure or permission to restart.
QQ examines the stored quadratic covariance-partition metadata at five
phases of these two paths (session1956,timeout20m). No closure, moment,
threshold, or running source is changed. Higher-resolution long evolution
is deferred until this localized time-control issue is understood better.


QQ completed both saved-path audits. Cell719 residual covariance fraction
changes smoothly .11735663→.11736261 across rejection17; neighbor750 stays
at the same saturated residual-map singular value1.4142e-6 throughout.
No sudden covariance-partition switch was found at these sampled phases.
However edge686–719 exists at s=0 with length2.627594e-8 and disappears
before s=.125; the adjacent polygons change6→5 and7→6 vertices.
Per-cell oriented face-normal closure remains within2.22e-16, with no false
boundary face in the inspected cells. This is evidence of a genuine local
mesh-topology event, not a demonstrated geometry closure bug. Quantifying
its contribution to the dense-path error is the next diagnostic; covariance
caps, solver tolerances, and accepted trajectories remain unchanged.


### 13:45 Topology-event trial and independent 3D audit continuation

QR finds the isolated weighted power-circle event for sites686/687/718/719
on PZ rejection17 at phase .04060305493, time .00109788894035.
The actual saved geometry confirms that face686–719 disappears across it.
The determinant utility has three targeted passing tests; roots are candidates
requiring adjacency/geometry verification, not a general production event solver.

QS independently recomputes a midpoint step from that rejected source with
dt=1.34180051102e-7, the predicted time to the event. Source RHS matches
exactly. Seventeen saved path nodes, with the unchanged full10period budget,
give mass ratio3.00229987e-5, quadrature ratio1.33956172e-5, scaled rate
6.02160761e-10, generator defect8.17679258e-14, positive minimum pressure
.69420346 and stage mass residual3.68802126e-15. This is promising single-step
evidence for event-aware proposals; it does not establish the next step across
the event, automatic event selection, or long-time physical accuracy.
PZ's accepted trajectory and all running sources remain unchanged.

At13:44–13:45, saved checkpoints: JU t=.17535210 (1219 steps),
MV t=.09394635 (2532), free-affine QB t=.03713651 (464), PQ32
t=.01120562 (150), full-horizon PZ t=.0017220511 (96), NX t=.015 (15),
OR t=.008 (4), QK t=.001 (1), auxiliary OT t=.01375. QD10 has no first
checkpoint yet but its actual Python process1936379 was verified live.
QK's process1947710 is also live and has completed its first step.
Canonical PK has integrated28096/36864 cells; independent QE has audited
15328/36864. These are partial counts, not official Kepler evolution.

OR's independent audit exited0 at13:35:37 and covers t=.008: absolute
relative mass1.20325057e-6, nested quadrature difference8.87151916e-10,
duration-scaled rate1.05332036e-5, generator defect1.48474610e-5.
NX's prior audit exited0 at13:14:43 through t=.012. After verifying no
audit process remained and verifying SHA256SUMS, resumed its existing
node caches at13:44 (session66768, timeout1h, execution_resume_1344.txt)
to audit newly accepted history. No fluid evolution was replayed.
All four full validation targets remain incomplete.

QT starts a diagnostic follow-up from QS's saved endpoint at twice its dt,
with the same17-node/full-horizon checks. It does not alter PZ. The first
launch exited1 before any RHS evaluation because the audit wrapper passed
the trial wrapper rather than its record to restore(). Corrected only that
loader; the failed log is retained. Retry session35606 uses execution_retry.txt
and run_retry.log with a20minute bound. Actual geometry crossing must still
be checked on the recomputed proposals: the old rejected polynomial's event
time need not equal that on a new midpoint polynomial.

QT retry exited0: dt2.683601022e-7, mass ratio1.03882009e-4,
quadrature ratio2.20509419e-5, rate2.09127090e-9, generator3.27460281e-13,
minimum pressure.69420346, stage mass3.68797213e-15. It satisfies the same
sampled per-step criteria. QU (session43861, timeout10m) now checks the
actual686–719 face at phases0,.5,1 on both recomputed dense paths.
Neither proposal is installed as a production checkpoint.

QU exited0 and confirms actual face length2.627594e-8→1.313797e-8→0
across QS, then0 at every inspected QT phase. The recomputed path really
passes the local topology event. This does not yet show a lower cost over
the original rejected interval: the first two small steps cover only part
of it. QV (session52454, timeout20m) therefore proposes the remaining
interval from QT to the exact original rejection17 endpoint. All budgets
remain unchanged; accepted production histories are untouched.

QV exited0 and passes the remaining dt2.902138555e-6: mass ratio.01235587,
quadrature ratio.0005909053, scaled rate2.47547582e-7, generator defect
3.83741927e-11, minimum pressure.69420347, stage mass3.68789421e-15.
QW verifies exact state/factor continuity across all three saved proposals,
checks the final time against rejection17, and sums absolute path and
quadrature estimates over the same original interval (session1875).
This is local evidence for an event-aware trial controller, not proof of
physical field convergence or of full10period vortex accuracy.

QW exited0: exact array continuity confirmed, final time.0011010594390064726
matches rejection17. Aggregate absolute relative mass3.42838245e-14,
allocated whole-interval mass budget ratio.010860468 and summed quadrature
budget ratio.0005191195. All three sampled step gates pass unchanged.
This verifies one local event split across an identical time interval.
Next: adjacency-qualified automatic proposals and additional real events
before production integration; retain all original full validation targets.
NX independent audit is processing step15 nodes at13:51, retaining earlier
validated node caches.


### 13:54 Automatic interior-event candidates and canonical native quadrature

QX adds an endpoint adjacency comparison to power_flip_time.py. For disappearing
interior edges with exactly two common neighbors it proposes isolated weighted
circle roots; ambiguous adjacency and persistent degeneracy are reported.
Five targeted tests pass. Boundary events and transient flips that restore the
original adjacency are outside this detector's coverage; it is not a complete
event-aware integrator and is not installed into running physics snapshots.

Actual PZ archived rejections7/8 yield edge491–522 at phases
.579374599/.105832536, both verified against actual geometry just before/after.
Rejection17 automatically recovers686–719 at phase.04060305488.
Rejection1 has no detected disappearing interior edge: topology handling alone
cannot explain or eliminate all early rejections. QY (session11336,20minute
bound) tests an unchanged-budget proposal from rejection8 to its automatically
predicted event; source RHS and original trial hashes are checked.

NX independent audit exited0 at13:53:30, covering the checkpoint captured
at launch rather than later fluid progress; its report is retained in JSON.
OR fluid reached t=.01; its earlier audit exited0 at13:35:37 through .008.
After checking no Python process occupied its audit directory and verifying
the source manifest, resumed existing caches (session52645,1h bound,
execution_resume_1353.txt). Fluid evolution was not repeated.

PK canonical IC integration reached32096/36864 cells and QE order16 audited
19488/36864 at the latest checkpoint observation. QZ now audits the same
original cached cells with density order6 (session51718,2h bound). It shares
only immutable physics source and original IC chunks, with a separate audit
cache/identity. This addresses the gap between higher-order diagnostic density
quadrature and native reconstruction; it neither replaces the physical IC
nor certifies topology-wide closure, potential enrichment, or100fluid orbits.

QY exited0: automatically identified event8 proposal dt1.932496893e-6,
source RHS exactly matches archived RHS; mass ratio.0050994115, quadrature
ratio.0002272658, rate1.02673264e-7, generator1.70251591e-11, positive
pressure.69420230, stage mass3.71502383e-15. RA continues this diagnostic
endpoint through the remaining1.632745378e-5 to the original rejected end:
mass ratio.390885215, quadrature ratio.018031286, rate7.82975555e-6,
generator1.21504262e-9, positive pressure.69420233. Both pass the same
sampled criteria. RB (session47303) checks actual recomputed face491–522
geometry and exact state continuity plus aggregate quadrature/mass over
this complete interval. No production states or controller settings changed.

RB exited0. QY/RA span exactly the original rejection8 interval with exact
state/factor continuity; aggregate absolute relative mass6.10565618e-12,
allocated mass budget ratio.350042893, summed quadrature budget ratio
.0161470414. Actual face491–522 length decreases3.772237e-7→1.886127e-7
→9.866766e-13 over QY, then disappears during RA. Thus the predicted time
is close but is not an exact root for the recomputed midpoint path. No
face clipping or forced geometry correction was introduced. A production
controller must treat the event time as a proposal and retain every gate.

JU's bounded execution ended13:58:21 with exit124; original PID1842321
was verified absent. Captured source hashes in resume_1358_receipt.json
at t=.17799820723,1254 steps, including pending trial/node caches.
Verified immutable source manifest and resumed existing run with a2h
bound (session52241, execution_resume_1358.txt/run_resume_1358.log).
This preserves the existing diagnostic continuation toward .2; its known
physical density and full10period accounting limitations still apply.


### 14:04 Canonical initialization complete, spatial 3D evidence, event retry control

PK exited0 at14:02:07 after integrating all36864 original canonical cells.
This completes only physical IC moments. QE order16 and QZ order6 per-cell
audits continue. RG (session9706,4h bound) now evaluates the retained
canonical initial state with whole-neighborhood quadratic projected velocity
and degree4 pressure, density quadrature16 and independent field integration
orders16/24. It reuses original moments and does not activate enriched
potential moments or evolve the disk. Official100fluid orbits remain unmet.

RC compares stored6³/8³ coarse endpoints at identical t=.01: order12
rho L1 .002317925→.001672361, velocity/amplitude L1 .040096018→.022879519,
pressure/amplitude² L1 .049811149→.023674951. All three improve spatially
at this time. Maximum sampled vertical velocity .001418047→.001820934
does not improve, and this is not a finalT=.02 convergence certificate.
The independently normalized mass-weighted velocity relative L1 is still
.2441105/.1392939; these differ from the fixed-amplitude acceptance norm.

RD adds actual-geometry verification to an optional rejected-step proposal.
Only isolated interior candidates contracting dt by at least10% and above
the existing minimum step are used. Otherwise ordinary halving remains.
Material-mass-lift geometry is explicitly unsupported. Retry proposal details
are archived with the rejected trial; every regular sampled mass, quadrature,
generator, pressure, factor, stage and accepted-state gate remains in force.
RD reproduces actual event proposals7/8/17 and falls back to halving for1.
Six targeted tests pass. RE immutable source snapshot passes505 tests with
one pre-existing Sod logarithm warning (exit0,14:02:37).

RF starts the first64² evolution from the saved QO original initial arrays,
with unchanged radial velocity fit and full10period budgets, and opts into
this retry proposal (session60970,4h bound). Execution endpoint .01 is the
first spatial/temporal physical gate; validation_end remains10.468597915605091.
Existing32² and accepted diagnostic histories are preserved. RF source is
the verified RE manifest, and configuration hashes bind the reused QO IC.

At14:04:52 actual Python processes RF1996150 and RG1995738 were verified
live with accumulating CPU time; neither had a terminal receipt. Empty
early progress logs are not treated as failures or grounds for restart.


### 14:08 Three-dimensional paths and shock generator mismatch

QD10 first actual step reached t=.002, and QK8 fine reached .002 with2steps.
Started independent order6/17-node history audits in NT/qd_n10_order6
(session75548) and NT/qk_n8fine_order6(session99357), each2h, using identical
library hashes and preserved actual path records. OR audit exited0 at14:01:07
through .01: absolute relative mass1.48493010e-6, quadrature8.87167992e-10,
scaled rate1.05332036e-5, generator1.48474610e-5. This is a partial interval.
RH (session79023) scores existing6³/8³/10³ coarse states at the same .002.

RI reads QB free-affine shock's actual rejections102–107 aroundt=.03731.
Their dt~5e-7, mass ratios.0008–.00164, and quadrature ratios.0021–.00495
are small; generator path defects1.25–1.77e-5 exceed the1e-5 gate and
dominate rejection. Stage mass residuals~1.9e-15 remain small. Worst mass
flux cell249 is not necessarily the worst generator-defect cell. RJ separately
locates the actual generator defect and affine/quadratic contributions on
rejection106 at five saved-path phases, with no modified path or thresholds.

RF first rejection archives a verified event proposal for face2196–2259
at phase.01454311178, retrydt3.635777944e-6. This confirms the optional
controller executed in the new64² run, not that its retry is accepted yet.

RH exited0 at14:07:36. At equal t=.002,6³/8³/10³ order12 field errors:
rho .000473037/.000343816/.000281365, velocity/amplitude
.038149375/.022065751/.013889374, pressure/amplitude²
.046082011/.020825988/.013078872. Spatial refinement improves all three
global norms at this early time; sampled vertical velocity maxima
.000307511/.000405582/.000368127 remain nonmonotonic. FullT=.02 and
independent history audits are still required.

RJ confirms actual worst generator mismatch at cell248, x-component,
from phase.25 onward, reaching1.77364978e-5. Its affine velocity changes
smoothly while the quadratic correction changes .00519827→.00519158.
Cell area~.00174174 and spatial covariance eigenvalues~2.58e-8/1.59e-7
remain positive. This is not evidence of a collapsed cell or failed mass root.
Source-phase mismatch is only5.93e-10. RK (session39379) observes the
unchanged batch projection inputs/outputs to distinguish unconstrained-fit
variation from covariance-budget effects; no model substitution or steering.

RK exited0. At fixed cell248 the quadratic Gram factor condition is only
5.47793→5.47799, and the dominant budget singular value remains saturated
at1−1e-12 across all five phases. The transverse curvature variation is
already present before budget projection: unconstrained x-velocity yy
coefficient3.53755e-5→1.90824e-5→1.99756e-5→3.59796e-5→6.76430e-5.
It is therefore not explained by a newly switching covariance cap. The
residual factor is highly anisotropic (~1.1567e-4 and3.629e-14 singular
values). Next diagnostic: conditioning and data sensitivity of the augmented
free-affine curvature fit, while preserving the accepted shock trajectory.

RF has now actually accepted three steps through t=1.409704574e-5,
with ratios.61149245/.08360135/.18815641. Its first automatically event-sized
retry passed the unchanged gates; this is initial progress, not a long-time
physical certificate.


### 14:16 Shock fit conditioning and independent stencil trials

RL actual cell248 augmented fit uses12neighbors and has condition
94112→94365, with smallest singular value~2.05e-5. Across the saved
rejected step, matrix/rhs maximum changes4.22e-6/6.53e-5. Holding either
matrix or rhs fixed still changes the weak transverse curvature; the
variation exists before covariance projection. All encoded fit inputs are
saved for further reproducibility and sensitivity analysis.

RM static alternatives retain the original moments and only change the
velocity fit's support: scalar support4 gives15neighbors, condition10.47,
and transverse x-velocity curvature~9.2e-10; bounding-cell width scaling
gives9neighbors, condition14.91, transverse curvature~1.52e-9. The
longitudinal curvature also changes (-.23033 baseline→-.13435/-.24199).
Thus these are alternative reconstructions, not an equivalent numerical
implementation. Their old-path generator defects are not acceptance tests.

RN independently recomputes source RHS, midpoint, factors and17dense nodes
from rejection106's retained source with dt1e-5, separately for support4
and cell_width (sessions25656/33150,20minute bounds). Pressure stencils
are explicitly left unchanged by restricting the interception to velocity
prepare. The full .15 shock budgets remain unchanged. These diagnostic
trials do not modify QB, overwrite its history, or certify full physical
accuracy of either new reconstruction.

NX fluid reached t=.019,19steps. Its prior independent audit exited0
through .015; after source verification and checking no Python process
occupied its audit directory, resumed existing node caches at14:16
(session34336,1h bound,execution_resume_1416.txt).

RN both trials completed. At dt1e-5 (19.26times rejection106's dt),
support4/cell_width path ratios are .06134479/.05001065 and generator
defects3.169575e-8/3.305882e-8, with positive pressure and factors. Stage
mass residuals1.914e-17/3.652e-15. These are genuine recomputed new-RHS
trials; source RHS differences.006343/.002583 explicitly confirm changed
reconstruction. All sampled numerical gates pass, but no long-interval
physical accuracy or convergence is inferred. Next: explicit optional
geometry-aware velocity stencil, manufactured compressed-cell tests,
and controlled full shock comparison before selecting any default.

MV previous bounded execution ended14:17:44,exit124; original PID1865974
was verified absent. Saved resume_1418_receipt.json at t=.10628814957,
2903steps with checkpoint/pending/configuration hashes, verified immutable
source manifest, and resumed the same run for2h(session22802). No accepted
step was replayed or erased. Original T=.15 remains the shock target.


### 14:21 Optional cell-width velocity stencil

Implemented explicit velocity_stencil=cell_width, default area unchanged.
Only tensor velocity fits use per-axis polygon bounding widths for neighbor
support. Polynomial basis, fixed means/moments and covariance budget remain
the original ones; pressure stencils are untouched. This is an axis-dependent
reconstruction option, not a rotation-invariance claim or geometry steering.
Environment/configuration identity checks and cached-source RHS fit/stencil
checks prevent importing a RHS from the other reconstruction.

Added tests for bitwise default support preservation, compressed-wall
manufactured planar curvature with a well-conditioned stencil, and invalid
configuration/widths. RO immutable snapshot is running regression and a
real saved-state comparison against RN's isolated cell-width implementation
(session56109). RP16/32 fullT=.15 controls are prepared from retained
original physical IC arrays; launch awaits those checks. These are distinct
reconstruction comparisons, not restarts of existing accepted trajectories.

RO exited0 at14:21:32:508 regression tests pass, with one existing Sod
initialization logarithm warning. Actual source RHS equals the RN cell-width
prototype exactly (maximum difference0), stage mass3.61955774e-15.
After manifest verification, launched RP16/32 fullT=.15 comparisons
at14:22 (sessions98950/99228,4h bounds). Each reuses original initial arrays
with explicit source hashes and velocity_stencil=cell_width; original
area-stencil runs continue separately. No event-controller change was added
to these controls, isolating the velocity stencil as the new option.

RQ compares n8 at exactlyt=.002 using the existing independently integrated
17-node/order6 caches and identical IC arrays. Halvingdt .002→.001
reduces absolute mass defect3.15484504e-7→7.82789952e-8, ratio4.03025745.
Generator defect1.48474610e-5→3.71289912e-6 also decreases approximately
fourfold. This supports temporal consistency at the shared early time,
not fullT=.02 path certification or final physical accuracy.

RP Python processes2019178/2019179 verified live at14:23:50. n16
has11accepted steps through .00037438851; n32 is reducing its initial
step under the unchanged gates. These rejections are not terminal failures.


### 14:27 Full6³ temporal run completed, pressure criterion still fails

NX dt=.001 reachedT=.02 with20stored actual paths. Final order12
rho/v/P=.004549458/.042642928/.055554603, versus coarse pressure
.055557249. Pressure remains above .05 for both time resolutions.
Mass arrays are exact, stage mass peak1.95008248e-13, energy ledger
error2.22044605e-16. This supports a spatial/reconstruction error limitation
at6³, not full physical validation. NX independent audit is still running
on the earlier captured prefix and must then include the final path.

NT10 first-step audit exited0 at14:22:52 through .002: absolute mass
6.51135394e-7, quadrature1.00709624e-9, scaled rate2.53220234e-5,
generator1.18101985e-5. Its mass defect exceeds8³ at the same dt.
RR10 time refinement (.001, target.02, session75520,4h bound) starts
from exact original QD10 arrays and score, with unchanged immutable source;
this separates temporal error from spatial physical-field improvement.
RS(session68750) independently scores the already saved8³ fine endpoint
at.002, reusing the coarse score only after matching its record hash.

PA32's2h control ended14:24:36 with exit124; original PID1875066
was verified absent. Captured checkpoint/pending hashes at t=.01859280656,
327steps, and resumed with unchanged source for1h(session71445).
This preserves the existing short comparison toward.02.

QE exited0 at14:26:39 and covers all36864 canonical original cells at
density order16. Maximum independently integrated relative mass error
6.43172734e-12,16/24 integration difference4.78931309e-14, minimum
decoded cell pressure2.77777771e-5. This completes that per-cell IC audit,
not native order6, global velocity/pressure closure, potential enrichment,
or100fluid orbital periods. RG and QZ still address separate requirements.

NX terminal receipt confirms exit0 at14:24:35. RS n8 same-time .002
field comparison exited0 at14:26:40: coarse/fine pressure
.02082598765/.02082560188, velocity .02206575087/.02206581906,
rho .00034381572/.00034383361. Time refinement changes these physical
errors very little while independent mass-path error improves fourfold.
Raw endpoint factor differences are not directly a covariance error:
RS supplementary audit(session79361) compares LLᵀ and signed column
permutations explicitly without modifying factors.

RS supplementary covariance audit completed: maxLLᵀ difference3.290922e-12;
after cellwise signed column permutations, maxfactor difference2.027717e-8.
The rawfactor-coordinate difference.00177958 therefore should not be
interpreted as a comparably large covariance error. No factor was changed.
RR10 fine Python PID2023529 verified live at14:29:49.

NX independent prefix audit exited0 at14:29:13 through.019,19steps,
absolute relative mass4.96954939e-7. After terminal/process/source checks,
resumed the retained cache to include finalstep20 (execution_resume_1430.txt,
1h bound). Final physical pressure failure remains separate from this
numerical path audit.


### 14:34 Native6 log incident, verified cache evidence, canonical initial fields

A status-refresh helper reused its documentation-path variable for the last
progress log and wrote documentation JSON into QZ/n192/run.log. This was
an assistant bookkeeping error. That log is contaminated and excluded from
progress evidence. RT preserves a hashed snapshot and incident record.
Python PID1984343 was verified live; no physics/cache files were rewritten
or numerical process restarted. RT independently verified source/geometry/
script identities, original chunk hashes and exact cell-index sequences for
18752 completed native6 cells. Their max mass error1.71612516e-5, quadrature
difference1.32244545e-14, minP2.77777771e-5. This is a partial prefix and
is markedly less accurate in independent mass integration than order16.
The documentation audit was refreshed using separate explicit output paths.

RG initial global order16 integral completed: relative densityL1 .0059349155,
mass-weighted velocityL1 1.52765263e-6, pressureL1 5.97968660e-10.
Minimum density/pressure are positive, but sampled max relative density
error36.96074 is large. Reference density jumps1↔.01 at radii.5 and2.
RU(session83845,1h bound) independently measures the density error in
polygons crossing those interfaces at orders16/24; it does not assume
that all global error is confined to them. RG order24 remains separate.
Canonical fluid100-orbit evolution has not started.

NX independent history audit exited0 at14:33:33 and covers all20steps
throughT=.02. Absolute relative mass5.17933943e-7, quadrature difference
1.15839977e-9, scaled rate3.87202949e-6, generator4.31441416e-6.
Coarse/fine full-interval absolute mass ratio=4.00842479103773.
Both physical pressure errors still exceed .05. This completes the6³
temporal path comparison, not the full three-dimensional accuracy goal.

RU exited0 at14:34:57:503polygons cross the actual r=.5 or2 density
interfaces. Their order16 densityL1 accounts for .999999697 of the global
order16 density error. Cut-cell normalized L1 is .00593491368(order16)
and .00592904903(order24). Worst sampled outer-interface density is
~.38 against reference .01, immediately outside r=2. This local smearing
is not a mass-integral failure, and its later transport remains unverified.

RV(session78661,4h bound) begins a separate canonical initial per-cell
potential-enriched density audit. It preserves the original m,q,Q and
reference potential moments, uses order16 base density, and independently
checks mass/potential integrals at16/24. Frozen cells retain the original
unenriched closure. Per-chunk identity caches and failure-cell receipts
are retained; no target moment is adjusted on failure. This is preparation
for a coupled canonical step, not a claim of fluid evolution or100orbits.

RV has produced352identity-bound initial-cell audits by14:38:09 and
its Python process is live. This remains an incomplete prefix.


### 14:42 Equal-budget shock field comparison and auxiliary angular ledger

RW exited0. At t=.02 with identical original IC and full .15 budgets,
area/cell_width16² rhoL1=.01330303/.01330646, vL1=.02960721/.02956488,
PL1=.02229074/.02229181. Both meet the .05 early physical criterion;
cell_width is not uniformly more accurate. Max meanvy1.0874e-11/8.5900e-12
retains planar symmetry. Dense accepted paths were sampled, not replayed.
The fullT=.15 and32² comparisons remain in progress.

OR fluid reached .014; prior independent audit exited0 through .01.
Verified terminal/process/source state and resumed its cached nodes
(session78662,1h bound,execution_resume_1440.txt) for the new paths.

RX(session31093) checks accepted auxiliary Kepler angular momentum using
Dxy−Dyx about the physical center, separating centroid orbital and intrinsic
parts and subtracting the frozen-reservoir angular/momentum ledger.
This avoids equating generator orbital angular momentum with total fluid
angular momentum. It does not establish the canonical100-orbit result.

RX completed on23accepted auxiliary steps from .01 to .01575. Fluid
angular change−1.47717216e-9 agrees with boundary exchange−1.47717139e-9;
residual−7.76569653e-16, relative5.59052179e-17. Centroid orbital and
intrinsic components change separately, reinforcing the need to retain
both. This is conservation evidence for the captured auxiliary interval,
not a physical-orbit accuracy or canonical-duration certificate.
