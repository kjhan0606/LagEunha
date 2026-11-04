# 면 기울기 수정과 실제 3D 연결

같은 W0 압력·면 기하 LSQ·중앙 점성·BE·KDK를 사용한 기본
1D/2D/3D 실제 진화를 검증했다. 전체 2D suite, 새 방법의 장기
Kepler, 다중 랭크/GPU 3D 및 GIZMO 대비 실행시간 2배 이하 조건은
미완료다. 기존 20pi 궤적과 이전 fit의 열화량은 새 방법의 결과가 아니다.
원고의 기존 결과도 실험 중인 공통 후보의 결과로 대체하지 않았다.

## 2D 축 shock의 원인과 수정

단계별 N=64 native Sod 재생에서 t=0 압력 RHS와 첫 pressure kick은
횡성분이 정확히 0이었다. 두 번째 step의 두 번째 점성 갱신에서 동일
x-column의 vx 산포가 약 0.00105로 커졌다. 열 산포 뒤 pressure RHS가
횡속도를 만들었다. 첫 원인은 초기 이미지 상태나 h=V^(1/D) 차이가 아니다.

기존 LSQ는 면적 약 0인 모서리 이웃에도 유한한 방향 비중을 주었다.
첫 fit의 가짜 Gxy는 3.04256, 다음 fit의 Gxx column 산포는 0.202772다.
같은 상태에서 면 기하 비중 A/r^(D-1)을 사용하면 Gxy는 2.53e-16,
column 산포는 0이다. 인과 자료는 `planar_symmetry_trace_v1` 및
`planar_symmetry_frozen_fit_v1`에 보존했다.

공유 `gfs_velocity_fit_add_face`를 1D, native CPU/CUDA 2D 및 CPU 3D에서
사용한다. 최소제곱 목적식은 다음과 같다.

\[
 \sum_f\frac{A_f}{r_f^{D+1}}
 \left|\Delta\boldsymbol v_f-G_i\boldsymbol r_f\right|^2.
\]

기존 방향 정규화에 A/r^(D-1)을 곱한 형태다. 1D의 단위 면적은
기존 fit과 정확히 같다. 비중은 정확한 polygon solid angle이 아니라
차원 없는 angular proxy다. Laguerre 가중치 W나 새로운 물리 자유도가
아니다. 새 조절 계수·cutoff·복원력은 없다. 소멸 면의 영향도 사라지며
affine 정확성·회전 공변성·Galilean invariance는 유지된다. 일반 LSQ
절단오차나 가까운 생성점의 큰 비중까지 없애는 방법은 아니다.

C/CUDA 세 차원 affine 감사가 통과했다. 소멸 면 대조의 가짜 횡기울기는
기존 비중에서 1/6, 새 비중에서 1.77e-17이다. 1D N=200의 9종 최종장은
이전과 byte-identical이다. 실제 2D 14개 계산도 지정 시각에 도달했다.
contact 여섯 경우는 저장 정밀도에서 평형과 이류를 보존했고 정상
사인 전단 N=32/64/128도 수렴했다. 축 Sod t=0.2의 결과는 아래와 같다.

| 한 변의 셀 수 | L1 rho | 최대 횡속도 | 상대 총에너지 변화 |
|---:|---:|---:|---:|
| 64 | 0.0314707 | 4.38e-16 | 4.65009e-4 |
| 128 | 0.0186933 | 9.80e-16 | 2.52680e-4 |
| 256 | 0.00997153 | 3.50e-12 | 1.29377e-4 |

이전 최대 횡속도 0.07087/0.04285는 새 N=64/128 진화에서 제거됐다.
충격파의 유한 해상도·시간간격 오차는 남는다. shock norm은 현재
생성점에서 표본화한 해석해에 대한 체적 가중 평균이다. 정확한 셀
평균 해에 대한 적분 norm이나 방향별 회전 불변성 증명이 아니다.

## 실제 native 3D 연결과 기하 실패

원래 `exam3d_gfs.c`의 linked-cell·Voro3D·OpenMP·LAG3DV1 IO를 재사용했다.
새 `gfs_unified3d.inc`가 초기 primitive부터 종료까지 공통 연산자를
적용한다. legacy Riemann·RK4·floor·entropy repair·jitter·유한 면적
cutoff는 섞이지 않는다. 실패하면 마지막 accepted 상태를 복원하고
비영 코드로 종료한다. 3D는 CPU/OpenMP 한 랭크이며 GPU 검증은 아니다.

셀 질량은 고정이고 rho=m/V, P=(gamma-1)U/V, h=V^(1/3)이다.
3차원 curl의 크기를 같은 압축/회전 비율에 쓴다. 주기 graph는 세 방향
image shift와 절단면 각운동량 flux를 포함한다. BE 점성 활성 구간의
시간 정확도는 1차다. Python은 실행 관리·분석에만 사용한다.

`unified_initial3d.c`는 실제 cubic 셀에 질량·내부에너지를 적분한다.
대각선 stripe 부피는 균일 변수 합의 누적 부피식으로 구한다. smooth
전단의 속도는 생성점 표본이며 미해상 속도 분산을 열로 바꾸지 않는다.
IO의 specific energy u와 내부 상태 U=m*u를 구분한다.

최초 3D 전단·대각선 Sod는 첫 drift 뒤 상호 면적 차이 약 5--6e-12로
중단됐다. 면적 반올림 bound의 약 167--216배다. 실제 실패 좌표와 면
배열을 `unified3d_geometry_capture_v1`에 보존했다.
기존 Voro3D는 교차점 inset EPS=1e-9와 내외 판정 상대 허용치 1e-6을 쓴다.
새 명시적 `Voro3D_FindVCStrict`는 long-double 중간 계산, 유한 inset 없는
교차점과 double 반올림 척도의 plane 판정을 사용한다. 기존 API는
유지하고 선택 상태는 thread-local이다. float geometry는 거부한다.
임의의 퇴화 배치를 보장하는 exact predicate는 아니다.

동일 실패 좌표의 재생에서 상호 면적 차이는 6.94--7.63e-17로 줄었다.
면적 가중 centroid moment 불일치도 반올림 bound의 약 0.001배다.
좌표 jitter나 reciprocal gate 완화는 없다. legacy API의 기존 세 차원
기하 감사도 통과했다.

Astra가 지적한 무조건적인 endpoint clamp도 보완했다. 최종 구현은
선택할 endpoint의 plane residual이 반올림 bound 안일 때만 clamp한다.
범위 밖·비유한 값·비유일한 평행 교차는 실패한다. 직접 fixture 5개와
float 거부 검사가 통과했다. 실제 진화의 clamp 수와 최대 residual/bound도
기록했고 모두 1 이하였다.

실제 adapter의 체적 유한차분은 cubic cospherical 상태에서 perturbation을
줄이며 검사했다. topology 변화가 있는 큰 perturbation에서 오차는
대략 1차로 감소한다. 최종 strict 감사의 최소 오차는 h^2로 정규화해
약 1.06e-7이다. 유한차분 허용치 1e-6은 유지했다. 일반 유체해의 공간
차수를 뜻하지 않는다. 각 accepted drift의 면 폐합·체적 identity·
법선 방향·전역 체적·pressure-work adjoint 감사도 통과했다.

## 3D 진화 및 시간간격 검사

최종 bounded-endpoint 실행 파일로 실제 15개 계산이 지정 시각에 도달했다.
정지 contact와 (1,0.5,0.25)의 이동 contact를 8^3/16^3에서 수행했다.
이동 시험은 축·대각선 면과 세 축 주기 이류를 포함한다. 최대 압력·
속도 오차는 약 2e-14 이하이며 질량은 정확히 유지됐다.

3D 정상 전단 vx=sin(pi*(y+z)), vy=vz=0과 축·대각선 Sod는
8^3/16^3/32^3에서 수행했다. 전단 t=0.25의 L1 속도 오차는
0.0192397, 0.00228680, 0.000736808로 감소했다. 축 Sod t=0.2의
L1rho는 0.130857, 0.0812088, 0.0547981이고 대각선 t=0.15에서는
0.115752, 0.0716808, 0.0472503이다. 축 shock의 횡속도는 저장
정밀도에서 0, 대각선은 약 3.2e-15 이하다. 방향별 시각·법선 해상도가
달라 norm의 단순 방향 비교는 하지 않는다.

32^3 축 Sod에서 공간해상도와 t=0.2를 고정한 시간간격 검사다.

| 시간간격 배율 | 상대 총에너지 변화 | L1 rho |
|---:|---:|---:|
| 1 | 8.76003e-4 | 0.0547981 |
| 1/2 | 2.83656e-4 | 0.0544210 |
| 1/4 | 7.87067e-5 | 0.0544172 |

에너지 오차는 감소하지만 유한-step 정확 총에너지 보존은 아니다.
에너지 오차 감소율만으로 BE를 포함한 유체해의 시간 차수를 높여
부르지 않는다. coarse 3D shock의 약 1e-3 에너지 오차와 공간 오차도 남는다.

## 위치와 다음 묶음

개발 checkout은 `/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1`,
계산 원본은 `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`이다.
`Research/FaceFit3D_20261104` manifest에 소스·자료·실행 파일의 해시와
큰 원본의 GPFS 위치를 기록한다. 최종 3D binary는
`unified3d_strict_build_v3/unified_native_general.exe`다. 변경 unit의
재컴파일과 부분 재링크이며 전체 clean rebuild 검증은 아니다.

다음 묶음은 Gresho와 고정된 물리 shear 폭의 KH다. 초기 pressure RHS와
이후 소산을 구분하고 mode 성장·횡방향 운동에너지·혼합층 폭을 함께 본다.
새 fit의 Kepler t0 RHS만 갱신하며 이전 fit의 열량을 재사용하지 않는다.
이후 Sedov/Noh와 별도 중력·경계 계약이 필요한 RT/cylinder로 진행한다.

성능 변경은 step당 2D 새 기하 생성 한 번을 목표로 한다. geometry와
pressure state의 불변성을 구분하고 마지막 kick 뒤의 CFL·fit·점성
계수·active sign은 다시 계산한다. 재시작·재배열·image 변경은 cache를
무효화해야 한다. contact·전단·shock 회귀와 단독 profiling 뒤 동등한
GIZMO 비교를 한다. 작은 시험 실행시간은 2배 목표의 증거가 아니다.
