# Astra: 단일 target 출판 종료 평가와 all-target registry 계약

2026-10-21. 계획·평가 Astra, 구현·실행 driver. SHARED_CONTEXT, missing-source exact 회계 장부, conditional radial 코드, potential P2 proof, missing-source 조립 및 독립 audit 코드를 정적으로 읽었다. v8 frozen radial header와 현재 파일의 diff도 확인했다. 원본 origin LagEunha/HEAD 8e25dfe 확인. 테스트·계산·폴링은 하지 않았다.

## 1. 종료 판정

**이번 target68526의 CPU reference transactional 출판을 통과로 인정한다.** 기존 원본12tau와 최종 R80 변환을 포함한 potential .712959619tau, kinetic .248198055tau, volume .102901761tau 등 모든 항이 통과했다. 질량·물리U 양성, 원본 state/heat SHA74062... 보존, 기존25 source 값/forward 및 다른11 remainder 불변이 독립 저장자료 회계로 확인됐다. 새로운 계산은 source301의2816점뿐이며 총268032 실제 점이다.

이는 mixed arithmetic CPU reference의 한 target 결과다. 독립 audit는 저장된 숫자·오차·직렬화 회계를 검증하며 PDE나 AD2 증명을 독립적으로 다시 계산하지 않는다. libm8ulp 등의 명시된 산술 가정은 그대로다. accepted RK4, production CUDA transfer, 일반 evolved heat 수송 및 장시간 Kepler는 미검증이다.

## 2. 구현·수학 점검

읽은 범위에서 출판 판정을 뒤집는 치명적 결함은 발견하지 못했다.

- v8 대비 radial 변경은 zero-containing s enclosure의 조건부 division-free recurrence와 계측에 한정된다. 분리된 기존 경로는 유지됐다. recurrence는 같은 binomial 함수와 미분을 나타내며 tail은 1-sh만으로 제한한다. th가 마지막 c_N*s^N의 전체 ball을 감싸므로 그 ball의 최대 크기를 사용한 무한 tail 상계는 zero-crossing에서도 유효하다. nominal 조기종료 검사만으로 정확도를 주장하지 않고 마지막 infinite-tail radius가 남는 것이 필수다.
- P2 proof는 full normalized potential integrand에서 density, map 및 measure를 함께 미분한다. rho0/J와 physical J를 대수적으로 소거하며, 원래 clipped graph의 implicit derivatives를 사용한다. 중심에서 정의한 고정 polynomial은 Gauss16에 정확하므로 저장 quadrature 값은 바뀌지 않는다. convex complex domain 전체 Hessian 차이를 적분한 상계와 계수 오차를 유지한다. center AD2 또는 표본의 차이를 analytic remainder로 대체하지 않았다.
- OMP panel reduction을 saved-moment serial tree와 비교하고 R80 n/e의 유효10byte payload를 확인했다. padding byte 일치와 혼동하지 않는다. point별 오류와 실제 합산오류를 모두 가진 결과가 최종 Q128 조립/직렬화 관문으로 연결된다.
- whole4 group은 한 번, clipped26 source는 각각 한 번 사용된다. P2에서는 유효한 old/new bound의 최소만 선택하며, 원래 값과 다른11 모멘트는 유지됐다.

추가 local certificate/정밀도/차수 조사는 이번 종료 후 반복하지 않는다. 다음은 단일 target을 더 다듬는 일이 아니라 유한 source patch 전체의 소유권·보존 전달 구조를 등록하는 단계다.

## 3. 다음 단일 묶음: all24 x384 registry/incidence/cost ONLY

전체 quadrature, weight Newton, compression, 시간진화를 시작하지 않는다. source384와 실제 target24의 frozen registry를 사용하며 네 먼 사이트도 삭제하거나 임의 외곽 filler로 바꾸지 않는다.

### A. 정의할 정확한 domain

    Omega_s = chi_s(T_s)
    Omega = union_s Omega_s
    C_j = {x: |x-q_j|^2-W_j <= |x-q_k|^2-W_k for all actual k}
    D_sj = Omega_s intersect C_j

모든 실제24 power site가 유효하면 C_j는 평면 전체를 분할한다. 그러므로 모든 source 물질은 어떤 target으로 분류된다. target C_j가 unbounded이어도 D_sj는 유한 source 이미지 안에 있으므로 유한하다. 임의 bounding box로 power cell을 자르지 않는다.

반면 C_j 밖이 아니라 **Omega 밖의 target 부분은 입력 물질이 없는 영역**이다. C_j 전체를 채운다는 뜻으로 밀도·속도·열을 연장하지 않는다. 출력 domain은 D_j=C_j intersect Omega이며, full_target_coverage와 partial_patch_coverage를 명확히 구별한다. bounded target이라는 사실만으로 C_j가 Omega에 포함됐다고 간주하지 않는다. 빈 target은 empty-domain record로 두고 0질량 셀의 EOS를 계산하지 않는다.

모든24 target을 등록하면 일반적으로 '어느 target에도 속하지 않은 source 질량'은 없어야 한다. external은 우선 **source support boundary의 기하 역할**이다. 나중에 물리 evolution 대상 target을 subset으로 정한다면 제외된 실제 target의 D_sj를 명시적 external recipient로 라우팅할 수 있지만, internal target과 external mass bin에 같은 조각을 중복 계상하지 않는다. 이번 registry에서 수신자 집합을 몰래 바꾸지 않는다.

### B. exact partition과 유일 소유권

- 기존 exact power-plane coefficients 및 canonical curved boundary를 사용한다. power tie는 (power value, stable global target ID)의 lexicographic rule로 한 번 소유한다. 곡선상의 부호가 미확정이면 ID 규칙으로 때우지 않고 unresolved로 남긴다. 정확 tie와 roundoff uncertainty를 구별한다.
- canonical source edge는 최대2 adjacent source와 고정 방향을 가진다. 한쪽만 있는 edge는 외부 support boundary다. 각 외부 edge 조각도 실제 target 소유자·방향·endpoint/event ID를 가진다. 네 먼 사이트로 생긴 unbounded cell과 source patch의 열린 외곽은 별개 개념이다.
- source 이미지들이 서로 겹쳐 물질을 이중 계상하지 않는다는 증거가 필요하다. 이웃은 공유 analytic trace·반대 방향·양의J를 확인하고, 비이웃은 기존 bbox pruning 및 필요한 실제 경계/containment 증거로 interior disjoint를 확인한다. local J>0만으로 global injectivity라고 선언하지 않는다. 이미 존재하는 증거는 그대로 재사용한다.
- source별 target 분할의 union이 원래 Omega_s이고 interior overlap이 없음을 oriented boundary cancellation/containment로 등록한다. 대략적인 면적합만으로 누락·중복을 상쇄하지 않는다. source→target CSR와 target→source CSR는 동일한 overlap ID의 transpose여야 한다.

### C. registry 행과 source/target 장부 계약

각 (epoch, sourceID, targetID, componentID)에 대해 최소한 다음을 저장한다:

1. source/target/state/heat/geometry SHA 및 원래 사이트·weight 값.
2. excluded/empty/whole/intersection/unresolved, bounded_target 및 full/partial patch coverage.
3. canonical boundary/vertex/event ID, 방향과 실제 support-exterior 역할.
4. 계산 owner rank와 immutable overlap ID, 향후 packet ABI 및 기존값 재사용 여부.
5. 증명 완료/미완료 범위, 예상 chart/panel/node 수와 storage 상계, target normalizer 출처/미확정 상태.

향후 Q_sj를 계산하면 동일 packet과 동일 absolute-error 객체를 source outgoing 및 target incoming에 연결한다. 두 숫자를 따로 적분해서 일치한다고 주장하지 않는다. 모든12 모멘트에 적용하며, 특히 K와 physicalU를 분리 유지하고 E=K+U 장부도 같은 두 항에서 만든다. registry 단계는 이 연결의 구조적 항등식만 검증한다. 아직 없는 숫자 적분의 물리보존 pass를 출력하지 않는다.

source68526 관련 출판 packet은 target 전체 aggregate다. 내부 source별 기여를 registry와 연결할 때 retained whole4 group이 합계 packet 하나라는 사실을 유지한다. aggregate를4번 사용하거나 각 source에 임의 분배하지 않는다. 이미 존재하는 whole-source별 적분 기록을 사용할 수 있으면 그 provenance와 group 합 일치를 확인한다. 없다면 grouped ownership을 명시하고 source별 수치 보존 검증은 pending으로 남긴다. target68526 최종값은 재적분하지 않는다.

### D. 비용과 실행 전 관문

9216 source-target 조합은 registry의 최대 조합 수이지9216개의 적분을 실행하라는 뜻이 아니다. exact/interval broad-phase pruning, canonical plane-pair 재사용, 전체 포함·공집합을 먼저 분리한다. 새 topology 증명은 실제 필요한 교차에만 수행한다. 목표는 source별 incidence와 외부 support 조각까지 완결하는 것이다.

비용표에는 target별 실제/예상 overlap·chart·panel·positive node 수, analytic proof 수, 출력 bytes, 재사용 절감량을 구분한다. 기존268032점은 이미 계산된 증거로 표시하고 다시 비용을 부과하지 않는다. 알려지지 않은 chart 비용을0으로 놓지 않는다. 기존 한 target의 실행시간을 GPU solver 속도나 전 target에 대한 확정 walltime으로 쓰지 않는다.

원래 target68526 tau/normalizer는 변경하지 않는다. 다른 target의 cold 기준은 동일한 기존 물리 계약에 따라 출처와 양의 U/m/domain 하한이 필요하다. 이번 registry에서 하한이 없으면 pending으로 명시한다. empty 또는 partial cell에 원래 full-cell U를 분모로 임의 대입하지 않는다. registry_complete와 all_target_quadrature_ready를 별도 flag로 둔다.

## 4. 이번 registry 묶음의 종료 조건

성공은 다음 네 가지다: (i) 실제24site/384source의 모든 조합이 증명된 배제 또는 유일한 overlap record로 귀결됨, (ii) source partition·shared/external boundary·inverse CSR ownership이 완결됨, (iii) target68526 증거가 정확히 한 번 재사용됨, (iv) 나머지 전체 적분의 비용과 미완료 normalizer/증명 요구가 숨김없이 나열됨.

unresolved topology, source-image overlap, registry에 없는 material recipient, 중복 ownership, 임의 filler 요구가 나오면 registry를 미완료로 저장하고 full quadrature로 넘어가지 않는다. 전체 적분 예산을 자동으로24배 늘리지 않는다. 정상 완료 후에는 이 비용/관문을 근거로 다음 whole-remap native 조립 묶음을 정한다.

이번 출판은 실제 진전이지만 cold CONTROL의 한 CPU target이다. 일반 evolved heat 및 checkpoint 계약, 실제 CUDA backend 전달, 시간법/충격 검증, t6/8/10과 최장 공전은 여전히 남아 있다. no steering, Re=infinity, 물리 U/분해된 회전 에너지 보존 및 조건부 W0 질량 환원은 유지한다.
