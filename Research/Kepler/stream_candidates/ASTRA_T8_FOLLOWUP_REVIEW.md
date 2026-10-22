# Astra 후속 묶음 평가 — 2026-10-20

운전자가 t8 해결 후보와 완료한 판별 실험을 `t8_followup_review_request.md`로 묶어 전달한 뒤 받은 평가. 코딩·실행은 운전자가 진행한다.

## 핵심 판정

최우선은 intrinsic spin 추가보다 **왜곡 셀의 선형 재구성·limiter 일관성**이다. 운전자의 입력/면→재구성→해상도/물질부피 순서에 동의한다. 중력 시간적분 단독 설명은 약해졌으며 초기 공간 연산자에서 직접 문제가 보인다. 이전 pressure alias 의혹은 철회한다(실제 쓰기는 후속 kernel B).

## 권고한 좁은 판별 묶음

1. affine 진단에서 rho/P 제한을 유지하고 velocity clipping만 해제한다. limiter 전 LSF gradient, clipping 계수, 양쪽의 동일한 실제 면점 재구성, 정상속도 jump, Pstar와 추가 열을 연결한다. 첫 force가 주 판정이며, 중력 없는 강체회전은 이후 강체회전을 유지하는 해가 아니다.
2. tiny-edge 병합을 별도의 geometry-only 대조로 분리한다. orphan, 면 벡터/작용점, constant-pressure 힘, pair work를 본다. float geometry의1e-8 초과 개수만으로 실패를 선언하지 않고 절대/flux 가중 잔차를 본다. 두 면 벡터를 단순 평균하여 셀 폐합을 깨는 수리는 피한다.
3. 그 뒤 불변 입력+halo+LSF+검증된 limiter+완전 final volume으로 cold 초기 RHS와 짧은 진화를 다시 평가한다. 이미 가열된 t6 재시작만으로 cold 안정성 개선을 판단하지 않는다.

## 수학적 설계 평가

이웃 중심 extrema 제한을 유지하면서 모든 왜곡 셀의 affine 정확성까지 요구할 수는 없다. 실제 면점이 그 표본 범위 밖이면 정확 선형장도 extrema를 초과한다. 따라서 affine 성분과 비선형 잔차를 구분하고, 좌표 회전 B→QBQᵀ 및 균일 boost에 불변인 판정을 설계한다. density/pressure 양성과 충격 처리는 별도 관문이다. 공통 velocity 제한계수 하나로 바꾸는 것만으로 충분하지 않다. 회전으로 보이면 무조건 limiter를 끄는 일반 규칙도 채택하지 않는다.

## spin / refinement / moment

기존 속도장을 정확히 면에 전달하지 못하는 오류를 새 spin 자유도로 보상하지 않는다. 기하와 work GCL이 맞아도 순간 Voronoi 셀이 material cell이 아닐 수 있다. 생성자 유체 이동, W0 Voronoi, 개별질량 고정, 임의 비압축성 affine 흐름에서 정확한 균일밀도를 모두 자동 만족한다고 가정하면 안 된다. 이 잔차가 지배적이면 mass/volume/work 또는 동적 Laguerre·moment 계약을 설계한다. spin 하나로 닫히지 않는다.

고정 eps에서 중심 refinement는 재구성·공유면 관문 뒤에도 곡률 잔차가 남을 때 비용 대비 좋은 다음 후보다. 같은 연속 IC의 실제 셀 적분 질량과 총질량, 전이부 이웃 완전성을 검사한다. 기존 기준선은 보존한다. 현재 결과는 장시간 성공의 증거가 아니다.

## 운전자의 후속 기록

평가 이후 예정했던 velocity-only/merge-only native 개입을 수행했다. 별도 진단 입력 v1의 KH trailer에서 이전 dt 슬롯에 Gamma를 잘못 기록한 운전자 오류를 발견했다. 기하·LSF·TVD 수식/첫 속도 gradient에는 영향이 없지만 pair restoring pressure의 dt cap에는 영향이 있으므로, v1 결과는 보존하고 `irregular_affine_probe_v2`에서 particle payload를 그대로 유지한 채 dt/dtold를 표준 cold 값1e-7로 교정해 다시 판정한다. 생산 baseline과 native restart reader는 변경하지 않았다. 후속 수치 결과의 주 근거는 v2를 사용한다.
