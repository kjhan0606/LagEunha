# Astra: 실제 signed ledger 이후 다음 묶음

2026-10-20. 완성된 장부·독립 Jacobian 검증과 과거10월9일 이력을 한 번 묶어 평가했다. 한국어 답신 요약이며 driver가 코딩한다. 요청 원문은 `astra_ledger_review_request.md`.

## 이번 판단

다음은 **압력비100 분기의 실제 작동과 U→0+ 경계 판별**이다. 호환 force-work 후보를 진화시키지 않는다. 냉각 위치는 특정됐으나, 열률 부호 전환이 불연속 분기인지 각 분기 내부의 공간 연산자인지 아직 분리되지 않았다.

확정 범위: 추가냉각은(Pi-Pstar)s_geom이고 pair/work-geometry 차이는 원인이 아니다. RK4 최종 조합이나 단열 팽창만의 결과도 아니다. 다만 Pstar≠Pi 또는 국소 entropy 감소만으로 모든Riemann 교환의 오류/엔트로피 부등식 위반을 단정하지 않는다. 허용 에너지·엔트로피 수송 계약 및 유효 열 영역의 경계를 확인한다.

## 제한 실행과 종료 기준

1. 기존 통과quarter4를 동일endpoint까지 한 번 계측 재생한다. 실제cell P, 분기조건/flag, reconstructed state/floor, HLLC/HLL의Pstar/vstar, 실제 선택값, geometric/work속도, force/work/Udot, ID·stage·입력hash를 잇는다. 기존수치와 일치하지 않으면 물리 해석을 중단한다.
2. 저장된 동일 입력에서 양 분기만 읽기 전용 평가한다. 실제 힘에는 적용하지 않는다. 압력변경과 work속도변경의 열률 jump 기여를 나눈다. threshold이동·hysteresis·분기강제 고정은 하지 않는다.
3. U→0+에서는 EOS·sound speed·재구성/분기를 함께 갱신한다. 고정gradient 평가는 부분 연산자로 표시한다. 실제 압력floor 포함 구현과 floor없는 한쪽 극한을 구분한다. 소수의 양의 값 외넓은scan은 하지 않는다.

- branch전환과jump가 일치하고 같은입력 양분기가재현하면 불연속이chatter/시간강직성에 기여했다고 판정한다. 분기제거를 해법으로 단정하지 않는다.
- 분기내 U→0+에서 유한 음의Udot가 남으면 공간양성결함이며dt축소로근본해결할수없다.
- 경계가안쪽을향하고RKstage만음수면admissible시간적분을 검토한다.
- floor/입력의존성 때문에 극한을 확정할수없으면unresolved로 종료하며 장시간/자동scan으로넘어가지 않는다.

## Adjoint 기준의 역할

실제체적Jacobian과pressureforce의power항등식은 연산자 장부의 기준이다. 단순adjoint/VPH복귀를 장시간 해결책으로 추천하지 않는다. 과거 보존적 가역 폐쇄의sliding-cell모드, weight-DAE유효성, staggered-remap한계를 보존한다. 후보가 현재 열고갈을 제거해도 별도sliding/shear선형응답과 충격·장시간관문이 남는다. 최종약1000공전요구를 축소하지 않는다.
