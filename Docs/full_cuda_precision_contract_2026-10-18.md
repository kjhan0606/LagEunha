# 전체 CUDA 경로의 현재 경계와 정밀도 계약

## 현재 코드에서 확인한 사실

- 라이브 환경은 GFS_TCF/FACE_Z/QBP, shared geometry, 고차 재구성, HG split 및 새 residual memory를 활성화한다.
- `Exam/exam_gpu_extract.c`의 getAccVoro2DBlend_GPU는 TCF/DYN_W2를 명시적으로 거부한다. 이 보호 조건을 제거하는 것만으로 같은 방정식의 GPU 구현이 되지 않는다.
- `Exam/exam_gpu.h`의 ParticleSoA는 위치·속도·질량을 double로 보유하지만 den/pressure/csound/volume/ie/w2/w2old/K 및 die_out/dK_out/dt_out은 float이다.
- 기존 tcf_full_face_probe의 CPU/CUDA 일치는 double 입력을 사용한 독립 전사식 검사다. 위 SoA 변환, 실제 QBP 재구성, 공유면 캐시, HG/RK4 호출은 검사하지 않았다.
- 새 residual memory의 CUDA 행/독립 조립 검사와 생산 OpenMP 검사를 전체 GPU hydro 검증으로 간주하지 않는다. 기존 CPU 물리 참조 궤도를 유지한다.

## 보존한 실제 면 자료의 압력 저장 정밀도 검사

`check_tcf_pressure_storage_precision.py`는 기존 t6 68,988개 면/11,498개 owner 자료를 재사용했다. double 계산식과 기하·속도·Q를 고정한 채 압력 두 값만 float32를 거쳐 저장하도록 바꿨다. 독립 FP64 계산은 기존 CPU 출력과 최대차0이었다.

| 항목 | 상대 L2 차이 | 최대 절대차 |
|---|---:|---:|
| 압력 저장 | 2.09838e-8 | 6.61116e-12 |
| 면 힘 | 1.56986e-8 | 1.76122e-13 |
| owner별 합산 힘 | 7.17727e-8 | 1.85207e-13 |
| owner별 합산 에너지 일 | 6.45891e-8 | 2.66175e-13 |

Q를 고정했으므로 열과 부피율은 그대로였다. 이 검사는 밀도/음속/재구성/Q 변화를 제외한 **압력 저장 하나의 영향**이며 전체 혼합정밀도 오차의 상한도 하한도 아니다. 측정한 차이가 이전 궤도 붕괴의 원인이라는 결론은 내리지 않는다.

증거: `canonical_tcf_pressure_storage_precision_2026-10-18.json`. 기존 계산/입력/바이너리/출력은 바꾸지 않았다.

## 다음 전체 GPU 구현의 필수 조건

1. 첫 같은 방정식 참조 구현은 현재 CPU TCF 입력/출력을 FP64로 운반한다. 기존 float SoA를 그대로 통과시킨 결과를 독립 double 커널의 정합성 결과와 혼동하지 않는다. 정밀도 축소는 이후 별도 수렴 시험 대상이다.
2. shared face geometry와 QBP/CHI/limiter/FACE_Z의 실제 단계 입력을 전달한다. generic HLLC 경로로 대체하지 않는다. 이미 끝난 독립 face algebra 검사는 재사용한다.
3. owner별 force/energy/volume/CFL 합산과 host writeback을 기존 실제 CPU RHS와 비교한다. 단순 면별 일치만으로 이 단계를 통과시키지 않는다.
4. 새 z/Bv/-B^Tz와 기존 열 갱신을 분리한다. GPU 출력 die와 host auxiliary acceleration 사이에 이중 열변환을 넣지 않는다.
5. 전체 RK4 stage/rejection/restart를 비영 상태에서 검증한 뒤 기존 체크포인트에서 연속 진화한다. MPI 소유권/ghost 교환은 별도 필수이며 기존 물리 우선 순서를 유지한다.
6. 전체 목적의 완료 조건은 여전히 역사적 t6–10과 이후 장시간 물리 정확도 및 Laguerre 확장 가능성이다. 위 구현 검사만으로 완료되지 않는다.

## 실제 CPU 호출부와 CUDA가 공유하는 함수 연결

`Exam/tcf_face_terms.h`에 FP64 정상/접선 힘, 에너지 일, 부피율, 열 계산을 묶었다. 재구성된 Q와 공통 면 기하는 호출자가 공급한다. 정상/접선 힘을 분리해 반환하므로 기존 합산 순서와 NOTAN 분기를 보존할 수 있다. `Exam/exam.c`는 `GFS_TCF_SHARED_ALGEBRA=1`일 때 이 함수를 사용하며 기본은 기존 계산식이다. 환경 변수는 스레드 루프 바깥에서 읽는다. 시험용 `tcf_full_face_probe.h`도 같은 함수를 호출하도록 바꿔 별도 전사식의 중복을 제거했다.

별도 빌드 `Codex_tcf_shared_build_20261018/eunha2.exe`의 SHA256은 `98c6f42060535d6c4149d625d24831df6335add999bd78e2e676f9fcd496c42c`이다. `build_memory_parallel.py --root ...`로 새 경로에만 빌드했다. 후속 빌드의 Exam staging에는 새 tcf_face_terms.h도 포함해야 한다.

`run_tcf_shared_control.py`는 보존한 t=.4 snapshot161/비영 z에서 새 옵션을 켜고 한 RK4 스텝을 수행했다. 최종t=.4025000000000003에서 기존 완료 old 출력과 reader11개 필드 및 z가 모두 bitwise 동일했다. 이전 대조 시뮬레이션은 다시 실행하지 않았다. 새 디렉토리에 원본/최종 합계 약50MB만 사용했다. 근거: `canonical_tcf_shared_control_2026-10-18.json`.

Ada에서 공용 함수를 사용해68,988개 기존 면 입력을 실행했다. OpenMP1/4는 bitwise 동일하고 새 CPU 출력은 이전 CPU probe 출력과 파일 바이트도 같았다. **GPU의 exact-byte 비교는 실패**했다. 따라서 GPU bitwise 일치라고 보고하지 않는다. 독립적으로 성분별 상대 L2 차이를 계산했다:

| 성분 | CPU/GPU 상대 L2 차이 | 최대 절대차 |
|---|---:|---:|
| Fx/Fy | 0 | 0 |
| 에너지 일 | 1.05925e-16 | 1.69407e-21 |
| 부피율 | 6.49860e-17 | 1.38778e-17 |
| 열 | 1.01675e-15 | 1.05879e-22 |

수치 오차 관문1e-12는 통과했으며 비트 동일 관문과 구별해 기록했다. 해당 차이의 정확한 코드 생성 원인을 분리 실험하지 않았으므로 특정 컴파일러 옵션 탓으로 단정하지 않는다. 근거: `canonical_tcf_shared_cuda_2026-10-18.json`, remote `Codex_tcf_shared_cuda_20261018/analysis_detailed.json`.

이 작업은 시험용 전사식을 실제 호출 가능한 공용 계산식으로 연결한 것이다. GPU TCF 전체 경로를 켜지 않았고 기존 unsupported guard는 유지한다. 실제 GPU 기하/재구성/합산/writeback/RK4 결합 및 장시간 물리 검증은 여전히 남아 있다. 실행 중인 장시간 바이너리는 교체하지 않았다.
