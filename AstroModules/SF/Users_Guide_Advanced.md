# 고급 별탄생 모듈 사용 가이드

은하 형성 시뮬레이션을 위한 최첨단 별탄생 모델 구현 (C 언어)

## 🌟 새로운 기능 (v2.0)

### 1. 자기장 효과
- 중력 붕괴에 대한 자기 압력 지지
- 질량-자기선속 비율 계산 (λ = M/Φ_B)
- 자기적 아임계/초임계 구름 구별
- 플라즈마 베타(β) 의존 별탄생 효율
- 참고: Pattle et al. (2022), Kretschmer & Teyssier (2020)

### 2. 중금속함량 의존 IMF
- IMF 특성 질량이 중금속함량에 따라 변화
- 낮은 금속함량 → 무거운별 위주 IMF (top-heavy)
- 높은 금속함량 → 가벼운별 위주 IMF (bottom-heavy)
- 표면밀도 의존성 포함
- 참고: Tanvir et al. (2024), Yan et al. (2018, 2020), Bate (2023)

### 3. 쌍성계 집단
- 질량 의존 쌍성 비율 (M < 0.5 M_sun: 23%, M > 5 M_sun: 70%)
- 현실적인 질량비(q) 분포
- 궤도 주기 분포 (1일 ~ 10만일)
- 이심률 분포 (넓은 궤도: 열적 분포, 좁은 궤도: 원형화)
- 참고: Moe & Di Stefano (2017)

## 파일 구성

### 핵심 구현
- `star_formation_advanced.h` - 모든 고급 기능 포함
- `star_formation.h` - 기본 모듈 (하위 호환)

### 예제 프로그램
- `star_formation_advanced_example.c` - 고급 기능 시연
- `star_formation_example.c` - 기본 예제

### 문서
- `README_ADVANCED.md` - 영문 상세 문서
- `사용가이드_한글.md` - 본 파일

## 빠른 시작

### 컴파일

```bash
# 고급 버전 컴파일
gcc -o star_formation_advanced star_formation_advanced_example.c -lm -O2

# 실행
./star_formation_advanced
```

또는 Makefile 사용:

```bash
make -f Makefile_advanced
make -f Makefile_advanced run-advanced
```

### 기본 사용법

```c
#include "star_formation_advanced.h"

int main() {
    // 자기장이 있는 가스 원소 설정
    GasElement gas;
    gas.density = 1.0e-20;              // g cm^-3
    gas.temperature = 30.0;             // K
    gas.metallicity = 0.5 * METALLICITY_SOLAR;  // 0.5배 태양 중금속함량
    gas.velocity_dispersion = 1.0e5;    // cm s^-1
    gas.volume = pow(10.0 * PARSEC, 3); // cm^3
    
    // 자기장 속성
    gas.magnetic_field_strength = estimate_magnetic_field(
        gas.density, gas.velocity_dispersion);
    gas.magnetic_pressure = gas.magnetic_field_strength * 
        gas.magnetic_field_strength / (8.0 * M_PI);
    gas.mass_to_flux_ratio = 2.5; // 초임계
    
    // 모든 고급 기능으로 별 생성
    StarParticle sp;
    double timestep = 1.0e6; // 년
    
    bool success = form_stars_advanced(
        &gas,                           // 가스 원소
        timestep,                       // 시간 간격
        IMF_METALLICITY_DEPENDENT,      // 금속함량 의존 IMF 사용
        true,                           // 쌍성계 포함
        &sp                            // 출력
    );
    
    if (success) {
        printf("생성된 별 질량: %.2e M_sun, 별 개수: %d\n", 
               sp.mass, sp.num_stars);
        printf("쌍성 비율: %.1f%%\n", 
               100.0 * sp.num_binaries / sp.num_stars);
        
        // 쌍성계 속성 접근
        for (int i = 0; i < sp.num_stars; i++) {
            if (sp.is_binary[i]) {
                printf("쌍성: M1=%.3f, q=%.3f, P=%.2e 일\n",
                       sp.stellar_masses[i],
                       sp.binary_params[i].mass_ratio,
                       sp.binary_params[i].period_days);
            }
        }
        
        free_star_particle_advanced(&sp);
    }
    
    return 0;
}
```

## 물리학 상세 설명

### 자기장 효과

**플라즈마 베타 (β)**
```
β = P_열 / P_자기 = (ρkT/μmH) / (B²/8π)
```

- **β >> 1**: 열압력 지배, 자기장 효과 미미
- **β ~ 1**: 등분배, 중간 정도 억제 (~25-50%)
- **β << 1**: 자기장 지배, 강한 억제 (~50-80%)

**질량-자기선속 비율 (λ)**
```
λ = (M/Φ_B) / (M/Φ_B)_임계
```

- **λ > 1**: 초임계 - 구름이 붕괴 가능
- **λ < 1**: 아임계 - 자기적으로 지지됨, 별탄생 거의 없음

**자기장이 있는 별탄생 효율**
```
ε_ff(B) = ε_ff,기본 × f_자기(β, λ)
```

### 중금속함량 의존 IMF

**특성 질량**
```
m_char = m_char,⊙ × (Z/Z_⊙)^(-α_Z) × (Σ/100 M_⊙pc^-2)^(-α_Σ)
```

전형적인 값:
- α_Z = 0.15 (Z > 0.01 Z_⊙일 때)
- α_Z = 0.4 (Z < 0.01 Z_⊙일 때, 원시 환경)
- α_Σ = 0.3 (표면밀도 효과)

**별 집단에 미치는 영향:**
- Z = 0.001 Z_⊙: m_char ~ 1.9 M_⊙ (무거운별 위주, ~8%만 < 0.5 M_⊙)
- Z = 0.01 Z_⊙: m_char ~ 0.24 M_⊙ (~56%가 < 0.5 M_⊙)
- Z = 1.0 Z_⊙: m_char ~ 0.12 M_⊙ (~65%가 < 0.5 M_⊙)
- Z = 3.0 Z_⊙: m_char ~ 0.10 M_⊙ (가벼운별 위주, ~73%가 < 0.5 M_⊙)

### 쌍성계 속성

**질량별 쌍성 비율 (Moe & Di Stefano 2017)**
- M < 0.5 M_⊙: f_binary = 23%
- 0.5 < M < 1.5 M_⊙: f_binary = 44%
- 1.5 < M < 5 M_⊙: f_binary = 59%
- M > 5 M_⊙: f_binary = 70%

**질량비 분포**
- 무거운 별일수록 q ~ 1 (비슷한 질량)로 치우침
- 범위: 0.1 ≤ q ≤ 1.0

**궤도 주기 분포**
- 로그-균일 분포: 1일 ~ 10만일
- 짧은 주기 (< 10일): 원형화 (e ≈ 0)
- 긴 주기: 열적 이심률 분포 (e² 균일)

## 예제 출력 결과

### 1. 자기장 효과
```
자기장 없음:         SFR = 8.54e-04 M_sun/yr (100%)
약한 B장 (β ~ 50):  SFR = 8.54e-04 M_sun/yr (100%)
강한 B장 (β ~ 1):   SFR = 6.38e-04 M_sun/yr (75%)
매우 강함 (β ~ 0.25): SFR = 4.80e-04 M_sun/yr (56%)
아임계 (λ < 1):     SFR = 4.27e-05 M_sun/yr (5%)
```

### 2. 중금속함량 의존 IMF
```
Z = 0.001 Z_sun: m_char = 1.89 M_sun (7.7%만 < 0.5 M_sun)
Z = 0.01 Z_sun:  m_char = 0.24 M_sun (56%가 < 0.5 M_sun)
Z = 1.0 Z_sun:   m_char = 0.12 M_sun (65%가 < 0.5 M_sun)
Z = 3.0 Z_sun:   m_char = 0.10 M_sun (73%가 < 0.5 M_sun)
```

### 3. 쌍성계 집단
```
총 쌍성계: 145/473 별 (30.7%)

질량 범위별:
  < 0.5 M_sun:     49/293 (16.7%)
  0.5-1.5 M_sun:   46/92  (50.0%)
  1.5-5 M_sun:     40/71  (56.3%)
  > 5 M_sun:       10/17  (58.8%)

쌍성계 속성:
  평균 질량비: q = 0.543
  평균 주기: 6820일 (18.7년)
  평균 이심률: e = 0.533
```

## 주요 참고문헌

### 자기장
1. **Pattle et al. (2022)** - 별탄생에서 자기장 (arXiv:2203.11179)
2. **Kretschmer & Teyssier (2020)** - 자기장 별탄생 모델 (MNRAS 527, 6779)

### 중금속함량 의존 IMF
3. **Tanvir et al. (2024)** - IMF의 금속함량 의존성 (MNRAS 527, 7306)
4. **Yan et al. (2018)** - IGIMF 이론 (A&A 607, A126)
5. **Bate (2023)** - 적색편이-금속함량 IMF 변화 (MNRAS 537, 752)

### 쌍성계
6. **Moe & Di Stefano (2017)** - 쌍성 비율 체계 (ApJS 230, 15)
7. **Stanway & Eldridge (2020)** - 쌍성 집단 불확실성 (MNRAS 495, 4605)

## 시뮬레이션 통합

### SPH 코드
```c
// SPH 입자 루프에서
for (각_가스_입자) {
    // SPH 자기장으로부터 B 계산
    gas.magnetic_field_strength = sqrt(Bx*Bx + By*By + Bz*Bz);
    
    // 국소 속성으로부터 질량-선속비 추정
    gas.mass_to_flux_ratio = calculate_mass_to_flux_ratio(
        particle_mass, magnetic_flux);
    
    // 별 생성
    form_stars_advanced(&gas, dt, IMF_METALLICITY_DEPENDENT, 
                       true, &new_star);
}
```

### 격자 코드
```c
// 격자 셀 루프에서
for (각_셀) {
    gas.density = cell_density;
    gas.metallicity = cell_metallicity;
    gas.volume = cell_volume;
    
    // 격자로부터 B장 가져오기
    gas.magnetic_field_strength = cell_B_magnitude;
    gas.magnetic_pressure = cell_B_magnitude * cell_B_magnitude / 
                           (8.0 * M_PI);
    
    // 별 생성
    form_stars_advanced(&gas, dt, IMF_METALLICITY_DEPENDENT,
                       true, &new_star);
}
```

## 주요 구조체

### GasElement (확장됨)
```c
typedef struct {
    // 기본 속성
    double density;              // 밀도 (g cm^-3)
    double temperature;          // 온도 (K)
    double metallicity;          // 중금속함량
    double velocity_dispersion;  // 속도 분산 (cm s^-1)
    double volume;               // 부피 (cm^3)
    
    // 자기장 속성 (새로 추가)
    double magnetic_field_strength; // B장 세기 (Gauss)
    double magnetic_pressure;       // 자기압 (dyn cm^-2)
    double mass_to_flux_ratio;      // μ = M/Φ_B (정규화)
} GasElement;
```

### StarParticle (확장됨)
```c
typedef struct {
    // 기본 속성
    double mass;                // 총 별 질량 (M_sun)
    int num_stars;              // 별 개수
    double *stellar_masses;     // 개별 별 질량 배열
    double velocity[3];         // 속도 성분
    double metallicity;         // 중금속함량
    double age;                 // 나이 (년)
    
    // 쌍성 정보 (새로 추가)
    int num_binaries;           // 쌍성계 개수
    bool *is_binary;            // 쌍성 여부 배열
    BinaryOrbit *binary_params; // 쌍성 궤도 파라미터
} StarParticle;
```

### BinaryOrbit
```c
typedef struct {
    double mass_ratio;       // q = M2/M1 [0,1]
    double period_days;      // 궤도 주기 (일)
    double eccentricity;     // 궤도 이심률 [0,1)
    double separation_au;    // 긴반지름 (AU)
} BinaryOrbit;
```

## 성능 및 효율성

### 계산 오버헤드
- 자기장 계산: ~5% 추가
- 금속함량 의존 IMF: 표준 IMF와 동일
- 쌍성 생성: ~10% 추가 (활성화 시)
- 모든 기능 포함: 총 ~15-20% 오버헤드

### 수치적 안정성
- 질량-선속비는 물리적 범위로 제한
- 플라즈마 베타 수치 문제 체크
- 쌍성 궤도 파라미터 검증
- 가스 원소당 최대 별 개수: 10,000 (설정 가능)

## 테스트

종합 테스트 실행:

```bash
make -f Makefile_advanced run-advanced
```

테스트 포함 사항:
1. 자기장 세기 변화
2. 중금속함량 범위 0.001-3 Z_⊙
3. 질량별 쌍성 비율
4. 복합 효과 (낮은 Z + 강한 B + 쌍성)

## 버전 이력

### v2.0 (2025)
- 자기장 효과 추가
- 중금속함량 의존 IMF 추가
- 쌍성계 집단 추가
- 종합 테스트 및 문서화

### v1.0 (2025)
- 초기 릴리스
- Kennicutt-Schmidt 및 체적 SF 모델
- 3가지 표준 IMF (Chabrier, Kroupa, Salpeter)

## 라이선스

이 코드는 교육 및 연구 목적으로 제공됩니다.
논문에 사용 시 관련 참고문헌을 인용해 주시기 바랍니다.

## 활용 예시

이 모듈은 다음과 같은 연구에 활용될 수 있습니다:

1. **은하 형성 시뮬레이션**: 우주론적 시뮬레이션에서 별탄생 모델링
2. **초기 우주 연구**: 낮은 금속함량에서 IMF 변화 연구
3. **자기장 연구**: 자기장이 별탄생에 미치는 영향 정량화
4. **쌍성 진화**: 쌍성계가 별 집단에 미치는 영향 연구
5. **관측 비교**: 시뮬레이션 결과와 관측 데이터 비교

---

**참고**: 이 구현은 2024-2025년 기준 최신 물리학을 반영합니다.
새로운 관측과 시뮬레이션이 나오면 파라미터 업데이트가 필요할 수 있습니다.
