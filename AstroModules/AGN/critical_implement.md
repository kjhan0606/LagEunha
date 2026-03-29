# Latest AGN/SMBH Research Trends and Code Update Recommendations (2024-2025)

## Executive Summary

Based on a comprehensive review of the latest research (2024-2025), several transformative discoveries-primarily driven by JWST-require significant updates to AGN/SMBH simulation codes. The most critical findings involve:

1. **Little Red Dots (LRDs)** - A new population of high-z AGN
2. **Super-Eddington physics** - Essential for early BH growth
3. **Gas-enshrouded AGN model** - Non-stellar Balmer break origin
4. **Advanced spin evolution** - BZ spindown and state-dependent models
5. **GW recoil dynamics** - Wandering/escaping BHs

---

## 1. Little Red Dots (LRDs) - The Most Surprising JWST Discovery

### What Are They?
- Compact (R < 30-100 pc), red optical sources discovered by JWST
- Show broad-line AGN signatures but lack X-ray emission
- Characteristic "V-shaped" SED with Balmer break
- Found at z ~ 4-9, with ~341 confirmed as of 2025
- M_BH ~ 10^6-10^7 M_sun with extreme M_BH/M_star ratios (up to 0.1-1!)

### Key Challenge
Traditional AGN models cannot explain:
- V-shaped SED with Balmer break AND broad emission lines
- X-ray weakness (L_bol/L_X > 100-1000)
- No dust reemission despite red colors
- Extreme compactness and high number density

### Leading Theories (Inayoshi & Ho 2025; Inayoshi & Maiolino 2025)

**Gas-Enshrouded AGN Model:**
- Dense gas (n ~ 10^9-10^11 cm^-3) surrounds the BH with high covering fraction
- Gas produces Balmer absorption AND broad emission simultaneously
- Red optical spectrum WITHOUT dust reddening (thermal gas emission at T_eff ~ 5000 K)
- Naturally explains X-ray weakness and flat IR continuum

### Code Implementation

```c
/* Dense gas enshrouded AGN parameters */
#define N_GAS_DENSE_MIN     1e9     /* Minimum dense gas density [cm^-3] */
#define COVERING_FRACTION   0.8     /* Gas covering fraction */
#define T_EFF_GAS           5000.0  /* Effective temperature [K] */

typedef struct {
    double covering_fraction;    /* Solid angle coverage */
    double gas_density;          /* Circumnuclear gas density */
    double column_density;       /* N_H column density */
    bool is_lrd_candidate;       /* LRD classification flag */
} GasEnvelopeProperties;

/* Check if system qualifies as LRD-like */
bool check_lrd_conditions(SinkParticle *sink, double r_inner) {
    double n_gas = sink->weighted_density * sim.scale_nH;
    double M_BH_Msun = sink->mass * sim.scale_m / MSUN;
    
    /* LRDs: compact, high-density, high M_BH/M_star */
    if (r_inner < 100.0 * PC_TO_CM / sim.scale_l &&    /* < 100 pc */
        n_gas > N_GAS_DENSE_MIN &&                      /* Dense gas */
        M_BH_Msun > 1e6) {                              /* Massive BH */
        return true;
    }
    return false;
}

/* Modified X-ray luminosity for gas-enshrouded AGN */
double compute_xray_luminosity(SinkParticle *sink, double L_bol) {
    if (sink->gas_envelope.covering_fraction > 0.5) {
        /* Heavy obscuration reduces X-ray by factor 10-100 */
        double tau = sink->gas_envelope.column_density / 1e24;  /* Compton-thick if >1 */
        return L_bol / (50.0 * (1.0 + tau * tau));
    }
    return L_bol / 20.0;  /* Standard bolometric correction */
}
```

### References
- **Inayoshi & Ho (2025)** - "A Critical Evaluation of the Physical Nature of the Little Red Dots"
- **Inayoshi & Maiolino (2025)** - "Extremely Dense Gas around Little Red Dots"
- **Greene et al. (2024)** - LRD population studies
- **de Graaff et al. (2025)** - "A remarkable ruby" (The Cliff analysis)

---

## 2. Enhanced Super-Eddington Accretion Model

### JWST Motivation
JWST has found SMBHs at z > 9 (e.g., CAPERS-LRD-z9 at z=9.288, UHZ-1 at z=10.1) that are far more massive than expected from Eddington-limited growth.

### Latest Models (Hu#ko et al. 2025 MNRAS)

The SWIFT/EAGLE simulations show:
- Super-Eddington bursts around z ~ 7.5 increase BH mass by **orders of magnitude**
- This triggers strong feedback that quenches subsequent growth
- Results in galaxies that are half as massive but 10x more extended

### Slim Disc + Jet Model

```c
/* Super-Eddington accretion with slim disc (Hu#ko et al. 2025) */
typedef struct {
    double f_edd;           /* Eddington ratio */
    double epsilon_rad;     /* Radiative efficiency */
    double epsilon_jet;     /* Jet efficiency */
    double epsilon_wind;    /* Wind efficiency */
    int accretion_state;    /* THIN, SLIM, THICK */
} AccretionDiscState;

#define STATE_THICK  0   /* f_Edd < 0.01: ADAF/RIAF */
#define STATE_THIN   1   /* 0.01 < f_Edd < 1: Standard thin disc */
#define STATE_SLIM   2   /* f_Edd > 1: Super-Eddington slim disc */

AccretionDiscState compute_disc_state(double f_edd, double spin_a) {
    AccretionDiscState state;
    state.f_edd = f_edd;
    
    if (f_edd < 0.01) {
        /* ADAF/RIAF regime - geometrically thick, radiatively inefficient */
        state.accretion_state = STATE_THICK;
        state.epsilon_rad = 0.1 * f_edd / 0.01;  /* Very low radiative efficiency */
        state.epsilon_jet = compute_mad_efficiency(spin_a);  /* MAD jets efficient */
        state.epsilon_wind = 0.01;
        
    } else if (f_edd < 1.0) {
        /* Standard thin disc (Novikov-Thorne) */
        state.accretion_state = STATE_THIN;
        state.epsilon_rad = compute_radiative_efficiency(spin_a);
        state.epsilon_jet = 0;  /* No jets in quasar mode (standard model) */
        state.epsilon_wind = 0.05;  /* Some wind driving */
        
    } else {
        /* Slim disc regime - photon trapping reduces efficiency */
        state.accretion_state = STATE_SLIM;
        double epsilon_thin = compute_radiative_efficiency(spin_a);
        
        /* Slim disc efficiency (Wang & Zhou 1999; Madau et al. 2014) */
        state.epsilon_rad = epsilon_thin / (1.0 + log(f_edd));
        
        /* CRITICAL: Jets CAN launch even at super-Eddington (Ricarte et al. 2023) */
        /* Slim disc becomes geometrically thick # MAD jets possible */
        if (f_edd > 3.0) {
            state.epsilon_jet = 0.5 * compute_mad_efficiency(spin_a);
        } else {
            state.epsilon_jet = 0.1 * compute_mad_efficiency(spin_a) * (f_edd - 1.0) / 2.0;
        }
        
        /* Strong winds in super-Eddington */
        state.epsilon_wind = 0.1 * (1.0 + log(f_edd));
    }
    
    return state;
}

/* Apply super-Eddington feedback (Massonneau et al. 2023) */
void apply_super_eddington_feedback(SinkParticle *sink, AccretionDiscState *state) {
    double dm = sink->dM_accreted * sim.scale_m;
    double c2 = CLIGHT * CLIGHT;
    
    /* Total feedback energy */
    double E_rad = state->epsilon_rad * dm * c2;
    double E_jet = state->epsilon_jet * dm * c2;
    double E_wind = state->epsilon_wind * dm * c2;
    
    /* In super-Eddington: jets + winds dominate over radiation */
    if (state->accretion_state == STATE_SLIM) {
        /* Bipolar jet injection */
        if (E_jet > 0) {
            inject_bipolar_jet(sink, E_jet);
        }
        /* Spherical wind injection */
        inject_isotropic_wind(sink, E_wind);
        /* Reduced thermal (most radiation trapped) */
        inject_thermal_feedback(sink, 0.1 * E_rad);
    }
}
```

### Key Physics Updates
1. **Photon trapping**: #_r decreases logarithmically above Eddington
2. **Jets at super-Eddington**: Slim disc can still launch MAD jets (Ricarte et al. 2023)
3. **Strong winds**: Up to 10-30% of accretion power in winds
4. **Self-regulation**: Super-Eddington bursts followed by feedback quenching

### References
- **Hu#ko et al. (2025)** - MNRAS, "Super-Eddington accretion and feedback"
- **Ricarte et al. (2023)** - GRRMHD slim disc simulations
- **Pacucci & Narayan (2024)** - Super-Eddington X-ray weakness
- **Lupi et al. (2024)** - Sustained super-Eddington accretion

---

## 3. BZ Jet Spindown - Critical for Spin Distribution

### NewHorizon Results (Beckmann et al. 2025 MNRAS)

Key finding: **Including BZ spindown dramatically increases spin scatter**, especially for M_BH = 10^5-10^7 M_sun.

Without BZ spindown:
- BHs in this mass range predicted to be ~universally maximally spinning
- Average radiative efficiency <#_r> ~ 0.30

With BZ spindown:
- Significant scatter in spins (0.4 < |a| < 0.998)
- Average radiative efficiency <#_r> ~ 0.19
- 3-8x more energy available for AGN feedback!

### Implementation

```c
/* BZ jet spindown (Tchekhovskoy et al. 2010; Beckmann et al. 2025) */
void apply_bz_spindown(SinkParticle *sink, double dt) {
    /* Only in jet/radio mode */
    if (sink->dM_bondi / sink->dM_eddington > X_FLOOR) return;
    
    double a = fabs(sink->spin_magnitude);
    double M_BH = sink->mass * sim.scale_m;
    
    /* BZ power extraction rate */
    /* P_BZ = (kappa/4pi) * Phi_BH^2 * Omega_H^2 * f(a) */
    /* Omega_H = a*c / (2*r_H) where r_H = (1 + sqrt(1-a^2)) * GM/c^2 */
    
    double r_H = (1.0 + sqrt(1.0 - a*a)) * GRAVITY_CGS * M_BH / (CLIGHT * CLIGHT);
    double Omega_H = a * CLIGHT / (2.0 * r_H);
    
    /* MAD magnetic flux (dimensionless) */
    double phi_BH = 50.0;  /* Saturated MAD value */
    
    /* BZ efficiency (Tchekhovskoy et al. 2011) */
    double kappa = 0.05;  /* Coupling constant */
    double eta_BZ = kappa * phi_BH * phi_BH * (a / (1.0 + sqrt(1.0 - a*a)));
    eta_BZ = eta_BZ * eta_BZ;  /* Quadratic in Omega_H */
    
    /* Angular momentum extraction rate */
    /* dJ/dt = P_BZ / Omega_H */
    double P_BZ = eta_BZ * sink->dM_bondi * sim.scale_m / sim.scale_t * CLIGHT * CLIGHT;
    double dJ_dt = P_BZ / Omega_H;
    
    /* Spin angular momentum: J = a * G * M^2 / c */
    double J_BH = a * GRAVITY_CGS * M_BH * M_BH / CLIGHT;
    
    /* Update spin */
    double dJ = dJ_dt * dt * sim.scale_t;
    double new_J = J_BH - dJ;
    
    /* New spin magnitude */
    double new_a = new_J * CLIGHT / (GRAVITY_CGS * M_BH * M_BH);
    new_a = fmax(0.0, fmin(new_a, MAX_SPIN));
    
    /* Apply spindown (maintain direction) */
    if (sink->spin_magnitude > 0) {
        sink->spin_magnitude = new_a;
    } else {
        sink->spin_magnitude = -new_a;
    }
}
```

### References
- **Beckmann et al. (2025)** - MNRAS 536, 1838, "BH spin evolution from NewHorizon"
- **Tchekhovskoy et al. (2010, 2011)** - BZ power calculations
- **Narayan et al. (2022)** - MAD jet efficiency

---

## 4. GW Recoil Kicks - Wandering and Escaping BHs

### Latest Results (Dong-Páez et al. 2025 A&A)

Using NewHorizon simulations with full GW recoil implementation:
- Kick magnitudes depend on mass ratio, spins, and orbital angular momentum
- "Dry" gas-poor mergers # random spin configurations # larger kicks
- "Wet" gas-rich mergers # aligned configurations # smaller kicks
- But: high vs low accretion rate is NOT sufficient to predict kick magnitude

### Key Finding
The distinction between wandering and escaping BHs matters for:
- Population statistics
- Second-generation mergers
- Nuclear BH occupation fraction

### Implementation

```c
/* GW recoil kick calculation (Varma et al. 2019 surfinBH; Dong-Páez et al. 2025) */
typedef struct {
    Vec3 v_kick;        /* Kick velocity vector */
    double v_mag;       /* Kick magnitude */
    bool will_escape;   /* True if v_kick > v_escape */
} GWRecoilResult;

GWRecoilResult compute_gw_recoil(SinkParticle *s1, SinkParticle *s2,
                                  Vec3 L_orbital, double v_escape) {
    GWRecoilResult result;
    
    /* Ensure M1 >= M2 */
    double M1 = fmax(s1->mass, s2->mass);
    double M2 = fmin(s1->mass, s2->mass);
    double q = M2 / M1;
    double eta = q / ((1.0 + q) * (1.0 + q));
    
    double a1 = fabs(s1->mass >= s2->mass ? s1->spin_magnitude : s2->spin_magnitude);
    double a2 = fabs(s1->mass >= s2->mass ? s2->spin_magnitude : s1->spin_magnitude);
    Vec3 chi1 = s1->mass >= s2->mass ? s1->spin_direction : s2->spin_direction;
    Vec3 chi2 = s1->mass >= s2->mass ? s2->spin_direction : s1->spin_direction;
    
    Vec3 L_hat = vec3_normalize(L_orbital);
    
    /* Spin components parallel and perpendicular to L */
    double chi1_para = a1 * vec3_dot(chi1, L_hat);
    double chi2_para = a2 * vec3_dot(chi2, L_hat);
    Vec3 chi1_perp = vec3_sub(vec3_scale(chi1, a1), vec3_scale(L_hat, chi1_para));
    Vec3 chi2_perp = vec3_sub(vec3_scale(chi2, a2), vec3_scale(L_hat, chi2_para));
    
    /* Fitting coefficients (Lousto & Zlochower 2013) */
    double A_m = 1.2e4;      /* km/s - mass asymmetry */
    double B_m = -0.93;
    double H = 6.9e3;        /* km/s - spin asymmetry (perpendicular) */
    double K_s = 6.0e4;      /* km/s - superkick (parallel) */
    double V_perp = 0.0;     /* Additional perpendicular contribution */
    
    /* Mass asymmetry kick (in orbital plane) */
    double v_m = A_m * eta * eta * (1.0 - q) / (1.0 + q) * (1.0 + B_m * eta);
    
    /* Spin asymmetry kick (perpendicular to orbital plane) */
    double Delta_para = (chi2_para - q * chi1_para) / (1.0 + q);
    double v_perp = H * eta * eta * Delta_para;
    
    /* Superkick (in orbital plane, perpendicular to separation) */
    Vec3 Delta_perp = vec3_scale(vec3_sub(chi2_perp, vec3_scale(chi1_perp, q)), 
                                  1.0 / (1.0 + q));
    double Delta_perp_mag = vec3_mag(Delta_perp);
    
    /* Phase angle (random if not tracked) */
    double phase = 2.0 * PI * ((double)rand() / RAND_MAX);
    double v_superkick = K_s * eta * eta * eta * Delta_perp_mag * cos(phase);
    
    /* Total kick magnitude */
    result.v_mag = sqrt(v_m * v_m + v_perp * v_perp + v_superkick * v_superkick);
    
    /* Convert to code units */
    result.v_mag *= 1e5 / sim.scale_v;  /* km/s to code units */
    
    /* Kick direction (simplified: along orbital plane + perpendicular) */
    Vec3 e_r = vec3_normalize(vec3_sub(s1->position, s2->position));
    Vec3 e_phi = vec3_cross(L_hat, e_r);
    
    result.v_kick = vec3_add(
        vec3_add(vec3_scale(e_r, v_m + v_superkick * cos(phase)),
                 vec3_scale(e_phi, v_superkick * sin(phase))),
        vec3_scale(L_hat, v_perp)
    );
    result.v_kick = vec3_scale(result.v_kick, 1e5 / sim.scale_v);
    
    /* Check if escaping */
    result.will_escape = (result.v_mag > v_escape);
    
    return result;
}

/* Apply kick after merger */
void apply_gw_kick(SinkParticle *remnant, GWRecoilResult *kick) {
    remnant->velocity = vec3_add(remnant->velocity, kick->v_kick);
    
    if (kick->will_escape) {
        printf("WARNING: BH %d escaping with v_kick = %.1f km/s\n",
               remnant->id, kick->v_mag * sim.scale_v / 1e5);
        remnant->is_wandering = true;
    }
}
```

### References
- **Dong-Páez et al. (2025)** - A&A 695, A231, "Wandering and escaping BHs"
- **Varma et al. (2019)** - surfinBH surrogate model
- **Lousto & Zlochower (2013)** - Nonlinear kick formula

---

## 5. Multi-Channel Heavy Seeding

### JWST-Driven Need

The "overmassive" BHs found by JWST (M_BH/M_star ~ 0.01-0.1 vs ~0.001 locally) require:
- Either super-Eddington growth from light seeds
- Or heavy seeds (10^4-10^6 M_sun) that grow normally
- Or primordial BHs as seeds

### Implementation

```c
/* Multi-channel BH seeding (Zhang et al. 2025; Inayoshi et al. 2020) */
typedef enum {
    SEED_POPIII = 0,      /* Pop III stellar remnant: ~100 Msun */
    SEED_DCBH = 1,        /* Direct collapse BH: ~10^5 Msun */
    SEED_NSC = 2,         /* Nuclear star cluster runaway: ~10^3 Msun */
    SEED_PBH = 3          /* Primordial BH: variable mass */
} SeedingChannel;

typedef struct {
    SeedingChannel channel;
    double M_seed;
    double prob;           /* Formation probability */
} SeedingModel;

SeedingModel determine_seeding(VoronoiCell *cell, double J_LW, double Z_gas) {
    SeedingModel model;
    
    double n_gas = cell->density * sim.scale_nH;
    double T_vir = cell->temperature;
    
    /* Critical values */
    double J_crit = 1000.0;      /* Lyman-Werner flux threshold [J_21] */
    double Z_crit = 1e-4;        /* Metallicity threshold [Zsun] */
    double T_atomic = 8000.0;    /* Atomic cooling threshold [K] */
    
    /* Direct Collapse BH channel (Latif et al. 2022) */
    if (J_LW > J_crit && Z_gas < Z_crit && T_vir > T_atomic) {
        model.channel = SEED_DCBH;
        model.M_seed = 1e5 * MSUN / sim.scale_m;
        model.prob = 0.01;  /* Rare but possible */
        return model;
    }
    
    /* Nuclear Star Cluster channel */
    if (cell->stellar_density > 1e6 &&  /* Very high stellar density */
        cell->velocity_dispersion > 50.0) {  /* High velocity dispersion */
        model.channel = SEED_NSC;
        model.M_seed = 1e3 * MSUN / sim.scale_m;
        model.prob = 0.1;
        return model;
    }
    
    /* Pop III remnant (default) */
    model.channel = SEED_POPIII;
    model.M_seed = 100.0 * MSUN / sim.scale_m;
    model.prob = 0.5;  /* Common in early universe */
    
    return model;
}
```

---

## 6. Summary: Priority Implementation Roadmap

### Highest Priority (Critical for JWST Interpretation)

| Feature | Impact | Difficulty | References |
|---------|--------|------------|------------|
| Super-Eddington + slim disc | Explains high-z massive BHs | Medium | Hu#ko+25, Lupi+24 |
| State-dependent accretion | Proper disc physics | Medium | Ricarte+23 |
| BZ jet spindown | Realistic spin distribution | Low | Beckmann+25 |

### High Priority (Major Physics Updates)

| Feature | Impact | Difficulty | References |
|---------|--------|------------|------------|
| GW recoil kicks | Wandering BH population | Medium | Dong-Páez+25 |
| Gas-enshrouded AGN | LRD interpretation | High | Inayoshi+25 |
| Jets at super-Eddington | Feedback at high-z | Medium | Ricarte+23 |

### Medium Priority (Refinements)

| Feature | Impact | Difficulty | References |
|---------|--------|------------|------------|
| Heavy seeding channels | BH mass function | Low | Zhang+25 |
| X-ray weakness modeling | Observable predictions | Low | Pacucci+24 |
| Spin-galaxy alignment | Jet orientation | Low | Peirani+24 |

---

## Key Takeaways

1. **LRDs are NOT dust-reddened AGN**: They require a new "gas-enshrouded" model with dense circumnuclear gas at high covering fraction.

2. **Super-Eddington is essential**: Without it, simulations cannot reproduce JWST high-z BH masses. Include slim disc + jet launching.

3. **BZ spindown matters**: Including jet spindown increases spin scatter significantly and affects feedback energy budget by factors of 3-8.

4. **GW kicks create wandering BHs**: Must be included for accurate BH population statistics and merger rates.

5. **State transitions are key**: Need smooth transitions between ADAF # thin # slim disc with appropriate efficiency and feedback modes for each.

---

## Complete Reference List

### Little Red Dots
- Inayoshi & Ho (2025) - arXiv:2512.03130
- Inayoshi & Maiolino (2025) - ApJL 980, L27
- Greene et al. (2024) - LRD population
- de Graaff et al. (2025) - A&A 701, A168

### Super-Eddington Accretion
- Hu#ko et al. (2025) - MNRAS 537, 2559
- Ricarte et al. (2023) - GRRMHD simulations
- Lupi et al. (2024) - Sustained super-Eddington
- Pacucci & Narayan (2024) - X-ray weakness

### BH Spin Evolution
- Beckmann et al. (2025) - MNRAS 536, 1838
- Sala et al. (2024) - A&A, OPENGADGET3
- Peirani et al. (2024) - A&A, spin-galaxy alignment

### GW Recoil
- Dong-Páez et al. (2025) - A&A 695, A231
- Varma et al. (2019) - surfinBH

### High-z JWST AGN
- Taylor et al. (2025) - ApJL 989, L7 (z=9.288 AGN)
- Goulding et al. (2024) - UHZ-1 at z=10.1
- Tripodi et al. (2025) - Nature Comm., LRD at z=8.6
