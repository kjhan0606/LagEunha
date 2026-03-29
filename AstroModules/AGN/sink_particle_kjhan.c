/*******************************************************************************
 * sink_particle_kjhan.c
 * 
 * AGN/SMBH Physics Module for Galaxy Formation Simulations
 * 
 * Converted from RAMSES Fortran code (sink_particle_kjhan.f90)
 * Uses Voronoi tessellation with Finite Volume Method (FVM)
 * 
 * Key Physics:
 *   - Black hole seeding and sink particle creation
 *   - Bondi-Hoyle accretion with boost factor
 *   - Spin evolution (warped disc model, self-gravity regime)
 *   - BH-BH mergers with Rezzolla spin formula
 *   - Dual-mode AGN feedback (thermal/quasar and kinetic/jet)
 *   - MAD jet efficiency model
 * 
 * References:
 *   - Dubois et al. (2012) - Dual AGN feedback
 *   - Booth & Schaye (2009) - Bondi accretion boost
 *   - Fanidakis et al. (2011) - Spin evolution
 *   - McKinney et al. (2012) - MAD jet efficiency
 *   - Rezzolla et al. (2008) - Merger spin formula
 * 
 * Author: Converted to C with Voronoi/FVM structure
 * Original: Romain Teyssier (2007), modified by kjhan
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

/*******************************************************************************
 * PHYSICAL CONSTANTS (CGS Units)
 ******************************************************************************/
#define CLIGHT          3.0e10          /* Speed of light [cm/s] */
#define GRAVITY_CGS     6.67e-8         /* Gravitational constant [cgs] */
#define MSUN            2.0e33          /* Solar mass [g] */
#define PROTON_MASS     1.66e-24        /* Proton mass [g] */
#define SIGMA_THOMSON   6.652e-25       /* Thomson cross section [cm^2] */
#define YEAR_TO_SEC     3.1557e7        /* Year in seconds */
#define PI              3.14159265358979
#define FOURPI          (4.0 * PI)

/*******************************************************************************
 * SIMULATION PARAMETERS
 ******************************************************************************/
#define NDIM            3               /* Number of dimensions */
#define MAX_SINKS       1000            /* Maximum number of sink particles */
#define MAX_CELLS       1000000         /* Maximum number of Voronoi cells */
#define MAX_NEIGHBORS   50              /* Maximum neighbors per cell */

/* Black hole seeding parameters */
#define N_SINK_THRESHOLD    1.0e4       /* Density threshold for BH seeding [H/cc] */
#define STELLAR_DENSITY_MIN 0.1         /* Minimum stellar density [H/cc] */
#define VELOCITY_DISP_MIN   20.0        /* Minimum velocity dispersion [km/s] */
#define SEED_MASS           1.0e5       /* Seed BH mass [Msun] */

/* Accretion parameters */
#define BOOST_ACC           2.0         /* Boost factor exponent */
#define N_STAR_THRESHOLD    0.1         /* Reference density for boost */
#define SIGMAV_MAX          10.0        /* Max velocity dispersion [km/s] */
#define F_BONDI             1.0         /* Bondi correction factor */

/* AGN feedback parameters */
#define X_FLOOR             0.01        /* Eddington ratio threshold for mode switch */
#define EAGN_THERMAL        0.15        /* Thermal feedback efficiency */
#define EAGN_KINETIC        1.0         /* Kinetic feedback efficiency */
#define F_EK_AGN            0.5         /* Kinetic energy fraction in jet mode */
#define T2_MAX_AGN          1.0e10      /* Maximum AGN heating temperature [K] */
#define R_AGN               50.0        /* AGN feedback radius [kpc] */

/* Spin parameters */
#define MAX_SPIN            0.998       /* Maximum allowed spin */
#define ALPHA_VISCOSITY     0.1         /* Shakura-Sunyaev alpha parameter */
#define NU_RATIO            0.1         /* Warp viscosity ratio */

/* Merger parameters */
#define MERGE_RADIUS_FACTOR 4.0         /* Merger distance in units of dx_min */

/* Rezzolla formula fitting constants */
#define S4_REZZOLLA         (-0.129)
#define S5_REZZOLLA         (-0.384)
#define T0_REZZOLLA         (-2.686)
#define T2_REZZOLLA         (-3.454)
#define T3_REZZOLLA         (+2.353)

/*******************************************************************************
 * DATA STRUCTURES
 ******************************************************************************/

/* 3D Vector structure */
typedef struct {
    double x, y, z;
} Vec3;

/* Voronoi cell structure (FVM approach) */
typedef struct {
    int id;                         /* Cell ID */
    Vec3 center;                    /* Cell centroid position */
    double volume;                  /* Cell volume */
    
    /* Hydrodynamic quantities (conserved variables) */
    double density;                 /* Gas density */
    Vec3 momentum;                  /* Momentum density */
    double energy;                  /* Total energy density */
    double metallicity;             /* Metal mass fraction */
    
    /* Derived quantities */
    Vec3 velocity;                  /* Gas velocity */
    double pressure;                /* Gas pressure */
    double temperature;             /* Gas temperature */
    double sound_speed;             /* Sound speed */
    
    /* Stellar properties */
    double stellar_density;         /* Local stellar density */
    double velocity_dispersion;     /* Stellar velocity dispersion */
    
    /* Neighbors for flux calculation */
    int n_neighbors;
    int neighbors[MAX_NEIGHBORS];
    double face_areas[MAX_NEIGHBORS];
    Vec3 face_normals[MAX_NEIGHBORS];
    
    /* Flags */
    bool is_leaf;                   /* True if not refined */
    bool has_sink;                  /* True if contains a sink */
} VoronoiCell;

/* Sink particle (Black Hole) structure */
typedef struct {
    int id;                         /* Unique sink ID */
    bool active;                    /* Is this sink active? */
    
    /* Position and velocity */
    Vec3 position;
    Vec3 velocity;
    
    /* Mass properties */
    double mass;                    /* Current BH mass */
    double dM_bondi;                /* Bondi accretion rate */
    double dM_eddington;            /* Eddington accretion rate */
    double dM_accreted;             /* Actually accreted mass this timestep */
    
    /* Spin properties */
    double spin_magnitude;          /* Signed spin parameter a (-1 to 1) */
    Vec3 spin_direction;            /* Spin axis unit vector */
    double radiative_efficiency;    /* Radiative efficiency epsilon_r */
    
    /* Angular momentum from gas */
    Vec3 gas_angular_momentum;
    
    /* AGN feedback */
    double E_save;                  /* Stored energy for delayed feedback */
    double E_AGN;                   /* AGN energy output */
    bool ok_blast;                  /* Ready for feedback blast */
    
    /* Weighted averages for accretion */
    double weighted_density;
    double weighted_c2;
    Vec3 weighted_velocity;
    double weighted_volume;
    
    /* Birth properties */
    double birth_time;
    int birth_level;
} SinkParticle;

/* Simulation state structure */
typedef struct {
    /* Grid */
    int n_cells;
    VoronoiCell *cells;
    
    /* Sink particles */
    int n_sinks;
    SinkParticle sinks[MAX_SINKS];
    
    /* Time stepping */
    double time;
    double dt;
    double dt_coarse;
    
    /* Unit conversions */
    double scale_l;                 /* Length scale [cm] */
    double scale_d;                 /* Density scale [g/cm^3] */
    double scale_t;                 /* Time scale [s] */
    double scale_v;                 /* Velocity scale [cm/s] */
    double scale_m;                 /* Mass scale [g] */
    double scale_T2;                /* Temperature scale */
    double scale_nH;                /* Number density scale */
    
    /* Physics flags */
    bool use_spin;
    bool use_mad_jet;
    bool use_self_gravity;
    
    /* Cosmology */
    double aexp;                    /* Expansion factor */
    double hubble;
    
    /* Box parameters */
    double box_size;
    double dx_min;                  /* Minimum cell size */
} SimulationState;

/* Global simulation state */
static SimulationState sim;

/*******************************************************************************
 * UTILITY FUNCTIONS
 ******************************************************************************/

/* Vector operations */
static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 vec3_scale(Vec3 v, double s) {
    return (Vec3){v.x * s, v.y * s, v.z * s};
}

static inline double vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static inline double vec3_mag(Vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static inline Vec3 vec3_normalize(Vec3 v) {
    double mag = vec3_mag(v);
    if (mag > 0) {
        return vec3_scale(v, 1.0 / mag);
    }
    return (Vec3){0, 0, 0};
}

/* Minimum/Maximum */
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }

/*******************************************************************************
 * ISCO AND RADIATIVE EFFICIENCY CALCULATIONS
 ******************************************************************************/

/*
 * Compute the Innermost Stable Circular Orbit (ISCO) radius
 * for a Kerr black hole with spin parameter a.
 * 
 * Based on Bardeen (1970) formula.
 * 
 * Input:  spin_a - signed spin parameter (-1 to 1)
 * Output: r_isco - ISCO radius in units of GM/c^2
 */
double compute_isco_radius(double spin_a) {
    double a = fabs(spin_a);
    double a2 = a * a;
    
    /* Auxiliary functions Z1 and Z2 (Bardeen 1970) */
    double Z1 = 1.0 + pow(1.0 - a2, 1.0/3.0) * 
                (pow(1.0 + a, 1.0/3.0) + pow(1.0 - a, 1.0/3.0));
    double Z2 = sqrt(3.0 * a2 + Z1 * Z1);
    
    /* ISCO radius depends on spin sign (prograde/retrograde) */
    double r_isco;
    if (spin_a >= 0) {
        /* Prograde orbit */
        r_isco = 3.0 + Z2 - sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));
    } else {
        /* Retrograde orbit */
        r_isco = 3.0 + Z2 + sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));
    }
    
    return r_isco;
}

/*
 * Compute radiative efficiency from ISCO radius.
 * Based on Novikov-Thorne thin disc model.
 * 
 * epsilon_r = 1 - sqrt(1 - 2/(3*r_isco))
 * 
 * Range: ~0.057 (a=0) to ~0.32 (a=0.998)
 */
double compute_radiative_efficiency(double spin_a) {
    double r_isco = compute_isco_radius(spin_a);
    return 1.0 - sqrt(1.0 - 2.0 / (3.0 * r_isco));
}

/*******************************************************************************
 * BONDI-HOYLE ACCRETION
 ******************************************************************************/

/*
 * Compute Bondi-Hoyle accretion rate for a sink particle.
 * 
 * Formula: dM/dt = alpha * 4*pi*G^2*M^2*rho / (c_s^2 + v_rel^2)^(3/2)
 * 
 * where:
 *   alpha = boost factor = max((rho/rho_star)^boost_acc, 1)
 *   rho = gas density (mass-weighted average)
 *   c_s = sound speed
 *   v_rel = relative velocity between BH and gas
 */
void compute_bondi_accretion_rate(SinkParticle *sink) {
    /* Skip if no volume sampled */
    if (sink->weighted_volume <= 0) {
        sink->dM_bondi = 0;
        sink->dM_eddington = 0;
        return;
    }
    
    /* Compute mass-weighted averages */
    double density = sink->weighted_density / sink->weighted_volume;
    double c2_mean = sink->weighted_c2 / sink->weighted_volume / density;
    Vec3 gas_vel = vec3_scale(sink->weighted_velocity, 
                              1.0 / (sink->weighted_volume * density));
    
    /* Relative velocity squared (capped at sigmav_max) */
    Vec3 v_rel = vec3_sub(gas_vel, sink->velocity);
    double v2_mean = vec3_dot(v_rel, v_rel);
    double sigmav2 = pow(SIGMAV_MAX * 1e5 / sim.scale_v, 2);
    v2_mean = dmin(v2_mean, sigmav2);
    
    /* Boost factor (Booth & Schaye 2009) */
    double d_star = N_STAR_THRESHOLD / sim.scale_nH;
    double alpha = dmax(pow(density / d_star, BOOST_ACC), 1.0);
    
    /* Bondi-Hoyle accretion rate */
    double factG = GRAVITY_CGS / (sim.scale_l * sim.scale_v * sim.scale_v);
    double bondi_rate = alpha * FOURPI * density * 
                        pow(factG * sink->mass, 2) / 
                        pow(c2_mean + v2_mean, 1.5);
    
    /* Update radiative efficiency if using spin */
    if (sim.use_spin) {
        sink->radiative_efficiency = compute_radiative_efficiency(sink->spin_magnitude);
    } else {
        sink->radiative_efficiency = 0.1;
    }
    
    /* Eddington accretion rate */
    double prefact = FOURPI * GRAVITY_CGS * sim.scale_m * PROTON_MASS / 
                     (SIGMA_THOMSON * CLIGHT) / (sim.scale_m / sim.scale_t);
    double eddington_rate = prefact * sink->mass / sink->radiative_efficiency;
    
    /* Apply sub-Eddington boost correction */
    if (bondi_rate / eddington_rate < X_FLOOR) {
        bondi_rate *= F_BONDI;
    }
    
    /* Store results */
    sink->dM_bondi = bondi_rate;
    sink->dM_eddington = eddington_rate;
    
    /* If energy is stored (delayed feedback), prevent new accretion */
    if (sink->E_save > 0) {
        sink->dM_bondi = 0;
    }
}

/*******************************************************************************
 * GAS WEIGHTING FOR SINK PARTICLES
 ******************************************************************************/

/*
 * Accumulate weighted gas properties around sink particle.
 * Uses Voronoi cell volumes for proper FVM integration.
 * 
 * This replaces the RAMSES kernel-weighted approach with
 * a volume-weighted Voronoi approach.
 */
void accumulate_gas_properties(SinkParticle *sink, double r_kernel) {
    /* Reset weighted quantities */
    sink->weighted_density = 0;
    sink->weighted_c2 = 0;
    sink->weighted_velocity = (Vec3){0, 0, 0};
    sink->weighted_volume = 0;
    sink->gas_angular_momentum = (Vec3){0, 0, 0};
    
    double r_kernel2 = r_kernel * r_kernel;
    
    /* Loop over all Voronoi cells */
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        
        if (!cell->is_leaf) continue;
        
        /* Distance from sink to cell center */
        Vec3 dr = vec3_sub(cell->center, sink->position);
        
        /* Apply periodic boundary conditions */
        if (dr.x > 0.5 * sim.box_size) dr.x -= sim.box_size;
        if (dr.x < -0.5 * sim.box_size) dr.x += sim.box_size;
        if (dr.y > 0.5 * sim.box_size) dr.y -= sim.box_size;
        if (dr.y < -0.5 * sim.box_size) dr.y += sim.box_size;
        if (dr.z > 0.5 * sim.box_size) dr.z -= sim.box_size;
        if (dr.z < -0.5 * sim.box_size) dr.z += sim.box_size;
        
        double r2 = vec3_dot(dr, dr);
        
        /* Only include cells within kernel radius */
        if (r2 > r_kernel2) continue;
        
        double vol = cell->volume;
        double d = cell->density;
        
        /* Mass-weighted accumulation */
        sink->weighted_density += d * vol;
        sink->weighted_c2 += cell->sound_speed * cell->sound_speed * d * vol;
        sink->weighted_velocity = vec3_add(sink->weighted_velocity,
                                          vec3_scale(cell->velocity, d * vol));
        sink->weighted_volume += vol;
        
        /* Angular momentum for spin evolution */
        Vec3 L = vec3_cross(dr, vec3_scale(cell->velocity, d * vol));
        sink->gas_angular_momentum = vec3_add(sink->gas_angular_momentum, L);
    }
}

/*******************************************************************************
 * MASS ACCRETION FROM GAS CELLS
 ******************************************************************************/

/*
 * Remove accreted mass from Voronoi cells (FVM conservative update).
 * Mass is removed proportionally to the cell's contribution to the
 * weighted density.
 */
void accrete_mass_from_cells(SinkParticle *sink, double dm_total, double r_kernel) {
    if (dm_total <= 0 || sink->weighted_density <= 0) return;
    
    double r_kernel2 = r_kernel * r_kernel;
    
    /* Loop over cells and remove mass proportionally */
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        
        if (!cell->is_leaf) continue;
        
        Vec3 dr = vec3_sub(cell->center, sink->position);
        
        /* Periodic BCs */
        if (dr.x > 0.5 * sim.box_size) dr.x -= sim.box_size;
        if (dr.x < -0.5 * sim.box_size) dr.x += sim.box_size;
        if (dr.y > 0.5 * sim.box_size) dr.y -= sim.box_size;
        if (dr.y < -0.5 * sim.box_size) dr.y += sim.box_size;
        if (dr.z > 0.5 * sim.box_size) dr.z -= sim.box_size;
        if (dr.z < -0.5 * sim.box_size) dr.z += sim.box_size;
        
        double r2 = vec3_dot(dr, dr);
        if (r2 > r_kernel2) continue;
        
        /* Fraction of mass to remove from this cell */
        double weight = cell->density * cell->volume / sink->weighted_density;
        double dm_cell = dm_total * weight;
        
        /* Ensure we don't remove more than available */
        double max_dm = 0.5 * cell->density * cell->volume;
        dm_cell = dmin(dm_cell, max_dm);
        
        /* Update conserved variables (FVM) */
        double d_old = cell->density;
        double d_new = d_old - dm_cell / cell->volume;
        
        if (d_new > 0) {
            /* Scale momentum to conserve velocity */
            double ratio = d_new / d_old;
            cell->density = d_new;
            cell->momentum = vec3_scale(cell->momentum, ratio);
            
            /* Compute new internal + kinetic energy */
            double ekin_old = 0.5 * vec3_dot(cell->momentum, cell->momentum) / d_old;
            double ekin_new = 0.5 * vec3_dot(cell->momentum, cell->momentum) / d_new;
            cell->energy = cell->energy - ekin_old + ekin_new;
            cell->energy *= ratio;
            
            /* Update metallicity */
            cell->metallicity *= ratio;
        }
    }
}

/*******************************************************************************
 * SPIN EVOLUTION - WARPED DISC MODEL
 ******************************************************************************/

/*
 * Evolve black hole spin due to gas accretion.
 * 
 * Based on Fanidakis et al. (2011) warped disc model:
 * - Compute warp radius and disc angular momentum
 * - Check for alignment or anti-alignment
 * - Update spin magnitude using ISCO-based formula
 * 
 * Also handles self-gravity regime (Dotti et al. 2013).
 */
void evolve_spin_accretion(SinkParticle *sink, double dm_coarse) {
    if (!sim.use_spin) return;
    if (dm_coarse <= 0) return;
    
    double mbh = sink->mass * sim.scale_m;  /* BH mass in grams */
    double chi = sink->dM_bondi / sink->dM_eddington;  /* Eddington ratio */
    
    if (chi <= 0) return;
    
    /* Gas angular momentum direction */
    Vec3 j_gas = vec3_normalize(sink->gas_angular_momentum);
    
    /* Minimum mass increment for iteration */
    double dm_ini = dm_coarse * sim.scale_m;
    double dm_remaining = dm_ini;
    
    while (dm_remaining > 0) {
        double amod = fabs(sink->spin_magnitude);
        
        if (amod > 0) {
            /* Compute ISCO and efficiency */
            double r_lso = compute_isco_radius(sink->spin_magnitude);
            double epsilon_r = compute_radiative_efficiency(sink->spin_magnitude);
            double chi_eps = chi / (epsilon_r / 0.1);
            
            /* BH mass in units of 10^8 Msun */
            double mass_8 = mbh / (1e8 * MSUN);
            
            /* Viscous timescale (Fanidakis et al. 2011, eq. 19) */
            double t_nu1 = 5.3e5 * pow(amod, 0.875) * pow(chi_eps, -0.75) * 
                          pow(NU_RATIO, -0.875) * pow(ALPHA_VISCOSITY / 0.1, -1.5) * 
                          pow(mass_8, 1.375);  /* in years */
            
            /* Warp radius (eq. 15) */
            double r_warp = 6.4e3 * pow(amod, 0.625) * pow(mass_8, 0.125) * 
                           pow(chi_eps, -0.25) * pow(NU_RATIO, -0.625) / 
                           sqrt(ALPHA_VISCOSITY / 0.1);
            
            /* Disc mass and effective radius */
            double m_disc, r_eff;
            
            if (sim.use_self_gravity) {
                /* Self-gravity radius (Dotti et al. 2013) */
                double r_sg = 5e2 * pow(ALPHA_VISCOSITY / 0.1, 28.0/45.0) * 
                             pow(mass_8, -52.0/45.0) * pow(chi_eps, -22.0/45.0);
                
                if (r_sg < r_warp) {
                    /* Self-gravity dominated */
                    double m_sg = 6e5 * pow(ALPHA_VISCOSITY / 0.1, -1.0/45.0) * 
                                 pow(mass_8, 34.0/45.0) * pow(chi_eps, 4.0/45.0) * MSUN;
                    m_disc = m_sg;
                    r_eff = r_sg;
                } else {
                    /* Warp dominated */
                    double dm_acc = sink->dM_bondi * sim.scale_m / sim.scale_t * chi;
                    m_disc = dm_acc * (t_nu1 * YEAR_TO_SEC);
                    r_eff = r_warp;
                }
            } else {
                double dm_acc = sink->dM_bondi * sim.scale_m / sim.scale_t * chi;
                m_disc = dm_acc * (t_nu1 * YEAR_TO_SEC);
                r_eff = r_warp;
            }
            
            /* BH spin direction */
            Vec3 j_bh = sink->spin_direction;
            
            /* Angle between spin and gas angular momentum */
            double cos_theta = vec3_dot(j_bh, j_gas);
            
            /* Mass to process in this iteration */
            double dm;
            if (sink->spin_magnitude >= MAX_SPIN && cos_theta >= 0) {
                dm = dm_remaining;  /* Already max spin, aligned */
            } else {
                dm = dmin(dmax(m_disc, dm_ini * 1e-3), dm_remaining);
            }
            
            /* Ratio of disc to BH angular momentum */
            double aac = dm / mbh * sqrt(r_eff) / amod;
            double jd = 2.0 * aac;  /* Disc AM in units of BH AM */
            
            /* Check for anti-alignment condition */
            if (-aac > cos_theta) {
                /* Anti-alignment: BH spin flips */
                Vec3 j_new = vec3_add(j_bh, vec3_scale(j_gas, jd));
                sink->spin_direction = vec3_normalize(j_new);
                if (sink->spin_magnitude > 0) {
                    sink->spin_magnitude = -sink->spin_magnitude;
                }
            } else {
                /* Alignment: spin aligns with total AM */
                Vec3 j_new = vec3_add(j_bh, vec3_scale(j_gas, jd));
                sink->spin_direction = vec3_normalize(j_new);
                if (sink->spin_magnitude < 0) {
                    sink->spin_magnitude = -sink->spin_magnitude;
                }
            }
            
            /* Update spin magnitude */
            if (sim.use_mad_jet && chi < X_FLOOR) {
                /* MAD jet spinup (McKinney et al. 2012) */
                double dm_deplete = dm;
                double mbh2 = mbh;
                while (dm_deplete > 0) {
                    double a = sink->spin_magnitude;
                    double spinup = 0.97166 - 12.0026 * a - 4.04337 * a * a + 
                                   5.81317 * a * a * a + 2.50482 * a * a * a * a;
                    double dm2 = dmin(0.01 * mbh2 / (fabs(spinup) + 0.01), dm_deplete);
                    sink->spin_magnitude += spinup * dm2 / mbh2;
                    mbh2 += dm2;
                    dm_deplete -= dm2;
                }
            } else {
                /* Standard thin disc spinup */
                double r_lso_new = compute_isco_radius(sink->spin_magnitude);
                double mass_ratio = mbh / (mbh + dm);
                
                if ((mbh + dm) / mbh <= sqrt(r_lso_new)) {
                    double new_spin = (1.0/3.0) * sqrt(r_lso_new) * mass_ratio * 
                                     (4.0 - sqrt(3.0 * r_lso_new * mass_ratio * mass_ratio - 2.0));
                    sink->spin_magnitude = dmin(new_spin, MAX_SPIN);
                } else {
                    sink->spin_magnitude = MAX_SPIN;
                }
            }
            
            /* Ensure spin stays in bounds */
            if (sink->spin_magnitude > MAX_SPIN) sink->spin_magnitude = MAX_SPIN;
            if (sink->spin_magnitude < -MAX_SPIN) sink->spin_magnitude = -MAX_SPIN;
            
            mbh += dm;
            dm_remaining -= dm;
            
        } else {
            double dm;
            /* Zero spin: acquire gas AM direction */
            dm = dm_remaining;
            sink->spin_direction = j_gas;
            dm_remaining = 0;
        }
    }
    
    /* Update radiative efficiency */
    sink->radiative_efficiency = compute_radiative_efficiency(sink->spin_magnitude);
}

/*******************************************************************************
 * BLACK HOLE MERGER
 ******************************************************************************/

/*
 * Check if two sinks should merge based on distance and velocity.
 */
bool should_merge(SinkParticle *s1, SinkParticle *s2) {
    Vec3 dr = vec3_sub(s1->position, s2->position);
    
    /* Periodic BCs */
    if (dr.x > 0.5 * sim.box_size) dr.x -= sim.box_size;
    if (dr.x < -0.5 * sim.box_size) dr.x += sim.box_size;
    if (dr.y > 0.5 * sim.box_size) dr.y -= sim.box_size;
    if (dr.y < -0.5 * sim.box_size) dr.y += sim.box_size;
    if (dr.z > 0.5 * sim.box_size) dr.z -= sim.box_size;
    if (dr.z < -0.5 * sim.box_size) dr.z += sim.box_size;
    
    double r = vec3_mag(dr);
    double r_merge = MERGE_RADIUS_FACTOR * sim.dx_min;
    
    if (r >= r_merge) return false;
    
    /* Optionally check if gravitationally bound */
    Vec3 dv = vec3_sub(s1->velocity, s2->velocity);
    double v2 = vec3_dot(dv, dv);
    double E_kin = 0.5 * (s1->mass * s2->mass / (s1->mass + s2->mass)) * v2;
    
    double factG = GRAVITY_CGS / (sim.scale_l * sim.scale_v * sim.scale_v);
    double E_grav = factG * s1->mass * s2->mass / r;
    
    return (E_kin < E_grav);
}

/*
 * Compute final spin after BH-BH merger using Rezzolla et al. (2008) formula.
 * 
 * Inputs: masses M1, M2 (M1 >= M2), spins a1, a2, spin directions, orbital AM
 * Output: final spin magnitude and direction
 */
void compute_merger_spin(SinkParticle *s1, SinkParticle *s2,
                         double *a_final, Vec3 *spin_final) {
    /* Ensure M1 >= M2 */
    SinkParticle *primary = (s1->mass >= s2->mass) ? s1 : s2;
    SinkParticle *secondary = (s1->mass >= s2->mass) ? s2 : s1;
    
    double M1 = primary->mass;
    double M2 = secondary->mass;
    double q = M2 / M1;  /* Mass ratio <= 1 */
    double mu = q / pow(1.0 + q, 2);  /* Reduced mass ratio */
    
    double a1 = fabs(primary->spin_magnitude);
    double a2 = fabs(secondary->spin_magnitude);
    Vec3 ax1 = primary->spin_direction;
    Vec3 ax2 = secondary->spin_direction;
    
    /* Compute orbital angular momentum */
    Vec3 x1 = vec3_sub(primary->position, 
                       vec3_scale(vec3_add(vec3_scale(primary->position, M1),
                                          vec3_scale(secondary->position, M2)),
                                 1.0 / (M1 + M2)));
    Vec3 x2 = vec3_sub(secondary->position,
                       vec3_scale(vec3_add(vec3_scale(primary->position, M1),
                                          vec3_scale(secondary->position, M2)),
                                 1.0 / (M1 + M2)));
    Vec3 v1 = vec3_sub(primary->velocity,
                       vec3_scale(vec3_add(vec3_scale(primary->velocity, M1),
                                          vec3_scale(secondary->velocity, M2)),
                                 1.0 / (M1 + M2)));
    Vec3 v2 = vec3_sub(secondary->velocity,
                       vec3_scale(vec3_add(vec3_scale(primary->velocity, M1),
                                          vec3_scale(secondary->velocity, M2)),
                                 1.0 / (M1 + M2)));
    
    Vec3 L1 = vec3_scale(vec3_cross(x1, v1), M1);
    Vec3 L2 = vec3_scale(vec3_cross(x2, v2), M2);
    Vec3 L = vec3_add(L1, L2);
    double Lmod = vec3_mag(L);
    Vec3 L_hat = (Lmod > 0) ? vec3_normalize(L) : (Vec3){0, 0, 1};
    
    /* Spin-orbit alignment angles */
    double a1L = (a1 > 0 && Lmod > 0) ? vec3_dot(ax1, L_hat) : 0;
    double a2L = (a2 > 0 && Lmod > 0) ? vec3_dot(ax2, L_hat) : 0;
    
    /* Spin-spin alignment */
    double a1a2 = vec3_dot(ax1, ax2);
    
    /* Rezzolla formula for orbital angular momentum parameter */
    double Lmodana = S4_REZZOLLA / pow(1.0 + q * q, 2) * 
                    (a1 * a1 + a2 * a2 * pow(q, 4) + 2.0 * a1 * a2 * q * q * a1a2) +
                    (S5_REZZOLLA * mu + T0_REZZOLLA + 2.0) / (1.0 + q * q) * 
                    (a1 * a1L + a2 * q * q * a2L) +
                    2.0 * sqrt(3.0) + T2_REZZOLLA * mu + T3_REZZOLLA * mu * mu;
    
    /* Final spin magnitude */
    double af = 1.0 / pow(1.0 + q, 2) * 
               sqrt(a1 * a1 + a2 * a2 * pow(q, 4) + 2.0 * a1 * a2 * q * q * a1a2 +
                    2.0 * (a1 * a1L + a2 * q * q * a2L) * Lmodana * q + 
                    Lmodana * Lmodana * q * q);
    
    /* Final spin direction */
    Vec3 sf;
    sf.x = 1.0 / pow(1.0 + q, 2) * (a1 * ax1.x + a2 * ax2.x * q * q + L_hat.x * Lmodana * q);
    sf.y = 1.0 / pow(1.0 + q, 2) * (a1 * ax1.y + a2 * ax2.y * q * q + L_hat.y * Lmodana * q);
    sf.z = 1.0 / pow(1.0 + q, 2) * (a1 * ax1.z + a2 * ax2.z * q * q + L_hat.z * Lmodana * q);
    
    /* Enforce spin limit */
    if (af > MAX_SPIN) af = MAX_SPIN;
    if (af < -MAX_SPIN) af = -MAX_SPIN;
    
    *a_final = af;
    *spin_final = vec3_normalize(sf);
}

/*
 * Merge two sink particles into one.
 */
void merge_sinks(SinkParticle *primary, SinkParticle *secondary) {
    double M1 = primary->mass;
    double M2 = secondary->mass;
    double M_total = M1 + M2;
    
    /* Mass-weighted position (with periodic BCs) */
    Vec3 dr = vec3_sub(secondary->position, primary->position);
    if (dr.x > 0.5 * sim.box_size) dr.x -= sim.box_size;
    if (dr.x < -0.5 * sim.box_size) dr.x += sim.box_size;
    if (dr.y > 0.5 * sim.box_size) dr.y -= sim.box_size;
    if (dr.y < -0.5 * sim.box_size) dr.y += sim.box_size;
    if (dr.z > 0.5 * sim.box_size) dr.z -= sim.box_size;
    if (dr.z < -0.5 * sim.box_size) dr.z += sim.box_size;
    
    primary->position = vec3_add(primary->position, vec3_scale(dr, M2 / M_total));
    
    /* Wrap to box */
    if (primary->position.x < 0) primary->position.x += sim.box_size;
    if (primary->position.x >= sim.box_size) primary->position.x -= sim.box_size;
    if (primary->position.y < 0) primary->position.y += sim.box_size;
    if (primary->position.y >= sim.box_size) primary->position.y -= sim.box_size;
    if (primary->position.z < 0) primary->position.z += sim.box_size;
    if (primary->position.z >= sim.box_size) primary->position.z -= sim.box_size;
    
    /* Momentum-conserving velocity */
    primary->velocity = vec3_scale(
        vec3_add(vec3_scale(primary->velocity, M1),
                 vec3_scale(secondary->velocity, M2)),
        1.0 / M_total);
    
    /* Compute final spin */
    double a_final;
    Vec3 spin_final;
    compute_merger_spin(primary, secondary, &a_final, &spin_final);
    
    primary->spin_magnitude = a_final;
    primary->spin_direction = spin_final;
    
    /* Update mass and other quantities */
    primary->mass = M_total;
    primary->dM_accreted += secondary->dM_accreted;
    primary->E_save += secondary->E_save;
    
    /* Deactivate secondary */
    secondary->active = false;
}

/*******************************************************************************
 * AGN FEEDBACK
 ******************************************************************************/

/*
 * Compute MAD jet efficiency based on McKinney et al. (2012).
 * 4th order polynomial fit to jet+wind power.
 */
double compute_mad_efficiency(double spin_a) {
    double a = fabs(spin_a);
    double eff = 4.10507 + 0.328712 * a + 76.0849 * a * a + 
                47.9235 * a * a * a + 3.86634 * a * a * a * a;
    return eff / 100.0;
}

/*
 * Apply AGN feedback to surrounding gas cells.
 * 
 * Two modes:
 * 1. Quasar mode (high Eddington ratio): Isotropic thermal feedback
 * 2. Radio/Jet mode (low Eddington ratio): Bipolar kinetic+thermal feedback
 */
void apply_agn_feedback(SinkParticle *sink) {
    if (!sink->ok_blast && sink->E_save <= 0) return;
    
    double dm_smbh = sink->dM_accreted * sim.scale_m;
    double epsilon_r = sink->radiative_efficiency;
    double X_radio = sink->dM_bondi / sink->dM_eddington;
    
    double E_AGN;
    double u_blast = 0;
    bool jet_mode = (X_radio < X_FLOOR);
    
    /* Compute feedback energy */
    if (sink->E_save > 0) {
        /* Use stored energy */
        E_AGN = sink->E_save;
        sink->E_save = 0;
    } else if (jet_mode) {
        /* Radio/Jet mode */
        if (sim.use_mad_jet) {
            double eff_mad = compute_mad_efficiency(sink->spin_magnitude);
            E_AGN = eff_mad * dm_smbh * CLIGHT * CLIGHT;
        } else {
            E_AGN = EAGN_KINETIC * epsilon_r * dm_smbh * CLIGHT * CLIGHT;
        }
    } else {
        /* Quasar/Thermal mode */
        E_AGN = EAGN_THERMAL * epsilon_r * dm_smbh * CLIGHT * CLIGHT;
    }
    
    /* Convert to code units */
    E_AGN /= (sim.scale_m * sim.scale_v * sim.scale_v);
    
    /* Feedback radius */
    double r_max = dmax(sim.dx_min * sim.scale_l / sim.aexp, R_AGN * 3.08e21);
    r_max /= sim.scale_l;
    double r_max2 = r_max * r_max;
    
    /* Jet direction (from gas angular momentum) */
    Vec3 j_hat = vec3_normalize(sink->gas_angular_momentum);
    
    /* Collect gas properties within feedback region */
    double m_gas = 0;
    double vol_gas = 0;
    
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        if (!cell->is_leaf) continue;
        
        Vec3 dr = vec3_sub(cell->center, sink->position);
        
        /* Periodic BCs */
        if (dr.x > 0.5 * sim.box_size) dr.x -= sim.box_size;
        if (dr.x < -0.5 * sim.box_size) dr.x += sim.box_size;
        if (dr.y > 0.5 * sim.box_size) dr.y -= sim.box_size;
        if (dr.y < -0.5 * sim.box_size) dr.y += sim.box_size;
        if (dr.z > 0.5 * sim.box_size) dr.z -= sim.box_size;
        if (dr.z < -0.5 * sim.box_size) dr.z += sim.box_size;
        
        double r2 = vec3_dot(dr, dr);
        
        if (jet_mode) {
            /* Check if in jet cylinder */
            double dz_jet = vec3_dot(dr, j_hat);
            double dr_jet = sqrt(r2 - dz_jet * dz_jet);
            
            if (dr_jet <= r_max && fabs(dz_jet) <= r_max) {
                m_gas += cell->density * cell->volume;
                vol_gas += cell->volume;
            }
        } else {
            /* Spherical region for thermal feedback */
            if (r2 <= r_max2) {
                m_gas += cell->density * cell->volume;
                vol_gas += cell->volume;
            }
        }
    }
    
    if (vol_gas <= 0) {
        /* Store energy for later */
        sink->E_save += E_AGN;
        return;
    }
    
    /* Compute blast velocity for jet mode */
    if (jet_mode && m_gas > 0) {
        u_blast = sqrt(F_EK_AGN * E_AGN / m_gas);
    }
    
    /* Pressure for thermal component */
    double p_gas = (jet_mode ? (1.0 - F_EK_AGN) : 1.0) * E_AGN / vol_gas;
    
    /* Apply feedback to cells */
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        if (!cell->is_leaf) continue;
        
        Vec3 dr = vec3_sub(cell->center, sink->position);
        
        /* Periodic BCs */
        if (dr.x > 0.5 * sim.box_size) dr.x -= sim.box_size;
        if (dr.x < -0.5 * sim.box_size) dr.x += sim.box_size;
        if (dr.y > 0.5 * sim.box_size) dr.y -= sim.box_size;
        if (dr.y < -0.5 * sim.box_size) dr.y += sim.box_size;
        if (dr.z > 0.5 * sim.box_size) dr.z -= sim.box_size;
        if (dr.z < -0.5 * sim.box_size) dr.z += sim.box_size;
        
        double r2 = vec3_dot(dr, dr);
        bool in_feedback_region = false;
        double dz_jet = 0;
        double psy = 1.0;  /* Weight function */
        
        if (jet_mode) {
            dz_jet = vec3_dot(dr, j_hat);
            double dr_jet = sqrt(r2 - dz_jet * dz_jet);
            
            if (dr_jet <= r_max && fabs(dz_jet) <= r_max) {
                in_feedback_region = true;
                /* Gaussian profile perpendicular to jet */
                psy = exp(-dr_jet * dr_jet / (2.0 * r_max * r_max));
            }
        } else {
            if (r2 <= r_max2) {
                in_feedback_region = true;
            }
        }
        
        if (!in_feedback_region) continue;
        
        double d = cell->density;
        
        /* Compute current internal energy */
        double ekin = 0.5 * vec3_dot(cell->momentum, cell->momentum) / d;
        double eint = cell->energy - ekin;
        double T2_old = (1.4 - 1.0) * eint / d * sim.scale_T2;
        
        if (jet_mode) {
            /* Jet mode: Add kinetic energy along jet axis */
            Vec3 v_jet;
            if (dz_jet < 0) {
                v_jet = vec3_scale(j_hat, -u_blast);
            } else {
                v_jet = vec3_scale(j_hat, u_blast);
            }
            
            /* Add mass from jet (optional) */
            double d_jet = m_gas / vol_gas * psy;
            cell->density += d_jet;
            
            /* Add momentum */
            cell->momentum = vec3_add(cell->momentum, vec3_scale(v_jet, d_jet));
            
            /* Add thermal energy (limited by T2_MAX_AGN) */
            double dE_thermal = p_gas * d * psy;
            double T2_new = (1.4 - 1.0) * (eint + dE_thermal) / cell->density * sim.scale_T2;
            
            if (T2_new <= T2_MAX_AGN && T2_old < T2_MAX_AGN) {
                cell->energy += dE_thermal;
            } else if (T2_old < T2_MAX_AGN) {
                /* Cap at maximum temperature */
                double dE_max = T2_MAX_AGN / sim.scale_T2 / (1.4 - 1.0) * cell->density - eint;
                cell->energy += dE_max;
                sink->E_save += (dE_thermal - dE_max) * cell->volume;
            }
            
            /* Update kinetic energy for momentum change */
            double ekin_new = 0.5 * vec3_dot(cell->momentum, cell->momentum) / cell->density;
            cell->energy += (ekin_new - ekin);
            
        } else {
            /* Thermal mode: Add only thermal energy */
            double dE_thermal = p_gas * d;
            double T2_new = (1.4 - 1.0) * (eint + dE_thermal) / d * sim.scale_T2;
            
            if (T2_new <= T2_MAX_AGN && T2_old < T2_MAX_AGN) {
                cell->energy += dE_thermal;
            } else if (T2_old < T2_MAX_AGN) {
                double dE_max = T2_MAX_AGN / sim.scale_T2 / (1.4 - 1.0) * d - eint;
                cell->energy += dE_max;
                sink->E_save += (dE_thermal - dE_max) * cell->volume;
            }
        }
    }
    
    sink->dM_accreted = 0;
}

/*******************************************************************************
 * SINK PARTICLE CREATION (BLACK HOLE SEEDING)
 ******************************************************************************/

/*
 * Check if a Voronoi cell satisfies conditions for BH seeding.
 * 
 * Conditions (Dubois et al. 2012):
 * - Gas density > n_sink threshold
 * - Stellar density > 0.1 H/cc
 * - Velocity dispersion > 20 km/s
 */
bool check_seeding_conditions(VoronoiCell *cell) {
    /* Density threshold */
    if (cell->density * sim.scale_nH < N_SINK_THRESHOLD) {
        return false;
    }
    
    /* Stellar density requirement */
    if (cell->stellar_density < STELLAR_DENSITY_MIN) {
        return false;
    }
    
    /* Velocity dispersion requirement */
    if (cell->velocity_dispersion < VELOCITY_DISP_MIN) {
        return false;
    }
    
    /* Check if already has a sink nearby */
    double r_excl = 4.0 * sim.dx_min;
    for (int i = 0; i < sim.n_sinks; i++) {
        if (!sim.sinks[i].active) continue;
        
        Vec3 dr = vec3_sub(cell->center, sim.sinks[i].position);
        if (vec3_mag(dr) < r_excl) {
            return false;
        }
    }
    
    return true;
}

/*
 * Create a new sink particle at the given cell.
 */
int create_sink(VoronoiCell *cell) {
    if (sim.n_sinks >= MAX_SINKS) {
        fprintf(stderr, "Warning: Maximum number of sinks reached\n");
        return -1;
    }
    
    int idx = sim.n_sinks;
    SinkParticle *sink = &sim.sinks[idx];
    
    /* Initialize sink properties */
    sink->id = idx + 1;
    sink->active = true;
    
    sink->position = cell->center;
    sink->velocity = cell->velocity;
    
    /* Seed mass */
    sink->mass = SEED_MASS * MSUN / sim.scale_m;
    sink->dM_bondi = 0;
    sink->dM_eddington = 0;
    sink->dM_accreted = 0;
    
    /* Initial spin (zero) */
    sink->spin_magnitude = 0;
    sink->spin_direction = (Vec3){0, 0, 1};
    sink->radiative_efficiency = 0.1;
    
    sink->gas_angular_momentum = (Vec3){0, 0, 0};
    
    sink->E_save = 0;
    sink->E_AGN = 0;
    sink->ok_blast = false;
    
    sink->weighted_density = 0;
    sink->weighted_c2 = 0;
    sink->weighted_velocity = (Vec3){0, 0, 0};
    sink->weighted_volume = 0;
    
    sink->birth_time = sim.time;
    sink->birth_level = 0;
    
    /* Remove seed mass from gas cell */
    double dm_seed = sink->mass;
    double ratio = (cell->density * cell->volume - dm_seed) / (cell->density * cell->volume);
    if (ratio > 0) {
        cell->density *= ratio;
        cell->momentum = vec3_scale(cell->momentum, ratio);
        cell->energy *= ratio;
    }
    
    cell->has_sink = true;
    sim.n_sinks++;
    
    printf("Created sink %d at (%.3f, %.3f, %.3f) with mass %.2e Msun\n",
           sink->id, sink->position.x, sink->position.y, sink->position.z,
           sink->mass * sim.scale_m / MSUN);
    
    return idx;
}

/*******************************************************************************
 * MAIN INTEGRATION ROUTINES
 ******************************************************************************/

/*
 * Main sink particle creation routine.
 * Loops over all Voronoi cells and creates sinks where conditions are met.
 */
void create_sink_particles(void) {
    printf("Checking for new sink particle formation...\n");
    
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        
        if (!cell->is_leaf) continue;
        if (cell->has_sink) continue;
        
        if (check_seeding_conditions(cell)) {
            create_sink(cell);
        }
    }
}

/*
 * Merge sink particles that are close together.
 * Uses FOF-like algorithm.
 */
void merge_sink_particles(void) {
    printf("Checking for sink particle mergers...\n");
    
    for (int i = 0; i < sim.n_sinks; i++) {
        if (!sim.sinks[i].active) continue;
        
        for (int j = i + 1; j < sim.n_sinks; j++) {
            if (!sim.sinks[j].active) continue;
            
            if (should_merge(&sim.sinks[i], &sim.sinks[j])) {
                printf("Merging sinks %d and %d\n", sim.sinks[i].id, sim.sinks[j].id);
                
                /* Merge into more massive one */
                if (sim.sinks[i].mass >= sim.sinks[j].mass) {
                    merge_sinks(&sim.sinks[i], &sim.sinks[j]);
                } else {
                    merge_sinks(&sim.sinks[j], &sim.sinks[i]);
                }
            }
        }
    }
    
    /* Compact the sink array */
    int n_active = 0;
    for (int i = 0; i < sim.n_sinks; i++) {
        if (sim.sinks[i].active) {
            if (i != n_active) {
                sim.sinks[n_active] = sim.sinks[i];
            }
            n_active++;
        }
    }
    sim.n_sinks = n_active;
}

/*
 * Main Bondi-Hoyle accretion routine.
 * Computes accretion rates and removes mass from gas.
 */
void bondi_hoyle_accretion(void) {
    printf("Computing Bondi-Hoyle accretion...\n");
    
    double r_kernel = 4.0 * sim.dx_min;
    
    for (int i = 0; i < sim.n_sinks; i++) {
        SinkParticle *sink = &sim.sinks[i];
        if (!sink->active) continue;
        
        /* Accumulate gas properties */
        accumulate_gas_properties(sink, r_kernel);
        
        /* Compute accretion rate */
        compute_bondi_accretion_rate(sink);
        
        /* Actual mass to accrete (limited by Eddington) */
        double dm = dmin(sink->dM_bondi, sink->dM_eddington) * sim.dt;
        
        /* Remove mass from gas cells */
        accrete_mass_from_cells(sink, dm, r_kernel);
        
        /* Update sink mass and spin */
        sink->dM_accreted = dm;
        sink->mass += dm;
        
        /* Evolve spin */
        evolve_spin_accretion(sink, dm);
        
        printf("  Sink %d: M=%.2e Msun, dM/dt=%.2e Msun/yr, spin=%.3f\n",
               sink->id, 
               sink->mass * sim.scale_m / MSUN,
               sink->dM_bondi * sim.scale_m / sim.scale_t * YEAR_TO_SEC / MSUN,
               sink->spin_magnitude);
    }
}

/*
 * Apply AGN feedback from all sink particles.
 */
void apply_all_agn_feedback(void) {
    printf("Applying AGN feedback...\n");
    
    for (int i = 0; i < sim.n_sinks; i++) {
        SinkParticle *sink = &sim.sinks[i];
        if (!sink->active) continue;
        
        /* Check if feedback should occur */
        double X_radio = sink->dM_bondi / sink->dM_eddington;
        sink->ok_blast = (sink->dM_accreted > 0);
        
        apply_agn_feedback(sink);
        
        if (sink->ok_blast) {
            printf("  Sink %d: AGN feedback E=%.2e erg, mode=%s\n",
                   sink->id,
                   sink->E_AGN * sim.scale_m * sim.scale_v * sim.scale_v,
                   (X_radio < X_FLOOR) ? "JET" : "THERMAL");
        }
    }
}

/*
 * Update sink particle positions and velocities.
 */
void update_sink_positions(void) {
    for (int i = 0; i < sim.n_sinks; i++) {
        SinkParticle *sink = &sim.sinks[i];
        if (!sink->active) continue;
        
        /* Simple drift (should include acceleration from gas/stars) */
        sink->position = vec3_add(sink->position, vec3_scale(sink->velocity, sim.dt));
        
        /* Periodic boundary conditions */
        if (sink->position.x < 0) sink->position.x += sim.box_size;
        if (sink->position.x >= sim.box_size) sink->position.x -= sim.box_size;
        if (sink->position.y < 0) sink->position.y += sim.box_size;
        if (sink->position.y >= sim.box_size) sink->position.y -= sim.box_size;
        if (sink->position.z < 0) sink->position.z += sim.box_size;
        if (sink->position.z >= sim.box_size) sink->position.z -= sim.box_size;
    }
}

/*******************************************************************************
 * VORONOI CELL UPDATES (FVM)
 ******************************************************************************/

/*
 * Update derived quantities in Voronoi cells from conserved variables.
 * This is called after each hydro step.
 */
void update_cell_primitives(void) {
    double gamma = 1.4;
    
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        
        if (!cell->is_leaf) continue;
        
        double d = cell->density;
        if (d <= 0) continue;
        
        /* Velocity from momentum */
        cell->velocity = vec3_scale(cell->momentum, 1.0 / d);
        
        /* Kinetic energy */
        double ekin = 0.5 * d * vec3_dot(cell->velocity, cell->velocity);
        
        /* Internal energy and pressure */
        double eint = cell->energy - ekin;
        if (eint < 0) eint = 1e-10 * cell->energy;  /* Floor */
        
        cell->pressure = (gamma - 1.0) * eint;
        
        /* Temperature (assuming ideal gas with mu=0.6) */
        double mu = 0.6;
        cell->temperature = cell->pressure / d * mu * PROTON_MASS / 1.38e-16 * 
                           sim.scale_v * sim.scale_v;
        
        /* Sound speed */
        cell->sound_speed = sqrt(gamma * cell->pressure / d);
    }
}

/*
 * Compute fluxes between Voronoi cells (FVM).
 * Uses simple HLL Riemann solver.
 */
void compute_fluxes(double dt) {
    /* This is a simplified flux computation.
     * A full implementation would use:
     * - Exact or approximate Riemann solvers (HLL, HLLC, Roe)
     * - Gradient reconstruction for 2nd order accuracy
     * - Slope limiters for TVD property
     */
    
    for (int i = 0; i < sim.n_cells; i++) {
        VoronoiCell *cell = &sim.cells[i];
        if (!cell->is_leaf) continue;
        
        for (int j = 0; j < cell->n_neighbors; j++) {
            int neighbor_idx = cell->neighbors[j];
            if (neighbor_idx < 0 || neighbor_idx >= sim.n_cells) continue;
            
            VoronoiCell *neighbor = &sim.cells[neighbor_idx];
            if (!neighbor->is_leaf) continue;
            
            /* Face properties */
            double A = cell->face_areas[j];
            Vec3 n = cell->face_normals[j];
            
            /* Left and right states */
            double dL = cell->density;
            double dR = neighbor->density;
            Vec3 vL = cell->velocity;
            Vec3 vR = neighbor->velocity;
            double pL = cell->pressure;
            double pR = neighbor->pressure;
            double csL = cell->sound_speed;
            double csR = neighbor->sound_speed;
            
            /* Normal velocities */
            double vnL = vec3_dot(vL, n);
            double vnR = vec3_dot(vR, n);
            
            /* Wave speeds (HLL) */
            double SL = dmin(vnL - csL, vnR - csR);
            double SR = dmax(vnL + csL, vnR + csR);
            
            /* HLL flux */
            double F_mass, F_energy;
            Vec3 F_momentum;
            
            if (SL >= 0) {
                /* Left state */
                F_mass = dL * vnL;
                F_momentum = vec3_add(vec3_scale(vL, dL * vnL), vec3_scale(n, pL));
                F_energy = (cell->energy + pL) * vnL;
            } else if (SR <= 0) {
                /* Right state */
                F_mass = dR * vnR;
                F_momentum = vec3_add(vec3_scale(vR, dR * vnR), vec3_scale(n, pR));
                F_energy = (neighbor->energy + pR) * vnR;
            } else {
                /* HLL average */
                double FL_mass = dL * vnL;
                double FR_mass = dR * vnR;
                Vec3 FL_mom = vec3_add(vec3_scale(vL, dL * vnL), vec3_scale(n, pL));
                Vec3 FR_mom = vec3_add(vec3_scale(vR, dR * vnR), vec3_scale(n, pR));
                double FL_en = (cell->energy + pL) * vnL;
                double FR_en = (neighbor->energy + pR) * vnR;
                
                double denom = SR - SL;
                F_mass = (SR * FL_mass - SL * FR_mass + SL * SR * (dR - dL)) / denom;
                F_momentum = vec3_scale(
                    vec3_add(vec3_sub(vec3_scale(FL_mom, SR), vec3_scale(FR_mom, SL)),
                            vec3_scale(vec3_sub(neighbor->momentum, cell->momentum), SL * SR)),
                    1.0 / denom);
                F_energy = (SR * FL_en - SL * FR_en + 
                           SL * SR * (neighbor->energy - cell->energy)) / denom;
            }
            
            /* Update conserved variables */
            double dV = A * dt;
            cell->density -= F_mass * dV / cell->volume;
            cell->momentum = vec3_sub(cell->momentum, vec3_scale(F_momentum, dV / cell->volume));
            cell->energy -= F_energy * dV / cell->volume;
        }
    }
}

/*******************************************************************************
 * MAIN SINK PARTICLE STEP
 ******************************************************************************/

/*
 * Main routine called once per coarse timestep.
 * Handles all sink particle physics.
 */
void sink_particle_step(void) {
    printf("\n=== Sink Particle Step at t=%.3e ===\n", sim.time);
    
    /* 1. Create new sink particles */
    create_sink_particles();
    
    /* 2. Update sink positions */
    update_sink_positions();
    
    /* 3. Merge close sinks */
    merge_sink_particles();
    
    /* 4. Bondi-Hoyle accretion */
    bondi_hoyle_accretion();
    
    /* 5. AGN feedback */
    apply_all_agn_feedback();
    
    /* 6. Update cell primitives */
    update_cell_primitives();
    
    printf("=== End Sink Particle Step ===\n\n");
}

/*******************************************************************************
 * INITIALIZATION AND CLEANUP
 ******************************************************************************/

/*
 * Initialize simulation state and allocate memory.
 */
void init_simulation(int n_cells, double box_size, double dx_min) {
    /* Allocate cells */
    sim.cells = (VoronoiCell *)calloc(n_cells, sizeof(VoronoiCell));
    if (!sim.cells) {
        fprintf(stderr, "Failed to allocate Voronoi cells\n");
        exit(1);
    }
    
    sim.n_cells = n_cells;
    sim.n_sinks = 0;
    
    sim.box_size = box_size;
    sim.dx_min = dx_min;
    
    sim.time = 0;
    sim.dt = 0;
    sim.dt_coarse = 0;
    
    /* Set unit conversions (example: 100 Mpc box) */
    sim.scale_l = box_size * 3.08e24;  /* cm */
    sim.scale_d = 1.88e-29;            /* g/cm^3 (cosmological) */
    sim.scale_t = 1.0 / (sqrt(GRAVITY_CGS * sim.scale_d));
    sim.scale_v = sim.scale_l / sim.scale_t;
    sim.scale_m = sim.scale_d * pow(sim.scale_l, 3);
    sim.scale_T2 = sim.scale_v * sim.scale_v * PROTON_MASS / 1.38e-16;
    sim.scale_nH = sim.scale_d / PROTON_MASS;
    
    /* Physics flags */
    sim.use_spin = true;
    sim.use_mad_jet = true;
    sim.use_self_gravity = true;
    
    /* Cosmology */
    sim.aexp = 1.0;
    sim.hubble = 0.7;
    
    printf("Simulation initialized:\n");
    printf("  Box size: %.2e cm = %.2f Mpc\n", sim.scale_l, box_size);
    printf("  n_cells: %d\n", n_cells);
    printf("  dx_min: %.2e\n", dx_min);
}

/*
 * Free all allocated memory.
 */
void cleanup_simulation(void) {
    if (sim.cells) {
        free(sim.cells);
        sim.cells = NULL;
    }
    sim.n_cells = 0;
    sim.n_sinks = 0;
}

/*******************************************************************************
 * MAIN FUNCTION (Example usage)
 ******************************************************************************/

int main(int argc, char **argv) {
    printf("RAMSES AGN/SMBH Module - C/Voronoi/FVM Version\n");
    printf("===============================================\n\n");
    
    /* Initialize with example parameters */
    int n_cells = 1000;
    double box_size = 100.0;  /* Mpc */
    double dx_min = box_size / 512.0;
    
    init_simulation(n_cells, box_size, dx_min);
    
    /* Example: Initialize some cells with density */
    for (int i = 0; i < n_cells; i++) {
        sim.cells[i].id = i;
        sim.cells[i].is_leaf = true;
        sim.cells[i].has_sink = false;
        
        /* Random position */
        sim.cells[i].center.x = (double)rand() / RAND_MAX * box_size;
        sim.cells[i].center.y = (double)rand() / RAND_MAX * box_size;
        sim.cells[i].center.z = (double)rand() / RAND_MAX * box_size;
        
        /* Random volume (should be from Voronoi tessellation) */
        sim.cells[i].volume = pow(dx_min, 3) * (0.5 + (double)rand() / RAND_MAX);
        
        /* Initial hydro state */
        sim.cells[i].density = 1e-25 / sim.scale_d;  /* ~10 H/cc */
        sim.cells[i].temperature = 1e4;
        sim.cells[i].velocity = (Vec3){0, 0, 0};
        sim.cells[i].momentum = (Vec3){0, 0, 0};
        sim.cells[i].sound_speed = 10e5 / sim.scale_v;  /* 10 km/s */
        sim.cells[i].pressure = sim.cells[i].density * 
                                sim.cells[i].sound_speed * sim.cells[i].sound_speed / 1.4;
        sim.cells[i].energy = sim.cells[i].pressure / 0.4 + 
                             0.5 * sim.cells[i].density * 
                             vec3_dot(sim.cells[i].velocity, sim.cells[i].velocity);
        
        /* Stellar properties (for seeding) */
        sim.cells[i].stellar_density = 0.05 + 0.1 * (double)rand() / RAND_MAX;
        sim.cells[i].velocity_dispersion = 15.0 + 10.0 * (double)rand() / RAND_MAX;
    }
    
    /* Add one high-density cell to trigger sink formation */
    sim.cells[0].density = 1e-20 / sim.scale_d;  /* Very high density */
    sim.cells[0].stellar_density = 0.2;
    sim.cells[0].velocity_dispersion = 25.0;
    
    /* Run a few timesteps */
    sim.dt = 1e6 * YEAR_TO_SEC / sim.scale_t;  /* 1 Myr */
    
    for (int step = 0; step < 5; step++) {
        sim.time += sim.dt;
        sink_particle_step();
    }
    
    /* Cleanup */
    cleanup_simulation();
    
    printf("Simulation complete.\n");
    return 0;
}
