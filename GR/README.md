# Z4c DG-FVM: Numerical Relativity Code

## Discontinuous Galerkin Finite Volume Method for Z4c Formulation of Einstein's Equations

### Overview

This code implements the Z4c formulation of Einstein's equations using a Discontinuous Galerkin Finite Volume Method (DG-FVM) on Voronoi meshes (currently using Cartesian approximation).

### Key References

Each function and equation in the code is documented with references to the following papers:

1. **[1] Baumgarte & Shapiro**, "Numerical Relativity: Solving Einstein's Equations on the Computer", Cambridge University Press (2010)
   - ADM 3+1 formalism: Chapters 2-3
   - Gauge conditions: Chapter 4
   - Constraint equations: Eqs. (2.125)-(2.126)

2. **[2] Dumbser et al.**, J. Comput. Phys. 348, 70-117 (2017)
   - "ADER discontinuous Galerkin schemes for general-relativistic ideal magnetohydrodynamics"
   - First-order reduction: Eqs. (15)-(18)
   - Numerical flux: Eq. (35)
   - DG weak formulation: Section 3

3. **[3] Hilditch et al.**, Phys. Rev. D 88, 084057 (2013)
   - "An introduction to well-posedness and free-evolution"
   - Z4c evolution equations: Eqs. (13)-(23)
   - Conformal metric constraint: Eq. (3.5)
   - Traceless constraint: Eq. (3.6)

4. **[4] Bernuzzi & Hilditch**, Phys. Rev. D 81, 084003 (2010)
   - "Constraint violation in free evolution schemes"
   - Constraint damping parameters κ₁, κ₂: Eq. (2.7)

5. **[12] Shu & Osher**, J. Comput. Phys. 77, 439 (1988)
   - SSP-RK3 time integration: Eq. (2.18)

6. **[14] Alcubierre et al.**, Class. Quantum Grav. 21, 589 (2004)
   - Gauge wave test problem: Eqs. (4.1)-(4.3)

7. **[16] Ruiz et al.**, Gen. Relativ. Gravit. 40, 1705 (2008)
   - Boundary conditions for numerical relativity

### Z4c Variables

The code evolves 22 primary variables plus 33 auxiliary variables:

#### Primary Variables (22)
| Variable | Description | Index |
|----------|-------------|-------|
| φ | Conformal factor: γᵢⱼ = e^{4φ} γ̃ᵢⱼ | 0 |
| γ̃ᵢⱼ | Conformal metric (6 components) | 1-6 |
| K | Trace of extrinsic curvature | 7 |
| Ãᵢⱼ | Traceless conformal extrinsic curvature (6) | 8-13 |
| Γ̃ⁱ | Conformal connection functions (3) | 14-16 |
| Θ | Z4 constraint variable | 17 |
| α | Lapse function | 18 |
| βⁱ | Shift vector (3) | 19-21 |

#### Auxiliary Variables (33) - First-Order Reduction
| Variable | Definition | Reference |
|----------|------------|-----------|
| Aᵢ | ∂ᵢ ln α | [2] Eq. (15) |
| Bⁱⱼ | ∂ⱼ βⁱ | [2] Eq. (16) |
| Dₖᵢⱼ | ½ ∂ₖ γ̃ᵢⱼ | [2] Eq. (17) |
| Pᵢ | ∂ᵢ φ | [2] Eq. (18) |

### Numerical Methods

#### Spatial Discretization
- **Finite Volume Method** on Cartesian/Voronoi meshes
- **Numerical Flux**: Local Lax-Friedrichs (Rusanov) - [2] Eq. (35)

#### Time Integration
- **SSP-RK3** (Strong Stability Preserving Runge-Kutta) - [12] Eq. (2.18)
```
U^(1) = U^n + Δt L(U^n)
U^(2) = ¾U^n + ¼U^(1) + ¼Δt L(U^(1))
U^(n+1) = ⅓U^n + ⅔U^(2) + ⅔Δt L(U^(2))
```

#### Gauge Conditions
- **Lapse**: 1+log slicing: ∂ₜα = -2αK - [1] Eq. (4.60)
- **Shift**: Gamma-driver: ∂ₜβⁱ = ¾α²Γ̃ⁱ - ηβⁱ - [1] Eq. (4.82)

#### Constraint Damping
Parameters from [4] Bernuzzi & Hilditch:
- κ₁ = 0.02 (constraint damping strength)
- κ₂ = 0.0 (typical for BBH simulations)
- κ₃ = 0.5 (Gamma constraint damping)

### Building and Running

```bash
# Compile
make

# Run tests
./bin/z4c_dg -t flat -n 10 -T 1.0
./bin/z4c_dg -t gauge_wave -n 10 -T 1.0 -A 0.01

# Options
-t TEST    Test problem: flat, gauge_wave
-n N       Grid cells per dimension (creates N×N×N mesh)
-T TIME    Final simulation time
-A AMP     Wave amplitude (for gauge_wave)
-o DIR     Output directory
-h         Show help
```

### Output

The code outputs VTK files for visualization with ParaView or VisIt:
- `z4c_XXXXX.vtk` - Contains α (lapse), φ (conformal factor), K (trace), Θ (constraint)

### Test Problems

#### 1. Flat Spacetime
Minkowski metric in Cartesian coordinates. Tests code stability.
- Expected: All constraints remain exactly zero

#### 2. Gauge Wave [14]
Exact solution representing a pure gauge (coordinate) wave:
```
ds² = -H(x-t)dt² + H(x-t)dx² + dy² + dz²
H(u) = 1 - A sin(2πu)
```
- Tests wave propagation and long-term stability
- Expected: Constraint violations ~O(10⁻⁷) for A=0.01

### Code Structure

```
z4c_dg_fvm/
├── include/
│   └── z4c_dg.h          # Header with all declarations
├── src/
│   └── z4c_combined.c    # Combined implementation
├── Makefile
└── README.md
```

### Key Functions

| Function | Reference | Description |
|----------|-----------|-------------|
| `z4c_compute_source()` | [3] Eqs. (13)-(23) | Z4c source terms |
| `z4c_compute_flux()` | [2] Section 2.2 | Advection flux |
| `z4c_characteristic_speeds()` | [2] Section 2.3 | Wave speeds |
| `riemann_llf()` | [2] Eq. (35) | Local Lax-Friedrichs flux |
| `time_step_rk3_ssp()` | [12] Eq. (2.18) | SSP-RK3 integration |
| `compute_hamiltonian()` | [1] Eq. (2.125) | Hamiltonian constraint |
| `enforce_constraints()` | [3] Eqs. (3.5)-(3.6) | Algebraic constraints |

### Constraints Monitored

1. **Hamiltonian Constraint** - [1] Eq. (2.125):
   ```
   H = R + K² - KᵢⱼKⁱʲ = 0
   ```

2. **Z4 Constraint (Theta)** - [3] Eq. (3.8):
   ```
   Θ ≈ H/2 (should decay due to constraint damping)
   ```

3. **Algebraic Constraints**:
   - det(γ̃ᵢⱼ) = 1 - [3] Eq. (3.5)
   - tr(Ãᵢⱼ) = 0 - [3] Eq. (3.6)

### License

MIT License

### Author

Generated by Claude (Anthropic)
