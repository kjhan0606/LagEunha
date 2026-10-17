# Code Review: sink_particle_kjhan.c (AGN/SMBH Physics Module)

## 1. Overview
This code is a simulation module that ports the RAMSES AGN/SMBH module (based on Dubois et al.) into a **Voronoi Tessellation-based Finite Volume Method (FVM)** structure.

Beyond a simple port, it is evaluated as a high-level prototype (Toy Model) where cutting-edge theories such as **Spin Evolution (Fanidakis et al.)**, **MAD Jet Efficiency (McKinney et al.)**, and **BH Merger Spin (Rezzolla et al.)** are sophisticatedly implemented.

---

## 2. Physics Consistency & Strengths

It faithfully follows standard recipes for modern galaxy formation simulations, with spin-related physics being implemented in exceptional detail.

### #  Sophisticated Spin Evolution Model
- **Based on Fanidakis et al. (2011):** Instead of merely increasing mass, it determines the alignment between the angular momentum of the accreted gas and the BH spin.
- **Regime Distinction:** The logic distinguishing between the Self-Gravity Regime and the Warped Disc Regime to calculate disc mass is detailed.
- **MAD Jet Implementation:** By applying the **McKinney et al. (2012)** model where jet efficiency changes according to the spin parameter ($a$), robust feedback is implemented for high-spin BHs.

### #  BH Merger Physics (Rezzolla et al. 2008)
- When two black holes merge, it calculates the final spin ($a_f$) by considering the **mass ratio ($q$)**, **spin vectors**, and **orbital angular momentum**, rather than just summing masses.
- This is essential for predicting remnant BH properties after Gravitational Wave (GW) emission.

### #  Dual Mode Feedback
- Follows the **Dubois et al. (2012)** standard recipe.
- It naturally switches between **Radio Mode (Kinetic/Jet)** and **Quasar Mode (Thermal)** based on the Eddington ratio ($\chi < 0.01$).

---

## 3. Critical Physics Issues

These are areas that must be addressed to ensure the physical accuracy of the simulation.

### ## 1. Violation of Angular Momentum Conservation
- **Observation:** `accumulate_gas_properties` calculates gas angular momentum to determine the change in BH Spin.
- **Issue:** The BH gains angular momentum, but the **gas cells do not lose it** (absence of back-reaction).
- **Consequence:** The rotational velocity of the gas does not decrease, potentially delaying accretion or failing to conserve the total angular momentum of the system.

### ## 2. Absence of Dynamical Friction
- **Observation:** There is no dynamical friction term in the sink particle motion equations other than gravity and fluid drag.
- **Issue:** Massive BHs may wander around instead of sinking to the galactic center (potential deepest point).
- **Correction:** A drag force based on the Chandrasekhar Dynamical Friction formula needs to be added.

### ## 3. "Action at a Distance" in FVM
- **Observation:** The `accrete_mass_from_cells` function immediately removes mass from surrounding cells.
- **Issue:** Unlike SPH, arbitrarily removing mass inside FVM cells disrupts pressure balance, potentially creating **unphysical cavities and shocks**. Approaches involving flux flow across faces should be considered.

### ## 4. Absence of Radiative Cooling
- **Observation:** Only the Ideal Gas EOS ($\gamma=1.4$) is used, and there is no cooling function.
- **Issue:** Gas heated by AGN feedback does not cool down. This leads to **over-suppression**, where hot bubbles remain permanently, blocking further gas accretion forever.

---

## 4. Implementation Issues

### #  Neighbor Search Efficiency
- **Current Status:** $O(N_{sinks} \times N_{cells})$
- Looping through all cells for every sink severely degrades performance as the cell count increases.
- **Improvement:** Use BFS search utilizing Voronoi Adjacency Lists or implement Octree/k-d trees.

### #  Simple Flux Calculation
- **Current Status:** Simple 1st-order HLL Solver.
- **Issue:** Excessive numerical viscosity may smear out rotating accretion disc structures around the BH too quickly (angular momentum transport occurs via numerical error rather than physical viscosity).

---

## 5. Upgrade Roadmap

### Phase 1: Basics (Essential Physics)
1.  **Add Radiative Cooling:** Implement `Cooling_Function(T, Z)` (e.g., Sutherland & Dopita).
2.  **Add Dynamical Friction:** Apply drag force where $F_{df} \propto - \frac{\mathbf{v}}{v^3}$.

### Phase 2: Advanced Physics (Accretion & Feedback)
1.  **Angular Momentum Transfer:** Decelerate the velocity vector of gas cells by the amount of angular momentum transferred to the BH (Apply Torque).
2.  **Viscous Time Delay:** Instead of accreting immediately at the Bondi rate, delay accretion by the viscous timescale ($\dot{M}_{acc} \neq \dot{M}_{Bondi}$).
3.  **Flux-based Accretion:** Consider generating mass flux toward the cell containing the sink instead of deleting mass directly.

### Phase 3: Optimization
1.  **Optimize Neighbor Search:** Implement Voronoi Graph search algorithms.
2.  **Sub-cycling:** Apply a shorter time step ($dt_{BH}$) only for the region around the BH ($R_{acc}$) to improve integration precision.

---

> **Overall Assessment:**
> This code attempts to faithfully reflect the latest AGN physics theories, going beyond a simple toy model. The **Spin-driven Jet** and **Merger** logic are particularly impressive. If the **Cooling** and **Dynamical Friction** issues mentioned above are addressed, it would serve as an excellent prototype for research-grade code.
