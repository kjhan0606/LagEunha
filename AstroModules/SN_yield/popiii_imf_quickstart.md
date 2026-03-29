# Pop III IMF Upgrade - Quick Start Guide

## #  What's New (Version 1.1)

Added **Pop III Top-Heavy IMF** for primordial star formation!

### Key Features

#  **New IMF Type**: `IMF_POPIII_TOPHEAVY`  
#  **Adjustable Mass Fraction**: `popiii_mass_fraction` (0-100%)  
#  **71× More PISN** compared to normal Kroupa IMF  
#  **Thread-Safe**: Full OpenMP support maintained  

---

## #  Quick Usage

### Minimal Example

```c
#include "cluster_supernova.h"

int main() {
    ClusterProperties cluster = {
        .total_mass = 1e6,
        .age = 3e6,
        .metallicity = 0.0,
        .imf_type = IMF_POPIII_TOPHEAVY,  // # New!
        .sf_mode = SF_INSTANTANEOUS,
        .star_formation_rate = 0.0,
        .popiii_mode = POPIII_ENABLED,
        .popiii_imf_max = 260.0,
        .popiii_mass_fraction = 0.10      // # Adjustable!
    };
    
    ClusterOutput output;
    calculate_cluster_supernovae(&cluster, &output);
    print_cluster_output(&output);
    
    return 0;
}
```

**Compile**:
```bash
gcc -o myprogram myprogram.c -L. -lcluster_supernova -lm
LD_LIBRARY_PATH=. ./myprogram
```

---

## #  Results Comparison

**Same cluster, different IMFs** (10^6 M#, Z=0, 3 Myr):

| IMF Type | PISN Rate | Energy | Fe Yield |
|----------|-----------|--------|----------|
| Kroupa (normal) | 5.26e-06 /yr | 2.30e+47 erg/yr | 1.54e-04 M#/yr |
| **Pop III** | **3.73e-04 /yr** | **1.56e+49 erg/yr** | **1.10e-02 M#/yr** |
| **Increase** | **71×** | **68×** | **71×** |

---

## #  When to Use

### Use Pop III IMF (`IMF_POPIII_TOPHEAVY`)

- #  Primordial clusters (Z = 0)
- #  Very metal-poor (Z < 10^-4)
- #  First stars simulations (z > 15)
- #  Metal-free minihalos

### Use Normal IMF (`IMF_KROUPA`)

- #  Milky Way-like galaxies
- #  Solar metallicity (Z ~ 0.02)
- #  Enriched gas (Z > 10^-3)
- #  Late-time star formation

---

## #  Build & Test

### Build Everything
```bash
make clean
make all
```

### Run Tests
```bash
# Basic Pop III test
LD_LIBRARY_PATH=. ./test_popiii_imf

# OpenMP test
OMP_NUM_THREADS=4 LD_LIBRARY_PATH=. ./test_popiii
```

### Expected Output
```
PISN rate: 3.734e-04 /yr  # 
Total energy: 1.563e+49 erg/yr  # 
Fe yield: 1.096e-02 M_sun/yr  # 
```

---

## #  Available IMF Types

```c
typedef enum {
    IMF_SALPETER,         // M^-2.35
    IMF_KROUPA,           // Broken power law
    IMF_CHABRIER,         // Log-normal + power law
    IMF_POPIII_TOPHEAVY   // M^-1.0 (NEW!)
} IMFType;
```

---

## ## Parameters

### Pop III Mass Fraction

```c
.popiii_mass_fraction = 0.10  // 10% (default)
```

**Recommendations**:
- **0.10** (10%): Conservative, mixed population
- **0.20-0.30** (20-30%): Pop III dominated
- **0.50** (50%): Pure Pop III (extreme)

---

## #  Files Included

### Source Code
- `cluster_supernova.h` - Updated header
- `cluster_supernova.c` - Updated implementation
- `libcluster_supernova.a` - Static library
- `libcluster_supernova.so` - Shared library

### Test Programs
- `test_popiii_imf.c` - Pop III IMF comparison test
- `test_popiii_imf` - Compiled binary
- `test_popiii.c` - OpenMP thread-safety test

### Documentation
- `POPIII_IMF_UPGRADE.md` - Complete documentation
- `POPIII_IMF_QUICKSTART.md` - This file

---

## #  Physical Basis

**Pop III Top-Heavy IMF**:
- Peak at ~100 M# (vs ~0.3 M# for normal IMF)
- No low-mass stars (M < 10 M#)
- Slope M^-1.0 (flatter than M^-2.3)

**Why?**
- Higher Jeans mass in primordial gas
- No metal/dust cooling
- Higher temperatures # more massive fragments

**References**:
- Hirano et al. (2014), ApJ, 781, 60
- Susa et al. (2014), ApJ, 792, 32

---

## #  Version Info

**Version**: 1.1  
**Date**: January 10, 2026  
**Status**: Production Ready #   
**Thread-Safe**: Yes #   
**OpenMP**: Compatible #   

---

## #  Tips

1. **Auto-select IMF based on metallicity**:
   ```c
   IMFType select_imf(double Z) {
       return (Z < 0.0001) ? IMF_POPIII_TOPHEAVY : IMF_KROUPA;
   }
   ```

2. **Scan mass fractions**:
   ```c
   for (double f = 0.01; f <= 0.50; f += 0.01) {
       cluster.popiii_mass_fraction = f;
       calculate_cluster_supernovae(&cluster, &output);
       // analyze results...
   }
   ```

3. **Time evolution study**:
   ```c
   for (double age = 1e6; age <= 10e6; age += 1e6) {
       cluster.age = age;
       calculate_cluster_supernovae(&cluster, &output);
       // track PISN evolution...
   }
   ```

---

For detailed documentation, see **POPIII_IMF_UPGRADE.md**
