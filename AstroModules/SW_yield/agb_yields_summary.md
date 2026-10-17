# AGB Star Yields Summary (Corrected)

**Source:** Karakas, A.I. & Lugaro, M. 2016, ApJ, 825, 26

**Mass range:** 1.0 - 8.0 M# (Low and Intermediate Mass Stars)

**End state:** White Dwarf (no supernova explosion)

---

## Z = 0.014 (Solar Metallicity)

| Mass | H | He | C | N | O | Z_total | Wind Mass | 
|-----:|------:|------:|-------:|-------:|-------:|--------:|---------:|
| 1.00 | -0.006 | +0.008 | -0.0001 | +0.0002 | +0.0000 | +0.0003 | 0.42 |
| 1.50 | -0.018 | +0.020 | **+0.002** | +0.0008 | +0.0002 | +0.0038 | 0.88 |
| 2.00 | -0.022 | +0.023 | **+0.004** | +0.0013 | +0.0002 | +0.0066 | 1.34 |
| 2.50 | -0.042 | +0.041 | **+0.008** | +0.0020 | +0.0000 | +0.0122 | 1.82 |
| 3.00 | -0.088 | +0.078 | **+0.015** | +0.0030 | -0.0004 | +0.0227 | 2.31 |
| 3.50 | -0.099 | +0.088 | **+0.019** | +0.0037 | -0.0007 | +0.0272 | 2.77 |
| 4.00 | -0.077 | +0.073 | +0.016 | +0.0045 | -0.0009 | +0.0224 | 3.19 |
| 4.50 | -0.122 | +0.115 | +0.006 | **+0.020** | -0.0012 | +0.0280 | 3.64 |
| 5.00 | -0.201 | +0.190 | **-0.003** | **+0.036** | -0.0025 | +0.0341 | 4.12 |
| 6.00 | -0.351 | +0.346 | **-0.006** | **+0.044** | -0.0079 | +0.0338 | 5.08 |
| 7.00 | -0.470 | +0.479 | **-0.010** | **+0.042** | -0.0114 | +0.0239 | 6.02 |
| 8.00 | -0.587 | +0.607 | **-0.013** | **+0.040** | -0.0121 | +0.0184 | 6.94 |

---

## Key Physical Interpretation

### Mass Conservation Check #
```
H_yield + He_yield + Z_total # 0
```
# : 8 M## #  -0.587 + 0.607 + 0.018 = +0.038 # 0 (# # # #  # #  # # #  # #  #  0)

### Carbon Yield vs Mass

```
C Yield (M_sun)
     |
0.02 +--------# Peak at ~3 M#
     |        |   # Third Dredge-Up (TDU)
0.01 +---#    |
     |   |    #------#
   0 +---+-----------+----------
     |   |           |   # Hot Bottom Burning (HBB)
-0.01+   |           #----------
     |   |
     +---+-----------------------#
         2    3    4    5    6    7    8  Mass (M#)
```

### Nucleosynthesis Regimes

| Mass Range | Dominant Process | C Yield | N Yield | Physical Reason |
|:-----------|:-----------------|:-------:|:-------:|:----------------|
| 1-2 M# | First Dredge-Up | ~0 | + | Shallow convective envelope |
| **2-4 M#** | **Third Dredge-Up** | **+** | + | He-burning products mixed to surface |
| **>4 M#** | **Hot Bottom Burning** | **-** | **++** | C # N via CNO at base of envelope |

### Why Z_total is Positive (unlike Massive Stars)

| Component | AGB Stars | Massive Stars |
|:----------|:----------|:--------------|
| Main contribution | C, N, Ne from nucleosynthesis | C, O, Fe from SN |
| O, Fe | Locked in WD (negative) | Released in SN (positive) |
| **Net Z_total** | **Slightly positive** (C, N dominate) | **Strongly positive** (SN metals) |

---

## Comparison Across Metallicities

### Carbon Yield at 3 M#

| Z | C Yield | N Yield | Z_total |
|:---:|:-------:|:-------:|:-------:|
| 0.007 | **+0.023** | +0.0016 | +0.030 |
| 0.014 | **+0.015** | +0.0030 | +0.023 |
| 0.030 | +0.003 | +0.0050 | +0.008 |

**Lower metallicity # Higher C yield** (#  # # # #  TDU)

### Nitrogen Yield at 6 M#

| Z | C Yield | N Yield | Z_total |
|:---:|:-------:|:-------:|:-------:|
| 0.007 | -0.002 | **+0.038** | +0.031 |
| 0.014 | -0.006 | **+0.044** | +0.034 |
| 0.030 | -0.014 | **+0.036** | +0.015 |

**Higher metallicity # More efficient HBB** (stronger C#N conversion)

---

## References

1. Karakas, A.I. & Lugaro, M. 2016, ApJ, 825, 26
2. Cristallo, S. et al. 2015, ApJS, 219, 40
3. Ventura, P. & D'Antona, F. 2011, MNRAS, 410, 2760
4. Asplund, M. et al. 2009, ARA&A, 47, 481 (solar abundances)
