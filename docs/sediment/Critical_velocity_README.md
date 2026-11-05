### Selection of the Critical Shear Velocity Exponent: 31/22 vs. 31/32

**Date:** November 5, 2026  
**Source file:** `main/MATERIAL/FAMILY_PARTICLE/MOD_MT_Sediment.F90`  
**Function:** `calc_critical_shear_vel_sq(diam)`  
**Related historical commit:** [`3bc14df4`](https://github.com/jellywzx/CoLM202X_SYSU/commit/3bc14df45f85b040f1a20d86386d79d837054c66)

**Decision: Retain the exponent `31/22` and do not restore the historical modification to `31/32`.**

#### 1. Literature Basis

Iwagaki (1956), in *Hydrodynamical Study on Critical Tractive Force*, presented a piecewise empirical relationship for the critical shear velocity (p. 15, Eq. 1.67).

For sediment particles with diameters in the range

\[
0.118 \leq d < 0.303\ \mathrm{cm}
\]

(equivalent to 1.18–3.03 mm), the relationship is expressed as

\[
u_{*c}^{2}=134.6d^{31/22}
\]

where \(d\) is the sediment particle diameter in centimeters and \(u_{*c}\) is the critical shear velocity in centimeters per second.

The exponent `31/22` is also adopted in the iRIC hydraulics textbook and the GETFLOWS Theory Manual.

This relationship is semi-empirical, combining theoretical considerations of sediment incipient motion with experimentally constrained coefficients.

#### 2. Numerical Consistency

Replacing the exponent `31/22` with `31/32` substantially increases the predicted critical shear velocity within the affected particle-size range.

For example, at \(d=2\ \mathrm{mm}\):

| Exponent | Critical Shear Velocity |
|---|---:|
| 31/22 | 0.03733 m/s |
| 31/32 | 0.05321 m/s |

The exponent `31/32` produces a critical shear velocity approximately **42.5% higher** than that obtained using `31/22`.

Furthermore, adopting `31/32` introduces substantial discontinuities at the particle-size boundaries of 1.18 mm and 3.03 mm, whereas `31/22` maintains approximate continuity between adjacent segments of the empirical relationship.

These results provide additional numerical support for retaining `31/22`.

#### 3. Implications for the Current Study and Implementation Decision

In CoLM, the sediment module first calculates the active-layer-weighted mean particle diameter (`dMean`) and uses this value to select the appropriate segment of the critical shear velocity relationship.

The Egiazaroff correction is subsequently applied to calculate the critical shear velocity for each sediment size class.

The current Amazon fine-sediment experiment uses particle diameters of 2, 20, and 100 μm. Since the active-layer-weighted mean diameter cannot exceed the largest prescribed size class (100 μm), it cannot enter the 1.18–3.03 mm interval affected by the exponent discrepancy.

Therefore, the difference between `31/22` and `31/32` does not affect the critical shear velocity calculations for this experimental configuration.

Based on the original empirical formulation, supporting literature, and numerical continuity across particle-size intervals, **the exponent `31/22` is retained in the current CoLM sediment implementation**.

The modification to `31/32` introduced in historical commit `3bc14df4` will not be restored.

Although CaMa-Flood v4 adopts `31/32` in its sediment implementation, a specific theoretical derivation or independent calibration justification for this choice has not been established. Consequently, the difference is documented as an implementation discrepancy rather than treated as an independently calibrated parameter.

For future experiments involving coarser sediment, particularly where the active-layer-weighted mean diameter falls within 1.18–3.03 mm, `31/32` may be evaluated as a separate sensitivity-test configuration rather than introduced as a correction to the baseline physical formulation.

#### 4. References

1. Iwagaki, Y. (1956). *Hydrodynamical Study on Critical Tractive Force*. Transactions of the Japan Society of Civil Engineers, 41, 1–21. See p. 15, Eq. (1.67). [DOI](https://doi.org/10.2208/jscej1949.1956.41_1)

2. iRIC. *Hydraulics for Practitioners*, Section 8.9, Eq. (8.15). [Online textbook](https://i-ric.org/yasu/nbook_e/06_Chapt06_e.html)

3. Geosphere Environmental Technology Corp. (2026). *GETFLOWS Theory Manual*, Version 2.0, Eq. (3.39). [Theory manual](https://www.getflows.getc.co.jp/docs/en/en-theory-manual.html)

4. CaMa-Flood v4. `src/sediment/cmf_calc_sedpar_mod.F90`, function `calc_criticalShearVelocity(diam)`. [GitHub source code](https://github.com/global-hydrodynamics/CaMa-Flood_v4/blob/master/src/sediment/cmf_calc_sedpar_mod.F90)
