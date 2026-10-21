# Experimental Hatono H1 sediment inventory integration

Status: **scaffolding only, NOT ready for production or enabled model runs**.

Base: `master` at `97433db8dc7cb74d9b9ae3f99622e5cc1ab31365`.
This branch currently introduces only a default-off `h1_experimental` field to the
`&nl_colm_sediment_parameter DEF_SEDIMENT` record and a restart mode guard.
**H1 numerical flux integration has not been implemented.** Never enable
`h1_experimental=.true.` until the work below is complete and tested.

## Required numerical contract

At each morphology substep, freeze donor bed bulk volumes `layer(k,i)`,
deep bulk volumes `seddep(k,l,i)`, current particle fractions, shear/critical shear,
and hydro carriers **before any bed removal or incoming deposition**.

All quantities used below for flux/stock are *solid volume*, not bulk volume:

- `Se(k,i)=(1-lambda)*layer0(k,i)`
- `Sb(k,i)=(1-lambda)*(layer0(k,i)+sum_l seddep0(k,l,i))`
- `P(k,i)=Es_raw(k,i)*dt + B_raw_donor(k,i)`, where `B_raw_donor`
  is the sum of **all true-donor** gross bedload demands in solid m3:
  ordinary forward, ordinary reverse, all BIF paths. Do **not** include
  incoming bedload, deposition, levee protected bed, or the suspended
  advective flow from the water column.
- `T(k,i)=Se(k,i)*(-expm1(-P(k,i)/Sb(k,i)))` for positive `Sb,P`,
  otherwise zero.
- `f(k,i)=T/P` if `P>0`; zero for an empty bed.
- Apply `f` **exactly once** to gross bedload outgoing from that donor
  and to gross suspended-bed entrainment. No independent second limiter
  should change these two fractions. Keep physical settling separate.

## Source accounting / ordering

1. Current `grid_sediment_calc` already computes RMS friction velocity.
2. Its `ordinary_sediment_donor_demand` gathers reverse-flow
   demands through `push_ups2ucat`.
3. `prepare_bif_sediment` gathers BIF reverse-path demands via
   `push_bif_influx` and **directly changes `layer`**. Modify its
   donor computation and debit atomically with the joint H1 fraction,
   or it will duplicate or erase physical stock.
4. `calc_sediment_advection` currently recomputes and limits directional
   bed fluxes. Use the precomputed gross demands and joint donor fractions;
   do not independently assign a second complete initial layer budget.
5. `calc_sediment_exchange` currently takes `Es-D` first and caps
   the net at available layer inventory. H1 requires **gross** `Es`
   to compete with gross bedload, followed by independently limited
   deposition. Diagnostic `Es_eff` must equal the actual gross extracted
   amount; `D_eff` is actual settled plus documented cap deposits.
6. Ordinary and BIF incoming credits, hillslope input and concentration-cap
   deposits are applied only after the donor H1 decision. They cannot be
   re-exported during the *same* morphology interval.
7. Keep `sedbed_protected` immobile. Preserve `sedsto_protected`
   and levee accounting. Do not change their physical closure.
8. Finally call existing `calc_layer_redistribution` to refill
   the active layer. H1 should have a finite continuous-time limit
   despite this operation.

## Critical implementation caveats

- The code's saved `layer`/`seddep` are **bulk** volumes; `Es`,
  `sedsto`, `bedout` are solid volumes.
- For nonempty `Sb` but zero `P`, no extracted mass: scaling can
  be 1 safely since all demand fluxes are zero.
- For empty donor or dry cells, *both* erosion and bedload must be zero.
- The input bed depth of 10 m is not independently validated; H1 is
  a numerical experiment only, not a calibrated physical solution.
- H1 should never create sediment from deep layers within the same step.
- Record the algorithm mode in restart metadata; refuse mode mismatch.
- The default-off path must remain bit-for-bit equivalent to baseline.

## Validation plan before enabling H1

- Pure single-grain exchange: analytical exponential, `dt=3600/1800/900/450`.
- Two interacting grain classes; zero stocks; high demand; mass balance.
- Common gross bed budget for bedload + suspension; sediment settling
  excluded from demand.
- Multi-outlet forward/reverse donor tests including split MPI workers.
- BIF and levee wet/dry transitions; same-step credit re-export prevention.
- Restart OFF/OFF and ON/ON; reject OFF/ON and ON/OFF.
- Multiple `deltime` windows and `dt_morph` values; reject near-`1/dt`
  dependence in the strong-supply-limited experiments.
- Preserve existing default-off checks and tests; run both USEMPI and serial,
  standard and `CoLMDEBUG` compiles.

This document is a **design and safety checklist**, not proof of a working H1
implementation.
