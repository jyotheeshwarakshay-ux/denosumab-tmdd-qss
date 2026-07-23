# Diagnostic trail: v1 → v2

Full record of the debugging process referenced in the main [README](README.md). Retained because the reasoning is the substance of the project, not an appendix to it.

---

## The symptom

The first fit (`denosumab_QSS_TMDD_nlmixr2.R`) converged cleanly. SAEM ran 700 iterations, parameters stabilised to three decimal places, no warnings. The diagnostics said otherwise:

- `mean(DV − IPRED) = +5.014` — systematic underprediction
- \>75% of residuals positive
- max IPRED 70.6 against an observed max of 120.8 — the model could not reach the peaks
- `add.err` estimated at 2.68 against a true 0.72
- `Q` at 499% RSE with 902% CV between-subject variability

A converged fit that underpredicts three-quarters of its observations is not a fit. The inflated `add.err` was the tell: SAEM was widening the noise term to absorb a misfit it could not otherwise accommodate.

---

## Hypothesis 1 — negative-value clipping biased the data upward

**Reasoning.** The simulation applied `DV = pmax(DV, 0)` to prevent negative concentrations. At troughs where the true concentration approaches zero, clipping deletes the lower half of the noise distribution, which would bias observations upward and make the model appear to underpredict.

**Test.** `test_error_model.py` — simulated the clipping directly across a realistic concentration distribution and measured the resulting mean residual.

**Result: rejected.** Clipping produced a mean residual of −0.02 overall and +0.13 at troughs. Explaining a +5.01 bias requires an effect two orders of magnitude larger. Real, but nowhere near sufficient.

---

## Hypothesis 2 — sampling too sparse to identify Vc, Q and Vp separately

**Reasoning.** `Vc` recovered at 3.43 against a true 1.58, `Q` at 1.14 against 0.20, `Vp` at 7.32 against 6.06 — all inflated together. Volume and inter-compartmental clearance parameters are classically correlated. If the sampling schedule cannot separate them, the estimator may settle in a wrong-but-self-consistent region.

**Test.** `test_identifiability.py` — simulated the true parameter set and the fitted parameter set, then compared the two concentration-time curves **at the actual sampled timepoints**. If the data cannot distinguish them, the curves should overlap there.

**Result: rejected.** The curves differ by 22–40% at every sampled timepoint, RMSE 6.15 nmol/L. The data contains ample information to reject the fitted values. A denser schedule would not have helped, and testing the proposed denser design gave RMSE 7.09 — no improvement in discriminating power.

---

## Hypothesis 3 — compartment (CMT) mis-mapping

**Reasoning.** The simulation dosed via `ev$add.dosing(dosing.to = "depot")`. The fit read `CMT = 1` for dosing and `CMT = 2` for observations from a data frame. A four-compartment model (depot, Ctot, Cp, Rtot) makes the numbering non-obvious; if observations were mapped to the wrong state, everything downstream would be wrong.

**Test.** Compared the fit's reported `IPRED` column against its `Ctot` column across all summary statistics.

**Result: rejected.** IPRED and Ctot were identical (median 8.181 vs 8.181, mean 11.719 vs 11.719, max 70.563 vs 70.563). The residual statement `Ctot ~ add + prop` bound correctly and the mapping worked as intended.

---

## Hypothesis 4 — between-subject variances collapsed

**Reasoning.** Comparing the estimated BSV against the true IIV revealed a pattern: `eta_Vc` estimated at 0.0156 against a true 0.3225 (20-fold collapse), `eta_R0` at 0.0254 against ~1.02 (40-fold). With random effects near zero, SAEM cannot let individuals differ, so it must find one compromise typical value across a population that genuinely varies threefold — and that compromise sits high, producing underprediction of peaks.

The likely cause: the fit model initialised those variances *below* their true values, and SAEM drove them toward zero rather than climbing.

**Test.** `diagnostic_test_v2.R` — re-fit on 20 subjects with short SAEM settings, every eta initialised **at** its true ω².

**Result: partly confirmed.** The collapse was real and the fix held — `eta_Vc` ratio 1.15× instead of 0.05×. But `mean(DV − IPRED)` came back at **+5.67**, essentially unchanged from +5.01.

So variance collapse was a genuine defect, worth correcting, and not the cause of the bias.

---

## Hypothesis 5 — simulation and fit models differ structurally

**Reasoning.** With the estimator exonerated, the remaining possibility was that the model generating the data and the model being fitted were not the same model. The fit model contained one line absent from the simulation: `Rtot(0) <- R0`, setting the initial target concentration inside the model block rather than externally via `rxSolve(inits = ...)`.

**Test.** `structural_test.R` — built both model objects, solved each at identical true parameters with identical event tables, and compared the resulting curves point by point.

**Result: rejected.** Maximum difference 1e-12 nmol/L — pure numerical noise. The two initial-condition mechanisms are equivalent and the models are identical.

---

## Hypothesis 6 — data handling in the NONMEM-style dataset

**Reasoning.** The simulation used an `eventTable` object; the fit read a data frame with EVID/AMT/CMT columns. Even with identical models, a dosing record that fails to deliver drug to the depot would produce systematically low predictions.

**Test.** `data_handling_test_v2.R` — solved the same model twice for subject 1, once feeding `pk_data` and once feeding an `eventTable`, and compared.

*(An earlier version of this test used `nlmixr2(est = "posthoc")` with every parameter fixed and failed with "no parameters to estimate" — the estimator requires at least one free parameter. Rewritten to use `rxSolve` directly.)*

**Result: rejected.** max Ctot 37.23 via `pk_data` against 37.24 via `eventTable`. The dosing records deliver correctly.

---

## Hypothesis 7 — washout sampling starved the likelihood

**Reasoning.** With model, data handling, and estimator all cleared, attention returned to the sampling schedule — not its density, but *where* it sampled. The v1 offsets included 3600 h and 4380 h post-dose. On the typical curve those sit at 1.41 and 0.095 nmol/L. Against an additive residual error of 0.72 nmol/L, noise at the 4380 h trough is roughly **7× the signal**.

SAEM weights observations by inverse variance. Hundreds of near-zero observations, carrying no information but weighted as though they did, would pull the fitted curve downward.

The source paper never faced this: their assay LLOQ was 20 ng/mL, and 26.49% of post-dose samples fell below it and were set to MDV=1, excluded from the likelihood entirely. This project simulated the assay but not its detection limit, generating observations the real dataset never contained.

**Test.** Redesigned the schedule — dropped 3600 h and 4380 h, added 72 h, 336 h and 1440 h in the informative absorption/distribution region. Minimum typical concentration rises from 0.095 to 8.55 nmol/L; worst-case noise-to-signal falls from 7.0 to 0.08. Re-fit (`denosumab_QSS_TMDD_v2.R`).

**Result: partly confirmed.**

Every structural parameter improved:

| | v1 error | v2 error |
|---|---|---|
| ka | +54% | +12.3% |
| Vc | +117% | +37.3% |
| Vp | +21% | +2.2% |
| CL | +32% | +17.7% |
| R0 | −26% | +2.7% |

Precision roughly halved across the board (`kint` %RSE 24.6 → 9.8; `R0` 20.2 → 11.4).

But `mean(DV − IPRED)` moved only from +5.01 to **+3.33**, with 75.4% of residuals still positive — and `add.err` got *worse*, 2.68 → 7.88.

---

## Where it stands

Seven hypotheses, five rejected outright, two partly confirmed. The two real defects — variance collapse and washout sampling — were both corrected, and both improved the fit substantially. Neither explains the residual bias.

**The cause of the systematic underprediction is unidentified.** The next things worth testing, in the order I would try them:

1. Whether the QSS closed-form solution behaves correctly when `Ctot` approaches `Rtot + Kss` (the discriminant nears zero and the square root becomes numerically delicate)
2. Whether `max(disc, 0)` and `max(C, 0)` in the model block behave identically in nlmixr2's SAEM path and in rxode2's simulation path — the structural test compared rxode2 to rxode2, not rxode2 to nlmixr2's internal compilation
3. An SSE across replicates, to establish whether the bias is systematic or an artifact of this particular simulated dataset

---

## Note on method

Every hypothesis above was tested by direct measurement rather than argued from plausibility, and each test script is in this repository. Several plausible-sounding hypotheses were wrong — clipping, sparse sampling, CMT mapping — and would have cost days of misdirected work if acted on without testing.

The reason for documenting the rejected hypotheses alongside the confirmed ones is that the rejections carry most of the information. Knowing the models are identical to 1e-12 and the dosing records deliver correctly narrows the search far more than any single positive result did.
