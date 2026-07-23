# QSS-Approximation TMDD Model — Denosumab (nlmixr2)

A simulation-based parameter-recovery study of a two-compartment target-mediated drug disposition (TMDD) model under the quasi-steady-state (QSS) approximation, implemented in R / nlmixr2 with SAEM estimation.

The project asks a specific question: **given a published, peer-reviewed TMDD-QSS parameter set, can those parameters be recovered by re-fitting a model to data simulated from them — and which parameters cannot?**

Author: Jyotheeshwar Akshay Ravikumar
Institution: Sri Ramachandra University (SRIHER), Chennai, India
Field: Quantitative Systems Pharmacology / Pharmacometrics
Tools: R 4.6.1, nlmixr2 (SAEM), rxode2, xpose.nlmixr2, ggplot2, vpc

---

## Methodology — read this first

**This is not an analysis of real patient data, and the repository does not claim to be.**

Patient-level denosumab concentration data are not publicly available. What *is* public is the **fitted parameter table** from Choi et al. (2025), estimated by those authors from 6,583 real serum concentrations in 615 real subjects across a Phase I and a Phase III study.

This project:

1. Takes those published parameter estimates as ground truth
2. Simulates a dataset from them (80 subjects, 24 samples each, with the published between-subject variability and residual error)
3. Re-fits the QSS-TMDD model to that simulated data
4. Compares recovered parameters against the known truth

This is **simulation-based parameter recovery**, a standard method for evaluating estimability. It is not a claim to have analysed clinical data.

---

## Ground truth source

Choi S, Park S, Jung J, Baek S, Lim H-S (2025). *Population pharmacokinetics/pharmacodynamics analysis confirming biosimilarity of SB16 to reference denosumab.* **Frontiers in Pharmacology** 16:1631034. doi:10.3389/fphar.2025.1631034

Open access under CC-BY. Table 3 (final population PK parameter estimates, postmenopausal osteoporosis population) supplies every parameter used here.

Model theory:
- Dua P, Hawkins E, van der Graaf PH (2015). *A Tutorial on Target-Mediated Drug Disposition (TMDD) Models.* CPT Pharmacometrics Syst Pharmacol 4(6):324–337.
- Gibiansky L, Gibiansky E, Kakkar T, Ma P (2008). *Approximations of the target-mediated drug disposition model and identifiability of model parameters.* J Pharmacokinet Pharmacodyn 35(5):573–591.

**Derived parameter:** `kdeg` is not reported separately in Table 3. It is derived from the steady-state assumption stated in the paper's Methods (their Eq. 7), `ksyn = kdeg × R0`, giving `kdeg = ksyn / R0 = 6.566e-4 /h`. This derivation is documented in the script and stated here because an undocumented derived parameter is a reproducibility gap.

---

## Model

Two-compartment TMDD with QSS approximation and first-order subcutaneous absorption.

```
depot --ka--> central (Ctot) <--Q--> peripheral (Cp)
                  |
                  +-- CL (linear elimination)
                  +-- kint * RC (target-mediated elimination)

target: Rtot, synthesised at ksyn, degraded at kdeg, internalised as complex at kint
```

QSS closed form for free drug:

```
C = ½ · [ (Ctot − Rtot − Kss) + √( (Ctot − Rtot − Kss)² + 4·Kss·Ctot ) ]
RC = Ctot − C
```

Dosing: 60 mg SC at months 0, 6, 12 (converted to 408.2 nmol using MW ≈ 147 kDa). Concentrations in nmol/L, volumes in L, time in hours.

Estimation: SAEM, 300 burn-in + 400 EM iterations.

---

## Results — parameter recovery

| Parameter | True (Table 3) | Recovered | % error | %RSE | Shrinkage |
|---|---|---|---|---|---|
| ka (1/h) | 0.0078 | 0.00876 | +12.3% | 2.6 | 18.9% |
| Vc (L) | 1.58 | 2.169 | +37.3% | 22.4 | 30.4% |
| Vp (L) | 6.06 | 6.190 | +2.2% | 3.1 | 87.2% |
| CL (L/h) | 0.006 | 0.00706 | +17.7% | 1.2 | 46.3% |
| Q (L/h) | 0.20 | *fixed* | — | — | — |
| kint (1/h) | 0.022 | 0.01680 | −23.7% | 9.8 | **94.2%** |
| Kss (nmol/L) | 1.56 | 1.263 | −19.1% | **660** | **92.1%** |
| ksyn (1/h) | 0.01 | 0.00765 | −23.5% | 6.6 | **90.8%** |
| R0 (nmol/L) | 15.23 | 15.638 | +2.7% | 11.4 | **86.6%** |

All eight estimated parameters recovered within ±37%; six within ±25%. Seven of eight have %RSE below 25%.

`Q` was fixed to its literature value. In the first fit it returned 499% RSE with 902% CV between-subject variability — it was not being estimated in any meaningful sense. Fixing an unidentifiable parameter to a known value is a standard remedy (the same approach used for Ka in this portfolio's Theophylline project).

---

## The main finding: TMDD parameters are precise-looking but uninformed

Point estimates for `kint`, `Kss`, `ksyn` and `R0` land within ~25% of truth. Taken alone that looks like successful recovery. It isn't.

**Shrinkage on those four parameters is 86.6% to 94.2%.**

Shrinkage above roughly 30% conventionally indicates that the data carry little individual-level information and the estimator is falling back on the population mean. At 94%, individual η estimates are essentially uninformative — the model cannot distinguish one subject's target dynamics from another's.

`Kss` makes the point most sharply: the point estimate is 19% from truth, but its %RSE is 660% and its between-subject variability is 182% CV. The estimate landed near the true value without the data meaningfully constraining it.

**This reproduces the source paper's own reported limitation.** Choi et al. observed high shrinkage in the TMDD-related parameters (R0, kint, KSS) and attributed it to the limited informativeness of PK-only data for receptor-mediated dynamics, noting that the absence of target concentration measurements is a known factor impairing TMDD parameter identifiability. They reported this from 615 real patients; this project reproduces it from 80 simulated ones.

The finding survived a complete redesign of the sampling schedule between v1 and v2, which is evidence that it reflects a property of the model–data combination rather than an artifact of one design choice.

---

## Development history: v1 → v2

The first fit converged but was diagnosably wrong. Documenting how it was found matters more than the fact that it happened.

### v1 symptoms

- `mean(DV − IPRED) = +5.014` — systematic underprediction
- \>75% of residuals positive
- max IPRED 70.6 against observed max 120.8 — model could not reach the peaks
- `add.err` estimated at 2.68 against a true 0.72
- `Q` at 499% RSE

### Hypotheses tested and eliminated

Each was tested by direct measurement rather than assumed. Test scripts are in this repository.

| # | Hypothesis | Test | Result |
|---|---|---|---|
| 1 | Negative-value clipping biased the data upward | Simulated the clipping directly | Rejected — produced −0.02 bias, needed +5 |
| 2 | Sampling too sparse to identify Vc/Q/Vp | Compared true vs fitted curves at sampled times | Rejected — curves differ 22–40%, RMSE 6.15; data *does* distinguish them |
| 3 | Compartment (CMT) mis-mapping | Compared IPRED against Ctot | Rejected — identical |
| 4 | Between-subject variances collapsed | Re-fit with variances initialised at true values | **Partly confirmed** — collapse was real (eta_Vc 20× too small) and fixing it held (ratio 1.15×), but bias persisted at +5.67 |
| 5 | Simulation and fit models differ structurally | Solved both at identical parameters | Rejected — agreement to ~1e-12 |
| 6 | Data-handling in the NONMEM-style dataset | Solved via pk_data and via eventTable | Rejected — max Ctot 37.23 vs 37.24 |
| 7 | Washout sampling starved the likelihood | Redesigned schedule, re-fit | **Partly confirmed** — see below |

### What v2 changed

**(a) Sampling design.** v1 sampled at 3600 h and 4380 h post-dose, where the typical curve sits at 1.41 and 0.095 nmol/L. Against an additive residual error of 0.72 nmol/L, noise at the 4380 h trough is ~7× the signal. SAEM weights observations by inverse variance, so these uninformative points carried disproportionate weight.

The source paper avoided this without needing to design around it: their assay LLOQ was 20 ng/mL and 26.49% of post-dose samples fell below it and were set to MDV=1, excluded from the likelihood. This project simulated the assay but not its detection limit, generating observations the real dataset never contained.

v2 offsets: 12, 24, 72, 168, 336, 720, 1440, 2160 h. Minimum typical concentration rises from 0.095 to 8.55 nmol/L; worst-case noise-to-signal falls from 7.0 to 0.08.

**(b) Initial variance estimates** set at the true ω² values rather than below them (v1 initialised eta_Vc at 0.20 against a true 0.3225 and it collapsed 20-fold).

**(c) Q fixed** at its literature value.

### v1 → v2 outcome

| Metric | v1 | v2 |
|---|---|---|
| mean(DV − IPRED) | +5.014 | **+3.327** |
| % positive residuals | >75% | 75.4% |
| max IPRED | 70.6 | 79.8 |
| ka % error | +54% | **+12.3%** |
| Vc % error | +117% | **+37.3%** |
| Vp % error | +21% | **+2.2%** |
| R0 % error | −26% | **+2.7%** |
| kint %RSE | 24.6 | **9.8** |
| R0 %RSE | 20.2 | **11.4** |
| Q %RSE | 499 | fixed |

Structural parameter recovery and precision improved substantially. The systematic bias reduced by a third but did not resolve.

---

## Known limitations

Stated plainly, because they are not resolved.

**1. Residual error is inflated and the cause is not identified.** `add.err` estimated at 7.88 against a true 0.72 (~11×), and `prop.err` at 0.298 against 0.07. In v1 these were 2.68 and 0.589. The v2 sampling redesign improved every structural parameter but made the residual error terms worse, which indicates the error model is absorbing a misspecification that the sampling change did not address.

**2. Systematic underprediction persists.** 75.4% of residuals are positive with a mean of +3.33 nmol/L. Seven hypotheses were tested and eliminated; the cause remains unexplained. A model that underpredicts systematically is not fit for inference, and no inference beyond parameter recovery is drawn here.

**3. Kss is not identifiable in this design.** 660% RSE. The point estimate should not be treated as meaningful.

**4. Simulation-based, not clinical.** See the methodology note above.

**5. Single replicate.** A proper estimability assessment would use a stochastic simulation and estimation (SSE) study across many replicate datasets. This is one dataset, one fit; the parameter errors reported carry no confidence interval across replicates.

---

## Repository contents

```
denosumab_QSS_TMDD_v2.R          primary script (simulation, fit, diagnostics)
denosumab_QSS_TMDD_nlmixr2.R     v1 script, retained for the v1→v2 comparison

diagnostic_test_v2.R             hypothesis 4 test (variance initialisation)
structural_test.R                hypothesis 5 test (sim vs fit model equivalence)
data_handling_test_v2.R          hypothesis 6 test (pk_data vs eventTable)
test_error_model.py              hypothesis 1 test (clipping bias)
test_identifiability.py          hypothesis 2 test (sampling design)
sanity_check_qss.py              pre-implementation ODE verification
  simulated_denosumab_data_v2.csv   simulated dataset (1920 rows, 80 subjects)
  fit_v2.rds                        fitted model object
  parameter_recovery_v2.csv         true vs recovered comparison
  v2_01_dv_vs_pred.png              goodness-of-fit: DV vs PRED
  v2_02_dv_vs_ipred.png             goodness-of-fit: DV vs IPRED
  v2_03_cwres_time.png              CWRES vs time
  v2_04_cwres_pred.png              CWRES vs PRED
  v2_05_vpc.png                     visual predictive check
  session_info_v2.txt               R session details
```

---

## Reproducing

```r
install.packages(c("nlmixr2", "rxode2", "ggplot2", "dplyr",
                   "xpose.nlmixr2", "vpc"))

source("denosumab_QSS_TMDD_v2.R")
```

Runtime approximately 45–60 minutes for the SAEM fit on Posit Cloud free tier. Developed on R 4.6.1, Linux; full session details in `outputs/session_info_v2.txt`.

Note: `eventTable()` is called without `amount.units`/`time.units` arguments. Supplying them requires the `units` package, which is not present on a default Posit Cloud image and caused a silent simulation failure during development.

---

## What this project demonstrates

- Implementing a QSS-approximation TMDD model as a four-state ODE system in rxode2/nlmixr2
- Simulation-based parameter recovery against a published, peer-reviewed parameter set
- Using shrinkage and %RSE to distinguish a parameter that is *estimated* from one that is *identified* — the central point of the main finding
- Systematic diagnosis: seven hypotheses, each tested by direct measurement, with the test scripts retained as evidence
- Independently reproducing a published identifiability limitation
- Documenting unresolved problems rather than presenting a clean result

---

Part of a pharmacometrics portfolio in preparation for doctoral research in Quantitative Systems Pharmacology.

Licensed MIT. Ground-truth parameters reproduced from Choi et al. (2025) under CC-BY with attribution.
