# =============================================================================
# PROJECT 3: QSS-Approximation TMDD Model — Denosumab
# Tool:       nlmixr2 / rxode2 in R
# Ground truth: Choi S, Park S, Jung J, Baek S, Lim H-S (2025). Population
#               pharmacokinetics/pharmacodynamics analysis confirming
#               biosimilarity of SB16 to reference denosumab. Front Pharmacol
#               16:1631034. doi:10.3389/fphar.2025.1631034 (CC-BY, open access)
#               -- Table 3 gives the full PK parameter set used below.
# Model theory: Dua P, Hawkins E, van der Graaf PH (2015). A Tutorial on
#               Target-Mediated Drug Disposition (TMDD) Models. CPT
#               Pharmacometrics Syst Pharmacol 4(6):324-337.
#               Gibiansky L, Gibiansky E, Kakkar T, Ma P (2008). Approximations
#               of the target-mediated drug disposition model and
#               identifiability of model parameters. J Pharmacokinet
#               Pharmacodyn 35(5):573-591.
# Author:     Jyotheeshwar Akshay Ravikumar
# Purpose:    Portfolio Project 3 — fixes Project 2's core gaps:
#               (1) real, published, fitted parameter set as ground truth
#                   (not made-up synthetic values)
#               (2) full QSS parameter estimation (not a simplified fit)
#               (3) practical identifiability analysis on the QSS-specific
#                   parameters (KSS, R0, kint) that the source paper itself
#                   flags as high-shrinkage / poorly identified from PK-only
#                   data — this project directly tests that claim
#
# METHODOLOGY NOTE (read this before the README):
#   Denosumab concentration-time DATA are not publicly available (no patient-
#   level dataset was released with the Choi et al. paper). What IS public and
#   real is their FITTED PARAMETER TABLE (Table 3), estimated by them from
#   6,583 real serum concentrations in 615 real subjects using NONMEM/Monolix.
#   This script simulates a dataset FROM those real fitted parameters (adding
#   realistic between-subject variability and residual error using their
#   reported IIV and error model), then re-fits the QSS model to recover them.
#   This is NOT the same as fitting real patient data -- be explicit about
#   this distinction in the README, exactly as this comment is explicit here.
#   It IS a meaningful upgrade from Project 2: the ground truth is a real,
#   published, peer-reviewed fitted model, not arbitrary invented parameters.
#
# UNITS NOTE: all concentrations in nmol/L, all volumes in L, time in hours,
# consistent with how Table 3 reports its parameters. Dose converted from mg
# to nmol using denosumab's approximate molecular weight (~147 kDa, IgG2 mAb).
# =============================================================================


# =============================================================================
# STAGE 1: INSTALL AND LOAD PACKAGES
# =============================================================================
# install.packages("nlmixr2")
# install.packages("rxode2")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("xpose.nlmixr2")
# install.packages("vpc")

library(nlmixr2)
library(rxode2)
library(ggplot2)
library(dplyr)
library(xpose.nlmixr2)

set.seed(42)

# =============================================================================
# STAGE 2: GROUND-TRUTH PARAMETERS (Choi et al. 2025, Table 3 — PMO population)
# =============================================================================
# These are the REAL fitted values from the paper. We treat them as the "true"
# population parameters we are trying to recover by simulating data and
# re-fitting. Units: rates in 1/h, volumes in L, concentrations in nmol/L.

true_params <- list(
  ka    = 0.0078,      # 1/h, absorption rate constant (PMO patients)
  Vc    = 1.58,         # L, central volume (apparent, Vc/F)
  Vp    = 6.06,         # L, peripheral volume (apparent, Vp/F)
  CL    = 0.006,        # L/h, clearance (Caucasian, apparent CL/F)
  Q     = 0.20,          # L/h, inter-compartmental clearance (PMO)
  kint  = 0.022,        # 1/h, drug-target complex internalization rate
  Kss   = 1.56,          # nmol/L, QSS equilibrium constant
  ksyn  = 0.01,          # 1/h, target (RANKL) synthesis rate
  R0    = 15.23          # nmol/L, baseline target concentration (PMO)
)
# kdeg is NOT separately reported in Table 3. It is derived from the
# steady-state assumption ksyn = kdeg * R0 (stated explicitly in the paper's
# Methods, their Eq. 7): kdeg = ksyn / R0. This derivation is documented here
# and MUST be documented identically in the README -- an undocumented derived
# parameter is exactly the kind of silent gap that hurt Project 2.
true_params$kdeg <- true_params$ksyn / true_params$R0
cat("Derived kdeg =", true_params$kdeg, "1/h  (= ksyn / R0)\n")

# Inter-individual variability (IIV, %CV from Table 3), converted to
# log-normal variance (omega^2) via: omega^2 = log(1 + CV^2)
# Table 3 IIVs: ka 56.57%, VC/F 61.69%, VP/F 15.55%, kint 7.88%, KSS 58.12%,
#               ksyn 22.46%, R0 158.3%, CL/F 26.39%, Q/F 295.99%
cv_to_omega2 <- function(cv_pct) log(1 + (cv_pct / 100)^2)

true_omega2 <- list(
  ka   = cv_to_omega2(56.57),
  Vc   = cv_to_omega2(61.69),
  Vp   = cv_to_omega2(15.55),
  kint = cv_to_omega2(7.88),
  Kss  = cv_to_omega2(58.12),
  ksyn = cv_to_omega2(22.46),
  R0   = cv_to_omega2(158.3),
  CL   = cv_to_omega2(26.39),
  Q    = cv_to_omega2(295.99)
)

# Residual error (Table 3): combined additive + proportional
#   sigma_add  = 0.72 nmol/L
#   sigma_prop = 0.07 (7%)
true_sigma_add  <- 0.72
true_sigma_prop <- 0.07

# Dose: 60 mg SC at months 0, 6, 12 (matching the paper's simulated regimen)
# Converted to nmol using denosumab MW ~147 kDa (IgG2 monoclonal antibody)
MW_kDa <- 147
dose_mg <- 60
dose_nmol <- (dose_mg * 1e6) / (MW_kDa * 1e3)  # mg -> ng -> nmol
cat("Dose in nmol:", round(dose_nmol, 1), "\n")


# =============================================================================
# STAGE 3: SIMULATE A REALISTIC PMO DATASET FROM THE TRUE PARAMETERS
# =============================================================================
# Sampling schedule loosely matches the paper's Phase III design (sparse,
# pre-dose + several post-dose timepoints per 6-month interval), scaled down
# to a portfolio-appropriate N. We simulate N subjects, 3 doses (months
# 0/6/12), with realistic sparse sampling -- NOT rich sampling -- because
# sparse sampling around a QSS-TMDD model is exactly the condition under
# which the paper reports high shrinkage on KSS, R0, and kint. Reproducing
# that data-sparsity condition is what makes the identifiability check in
# Stage 6 meaningful rather than trivial.

n_subjects <- 80          # portfolio-scale; paper used 456 PMO patients
month_h <- 30.44 * 24     # hours per month

# Sparse sampling times within each 6-month interval (hours since that dose)
sample_offsets_h <- c(12, 24, 168, 720, 2160, 3600, 4380)  # ~0.5d,1d,1wk,1mo,3mo,5mo,pre-next-dose
dose_times_h <- c(0, 6, 12) * month_h

qss_model_rxode <- rxode2({
  ka    <- exp(lka + eta_ka)
  Vc    <- exp(lVc + eta_Vc)
  Vp    <- exp(lVp + eta_Vp)
  CL    <- exp(lCL + eta_CL)
  Q     <- exp(lQ  + eta_Q)
  kint  <- exp(lkint + eta_kint)
  Kss   <- exp(lKss + eta_Kss)
  ksyn  <- exp(lksyn + eta_ksyn)
  R0    <- exp(lR0 + eta_R0)
  kdeg  <- ksyn / R0

  # QSS closed-form free-drug concentration
  disc  <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot
  discP <- max(disc, 0)
  C     <- 0.5 * ((Ctot - Rtot - Kss) + sqrt(discP))
  Cfree <- max(C, 0)
  RC    <- Ctot - Cfree            # bound complex concentration

  d/dt(depot) <- -ka * depot
  d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * Cfree - (Q / Vc) * Cfree +
                 (Q / Vp) * Cp - kint * RC
  d/dt(Cp)    <- (Q / Vc) * Cfree - (Q / Vp) * Cp
  d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC
})

# Build subject-level parameter sets (log-normal IIV) and initial conditions
sim_params <- data.frame(
  id = 1:n_subjects,
  eta_ka   = rnorm(n_subjects, 0, sqrt(true_omega2$ka)),
  eta_Vc   = rnorm(n_subjects, 0, sqrt(true_omega2$Vc)),
  eta_Vp   = rnorm(n_subjects, 0, sqrt(true_omega2$Vp)),
  eta_CL   = rnorm(n_subjects, 0, sqrt(true_omega2$CL)),
  eta_Q    = rnorm(n_subjects, 0, sqrt(true_omega2$Q)),
  eta_kint = rnorm(n_subjects, 0, sqrt(true_omega2$kint)),
  eta_Kss  = rnorm(n_subjects, 0, sqrt(true_omega2$Kss)),
  eta_ksyn = rnorm(n_subjects, 0, sqrt(true_omega2$ksyn)),
  eta_R0   = rnorm(n_subjects, 0, sqrt(true_omega2$R0))
)

# NOTE: R0 sets the baseline Rtot(0) for each subject -- this must be wired
# into the event table's initial condition per-subject, not left at a single
# population value, or the simulation silently loses the between-subject
# baseline-target variability that Table 3 explicitly reports (CV 158.3%,
# the single largest IIV in the whole parameter set).

theta <- log(c(
  lka = true_params$ka, lVc = true_params$Vc, lVp = true_params$Vp,
  lCL = true_params$CL, lQ = true_params$Q, lkint = true_params$kint,
  lKss = true_params$Kss, lksyn = true_params$ksyn, lR0 = true_params$R0
))
names(theta) <- c("lka","lVc","lVp","lCL","lQ","lkint","lKss","lksyn","lR0")

sim_results <- list()

for (i in 1:n_subjects) {
  p <- sim_params[i, ]
  R0_i <- exp(theta["lR0"] + p$eta_R0)

  ev <- eventTable(amount.units = "nmol", time.units = "hours")
  for (dt in dose_times_h) {
    ev$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
  }
  samp_times <- sort(unique(as.vector(outer(dose_times_h, sample_offsets_h, "+"))))
  samp_times <- samp_times[samp_times <= max(dose_times_h) + max(sample_offsets_h)]
  ev$add.sampling(samp_times)

  inits <- c(depot = 0, Ctot = 0, Cp = 0, Rtot = as.numeric(R0_i))

  subj_params <- c(theta,
                    eta_ka = p$eta_ka, eta_Vc = p$eta_Vc, eta_Vp = p$eta_Vp,
                    eta_CL = p$eta_CL, eta_Q = p$eta_Q, eta_kint = p$eta_kint,
                    eta_Kss = p$eta_Kss, eta_ksyn = p$eta_ksyn, eta_R0 = p$eta_R0)

  sol <- rxSolve(qss_model_rxode, subj_params, ev, inits = inits)
  sol_df <- as.data.frame(sol)
  sol_df$ID <- i
  sim_results[[i]] <- sol_df
}

sim_all <- bind_rows(sim_results)

# Add combined residual error to Ctot -> observed DV
sim_all <- sim_all %>%
  mutate(
    IPRED = Ctot,
    err_add  = rnorm(n(), 0, true_sigma_add),
    err_prop = rnorm(n(), 0, true_sigma_prop),
    DV = IPRED * (1 + err_prop) + err_add,
    DV = pmax(DV, 0)   # concentrations can't be negative
  )

cat("Simulated dataset: ", nrow(sim_all), " total records across ", n_subjects, " subjects\n")


# =============================================================================
# STAGE 4: FORMAT AS NLMIXR2 (NONMEM-STYLE) DATASET
# =============================================================================
obs_rows <- sim_all %>%
  filter(time %in% unlist(lapply(dose_times_h, function(dt) dt + sample_offsets_h))) %>%
  transmute(ID, TIME = time, AMT = 0, DV = DV, EVID = 0, CMT = 2)

dose_rows <- expand.grid(ID = 1:n_subjects, TIME = dose_times_h) %>%
  transmute(ID, TIME, AMT = dose_nmol, DV = 0, EVID = 1, CMT = 1)

pk_data <- bind_rows(obs_rows, dose_rows) %>%
  arrange(ID, TIME, desc(EVID))

write.csv(pk_data, "outputs/simulated_denosumab_data.csv", row.names = FALSE)
cat("Formatted dataset saved. Rows:", nrow(pk_data), "\n")


# =============================================================================
# STAGE 5: DEFINE AND FIT THE QSS-TMDD MODEL IN NLMIXR2
# =============================================================================
# This is the model we FIT to the simulated data -- structurally identical to
# the simulation model above, but now with unknown parameters to be estimated
# starting from reasonable (not identical-to-truth) initial guesses.

qss_fit_model <- function() {
  ini({
    lka   <- log(0.01)      # initial guess, not equal to true value
    lVc   <- log(2)
    lVp   <- log(5)
    lCL   <- log(0.008)
    lQ    <- log(0.15)
    lkint <- log(0.02)
    lKss  <- log(1.5)
    lksyn <- log(0.01)
    lR0   <- log(12)

    eta_ka   ~ 0.2
    eta_Vc   ~ 0.2
    eta_Vp   ~ 0.1
    eta_CL   ~ 0.1
    eta_Q    ~ 0.3
    eta_kint ~ 0.05
    eta_Kss  ~ 0.2
    eta_ksyn ~ 0.05
    eta_R0   ~ 0.5

    add.err  <- 0.5
    prop.err <- 0.1
  })

  model({
    ka    <- exp(lka + eta_ka)
    Vc    <- exp(lVc + eta_Vc)
    Vp    <- exp(lVp + eta_Vp)
    CL    <- exp(lCL + eta_CL)
    Q     <- exp(lQ  + eta_Q)
    kint  <- exp(lkint + eta_kint)
    Kss   <- exp(lKss + eta_Kss)
    ksyn  <- exp(lksyn + eta_ksyn)
    R0    <- exp(lR0 + eta_R0)
    kdeg  <- ksyn / R0

    disc  <- (Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot
    discP <- max(disc, 0)
    C     <- 0.5 * ((Ctot - Rtot - Kss) + sqrt(discP))
    Cfree <- max(C, 0)
    RC    <- Ctot - Cfree

    d/dt(depot) <- -ka * depot
    d/dt(Ctot)  <- (ka * depot) / Vc - (CL / Vc) * Cfree - (Q / Vc) * Cfree +
                   (Q / Vp) * Cp - kint * RC
    d/dt(Cp)    <- (Q / Vc) * Cfree - (Q / Vp) * Cp
    d/dt(Rtot)  <- ksyn - kdeg * (Rtot - RC) - kint * RC

    Rtot(0) <- R0   # baseline target at steady state, per-subject via eta_R0

    Ctot ~ add(add.err) + prop(prop.err)
  })
}

fit <- nlmixr2(
  object  = qss_fit_model,
  data    = pk_data,
  est     = "saem",
  control = saemControl(seed = 42, print = 50, nBurn = 300, nEm = 400)
)

print(fit)
summary(fit)


# =============================================================================
# STAGE 6: COMPARE RECOVERED PARAMETERS AGAINST GROUND TRUTH
# =============================================================================
# THIS TABLE IS THE HEADLINE RESULT OF THE PROJECT.

fixed_effects <- fixef(fit)

comparison <- data.frame(
  parameter = c("ka","Vc","Vp","CL","Q","kint","Kss","ksyn","R0"),
  true_value = c(true_params$ka, true_params$Vc, true_params$Vp, true_params$CL,
                 true_params$Q, true_params$kint, true_params$Kss,
                 true_params$ksyn, true_params$R0),
  recovered_value = exp(c(fixed_effects["lka"], fixed_effects["lVc"],
                          fixed_effects["lVp"], fixed_effects["lCL"],
                          fixed_effects["lQ"], fixed_effects["lkint"],
                          fixed_effects["lKss"], fixed_effects["lksyn"],
                          fixed_effects["lR0"]))
)
comparison$pct_error <- 100 * (comparison$recovered_value - comparison$true_value) /
                         comparison$true_value

print(comparison)
write.csv(comparison, "outputs/parameter_recovery_comparison.csv", row.names = FALSE)

cat("\n=== INTERPRETATION GUIDE (fill in with actual numbers after running) ===\n")
cat("Expect CL, Vc, Vp, ka to recover well (rich linear-PK information in sparse\n")
cat("PK-only data). Expect KSS, R0, kint to show LARGER errors and/or high\n")
cat("shrinkage -- this matches the source paper's own reported finding that\n")
cat("these three parameters are poorly identified from PK-only data (no\n")
cat("target/RANKL concentration measurements). If your fit reproduces this\n")
cat("pattern, that IS the identifiability result -- document it as such.\n")


# =============================================================================
# STAGE 7: GOODNESS-OF-FIT DIAGNOSTICS
# =============================================================================
xpdb <- xpose_data_nlmixr2(fit)

dv_vs_pred(xpdb) +
  labs(title = "DV vs Population Predicted (PRED)",
       subtitle = "Denosumab QSS-TMDD model — nlmixr2 SAEM")

dv_vs_ipred(xpdb) +
  labs(title = "DV vs Individual Predicted (IPRED)")

res_vs_idv(xpdb, res = "CWRES") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "CWRES vs Time")

res_vs_pred(xpdb, res = "CWRES") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "CWRES vs Population Predicted")


# =============================================================================
# STAGE 8: VISUAL PREDICTIVE CHECK
# =============================================================================
vpc_fit <- vpcPlot(
  fit,
  n    = 300,
  bins = 8,
  pi   = c(0.05, 0.95),
  ci   = c(0.025, 0.975)
)
print(vpc_fit)


# =============================================================================
# STAGE 9: PRACTICAL IDENTIFIABILITY ANALYSIS ON QSS PARAMETERS
# =============================================================================
# Profile likelihood is the gold-standard approach; nlmixr2 supports this via
# repeated re-fits with a parameter fixed across a grid of values around its
# estimate, tracking the change in objective function value (OFV).
# We focus on KSS, R0, and kint -- the three parameters the source paper
# itself flags as high-shrinkage / poorly identified from PK-only data.
#
# This is a computationally expensive step (multiple re-fits per parameter).
# Run this ONLY after Stage 5-6 confirm the base fit is stable and sensible.

profile_identifiability <- function(fit, param_name, param_col, grid_frac = 0.3, n_points = 7) {
  # param_col: the theta name in the model, e.g. "lKss"
  center <- fixef(fit)[param_col]
  grid_vals <- seq(center * (1 - grid_frac), center * (1 + grid_frac), length.out = n_points)
  ofv_profile <- numeric(n_points)

  cat("\nProfiling", param_name, "across", n_points, "fixed values...\n")
  for (i in seq_along(grid_vals)) {
    # Refit with this theta FIXED at grid_vals[i] -- requires rebuilding the
    # model function with that parameter as fixed() instead of estimated.
    # (Implementation detail: nlmixr2's `fixed()` ini() syntax, as used for
    # Ka in Project 1 v2, is the mechanism here -- adapt qss_fit_model() to
    # accept a fixed value for the target parameter and re-run nlmixr2().)
    cat("  Grid point", i, ":", param_name, "=", exp(grid_vals[i]), "\n")
    # Placeholder -- fill in actual refit call once base model is confirmed
    # working on Posit Cloud; this loop structure is the scaffold.
  }
  return(data.frame(grid_vals = grid_vals, ofv = ofv_profile))
}

cat("\n=== IDENTIFIABILITY ANALYSIS SCAFFOLD ===\n")
cat("Function profile_identifiability() above is a SCAFFOLD, not a finished\n")
cat("analysis -- it needs the per-parameter fixed() refit loop filled in\n")
cat("once Stage 5's base fit is confirmed stable. Do not claim a completed\n")
cat("identifiability analysis in the README until this loop actually runs\n")
cat("and produces real OFV profiles.\n")


# =============================================================================
# STAGE 10: SAVE OUTPUTS
# =============================================================================
dir.create("outputs", showWarnings = FALSE)

write.csv(as.data.frame(fixed_effects), "outputs/population_parameters.csv", row.names = TRUE)

ggsave("outputs/01_dv_vs_pred.png", dv_vs_pred(xpdb), width = 6, height = 5)
ggsave("outputs/02_dv_vs_ipred.png", dv_vs_ipred(xpdb), width = 6, height = 5)
ggsave("outputs/03_cwres_time.png", res_vs_idv(xpdb, res = "CWRES"), width = 6, height = 5)
ggsave("outputs/04_cwres_pred.png", res_vs_pred(xpdb, res = "CWRES"), width = 6, height = 5)
ggsave("outputs/05_vpc.png", vpc_fit, width = 7, height = 5)

writeLines(capture.output(sessionInfo()), "outputs/session_info.txt")

cat("\n=== All outputs saved to /outputs/ ===\n")
cat("NEXT: fill in Stage 9's identifiability refit loop, then write the\n")
cat("README documenting: (1) ground-truth source and its CC-BY license,\n")
cat("(2) simulated-from-real-parameters methodology, explicitly labeled as\n")
cat("such, (3) the parameter recovery comparison table, (4) the\n")
cat("identifiability findings on KSS/R0/kint versus the paper's own claim.\n")
