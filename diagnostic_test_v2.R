# =============================================================================
# PROJECT 3 — DIAGNOSTIC TEST (not the final analysis)
#
# PURPOSE: test ONE hypothesis cheaply before committing to a full re-fit.
#
# HYPOTHESIS: in the first fit, the between-subject variances (etas) collapsed
#   toward zero because I initialised them BELOW their true values:
#       eta_Vc initialised 0.20  vs true 0.32   -> estimated 0.0156 (20x collapse)
#       eta_R0 initialised 0.50  vs true ~1.00  -> estimated 0.0254 (40x collapse)
#   With the etas collapsed, SAEM cannot let individuals differ, so it must find
#   one compromise typical Vc. That compromise sat high (3.43 vs true 1.58),
#   which produced systematic underprediction of peaks (mean IRES +5.0).
#
# TEST: re-fit with variances initialised AT the true values, on a small
#   dataset with short SAEM settings. ~5-10 min instead of ~50.
#
#   PASS  = eta_Vc stays near 0.32 and lVc lands near log(1.58) = 0.457
#   FAIL  = eta_Vc collapses again -> hypothesis wrong, look elsewhere
#
# This script REUSES objects already in memory from the main script:
#   true_params, true_omega2, pk_data, dose_times_h, sample_offsets_h,
#   dose_nmol, month_h, theta, cv_to_omega2
# Run the main script's Stages 1-4 first (or restore the session) before this.
# =============================================================================

library(nlmixr2)
library(rxode2)
library(dplyr)

set.seed(42)

# ---- guard: make sure prerequisites exist -----------------------------------
needed <- c("true_params", "true_omega2", "pk_data", "dose_times_h",
            "sample_offsets_h", "dose_nmol", "month_h")
missing <- needed[!sapply(needed, exists)]
if (length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "),
       "\nRun Stages 1-4 of the main script first.")
}

# =============================================================================
# STEP 1: SHRINK THE DATASET TO 20 SUBJECTS (speed)
# =============================================================================
pk_small <- pk_data %>% filter(ID <= 20)
cat("Diagnostic dataset:", nrow(pk_small), "rows,",
    length(unique(pk_small$ID)), "subjects\n")

# =============================================================================
# STEP 2: MODEL WITH CORRECTED INITIAL VARIANCES
# =============================================================================
# Two changes vs the original fit model:
#   (a) every eta initialised AT its true omega^2 (from true_omega2)
#   (b) lQ FIXED at its true value -- it had 499% RSE in the first fit, i.e.
#       it was not being estimated in any meaningful sense. Fixing an
#       unidentifiable parameter to a known value is the same standard,
#       defensible move used for Ka in Project 1 v2.
#
# Everything else (structure, equations, error model) is UNCHANGED, so this
# isolates the initialisation as the single variable being tested.

cat("\nInitial variances being used (= true values):\n")
for (nm in names(true_omega2)) {
  cat(sprintf("  eta_%-5s = %.4f\n", nm, true_omega2[[nm]]))
}
cat(sprintf("\nlQ FIXED at log(%.2f) = %.4f\n\n",
            true_params$Q, log(true_params$Q)))

qss_diag_model <- function() {
  ini({
    lka   <- log(0.01)
    lVc   <- log(2)
    lVp   <- log(5)
    lCL   <- log(0.008)
    lQ    <- fixed(log(0.20))    # FIXED: unidentifiable in the first fit
    lkint <- log(0.02)
    lKss  <- log(1.5)
    lksyn <- log(0.01)
    lR0   <- log(12)

    # Initialised AT the true omega^2 values (this is the variable under test)
    eta_ka   ~ 0.3225   # true: cv_to_omega2(56.57)
    eta_Vc   ~ 0.3225   # true: cv_to_omega2(61.69)  -- was 0.20, collapsed
    eta_Vp   ~ 0.0239   # true: cv_to_omega2(15.55)
    eta_CL   ~ 0.0680   # true: cv_to_omega2(26.39)
    eta_kint ~ 0.0062   # true: cv_to_omega2(7.88)
    eta_Kss  ~ 0.2870   # true: cv_to_omega2(58.12)
    eta_ksyn ~ 0.0492   # true: cv_to_omega2(22.46)
    eta_R0   ~ 1.0195   # true: cv_to_omega2(158.3)  -- was 0.50, collapsed

    add.err  <- 0.72    # start at the TRUE value rather than 0.5
    prop.err <- 0.07    # start at the TRUE value rather than 0.1
  })

  model({
    ka    <- exp(lka + eta_ka)
    Vc    <- exp(lVc + eta_Vc)
    Vp    <- exp(lVp + eta_Vp)
    CL    <- exp(lCL + eta_CL)
    Q     <- exp(lQ)              # no eta: fixed parameter
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

    Rtot(0) <- R0

    Ctot ~ add(add.err) + prop(prop.err)
  })
}

# =============================================================================
# STEP 3: SHORT FIT
# =============================================================================
# nBurn 150 + nEm 200 instead of 300 + 400. Enough to see whether the etas
# hold or collapse, which is all this test needs to determine.

cat("Starting diagnostic fit (expect ~5-10 min)...\n\n")

fit_diag <- nlmixr2(
  object  = qss_diag_model,
  data    = pk_small,
  est     = "saem",
  control = saemControl(seed = 42, print = 25, nBurn = 150, nEm = 200)
)

# =============================================================================
# STEP 4: VERDICT
# =============================================================================
cat("\n\n")
cat("======================================================================\n")
cat("DIAGNOSTIC VERDICT\n")
cat("======================================================================\n\n")

est_omega <- diag(fit_diag$omega)
fe <- fixef(fit_diag)

cat("BETWEEN-SUBJECT VARIANCES (the thing under test):\n")
cat(sprintf("%-10s %12s %12s %10s\n", "eta", "true", "estimated", "ratio"))
truth_map <- c(eta_ka = 0.3225, eta_Vc = 0.3225, eta_Vp = 0.0239,
               eta_CL = 0.0680, eta_kint = 0.0062, eta_Kss = 0.2870,
               eta_ksyn = 0.0492, eta_R0 = 1.0195)
for (nm in names(truth_map)) {
  if (nm %in% names(est_omega)) {
    tv <- truth_map[[nm]]; ev <- est_omega[[nm]]
    cat(sprintf("%-10s %12.4f %12.4f %9.2fx\n", nm, tv, ev, ev / tv))
  }
}

cat("\nKEY STRUCTURAL PARAMETERS:\n")
cat(sprintf("%-8s %12s %12s %10s\n", "param", "true", "recovered", "%error"))
show_par <- function(label, lname, truth) {
  if (lname %in% names(fe)) {
    rec <- exp(fe[[lname]])
    cat(sprintf("%-8s %12.5f %12.5f %9.1f%%\n",
                label, truth, rec, 100 * (rec - truth) / truth))
  }
}
show_par("ka",   "lka",   true_params$ka)
show_par("Vc",   "lVc",   true_params$Vc)
show_par("Vp",   "lVp",   true_params$Vp)
show_par("CL",   "lCL",   true_params$CL)
show_par("kint", "lkint", true_params$kint)
show_par("Kss",  "lKss",  true_params$Kss)
show_par("ksyn", "lksyn", true_params$ksyn)
show_par("R0",   "lR0",   true_params$R0)

cat("\nRESIDUAL ERROR:\n")
cat(sprintf("  add.err  true 0.72  ->  estimated %.4f\n", fe[["add.err"]]))
cat(sprintf("  prop.err true 0.07  ->  estimated %.4f\n", fe[["prop.err"]]))

cat("\nSYSTEMATIC BIAS CHECK (this is what was broken):\n")
ires <- fit_diag$DV - fit_diag$IPRED
cat(sprintf("  mean(DV - IPRED) = %+.4f   (first fit: +5.014; want near 0)\n",
            mean(ires, na.rm = TRUE)))
cat(sprintf("  frac positive    = %.1f%%      (first fit: >75%%; want ~50%%)\n",
            100 * mean(ires > 0, na.rm = TRUE)))
cat(sprintf("  max IPRED        = %.1f       (first fit: 70.6; data max ~120)\n",
            max(fit_diag$IPRED, na.rm = TRUE)))

cat("\n----------------------------------------------------------------------\n")
vc_ratio <- if ("eta_Vc" %in% names(est_omega)) est_omega[["eta_Vc"]] / 0.3225 else NA
mean_ires <- mean(ires, na.rm = TRUE)
if (!is.na(vc_ratio) && vc_ratio > 0.5 && abs(mean_ires) < 2) {
  cat("PASS: variances held and systematic bias is much reduced.\n")
  cat("      -> Hypothesis CONFIRMED. Scale back up to 80 subjects,\n")
  cat("         nBurn 300 / nEm 400, and this becomes the final v2 fit.\n")
} else {
  cat("FAIL: the problem persists.\n")
  cat(sprintf("      eta_Vc ratio = %.2fx (want > 0.5)\n", vc_ratio))
  cat(sprintf("      mean IRES    = %+.2f (want |value| < 2)\n", mean_ires))
  cat("      -> Hypothesis REJECTED. Do NOT scale up. Report back and\n")
  cat("         we look elsewhere (likely the QSS algebra or the\n")
  cat("         simulation/fit structural mismatch).\n")
}
cat("----------------------------------------------------------------------\n")

saveRDS(fit_diag, "outputs/fit_diagnostic_v2.rds")
cat("\nSaved to outputs/fit_diagnostic_v2.rds\n")
