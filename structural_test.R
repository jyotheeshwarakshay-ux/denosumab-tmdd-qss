# =============================================================================
# DEFINITIVE STRUCTURAL TEST
#
# Question: do the SIMULATION model and the FIT model produce the same
# concentration-time curve when given IDENTICAL parameters?
#
# If YES -> the models match; the bias comes from somewhere else entirely.
# If NO  -> found it. The simulation and fit are different models, which
#           fully explains why the fit can never reproduce the data.
#
# This runs in seconds. No SAEM, no estimation -- pure forward simulation
# of both model objects at the true population parameters (all etas = 0).
# =============================================================================

library(rxode2)
library(nlmixr2)

# --- true population parameters, etas all zero -------------------------------
p_true <- c(
  lka   = log(0.0078),
  lVc   = log(1.58),
  lVp   = log(6.06),
  lCL   = log(0.006),
  lQ    = log(0.20),
  lkint = log(0.022),
  lKss  = log(1.56),
  lksyn = log(0.01),
  lR0   = log(15.23),
  eta_ka = 0, eta_Vc = 0, eta_Vp = 0, eta_CL = 0, eta_Q = 0,
  eta_kint = 0, eta_Kss = 0, eta_ksyn = 0, eta_R0 = 0
)
R0_true <- 15.23
dose_nmol <- 408.2
month_h   <- 30.44 * 24
dose_times_h <- c(0, 6, 12) * month_h

# --- MODEL A: exactly as used for SIMULATION (initial cond. passed externally)
model_sim <- rxode2({
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
})

# --- MODEL B: exactly as used for FITTING (initial cond. INSIDE the block) ---
model_fit <- rxode2({
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

  Rtot(0) <- R0
})

# --- identical event table for both -----------------------------------------
make_ev <- function() {
  ev <- eventTable()
  for (dt in dose_times_h) {
    ev$add.dosing(dose = dose_nmol, start.time = dt, dosing.to = "depot")
  }
  ev$add.sampling(seq(0, 18 * month_h, length.out = 400))
  ev
}

cat("Solving MODEL A (simulation style: inits passed to rxSolve)...\n")
solA <- rxSolve(model_sim, p_true, make_ev(),
                inits = c(depot = 0, Ctot = 0, Cp = 0, Rtot = R0_true))
dfA <- as.data.frame(solA)

cat("Solving MODEL B (fit style: Rtot(0) <- R0 inside model)...\n")
solB <- rxSolve(model_fit, p_true, make_ev())
dfB <- as.data.frame(solB)

# --- compare -----------------------------------------------------------------
cat("\n")
cat("======================================================================\n")
cat("COMPARISON AT IDENTICAL PARAMETERS\n")
cat("======================================================================\n\n")

cat(sprintf("%-22s %14s %14s\n", "quantity", "MODEL A (sim)", "MODEL B (fit)"))
cat(sprintf("%-22s %14.4f %14.4f\n", "Rtot at t=0",
            dfA$Rtot[1], dfB$Rtot[1]))
cat(sprintf("%-22s %14.4f %14.4f\n", "max Ctot",
            max(dfA$Ctot), max(dfB$Ctot)))
cat(sprintf("%-22s %14.1f %14.1f\n", "time of max Ctot (h)",
            dfA$time[which.max(dfA$Ctot)], dfB$time[which.max(dfB$Ctot)]))
cat(sprintf("%-22s %14.4f %14.4f\n", "mean Ctot",
            mean(dfA$Ctot), mean(dfB$Ctot)))
cat(sprintf("%-22s %14.4f %14.4f\n", "min Rtot",
            min(dfA$Rtot), min(dfB$Rtot)))

d <- dfA$Ctot - dfB$Ctot
cat(sprintf("\n  max |difference| in Ctot : %.6f nmol/L\n", max(abs(d))))
cat(sprintf("  RMSE between curves      : %.6f nmol/L\n", sqrt(mean(d^2))))

cat("\n----------------------------------------------------------------------\n")
if (max(abs(d)) < 1e-4) {
  cat("MODELS ARE IDENTICAL.\n")
  cat("  -> The simulation and fit models are the same. Structural mismatch\n")
  cat("     is RULED OUT. The bias must come from the data-handling layer\n")
  cat("     (Stage 4 formatting) or from the observation/dosing records.\n")
} else {
  cat("*** MODELS DIFFER ***\n")
  cat(sprintf("  -> Max discrepancy %.4f nmol/L. This is the bug.\n", max(abs(d))))
  cat("     The fit was estimating a DIFFERENT model than the one that\n")
  cat("     generated the data, so it could never recover the parameters.\n")
}
cat("----------------------------------------------------------------------\n")

# --- also check what the SIMULATED data actually looks like vs model A -------
if (exists("sim_all")) {
  cat("\nCross-check against the actual simulated dataset:\n")
  cat(sprintf("  sim_all$Ctot  max = %.3f   (population spread, includes IIV)\n",
              max(sim_all$Ctot, na.rm = TRUE)))
  cat(sprintf("  MODEL A typical max = %.3f  (no IIV, etas = 0)\n", max(dfA$Ctot)))
}
