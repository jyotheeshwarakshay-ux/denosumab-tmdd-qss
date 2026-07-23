# =============================================================================
# DATA-HANDLING TEST  (v2 -- fixed)
#
# Previous version failed: nlmixr2 refuses est="posthoc" when every parameter
# is fixed ("no parameters to estimate"). My error.
#
# This version uses rxSolve directly with pk_data as the event table, which
# tests the same question without needing an estimator:
#
#   Does the pk_data dosing record (AMT=408.2, EVID=1, CMT=1) deliver drug
#   to the depot compartment the same way ev$add.dosing(dosing.to="depot")
#   did during simulation?
#
# Expected typical max Ctot (verified by the structural test): ~40.2 nmol/L
# =============================================================================

library(rxode2)
library(dplyr)

if (!exists("pk_data")) stop("pk_data not found -- run Stages 1-4 first.")

# --- the model (identical to both sim and fit; already verified equivalent) --
m <- rxode2({
  ka    <- exp(lka)
  Vc    <- exp(lVc)
  Vp    <- exp(lVp)
  CL    <- exp(lCL)
  Q     <- exp(lQ)
  kint  <- exp(lkint)
  Kss   <- exp(lKss)
  ksyn  <- exp(lksyn)
  R0    <- exp(lR0)
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

p_true <- c(lka = log(0.0078), lVc = log(1.58), lVp = log(6.06),
            lCL = log(0.006),  lQ  = log(0.20), lkint = log(0.022),
            lKss = log(1.56),  lksyn = log(0.01), lR0 = log(15.23))

# --- ROUTE 1: solve using pk_data as the event table (the FIT pathway) -------
# Take subject 1 only, in the exact format nlmixr2 received.
d1 <- pk_data %>% filter(ID == 1) %>%
  transmute(id = ID, time = TIME, amt = AMT, evid = EVID, cmt = CMT, dv = DV)

cat("Subject 1 records passed to the fit (first 6 rows):\n")
print(head(d1, 6))
cat(sprintf("\n  dosing rows: %d   observation rows: %d\n",
            sum(d1$evid == 1), sum(d1$evid == 0)))
cat(sprintf("  dose amount on dosing rows: %.1f\n\n",
            unique(d1$amt[d1$evid == 1])))

cat("Solving via pk_data event table (FIT pathway)...\n")
sol_data <- rxSolve(m, p_true, d1,
                    inits = c(depot = 0, Ctot = 0, Cp = 0, Rtot = 15.23))
df_data <- as.data.frame(sol_data)

# --- ROUTE 2: solve using eventTable (the SIMULATION pathway) ---------------
month_h <- 30.44 * 24
dose_times_h <- c(0, 6, 12) * month_h
sample_offsets_h <- c(12, 24, 168, 720, 2160, 3600, 4380)

ev <- eventTable()
for (dt in dose_times_h) {
  ev$add.dosing(dose = 408.2, start.time = dt, dosing.to = "depot")
}
st <- sort(unique(as.vector(outer(dose_times_h, sample_offsets_h, "+"))))
st <- st[st <= max(dose_times_h) + max(sample_offsets_h)]
ev$add.sampling(st)

cat("Solving via eventTable (SIMULATION pathway)...\n\n")
sol_ev <- rxSolve(m, p_true, ev,
                  inits = c(depot = 0, Ctot = 0, Cp = 0, Rtot = 15.23))
df_ev <- as.data.frame(sol_ev)

# --- compare ----------------------------------------------------------------
cat("======================================================================\n")
cat("COMPARISON: pk_data pathway vs eventTable pathway\n")
cat("======================================================================\n\n")

cat(sprintf("%-26s %14s %14s\n", "quantity", "pk_data", "eventTable"))
cat(sprintf("%-26s %14d %14d\n", "n output rows",
            nrow(df_data), nrow(df_ev)))
cat(sprintf("%-26s %14.4f %14.4f\n", "max Ctot",
            max(df_data$Ctot, na.rm = TRUE), max(df_ev$Ctot, na.rm = TRUE)))
cat(sprintf("%-26s %14.4f %14.4f\n", "mean Ctot",
            mean(df_data$Ctot, na.rm = TRUE), mean(df_ev$Ctot, na.rm = TRUE)))
cat(sprintf("%-26s %14.4f %14.4f\n", "max depot",
            max(df_data$depot, na.rm = TRUE), max(df_ev$depot, na.rm = TRUE)))
cat(sprintf("%-26s %14.4f %14.4f\n", "Rtot at first row",
            df_data$Rtot[1], df_ev$Rtot[1]))

cat("\n  Expected typical max Ctot (from verified structural test): ~40.2\n")

cat("\n----------------------------------------------------------------------\n")
mx_data <- max(df_data$Ctot, na.rm = TRUE)
mx_ev   <- max(df_ev$Ctot, na.rm = TRUE)

if (max(df_data$depot, na.rm = TRUE) < 1) {
  cat("*** DOSING NOT REACHING DEPOT VIA pk_data ***\n")
  cat("  max depot is ~0 -> the CMT=1 records are not delivering drug.\n")
  cat("  THIS IS THE BUG.\n")
} else if (abs(mx_data - mx_ev) / mx_ev < 0.05) {
  cat("BOTH PATHWAYS AGREE.\n")
  cat(sprintf("  pk_data max %.2f vs eventTable max %.2f (within 5%%).\n",
              mx_data, mx_ev))
  cat("  Data handling is correct. The bug is in the ESTIMATION step:\n")
  cat("  most likely the residual error model specification.\n")
} else {
  cat("*** PATHWAYS DISAGREE ***\n")
  cat(sprintf("  pk_data max %.2f vs eventTable max %.2f (%.0f%% apart).\n",
              mx_data, mx_ev, 100 * (mx_data - mx_ev) / mx_ev))
  cat("  The dataset is not reproducing the simulation. THIS IS THE BUG.\n")
}
cat("----------------------------------------------------------------------\n")
