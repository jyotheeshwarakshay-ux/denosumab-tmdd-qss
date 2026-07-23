"""
Sanity check for the QSS-TMDD model before building the full nlmixr2 pipeline.
Goal: confirm that simulating with the Table 3 parameters (Choi et al. 2025,
Front Pharmacol) produces a concentration-time profile with Cmax/Tmax in the
same ballpark as the paper's own reported simulation results (Table 5):

    Cmax (nmol/L): SB16 40.66, DEN 40.29 (median, 3x 60mg SC q6mo dosing)
    Tmax (h): SB16 259, DEN 267

If this sanity check doesn't land near those numbers, something is wrong with
my equations or units BEFORE any R/nlmixr2 time is spent chasing it.
"""

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Table 3 parameters (PMO / postmenopausal osteoporosis population, typical
# values -- no IIV, no residual error -- this is the population-mean profile)
# ---------------------------------------------------------------------------
ka    = 0.0078      # 1/h, absorption rate constant (PMO)
Vc    = 1.58         # L, central volume (population-typical, i.e. at median weight)
Vp    = 6.06         # L, peripheral volume
CL    = 0.006        # L/h, clearance (Caucasian, population-typical)
Q     = 0.20          # L/h, inter-compartmental clearance (PMO)
kint  = 0.022        # 1/h, complex internalization rate
Kss   = 1.56          # nmol/L, QSS constant
ksyn  = 0.01          # 1/h, target synthesis rate
R0    = 15.23         # nmol/L, baseline target (PMO)
kdeg  = ksyn / R0      # 1/h, derived (not directly reported -- see chat notes)

F     = 1.0            # bioavailability not separately reported in Table 3;
                        # CL/F, Vc/F etc are apparent (F folded in) -- use F=1
                        # and treat all volumes/CL as "apparent" per the paper

dose_mg = 60           # mg, subcutaneous
MW_kDa  = 147           # denosumab MW ~ 147 kDa (IgG2 mAb) -- for unit conversion
dose_nmol = (dose_mg * 1e6) / (MW_kDa * 1e3)  # mg -> ng -> nmol; check magnitude

print(f"Dose in nmol: {dose_nmol:.1f} nmol (sanity: should be roughly 400 nmol for 60mg IgG)")

# ---------------------------------------------------------------------------
# QSS-TMDD ODEs
# State vector: [Adepot, Ctot, Cp, Rtot]
#   Adepot = amount in depot (nmol)
#   Ctot   = total drug conc in central compartment (free + bound), nmol/L
#   Cp     = free drug conc in peripheral compartment, nmol/L
#   Rtot   = total target conc (free + bound), nmol/L
# ---------------------------------------------------------------------------
def qss_rhs(t, y):
    Adepot, Ctot, Cp, Rtot = y

    # QSS closed-form solution for free drug concentration C
    # C = 1/2 * [ (Ctot - Rtot - Kss) + sqrt((Ctot - Rtot - Kss)^2 + 4*Kss*Ctot) ]
    disc = (Ctot - Rtot - Kss)**2 + 4 * Kss * Ctot
    disc = max(disc, 0.0)  # numerical guard
    C = 0.5 * ((Ctot - Rtot - Kss) + np.sqrt(disc))
    C = max(C, 0.0)

    RC = Rtot - (Ctot - C)  # bound complex = Rtot - free target;
    # equivalently RC = Ctot - C (bound drug = total - free), consistent check below
    RC_check = Ctot - C
    # (RC and RC_check should match; both represent bound complex conc)

    dAdepot = -ka * Adepot
    dCtot   = (ka * Adepot) / Vc - (CL / Vc) * C - (Q / Vc) * C + (Q / Vp) * Cp - kint * RC_check
    dCp     = (Q / Vc) * C - (Q / Vp) * Cp   # FIXED: flux out of central scales by Vc, not Vp
    dRtot   = ksyn - kdeg * (Rtot - RC_check) - kint * RC_check

    return [dAdepot, dCtot, dCp, dRtot]

# ---------------------------------------------------------------------------
# Dosing: 60 mg SC at t=0, t=6 months, t=12 months (matching paper's regimen)
# ---------------------------------------------------------------------------
month_h = 30.44 * 24  # approx hours per month
dose_times = [0, 6 * month_h, 12 * month_h]
t_end = 18 * month_h  # paper simulates out to 18 months

y0 = [0.0, 0.0, 0.0, R0]  # baseline: no drug, target at steady state

all_t = []
all_Ctot = []

state = y0
t_prev = 0
for i, dt in enumerate(dose_times):
    if dt > t_prev:
        sol = solve_ivp(qss_rhs, [t_prev, dt], state, method='LSODA',
                         max_step=1.0, rtol=1e-8, atol=1e-10, dense_output=True)
        t_eval = np.linspace(t_prev, dt, 200)
        y_eval = sol.sol(t_eval)
        all_t.extend(t_eval)
        all_Ctot.extend(y_eval[1])
        state = sol.y[:, -1]
    # add dose to depot
    state = list(state)
    state[0] += dose_nmol
    t_prev = dt

# final segment to t_end
sol = solve_ivp(qss_rhs, [t_prev, t_end], state, method='LSODA',
                 max_step=1.0, rtol=1e-8, atol=1e-10, dense_output=True)
t_eval = np.linspace(t_prev, t_end, 400)
y_eval = sol.sol(t_eval)
all_t.extend(t_eval)
all_Ctot.extend(y_eval[1])

all_t = np.array(all_t)
all_Ctot = np.array(all_Ctot)

# Look at first dosing interval only (0 to 6 months) for Cmax/Tmax comparison
mask = all_t <= 6 * month_h
t1 = all_t[mask]
c1 = all_Ctot[mask]
cmax_idx = np.argmax(c1)

print(f"\n--- First-dose interval (0-6mo) ---")
print(f"Simulated Cmax: {c1[cmax_idx]:.2f} nmol/L  (paper reports ~40 nmol/L median)")
print(f"Simulated Tmax: {t1[cmax_idx]:.1f} h        (paper reports ~260 h median)")
print(f"Concentration at t=6mo (trough before 2nd dose): {c1[-1]:.3f} nmol/L")
