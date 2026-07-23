"""
HYPOTHESIS TO TEST
------------------
The fit underpredicted systematically (mean IRES +5.0, IPRED max 70.6 vs DV max
120.8) and Vc came back at 3.43 vs true 1.58. My hypothesis: the sampling design
is too sparse in the absorption/distribution phase to identify Vc, Q, Vp
separately, so they trade off and land in a wrong-but-self-consistent place.

TEST
----
1. Simulate the TRUE model at the CURRENT sampling times.
2. Simulate a DELIBERATELY WRONG parameter set (Vc inflated ~2.2x, Q inflated,
   Vp inflated -- matching what the fit actually returned).
3. Compare the two curves AT THE SAMPLED TIMEPOINTS ONLY.

If the two are nearly indistinguishable at those timepoints, the data genuinely
cannot tell them apart -> hypothesis CONFIRMED, sampling is the problem.
If they are clearly different, the data DOES contain the information and the
problem lies elsewhere (estimator, initial values, model coding) -> hypothesis
REJECTED, and re-running with denser sampling would NOT help.
"""

import numpy as np
from scipy.integrate import solve_ivp

month_h = 30.44 * 24
dose_nmol = 408.2
dose_times = [0, 6 * month_h, 12 * month_h]

# ---- TRUE parameters (Choi et al. Table 3, PMO) ----
TRUE = dict(ka=0.0078, Vc=1.58, Vp=6.06, CL=0.006, Q=0.20,
            kint=0.022, Kss=1.56, ksyn=0.01, R0=15.23)

# ---- What the fit actually returned ----
FITTED = dict(ka=0.01203, Vc=3.4344, Vp=7.3206, CL=0.007928, Q=1.1408,
              kint=0.023621, Kss=1.92261, ksyn=0.0083029, R0=11.328)


def simulate(p, t_eval_max, n_grid=4000):
    kdeg = p['ksyn'] / p['R0']

    def rhs(t, y):
        Ad, Ctot, Cp, Rtot = y
        disc = (Ctot - Rtot - p['Kss'])**2 + 4 * p['Kss'] * Ctot
        disc = max(disc, 0.0)
        C = 0.5 * ((Ctot - Rtot - p['Kss']) + np.sqrt(disc))
        C = max(C, 0.0)
        RC = Ctot - C
        dAd = -p['ka'] * Ad
        dCtot = (p['ka'] * Ad) / p['Vc'] - (p['CL'] / p['Vc']) * C \
                - (p['Q'] / p['Vc']) * C + (p['Q'] / p['Vp']) * Cp - p['kint'] * RC
        dCp = (p['Q'] / p['Vc']) * C - (p['Q'] / p['Vp']) * Cp
        dRtot = p['ksyn'] - kdeg * (Rtot - RC) - p['kint'] * RC
        return [dAd, dCtot, dCp, dRtot]

    state = [0.0, 0.0, 0.0, p['R0']]
    ts, cs = [], []
    t_prev = 0.0
    for dt in dose_times:
        if dt > t_prev:
            sol = solve_ivp(rhs, [t_prev, dt], state, method='LSODA',
                            rtol=1e-9, atol=1e-11, dense_output=True)
            te = np.linspace(t_prev, dt, n_grid // 3)
            ts.extend(te); cs.extend(sol.sol(te)[1])
            state = sol.y[:, -1]
        state = list(state); state[0] += dose_nmol
        t_prev = dt
    sol = solve_ivp(rhs, [t_prev, t_eval_max], state, method='LSODA',
                    rtol=1e-9, atol=1e-11, dense_output=True)
    te = np.linspace(t_prev, t_eval_max, n_grid // 3)
    ts.extend(te); cs.extend(sol.sol(te)[1])
    return np.array(ts), np.array(cs)


def interp_at(ts, cs, targets):
    return np.interp(targets, ts, cs)


# ---- CURRENT sampling design ----
current_offsets = np.array([12, 24, 168, 720, 2160, 3600, 4380])
current_times = np.unique(np.concatenate(
    [np.array(dose_times)[:, None] + current_offsets[None, :]]).ravel())
t_max = max(dose_times) + current_offsets.max()
current_times = current_times[current_times <= t_max]

ts_true, cs_true = simulate(TRUE, t_max)
ts_fit, cs_fit = simulate(FITTED, t_max)

c_true_cur = interp_at(ts_true, cs_true, current_times)
c_fit_cur = interp_at(ts_fit, cs_fit, current_times)

print("=" * 72)
print("CURRENT SAMPLING DESIGN  (7 offsets per dose interval)")
print("=" * 72)
print(f"{'time(h)':>9} {'TRUE':>10} {'FITTED':>10} {'diff':>10} {'%diff':>9}")
for t, a, b in zip(current_times, c_true_cur, c_fit_cur):
    pct = 100 * (b - a) / a if a > 1e-6 else np.nan
    print(f"{t:9.0f} {a:10.3f} {b:10.3f} {b-a:10.3f} {pct:8.1f}%")

resid = c_fit_cur - c_true_cur
rmse_cur = np.sqrt(np.mean(resid**2))
print(f"\n  RMSE between TRUE and FITTED at sampled times: {rmse_cur:.3f} nmol/L")
print(f"  Mean |%diff| (where true>1): "
      f"{np.nanmean(np.abs(100*(c_fit_cur-c_true_cur)/np.where(c_true_cur>1, c_true_cur, np.nan))):.1f}%")

# ---- PROPOSED denser sampling design ----
# Add points through absorption/distribution: Tmax is ~300h, so sample around it.
proposed_offsets = np.array([12, 24, 72, 168, 336, 504, 720, 1080,
                              1440, 2160, 2880, 3600, 4380])
proposed_times = np.unique(np.concatenate(
    [np.array(dose_times)[:, None] + proposed_offsets[None, :]]).ravel())
t_max2 = max(dose_times) + proposed_offsets.max()
proposed_times = proposed_times[proposed_times <= t_max2]

ts_true2, cs_true2 = simulate(TRUE, t_max2)
ts_fit2, cs_fit2 = simulate(FITTED, t_max2)
c_true_prop = interp_at(ts_true2, cs_true2, proposed_times)
c_fit_prop = interp_at(ts_fit2, cs_fit2, proposed_times)

resid2 = c_fit_prop - c_true_prop
rmse_prop = np.sqrt(np.mean(resid2**2))

print()
print("=" * 72)
print("PROPOSED SAMPLING DESIGN  (13 offsets per dose interval)")
print("=" * 72)
print(f"  RMSE between TRUE and FITTED at sampled times: {rmse_prop:.3f} nmol/L")
print(f"  n timepoints: current {len(current_times)}  ->  proposed {len(proposed_times)}")

print()
print("=" * 72)
print("VERDICT")
print("=" * 72)
if rmse_cur < 2.0:
    print("  Curves are NEARLY IDENTICAL at current sampling times.")
    print("  -> Data cannot distinguish TRUE from FITTED. HYPOTHESIS CONFIRMED.")
    print("  -> Denser sampling should help.")
else:
    print(f"  Curves DIFFER substantially at current sampling times "
          f"(RMSE {rmse_cur:.2f}).")
    print("  -> The data DOES contain information to reject the fitted values.")
    print("  -> HYPOTHESIS REJECTED. Sparse sampling is NOT the main cause.")
    print("  -> Denser sampling alone would likely NOT fix this.")
