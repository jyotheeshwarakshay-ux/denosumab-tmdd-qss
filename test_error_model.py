"""
Before rewriting the R script, verify the diagnosis of the error-model problem.

DIAGNOSIS TO TEST:
  Original code did:  DV = IPRED*(1+eps_prop) + eps_add,  then DV = max(DV, 0)
  At troughs IPRED ~ 0, so DV ~ eps_add ~ N(0, 0.72).
  Clipping at 0 deletes the negative half  -> upward bias.

  Expected consequence: mean(DV - IPRED) > 0, i.e. residuals systematically
  positive, which is what the fit reported (IWRES mean = 0.58).

FIX TO TEST:
  Apply an LLOQ (like the paper's 20 ng/mL) and DROP sub-LLOQ points as
  missing, instead of clipping them to zero.
"""

import numpy as np

rng = np.random.default_rng(42)

sigma_add = 0.72
sigma_prop = 0.07

# Simulate a realistic spread of true concentrations across a dosing interval:
# many near-zero troughs, some mid, some peaks -- matching the observed
# quartiles from the actual run (median ~12, min 0, max ~120).
ipred = np.concatenate([
    rng.uniform(0.0, 0.5, 600),    # deep troughs
    rng.uniform(0.5, 5.0, 400),    # low
    rng.uniform(5.0, 30.0, 500),   # mid
    rng.uniform(30.0, 120.0, 180), # peaks
])

eps_add = rng.normal(0, sigma_add, ipred.size)
eps_prop = rng.normal(0, sigma_prop, ipred.size)
dv_raw = ipred * (1 + eps_prop) + eps_add

print("=== ORIGINAL APPROACH: clip negatives to zero ===")
dv_clipped = np.maximum(dv_raw, 0.0)
resid_clipped = dv_clipped - ipred
print(f"  mean residual : {resid_clipped.mean():+.4f}   (should be ~0)")
print(f"  frac of points that got clipped: {(dv_raw < 0).mean():.1%}")
print(f"  mean residual among troughs (ipred<0.5): "
      f"{(dv_clipped[:600] - ipred[:600]).mean():+.4f}")

print()
print("=== FIX A: LLOQ, drop sub-LLOQ as missing ===")
# Paper's LLOQ was 20 ng/mL. Convert to nmol/L for denosumab (147 kDa):
#   20 ng/mL = 20 ug/L; nmol/L = (ug/L * 1000) / MW_kDa / 1000 * ... let's be careful:
#   conc_nmol_per_L = (conc_ng_per_mL * 1e-9 g/mL * 1e3 mL/L) / (147000 g/mol) * 1e9 nmol/mol
MW = 147000.0  # g/mol
lloq_ng_per_mL = 20.0
lloq_nmol_L = (lloq_ng_per_mL * 1e-6) / MW * 1e9  # ng/mL -> g/L -> mol/L -> nmol/L
print(f"  LLOQ = {lloq_ng_per_mL} ng/mL = {lloq_nmol_L:.4f} nmol/L")

keep = dv_raw >= lloq_nmol_L
dv_lloq = dv_raw[keep]
ipred_lloq = ipred[keep]
resid_lloq = dv_lloq - ipred_lloq
print(f"  frac retained : {keep.mean():.1%}  (paper retained ~73.5%)")
print(f"  mean residual : {resid_lloq.mean():+.4f}   (should be ~0)")

print()
print("=== FIX B: smaller additive error, no clipping, no LLOQ ===")
# Alternative: the additive term may simply be too large relative to the
# trough concentrations. Try proportional-dominant error.
sigma_add_small = 0.05
dv_b = ipred * (1 + rng.normal(0, sigma_prop, ipred.size)) + \
       rng.normal(0, sigma_add_small, ipred.size)
resid_b = dv_b - ipred
print(f"  sigma_add reduced to {sigma_add_small}")
print(f"  frac negative : {(dv_b < 0).mean():.2%}")
print(f"  mean residual : {resid_b.mean():+.4f}   (should be ~0)")
