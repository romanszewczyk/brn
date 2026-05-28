# Script 2 of 2: run_validation.py
#
# Empirical validation of the bias-corrected standard deviation estimators
# on real measurement data, using the empirical-population bootstrap.
#
# Input  : reference_population.csv  (produced by derive_reference_population.py)
# Output : bias_normal.csv, bias_truncated.csv,
#          fig_validation_normal.png, fig_validation_truncated.png
#
# Method
#   The reference population is treated as the population. Its population
#   standard deviation sigma_ref is the ground truth. Samples of size n are
#   drawn with replacement, every estimator is applied, and the relative bias
#   is measured against sigma_ref.
#
#   Part 1: normal-distribution estimator. Estimators of the paper:
#           Eq.1 biased, Eq.2 Bessel, Eq.3 the n-1.5 heuristic,
#           Eq.4 Gurland, Eq.8 the proposed LFF, and the proposed tabulated
#           correction which equals the exact c4 factor.
#   Part 2: truncated-normal estimator. Eq.11 uncorrected and Eq.13 proposed.
#           The factor b(n,r) is calibrated by the same data-driven
#           procedure used in the paper.
#

import numpy as np
from scipy import stats, optimize, special
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20240523)
POOL_CSV = "reference_population.csv"
RES_DIR = "."

# ----- load the reference population -----------------------------------------
pool = np.loadtxt(POOL_CSV, skiprows=1)
N = pool.size
SIGMA_REF = pool.std(ddof=0)
print("reference population N = {}, sigma_ref = {:.5f}".format(N, SIGMA_REF))
print("skewness = {:+.4f}, excess kurtosis = {:+.4f}".format(
    stats.skew(pool), stats.kurtosis(pool)))

# ----- estimator definitions -------------------------------------------------
def c4(n):
    return np.sqrt(2.0/(n-1.0)) * np.exp(special.gammaln(n/2.0)
                                         - special.gammaln((n-1.0)/2.0))

def b_lff(n):
    # proposed linear fractional function, paper Eq.8
    return (1.522*n - 1.6417) / (1.015*n - 1.0)

def normal_denominators(n):
    # name -> denominator d such that estimator = sqrt(SSD / d)
    return {
        "Biased_Eq1":       n,
        "Bessel_Eq2":       n - 1.0,
        "n_minus_1.5_Eq3":  n - 1.5,
        "Gurland_Eq4":      n - 1.5 + 1.0/(8.0*(n-1.0)),
        "Proposed_tab_c4":  (n-1.0)*c4(n)**2,
        "Proposed_LFF_Eq8": n - b_lff(n),
    }

EST = list(normal_denominators(10).keys())

def bootstrap_normal(n, M, chunk=20000):
    dens = normal_denominators(n)
    s1 = {k: 0.0 for k in EST}
    s2 = {k: 0.0 for k in EST}
    done = 0
    while done < M:
        m = min(chunk, M - done)
        x = pool[RNG.integers(0, N, size=(m, n))]
        ssd = np.sum((x - x.mean(axis=1, keepdims=True))**2, axis=1)
        for k in EST:
            s = np.sqrt(ssd / dens[k])
            s1[k] += s.sum()
            s2[k] += (s*s).sum()
        done += m
    out = {}
    for k in EST:
        mean = s1[k]/M
        var = max(s2[k]/M - mean*mean, 0.0)
        out[k] = (mean, var)
    return out

# ----- Part 1: normal-distribution estimator ---------------------------------
N_GRID = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100]
M_NORMAL = 400000

print()
print("PART 1: normal-distribution estimator (relative bias, percent)")
bias = {k: [] for k in EST}
mse = {k: [] for k in EST}
for n in N_GRID:
    res = bootstrap_normal(n, M_NORMAL)
    for k in EST:
        mean, var = res[k]
        bias[k].append(100.0*(mean - SIGMA_REF)/SIGMA_REF)
        mse[k].append(var + (mean - SIGMA_REF)**2)

with open(RES_DIR + "/bias_normal.csv", "w") as f:
    f.write("n," + ",".join(EST) + "," + ",".join(k+"_MSEratio" for k in EST) + "\n")
    for i, n in enumerate(N_GRID):
        base = mse["Bessel_Eq2"][i]
        f.write("{}".format(n))
        for k in EST:
            f.write(",{:.4f}".format(bias[k][i]))
        for k in EST:
            f.write(",{:.4f}".format(mse[k][i]/base))
        f.write("\n")
print("written bias_normal.csv")
for i, n in enumerate(N_GRID):
    print("  n={:>3} : ".format(n) + "  ".join(
        "{}={:+7.3f}".format(k, bias[k][i]) for k in EST))

# ----- Part 2: truncated-normal estimator ------------------------------------
def alpha_trunc(r):
    # variance scaling factor of a standard normal truncated to [-r, r], Eq.12
    return 1.0 - 2.0*r*stats.norm.pdf(r) / (2.0*stats.norm.cdf(r) - 1.0)

def truncated_normal_samples(n, M, r):
    lo, hi = stats.norm.cdf(-r), stats.norm.cdf(r)
    return stats.norm.ppf(RNG.uniform(lo, hi, size=(M, n)))

def calibrate_b(n, r, M_cal=300000):
    # paper's data-driven calibration: find b(n,r) so that
    # E[ sqrt(SSD / ((n-b)*alpha(r))) ] = 1 for a truncated-normal parent
    a = alpha_trunc(r)
    x = truncated_normal_samples(n, M_cal, r)
    ssd = np.sum((x - x.mean(axis=1, keepdims=True))**2, axis=1)
    def g(b):
        return np.mean(np.sqrt(ssd / ((n - b)*a))) - 1.0
    b = optimize.brentq(g, -5.0, n - 1e-6, xtol=1e-6)
    return b, g(b)

R_GRID = [0.3, 0.5, 1.0, 2.0]
NT_GRID = [2, 3, 5, 10, 20, 50, 100]
M_TRUNC = 400000
mu = pool.mean()

print()
print("PART 2: truncated-normal estimator (relative bias, percent)")
trunc = {}
with open(RES_DIR + "/bias_truncated.csv", "w") as f:
    f.write("r,n,b_nr,alpha_theory,alpha_empirical,"
            "uncorrected_bias_pct,proposed_bias_pct\n")
    for r in R_GRID:
        a = alpha_trunc(r)
        sub = pool[np.abs(pool - mu) <= r*SIGMA_REF]
        a_emp = (sub.std(ddof=0)/SIGMA_REF)**2
        rows = []
        for n in NT_GRID:
            b, _ = calibrate_b(n, r)
            s1u = s1p = 0.0
            done = 0
            while done < M_TRUNC:
                m = min(20000, M_TRUNC - done)
                x = sub[RNG.integers(0, sub.size, size=(m, n))]
                ssd = np.sum((x - x.mean(axis=1, keepdims=True))**2, axis=1)
                s1u += np.sqrt(ssd/((n-1.0)*a)).sum()
                s1p += np.sqrt(ssd/((n-b)*a)).sum()
                done += m
            bu = 100.0*(s1u/M_TRUNC - SIGMA_REF)/SIGMA_REF
            bp = 100.0*(s1p/M_TRUNC - SIGMA_REF)/SIGMA_REF
            rows.append((n, b, bu, bp))
            f.write("{},{},{:.4f},{:.5f},{:.5f},{:.4f},{:.4f}\n".format(
                r, n, b, a, a_emp, bu, bp))
        trunc[r] = rows
        print("  r={}: alpha_theory={:.5f} alpha_emp={:.5f}".format(r, a, a_emp))
        for n, b, bu, bp in rows:
            print("    n={:>3} b={:.3f} uncorrected={:+7.3f} proposed={:+7.3f}".format(
                n, b, bu, bp))
print("written bias_truncated.csv")

# ----- figures ---------------------------------------------------------------
styles = {
    "Biased_Eq1":       ("o-", "#c0392b", "Biased (Eq.1)"),
    "Bessel_Eq2":       ("s-", "#e67e22", "Bessel (Eq.2)"),
    "n_minus_1.5_Eq3":  ("^-", "#8e44ad", "n-1.5 heuristic (Eq.3)"),
    "Gurland_Eq4":      ("d-", "#16a085", "Gurland (Eq.4)"),
    "Proposed_tab_c4":  ("v-", "#2c3e50", "Proposed tabulated (= c4)"),
    "Proposed_LFF_Eq8": ("*-", "#2980b9", "Proposed LFF (Eq.8)"),
}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for k in EST:
    st, col, lab = styles[k]
    ax[0].plot(N_GRID, bias[k], st, color=col, label=lab, markersize=5)
ax[0].axhline(0, color="k", lw=0.8)
ax[0].set_xscale("log")
ax[0].set_xlabel("sample size n")
ax[0].set_ylabel("relative bias [%]")
ax[0].set_title("Normal case: relative bias")
ax[0].legend(fontsize=7)
ax[0].grid(alpha=0.3)
for k in EST:
    st, col, lab = styles[k]
    ax[1].plot(N_GRID, np.abs(bias[k]) + 1e-4, st, color=col, markersize=5)
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_xlabel("sample size n")
ax[1].set_ylabel("absolute relative bias [%]")
ax[1].set_title("Normal case: magnitude of bias")
ax[1].grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(RES_DIR + "/fig_validation_normal.png", dpi=150)

fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
for j, r in enumerate(R_GRID):
    rows = trunc[r]
    ns = [x[0] for x in rows]
    axes[j].plot(ns, [x[2] for x in rows], "s-", color="#c0392b",
                 label="uncorrected (Eq.11)")
    axes[j].plot(ns, [x[3] for x in rows], "o-", color="#2980b9",
                 label="proposed (Eq.13)")
    axes[j].axhline(0, color="k", lw=0.8)
    axes[j].set_xscale("log")
    axes[j].set_xlabel("sample size n")
    axes[j].set_title("truncation r = {}".format(r))
    axes[j].grid(alpha=0.3)
    if j == 0:
        axes[j].set_ylabel("relative bias [%]")
        axes[j].legend(fontsize=8)
fig.suptitle("Truncated-normal case: relative bias on real data", y=1.04)
fig.tight_layout()
fig.savefig(RES_DIR + "/fig_validation_truncated.png", dpi=150,
            bbox_inches="tight")
print("written fig_validation_normal.png and fig_validation_truncated.png")
print("DONE")
