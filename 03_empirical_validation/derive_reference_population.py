# Script 1 of 2: derive_reference_population.py
#
# Builds the reference measurement population used to validate the
# bias-corrected standard deviation estimators.
#
# Input  : household_power_consumption.txt  (raw UCI dataset, see README)
# Output : reference_population.csv         (one column, detrended Voltage
#                                            residuals of a stationary pool)
#
# Data source:
#   UCI Machine Learning Repository, dataset 235,
#   "Individual Household Electric Power Consumption",
#   Hebrail G. and Berard A., DOI 10.24432/C58K54, license CC BY 4.0.
#

import numpy as np
from scipy import stats

RAW = "household_power_consumption.txt"
OUT = "reference_population.csv"
DAY = 1440  # minutes per day; sampling rate is one record per minute

# ----- load the Voltage channel (column index 4); missing values are '?' -----
volt = []
with open(RAW, "r") as f:
    f.readline()  # header
    for line in f:
        p = line.rstrip("\r\n").split(";")
        if len(p) < 5:
            continue
        volt.append(np.nan if p[4] in ("?", "") else float(p[4]))
v = np.array(volt, dtype=float)
v = v[~np.isnan(v)]
print("valid Voltage readings : {}".format(v.size))

# ----- remove slow daily and seasonal drift with a one-day moving average ----
def moving_average(x, w):
    c = np.cumsum(np.insert(x, 0, 0.0))
    m = (c[w:] - c[:-w]) / w
    pl = w // 2
    return np.concatenate([np.full(pl, m[0]), m, np.full(w - 1 - pl, m[-1])])

resid = v - moving_average(v, DAY)

# ----- select day-blocks with homogeneous variance (stationary pool) ---------
# Pooling days of differing variance inflates kurtosis through a variance
# mixture. Keeping same-variance blocks yields a near-normal stationary
# population whose population standard deviation is a well-defined ground
# truth for the estimator validation.
nb = resid.size // DAY
bstd = np.array([resid[i*DAY:(i+1)*DAY].std(ddof=0) for i in range(nb)])
med = np.median(bstd)
keep = np.abs(bstd - med) / med < 0.05
pool = np.concatenate([resid[i*DAY:(i+1)*DAY] for i in range(nb) if keep[i]])

print("stationary day-blocks kept : {} of {}".format(int(keep.sum()), nb))
print("reference population size  : {}".format(pool.size))
print("population sigma (ddof=0)  : {:.5f} V".format(pool.std(ddof=0)))
print("skewness                   : {:+.4f}".format(stats.skew(pool)))
print("excess kurtosis            : {:+.4f}".format(stats.kurtosis(pool)))

# ----- write the reference population ----------------------------------------
np.savetxt(OUT, pool, fmt="%.6f", header="voltage_residual_V", comments="")
print("written : {}".format(OUT))
