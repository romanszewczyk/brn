# 03_empirical_validation -- Real measurement data

Empirical validation on the UCI Individual Household Electric Power
Consumption dataset.  The scripts here implement the bootstrap protocol
reported in Section 5 of the manuscript.

## Data source

  Name     : Individual Household Electric Power Consumption
  Provider : UCI Machine Learning Repository, dataset 235
  Authors  : Hebrail G. and Berard A.
  DOI      : 10.24432/C58K54
  URL      : https://archive.ics.uci.edu/dataset/235
  License  : Creative Commons Attribution 4.0 International

The raw file `household_power_consumption.txt` (about 127 MB) is NOT
included in the repository.  Download it from the UCI link above, unzip if
necessary, and place the `.txt` file in this directory.

## Contents

  derive_reference_population.py    raw dataset -> reference_population.csv
  run_validation.py                 reference_population.csv -> bias tables
  bias_normal.csv                   archived normal-case results
  bias_truncated.csv                archived truncated-normal results

## How the reference population was derived

The Voltage channel (column 4) is loaded and detrended with a one-day
moving average to remove slow daily and seasonal drift.  The detrended
series is split into one-day blocks, and only blocks whose standard
deviation lies within five percent of the median block standard deviation
are retained and pooled.  This selection removes the variance-mixture
effect that would otherwise inflate kurtosis, and yields a stationary,
near-normal population with a well-defined population standard deviation.

Resulting reference population:

  size N                       240,480 points
  population standard deviation 2.30193 V  (ground truth, divisor N)
  skewness                      -0.315
  excess kurtosis               -0.021
  relative uncertainty of sigma 0.144 percent  (equal to 1/sqrt(2N))

## Running the validation

Requirements: Python 3, numpy, scipy, matplotlib.

Reproduction from the raw source:

  # place household_power_consumption.txt in this directory
  python derive_reference_population.py    # writes reference_population.csv
  python run_validation.py                 # reads reference_population.csv

The raw dataset (~127 MB) is not redistributed here, so the reference
population is rebuilt by the first script; run_validation.py then consumes the
reference_population.csv it produces.  Both scripts use a fixed random seed, so
the numerical results are exactly reproducible.  The validation draws 400,000
bootstrap samples per sample size.  The archived bias_normal.csv and
bias_truncated.csv hold the exact tables reported in the paper.

## Validation method

The reference population is treated as the statistical population.  Its
population standard deviation is the ground truth.  Samples of size n are
drawn with replacement, every estimator is applied, and the relative bias
is measured against the ground truth.

Normal-distribution estimators compared (paper equation numbers): the
biased estimator (Eq. 1), the Bessel estimator (Eq. 2), the n-1.5 heuristic
(Eq. 3), the Gurland estimator (Eq. 4), the proposed linear fractional
function (Eq. 8), and the proposed tabulated correction, which equals the
exact c4 factor.

Truncated-normal estimators compared: the uncorrected estimator (Eq. 11
with the variance scaling factor of Eq. 12) and the proposed estimator
(Eq. 13).  The correction factor b(n,r) is calibrated by the same data-driven
procedure used in the paper, since the manuscript reports the coefficient
functions only graphically.

## Result files

bias_normal.csv columns: n, then the relative bias in percent for each of
the six estimators, then the MSE ratio of each estimator relative to the
Bessel estimator.

bias_truncated.csv columns: r, n, the calibrated b(n,r), the theoretical
and empirical variance scaling factors, and the relative bias in percent
of the uncorrected and the proposed estimator.
