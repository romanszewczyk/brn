# Bias-Corrected Standard Deviation Estimators for Normal and Truncated-Normal Data

MIT License. See [LICENSE](LICENSE).

This repository holds the Monte Carlo generators, identification scripts, validation suites, and tabulated results used to develop closed-form bias corrections for the sample standard deviation under normal and symmetrically truncated-normal sampling.

Two estimators are provided:

- **Normal distribution** -- A linear fractional correction `b(n)` obtained from `k = 2 x 10^8` Monte Carlo trials per sample size, validated against the exact chi-distribution `c4` benchmark.
- **Truncated normal** -- A second-degree rational correction `b(n,r)` obtained from `k = 2 x 10^8` trials per `(n,r)` configuration, covering truncation levels `r = 0.2..10.0` and sample sizes `n = 2..100`.

Both corrections require only elementary arithmetic and are suitable for embedded controllers and PLC-based statistical process control systems.

## Repository Layout

```
.
|-- LICENSE                               # MIT License
|-- README.md                             # This file
|-- manuscript.pdf                        # Associated publication
|
|-- 01_normal_bn/                         # Standard normal N(0,1)
|   |-- 01_monte_carlo/                   # Authoritative b(n) Monte Carlo (500 reps, k=2e8)
|   |-- 02_identification/                # Rational-model fit & tabulated lookup
|   |-- 03_first_principle_check/         # Validation vs. exact c4 benchmark
|   |-- 04_symbolic_regression/           # PySR experiments
|   `-- requirements.txt
|
|-- 02_truncated_normal_bnr/              # Symmetrically truncated normal
|   |-- 01_monte_carlo/                   # b(n,r) Monte Carlo (25 reps, k=2e8, inverse CDF)
|   |-- 02_surface_identification/        # differential-evolution fit of the rational surface
|   |-- 03_lookup_and_scripts/            # Lookup function std_bnr_opti.m
|   |-- 04_first_principle_check/         # Small-scale (Octave) + Large-scale (Python/JAX) validation
|   `-- requirements.txt
|
`-- 03_empirical_validation/              # Empirical validation on real measurement data (UCI household power)
    |-- derive_reference_population.py
    |-- run_validation.py
    |-- bias_normal.csv
    |-- bias_truncated.csv
    `-- README.md
```

## Running the Code

### Standard Normal `b(n)` -- Validate against exact theory

```bash
cd 01_normal_bn/03_first_principle_check
octave validate_vs_c4.m
```

Loads the pre-computed `optimal_bn_500times.mat` and compares the Monte Carlo `b(n)` against the exact chi-distribution `c4` benchmark. Deviations are below `10^-4`.

### Truncated Normal `b(r,n)` -- Run a first-principle check

```bash
cd 02_truncated_normal_bnr/04_first_principle_check/Large_scale
pip install -r requirements.txt
python test_installation.py
python monte_carlo_validation.py
python analyze_results.py
```

### Symbolic Regression -- Discover the LFF form

```bash
cd 01_normal_bn/04_symbolic_regression
pip install -r requirements.txt
python pysr_symbolic_regression.py   # Experiment A -- general search
python pysr_lff_search.py            # Experiment B -- LFF-constrained search
```

Full PySR runs take 30-120 minutes depending on hardware.

## Technology Stack

| Component | Technology | Notes |
|-----------|-----------|-------|
| Primary numeric | MATLAB / GNU Octave | `.m` scripts; Octave compatibility maintained (`pkg load optim`) |
| Accelerated compute | Python + JAX (CUDA) | 30-70x speed-up for Monte Carlo; `jax_enable_x64=true` |
| Symbolic regression | Python + PySR | Julia back-end; pinned to v1.5.8 in `requirements.txt` |
| Data exchange | `.mat` (v7 format) | Read/write via `scipy.io` |

## Hardware & Runtime Notes

| Script | Hardware | Approx. Runtime |
|--------|----------|-----------------|
| `ident_bn_batched.py` (full) | GPU (CUDA), 64 GB RAM | ~9-10 hours |
| `ident_brn.py` (full) | GPU (CUDA), 16 GB RAM | ~100-300 hours |
| `validate_vs_c4.m` | Any CPU | < 1 minute |
| `monte_carlo_validation.py` (k=1e6) | GPU recommended | ~1 minute |
| `params_ident_brn.m` | CPU, Octave + optim pkg | ~2-6 hours |

Pre-computed `.mat` result files are included, so validation and plotting scripts run without the long Monte Carlo generation.

## Dependencies

### GNU Octave
- Octave 6.0+ with the `optim` package (`pkg load optim`)
- Required for: `.m` identification scripts and small-scale first-principle checks

### Python
Each module has its own `requirements.txt`:
- `01_normal_bn/requirements.txt` -- JAX, NumPy, SciPy, Matplotlib
- `02_truncated_normal_bnr/04_first_principle_check/Large_scale/requirements.txt` -- JAX, NumPy, SciPy, Matplotlib, Seaborn
- `01_normal_bn/04_symbolic_regression/requirements.txt` -- PySR v1.5.8 + above

Install JAX with CUDA 12 support:
```bash
pip install jax[cuda12]
```

## Estimators

### Normal Distribution

Compact linear fractional correction:

```
b_LFF(n) = (1.522 n - 1.6417) / (1.015 n - 1)
s = sqrt( sum(xi - xbar)^2 / (n - b_LFF(n)) )
```

Residual bias < 1.5e-3 at n=2, falling below 1e-4 for n >= 9.

### Truncated Normal Distribution

Second-degree rational correction for truncation level `r`:

```
          P(r) n^2 + A(r) n + B(r)
b(n,r) = ---------------------------
          Q(r) n^2 + C(r) n + 1

sr(n) = sqrt( sum(xi - xbar)^2 / ((n - b(n,r)) * alpha(r)) )
```

where `alpha(r)` is the exact population variance scaling factor of the truncated normal.

## Citation

If you use this code or data, please cite:

```bibtex
@article{szewczyk2025bias,
  title={Data-driven bias correction of standard deviation estimators via symbolic regression:
         compact formulas for normal and symmetrically truncated distributions
         in industrial measurement applications},
  author={Szewczyk, Roman and Szumiata, Tadeusz and Nowicki, Michal},
  journal={Meas. Sci. Technol.},
  year={2025}
}
```

## Empirical Validation

```bash
cd 03_empirical_validation
# Download household_power_consumption.txt from UCI (see README.md)
python derive_reference_population.py
python run_validation.py
```

The archived `bias_normal.csv` and `bias_truncated.csv` contain the exact
results reported in Section 5 of the manuscript.
