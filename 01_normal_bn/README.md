# 01_normal_bn -- Standard Normal Distribution b(n)

Monte Carlo generators, curve-fitting scripts, and validation suites for the optimal correction factor `b(n)` that makes the sample standard deviation unbiased under `N(0,1)`.

## Directory Layout

### 01_monte_carlo/
Authoritative Monte Carlo computation of `b(n)` for `n = 2 .. 100`.

| File | Purpose |
|------|---------|
| `ident_bn_batched.py` | Main generator. 500 independent repetitions, `k = 2e8` samples per `n`. Batched to fit 64 GB RAM. Runtime ~9-10 hours on GPU. |
| `ident_bn_view.m` | Octave/MATLAB viewer for the results |
| `optimal_bn_500times.mat` | Stored result. Contains `bn_mean`, `bn_std`, `verification`, etc. |
| `Results_ident_bn.mat` | Single-run reference |

You do not need to re-run `ident_bn_batched.py`. The `.mat` file is the final data.

### 02_identification/
Curve fitting and tabulated lookup function.

| File | Purpose |
|------|---------|
| `params_ident_bn_calc_FRACT_b.m` | Bootstrap-refined rational-model fit |
| `prepare_opti.m` | Pre-processing helper |
| `std_opti.m` | Production lookup function. Interpolated tabulated `b(n)` for any sample size `n = 2..100`. |
| `Results_ident_bn.mat` | Tabulated values and statistics |
| `Results_rational_model.mat` | Fitted rational-model parameters |

### 03_first_principle_check/
Independent validation suite.

| File | Purpose |
|------|---------|
| `check_normal.m` | Small-scale FPC (`k = 5e6`, `n = 2..20`). Compares biased, Bessel, n-1.5, Gurland, LFF, and tabulated `b(n)` estimators. |
| `check_normal_calc.m` / `check_normal_calc_fast.m` | Variants with different performance trade-offs |
| `check_normal_view.m` | Viewer for stored FPC results |
| `check_normal_234.m` | Focused check for `n = 2, 3, 4` |
| `validate_vs_c4.m` | Validation of Monte Carlo `b(n)` against the exact chi-distribution `c4` benchmark. |
| `compute_mse_table.m` | Exact bias-variance-MSE table for standard-deviation estimators, computed analytically from chi-distribution theory. |
| `Res_check_normal.mat` | Stored FPC results |
| `std_opti.m` | Self-contained copy of the lookup function |

### 04_symbolic_regression/
PySR experiments that discovered the LFF form.

| File | Purpose |
|------|---------|
| `pysr_symbolic_regression.py` | Experiment A -- general symbolic-regression search |
| `pysr_lff_search.py` | Experiment B -- linear-fractional-function constrained search |
| `visualize_results.py` | Plotting helper |
| `optimal_bn_500times.mat` | Input data (copy from `01_monte_carlo/`) |
| `pysr_hyperparameters.md` | PySR configuration and hyperparameters |
| `requirements.txt` | PySR dependencies |

## How to Run

Validate against exact theory:
```bash
cd 03_first_principle_check
octave validate_vs_c4.m
```

Compute exact bias-variance-MSE table:
```bash
cd 03_first_principle_check
octave compute_mse_table.m
```

Run the Monte Carlo generator (long run, only needed if the .mat file is missing):
```bash
cd 01_monte_carlo
pip install -r ../requirements.txt
python ident_bn_batched.py
```

Re-run symbolic regression:
```bash
cd 04_symbolic_regression
pip install -r requirements.txt
python pysr_symbolic_regression.py
python pysr_lff_search.py
```

## Notes

- Octave scripts assume `pkg load optim` has been executed if the `optim` package is needed.
- The `std_opti.m` lookup function uses `interp1(..., "nearest")` for sample sizes outside the tabulated range. For integer `n` in `2..100` it returns the exact tabulated value.
