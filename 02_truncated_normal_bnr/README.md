# 02_truncated_normal_bnr -- Truncated normal distribution b(n,r)

Monte Carlo generators, surface-fitting scripts, and validation suites for the
correction factor `b(n,r)` that makes the sample standard deviation unbiased
under a symmetrically truncated normal distribution.

## Directory layout

### 01_monte_carlo/
Monte Carlo computation of `b(n,r)` for `n = 2..100` and `r = 0.2..10.0`.

| File | Purpose |
|------|---------|
| `ident_brn.py` | Main generator. 25 repetitions, `k = 2e8` samples per `(n,r)`, inverse-CDF sampling, per-`n` processing (peak host RAM ~2 GB), consolidated result rewritten after each repetition. Runtime ~100-300 h on GPU. |
| `optimal_brn_25reps.mat` | Stored result: `b(n,r)` over all repetitions. |
| `rbn_MC_data_uniform.m` | Averages the repetitions and pchip-interpolates to the uniform `(n,r)` grid used for fitting; writes `rbn_MC_data.mat`. |
| `rbn_MC_data.mat` | Uniform-grid `b(n,r)` surface (input to the surface fit). |
| `brn_view.m` | Surface viewer. |

You do not need to re-run the generator; the `.mat` files are the stored data.

### 02_surface_identification/
Fit of the second-degree rational correction (Eq. 15) to the Monte Carlo surface.

| File | Purpose |
|------|---------|
| `params_ident_brn.m` | Fits the five coefficient functions `P(r), A(r), B(r), Q(r), C(r)` by constrained differential evolution (`de_min`, Octave `optim`). |
| `MATHEMATICAL_BACKGROUND.md` | Model and regularization specification. |
| `analyze_results_PLOTS.m` | Coefficient-function and fit/error figures. |
| `rbn_MC_data.mat` | Input surface (copy from `01_monte_carlo/`). |
| `Results_surface_fit.mat` | Fitted coefficients and the realised surface. |

### 03_lookup_and_scripts/
| File | Purpose |
|------|---------|
| `std_bnr_opti.m` | Tabulated `b(n,r)` lookup with bilinear interpolation; accepts vector `n` and scalar `r`. |

### 04_first_principle_check/
Independent validation of the truncated-normal estimator.

#### Small_scale/ -- Octave reference
Quick sanity checks in Octave.

| File | Purpose |
|------|---------|
| `std_bnr_opti.m` | Lookup function |
| `std_bnr.m` | Baseline estimator (`alpha(r)` only) |
| `trandn.m` | Truncated normal sampler |
| `prepare_bnr_opti.m` | Pre-processing |
| `std_bnr_opti_test03_calc.m` / `_view.m` | Test scripts |
| `analyze_results_octave.m` | Plotting |
| `Results_std_brn.mat` | Stored results |

#### Large_scale/ -- Python/JAX validation
GPU-accelerated validation with chunking for large `k`. This check uses an
independent sampler (Botev 2016) from the inverse-CDF method of the main
identification, so agreement is a cross-method confirmation.

| File | Purpose |
|------|---------|
| `truncated_normal_estimator.py` | Optimal and baseline estimators (JIT) |
| `truncated_normal_generator.py` | Truncated-normal sampler |
| `monte_carlo_validation.py` | Chunked Monte Carlo framework |
| `analyze_results.py` | Bias, RMSE, and surface plots |
| `test_installation.py` | Installation check |
| `convert_to_mat.py` | Convert `.npz` results to `.mat` |
| `results_std_bnr.mat` / `.npz` | Stored validation results |

## How to run

Inspect the fitted surface:
```bash
cd 02_surface_identification
octave analyze_results_PLOTS.m
```

Small-scale Octave check:
```bash
cd 04_first_principle_check/Small_scale
octave std_bnr_opti_test03_calc.m
octave std_bnr_opti_test03_view.m
```

Large-scale Python/JAX check:
```bash
cd 04_first_principle_check/Large_scale
pip install -r requirements.txt
python test_installation.py
python monte_carlo_validation.py   # edit config for k, n_values, r_values
python analyze_results.py
```

Re-run the surface fit (long run, only needed if `Results_surface_fit.mat` is missing):
```bash
cd 02_surface_identification
octave
>> pkg load optim
>> params_ident_brn
```

## Notes

- `std_bnr_opti.m` reverts to the normal `b(n)` form for `r >= 4.0`, where
  truncation effects are negligible.
- The supported truncation range is `r >= 0.2`; for `r < 0.3` the generator can
  show numerical instability.
- The Large_scale Python scripts set `XLA_PYTHON_CLIENT_PREALLOCATE=false` to
  avoid reserving the whole GPU at startup.
