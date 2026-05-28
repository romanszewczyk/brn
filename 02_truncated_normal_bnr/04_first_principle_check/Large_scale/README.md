# Large-Scale Truncated-Normal First-Principle Check (Python/JAX)

GPU-accelerated validation suite for the `b(n,r)` truncated-normal estimator. Python/JAX counterpart to the Octave small-scale reference. Chunked Monte Carlo for `k = 1e9+` trials.

## Files

| File | Purpose |
|------|---------|
| `truncated_normal_estimator.py` | Optimal and baseline estimators, pre-computed `b(n,r)` lookup table, JIT-compiled functions |
| `truncated_normal_generator.py` | Truncated-normal sampler using the Botev (2016) algorithm |
| `monte_carlo_validation.py` | Chunked Monte Carlo framework with online statistics accumulation |
| `analyze_results.py` | Bias comparison, RMSE, 3D surfaces, heat-maps |
| `test_installation.py` | Verify installation and GPU setup |
| `convert_to_mat.py` | Convert `.npz` results to `.mat` for Octave compatibility |

## Quick Start

```bash
pip install -r requirements.txt
python test_installation.py
python monte_carlo_validation.py
python analyze_results.py
```

Edit `monte_carlo_validation.py` to set `n_values`, `r_values`, `k_total`, and `chunk_size` before running.

## Chunking

Memory per chunk is approximately `3 * n * chunk_size * 8` bytes. For an RTX 4090 (24 GB), practical chunk sizes are 10k-100k. The chunking formula used is:

```python
chunk_size = int((gpu_memory_gb * 0.5 * 1e9) / (n_max * 8))
```

## Performance (RTX 4090)

| k_total | chunk_size | Time |
|---------|------------|------|
| 1e6 | 10,000 | ~1 min |
| 1e7 | 50,000 | ~5 min |
| 1e8 | 100,000 | ~30 min |
| 1e9 | 100,000 | ~2-3 h |

## Estimators

**Optimal:** `sigma_hat = sqrt(SS / ((n - b(n,r)) * alpha(r)))`

**Baseline:** `sigma_hat = sqrt(SS / ((n - 1) * alpha(r)))`

where `alpha(r) = 1 - 2*r*phi(r)/(2*Phi(r) - 1)`.

## Notes

- All computations use float64.
- Scripts set `XLA_PYTHON_CLIENT_PREALLOCATE=false` to avoid grabbing the entire GPU.
- Random seeds are fixed for reproducibility.
- The `b(n,r)` table was extracted from the Octave `std_bnr_opti.m` with full 5-digit precision.
