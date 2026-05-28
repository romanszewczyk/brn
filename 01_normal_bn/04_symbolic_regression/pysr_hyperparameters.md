# PySR Configuration and Reproducibility

PySR version 1.5.8 was used for all symbolic-regression runs. Each run used a single independent execution (`n_runs = 1`) with fixed random seed 42. The 500 trials in `optimal_bn_500times.mat` refer to Monte-Carlo averaging of the `b_n` time series, not to repeated symbolic-regression runs.

## Experiment A -- General Symbolic-Regression Search

*Source: `pysr_symbolic_regression.py`*

| Hyperparameter | Value |
|----------------|-------|
| Operator set | Binary: `+`, `-`, `*`, `/` <br> Unary: `inv(x) = 1/x`, `sqr(x) = x^2`, `safe_sqrt(x)` |
| Population size | 300 |
| Number of populations | 100 |
| Number of iterations | 1,000 |
| Loss function | Mean-squared error: `loss(x, y) = (x - y)^2` |
| Parsimony / complexity | `parsimony = 0.01`; `maxsize = 30`; `maxdepth = 6` |
| Independent runs | 1 |
| Random seed | 42 |
| Parallel workers | 4 (`procs=4`) |

## Experiment B -- Linear-Fractional-Function (LFF) Constrained Search

*Source: `pysr_lff_search.py`*

| Hyperparameter | Value |
|----------------|-------|
| Operator set | Binary: `+`, `-`, `*`, `/` <br> Unary: `inv(x) = 1/x` |
| Population size | 100 |
| Number of populations | 40 |
| Number of iterations | 200 |
| Loss function | Mean-squared error: `loss(x, y) = (x - y)^2` |
| Parsimony / complexity | `parsimony = 0.05`; `maxsize = 15`; `maxdepth = 4` <br> Operator complexities: `{'/': 2, '+': 1, '-': 1, '*': 1}` <br> Constant complexity: 1 |
| Structural constraints | `nested_constraints = {'/': {'/': 0}}` (forbids nested division) |
| Annealing parameter | `alpha = 0.1` |
| Independent runs | 1 |
| Random seed | 42 |
| Parallel workers | 4 (`procs=4`) |

## Reproducibility

Both scripts load the averaged data vector `b_n` (99 points, `n` from the `.mat` file) and call `model.fit(X, y)` once. Because the random seed is fixed to 42 and the PySR version is pinned in `requirements.txt`, re-running the scripts on the same hardware recovers the same discovered equations.

The Pareto-optimal LFF discovered in Experiment A is:

```
b(n) = (1.500 n - 1.637) / (n - 0.999)
```

which was subsequently bootstrap-refined to the final LFF form.
