# Model and constraints for the b(n,r) surface fit

This note specifies the model and the regularization used by `params_ident_brn.m`
to fit the bias-correction factor `b(n,r)` of the symmetrically truncated normal
distribution.

## Model

The correction is the second-degree rational function (Eq. 15 of the paper):

```
         P(r) n^2 + A(r) n + B(r)
b(n,r) = ------------------------
         Q(r) n^2 + C(r) n + 1
```

with `n` the sample size (2 to 100) and `r` the truncation level (0.2 to 10.0).
The denominator constant is fixed to 1 to remove the overall scaling ambiguity
of the ratio. The five coefficient functions `P(r), A(r), B(r), Q(r), C(r)` are
represented at discrete control points in `r` and interpolated with monotone
cubic (pchip) splines:

```
r_control = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8,
             2.0, 2.3, 2.6, 2.9, 3.0, 3.3, 3.6, 4.0,
             4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0]   (26 points)
```

giving 26 x 5 = 130 free parameters.

## Regional behaviour

| Region                | r range      | Behaviour                                        |
|-----------------------|--------------|--------------------------------------------------|
| Truncated             | r < 4        | Full quadratic numerator and denominator         |
| Boundary              | r = 4        | Problem structure changes; correction reverts to the normal `b(n)` form for r >= 4 |
| Transition            | 4 <= r < 5.5 | `P, Q` decay toward zero                          |
| Asymptotic            | r >= 5.5     | `P ~ 0`, `Q ~ 0`; model reduces to `(A n + B)/(C n + 1)` |

## Objective

The fit minimises a weighted sum of squared residuals plus regularization
penalties:

```
Objective = sum_{n,r} w(n,r) (b_model(n,r) - b_data(n,r))^2 + penalties
```

The weights `w(n,r)` emphasise the small-`n`, low-`r` region where the
correction is largest and most sensitive, and de-emphasise the saturated region
`r >= 5.5` where `b(n,r)` is close to the normal case.

The regularization penalties encode the physical and asymptotic expectations of
the coefficient functions:

1. monotonicity of `P(r)` and `Q(r)`;
2. asymptotic decay `P, Q -> 0` and `A, B, C -> ` their measured asymptotic
   values for `r >= 5.5`;
3. high-`n` asymptotic guidance for `A(r)` and `C(r)` at `r > 5.5`;
4. second-derivative smoothness of all coefficient functions;
5. low-`r` stability of the first control points;
6. extra smoothness across the `4 <= r <= 5.5` transition;
7. continuity of the transition zone `5.0 <= r <= 5.5`.

The asymptotic targets `A_inf, B_inf, C_inf` (and the high-`n` targets) are
measured directly from the Monte Carlo data at large `r` rather than imposed,
and are stored in the result file.

## Optimisation

The 130 parameters are fitted by differential evolution (`de_min` from the
Octave `optim` package) in three stages with progressively tighter bounds around
the incumbent solution: stage 1 fits the weighted data alone, stages 2 and 3 add
the regularization penalties and refine. Fit quality and the realised
coefficient functions are written to `Results_surface_fit.mat` and can be
inspected with `analyze_results_PLOTS.m`.
