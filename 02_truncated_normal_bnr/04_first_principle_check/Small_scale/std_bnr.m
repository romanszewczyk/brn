function s = std_bnr(x,r)

% r must be the scalar

  % If x is a vector, treat it as n x 1
  if isvector(x)
    x = x(:);
  end

  % Now x is n x m
  [n, m] = size(x);

  % Column means: 1 x m
  mu = sum(x, 1) / n;

  % Center data: subtract column means
  xc = x - ones(n, 1) * mu;

  % Baseline (uncorrected) estimator: divisor (n-1)*alpha(r), Eq. 11
  bn = 1;

  % Variance scaling factor alpha(r) of the truncated normal (Eq. 12).
  % Standard normal pdf and cdf written out so no statistics package is needed.
  phi_r = exp(-r.^2 / 2) / sqrt(2*pi);
  Phi_r = 0.5 * (1 + erf(r / sqrt(2)));
  ar = 1 - 2*r*phi_r/(2*Phi_r - 1);

  % Adjusted variance and std
  s2 = sum(xc.^2, 1) / ((n - bn).*ar);


  s = sqrt(s2);

  if m == 1
    s = s(1);
  end

end

