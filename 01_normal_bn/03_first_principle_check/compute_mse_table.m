% compute_mse_table.m
% Exact bias, variance, and MSE for the Bessel, Proposed (unbiased), and
% Gurland standard-deviation estimators under N(0,1).
%
% The values are computed analytically from the chi-distribution c4 factor,
% not from Monte Carlo simulation, ensuring full numerical precision.
%
% Usage:
%   compute_mse_table

clear all
clc
more off

% -------------------------------------------------------------------------
% Sample sizes to report
% -------------------------------------------------------------------------
nT = [2, 3, 5, 10, 20, 30, 50, 100];

% -------------------------------------------------------------------------
% Exact c4(n) and related quantities
% -------------------------------------------------------------------------
c4 = @(n) sqrt(2./(n-1)) .* gamma(n./2) ./ gamma((n-1)./2);

% Bessel-corrected std: E[su] = c4(n)*sigma, sigma = 1
% Bias = c4(n) - 1
% Variance of su = 1 - c4(n)^2  (for sigma = 1)
% MSE = Bias^2 + Var

% Proposed unbiased estimator: E[s_opt] = sigma (bias = 0 by construction)
% Var(s_opt) = sigma^2 * (1/c4(n)^2 - 1)
% MSE = Var (since bias = 0)

% Gurland estimator denominator: n - 1.5 + 1/(8*(n-1))
% E[s_g] = c4(n) * sqrt((n-1) / (n - 1.5 + 1/(8*(n-1))))
% Var(s_g) = (n-1)/(n - 1.5 + 1/(8*(n-1))) * (1 - c4(n)^2)
%   (scale factor multiplies the Bessel variance)
% Bias = E[s_g] - 1
% MSE = Bias^2 + Var

fprintf('Table 2. Bias, Variance, and MSE for standard-deviation estimators (sigma = 1)\n');
fprintf('-------------------------------------------------------------------------------\n');
fprintf('%4s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
        'n', 'Bias_Bes', 'Var_Bes', 'MSE_Bes', 'Bias_Prop', 'Var_Prop', 'MSE_Prop');
fprintf('-------------------------------------------------------------------------------\n');

for i = 1:numel(nT)
    n = nT(i);
    c4n = c4(n);

    % Bessel
    bias_bes = c4n - 1;
    var_bes  = 1 - c4n^2;
    mse_bes  = bias_bes^2 + var_bes;

    % Proposed (unbiased)
    bias_prop = 0;
    var_prop  = 1/c4n^2 - 1;
    mse_prop  = var_prop;

    % Gurland
    denom_g = n - 1.5 + 1/(8*(n-1));
    scale_g = (n-1) / denom_g;
    bias_g  = c4n * sqrt(scale_g) - 1;
    var_g   = scale_g * (1 - c4n^2);
    mse_g   = bias_g^2 + var_g;

    fprintf('%4d  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f\n', ...
            n, bias_bes, var_bes, mse_bes, bias_prop, var_prop, mse_prop);
end

fprintf('-------------------------------------------------------------------------------\n');
fprintf('\n');

fprintf('Gurland estimator (for comparison):\n');
fprintf('-------------------------------------------------------------------------------\n');
fprintf('%4s  %12s  %12s  %12s  %12s\n', 'n', 'Bias_Gur', 'Var_Gur', 'MSE_Gur', 'MSEratio_Gur/Bes');
fprintf('-------------------------------------------------------------------------------\n');

for i = 1:numel(nT)
    n = nT(i);
    c4n = c4(n);

    % Bessel (reference)
    bias_bes = c4n - 1;
    var_bes  = 1 - c4n^2;
    mse_bes  = bias_bes^2 + var_bes;

    % Gurland
    denom_g = n - 1.5 + 1/(8*(n-1));
    scale_g = (n-1) / denom_g;
    bias_g  = c4n * sqrt(scale_g) - 1;
    var_g   = scale_g * (1 - c4n^2);
    mse_g   = bias_g^2 + var_g;

    fprintf('%4d  %12.4f  %12.4f  %12.4f  %12.4f\n', ...
            n, bias_g, var_g, mse_g, mse_g/mse_bes);
end

fprintf('-------------------------------------------------------------------------------\n');
