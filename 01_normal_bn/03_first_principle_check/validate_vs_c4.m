% validate_vs_c4.m
% Validation of Monte-Carlo b(n) against exact chi-distribution c4 benchmark.
%
% Usage:
%   validate_vs_c4
%
% Loads optimal_bn_500times.mat, computes the exact theoretical correction
% b_theory(n) from the gamma function, and compares it with the tabulated
% Monte-Carlo b(n).

clear all
clc
more off

% -------------------------------------------------------------------------
% Load Monte-Carlo results
% -------------------------------------------------------------------------
data = load('../01_monte_carlo/optimal_bn_500times.mat');
n_vals = double(data.n(1,:));   % 1 x 99, n = 2 .. 100 (cast: stored as int64)
bn_mc  = data.bn_mean(1,:);     % 1 x 99, mean b(n) over 500 repetitions

% -------------------------------------------------------------------------
% Exact c4(n) from chi distribution and theoretical b_theory(n)
% -------------------------------------------------------------------------
% c4(n) = sqrt(2/(n-1)) * gamma(n/2) / gamma((n-1)/2)
% b_theory(n) = n - (n-1) * c4(n)^2
% -------------------------------------------------------------------------

c4 = @(n) sqrt(2./(n-1)) .* gamma(n./2) ./ gamma((n-1)./2);
b_theory = @(n) n - (n-1) .* c4(n).^2;

bn_theory = b_theory(n_vals);

% -------------------------------------------------------------------------
% Deviations
% -------------------------------------------------------------------------
abs_dev = abs(bn_mc - bn_theory);
rel_dev = abs_dev ./ bn_theory;

max_abs_dev = max(abs_dev);
mean_abs_dev = mean(abs_dev);
max_rel_dev = max(rel_dev);

fprintf('Validation of tabulated b(n) against exact c4 benchmark\n');
fprintf('--------------------------------------------------------\n');
fprintf('Max absolute deviation  : %1.4e\n', max_abs_dev);
fprintf('Mean absolute deviation : %1.4e\n', mean_abs_dev);
fprintf('Max relative deviation  : %1.4e (%.3f%%)\n', max_rel_dev, max_rel_dev*100);
fprintf('\n');

% Coverage check: how many bn_mc fall inside 95% CI of theoretical prediction
% The MC standard error is std(bn)/sqrt(500)
se_bn = data.bn_std(1,:) / sqrt(500);
% For df = 499, the two-sided 95% t-critical value is 1.9647.
% Hard-coded to avoid dependency on the Octave statistics package.
t_crit = 1.9647;
ci_half = t_crit * se_bn;
inside = abs(bn_mc - bn_theory) <= ci_half;
coverage = sum(inside) / numel(inside);
fprintf('Coverage inside 95%% CI  : %d / %d (%.1f%%)\n', sum(inside), numel(inside), coverage*100);

% -------------------------------------------------------------------------
% Plot comparison and deviation
% -------------------------------------------------------------------------
figure('Position', [100 100 1200 400]);

subplot(1, 2, 1);
plot(n_vals, bn_mc,   'o', 'MarkerSize', 4, 'DisplayName', 'Monte Carlo b(n)');
hold on;
plot(n_vals, bn_theory, '-r', 'LineWidth', 1.5, 'DisplayName', 'Exact c_4 benchmark');
hold off;
xlim([2 100]);
xlabel('\it{n}');
ylabel('\it{b(n)}');
legend('Location', 'southeast');
grid on;
title('(a) Comparison of b(n) values');

subplot(1, 2, 2);
semilogy(n_vals, abs_dev, '-b', 'LineWidth', 1.5);
xlim([2 100]);
xlabel('\it{n}');
ylabel('|b_{MC}(n) - b_{theory}(n)|');
grid on;
title('(b) Absolute deviation (log scale)');

fprintf('\nFigure generated.\n');
