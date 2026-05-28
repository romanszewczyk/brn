% analyze_results_PLOTS.m
% Coefficient-function and error figures for the b(n,r) surface fit.
% Reproduces the P(r), A(r), B(r), Q(r), C(r) dependencies and the fit/error
% surfaces. Octave compatible.

clear all
close all
more off;

if ~exist('Results_surface_fit.mat', 'file')
    error('Results_surface_fit.mat not found. Run params_ident_brn.m first.');
end
load('Results_surface_fit.mat');
fprintf('Loaded Results_surface_fit.mat\n\n');

% Extract dimensions
[n_r, n_n] = size(rbn_mean);
n_vals = nn(1,:);
r_vals = rr(:,1);
err_all = rbn_calc - rbn_mean;
n_control = length(r_control_points);

% Extract coefficients
P_vals = a_res(1:n_control);
A_vals = a_res(n_control+1:2*n_control);
B_vals = a_res(2*n_control+1:3*n_control);
Q_vals = a_res(3*n_control+1:4*n_control);
C_vals = a_res(4*n_control+1:5*n_control);

% Check for variables that may not exist
if ~exist('A_high_n', 'var')
    has_high_n = false;
    A_high_n = -1.5;
    C_high_n = -1.0;
else
    has_high_n = true;
end

if ~exist('rmse_all', 'var')
    rmse_all = sqrt(mean(err_all(:).^2));
end

if ~exist('rmse_crit', 'var')
    idx_crit = (rr < 2.0) & (nn <= 4);
    rmse_crit = sqrt(mean(err_all(idx_crit).^2));
end

if ~exist('R2', 'var')
    R2 = 1 - sum(err_all(:).^2) / sum((rbn_mean(:) - mean(rbn_mean(:))).^2);
end

if ~exist('A_inf', 'var')
    A_inf = -1.0;
    B_inf = 1.5;
    C_inf = -0.7;
end

fprintf('Creating visualizations...\n\n');

%==========================================================================
% FIGURE 1: COEFFICIENT EVOLUTION WITH CONTROL POINTS
%==========================================================================

fprintf('Figure 1: Coefficient evolution with control points...\n');

figure(1);
set(gcf, 'Position', [100, 100, 1400, 900]);

r_plot = linspace(min(r_vals), max(r_vals), 300);

% Use pchip (monotonic) - same as optimization
P_plot = pchip(r_control_points, P_vals, r_plot);
A_plot = pchip(r_control_points, A_vals, r_plot);
B_plot = pchip(r_control_points, B_vals, r_plot);
Q_plot = pchip(r_control_points, Q_vals, r_plot);
C_plot = pchip(r_control_points, C_vals, r_plot);

% P coefficient
subplot(2,3,1);
plot(r_plot, P_plot, 'b-', 'LineWidth', 2.5); hold on;
plot(r_control_points, P_vals, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot([4 5], [0 0], 'g--', 'LineWidth', 2);
xlabel('r', 'FontSize', 11);
ylabel('P(r)', 'FontSize', 11);
title('P Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Interpolation', 'Control points', 'Problem region 4-5', 'Location', 'northeast');

% A coefficient
subplot(2,3,2);
plot(r_plot, A_plot, 'b-', 'LineWidth', 2.5); hold on;
plot(r_control_points, A_vals, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
if has_high_n
    plot([5 max(r_vals)], [A_high_n A_high_n], 'g--', 'LineWidth', 2);
end
xlabel('r', 'FontSize', 11);
ylabel('A(r)', 'FontSize', 11);
title('A Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
if has_high_n
    legend('Interpolation', 'Control points', 'High-n target', 'Location', 'northeast');
else
    legend('Interpolation', 'Control points', 'Location', 'northeast');
end

% B coefficient
subplot(2,3,3);
plot(r_plot, B_plot, 'b-', 'LineWidth', 2.5); hold on;
plot(r_control_points, B_vals, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('r', 'FontSize', 11);
ylabel('B(r)', 'FontSize', 11);
title('B Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Interpolation', 'Control points', 'Location', 'northeast');

% Q coefficient
subplot(2,3,4);
plot(r_plot, Q_plot, 'b-', 'LineWidth', 2.5); hold on;
plot(r_control_points, Q_vals, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot([4 5], [0 0], 'g--', 'LineWidth', 2);
xlabel('r', 'FontSize', 11);
ylabel('Q(r)', 'FontSize', 11);
title('Q Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Interpolation', 'Control points', 'Problem region 4-5', 'Location', 'northeast');

% C coefficient
subplot(2,3,5);
plot(r_plot, C_plot, 'b-', 'LineWidth', 2.5); hold on;
plot(r_control_points, C_vals, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
if has_high_n
    plot([5 max(r_vals)], [C_high_n C_high_n], 'g--', 'LineWidth', 2);
end
xlabel('r', 'FontSize', 11);
ylabel('C(r)', 'FontSize', 11);
title('C Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
if has_high_n
    legend('Interpolation', 'Control points', 'High-n target', 'Location', 'northeast');
else
    legend('Interpolation', 'Control points', 'Location', 'northeast');
end

% Summary panel
subplot(2,3,6);
axis off;
text(0.1, 0.90, 'SUMMARY', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.78, sprintf('Overall RMSE: %.6f', rmse_all), 'FontSize', 11);
text(0.1, 0.68, sprintf('R-squared: %.6f', R2), 'FontSize', 11);
text(0.1, 0.58, sprintf('Control points: %d', n_control), 'FontSize', 10);
text(0.1, 0.48, sprintf('Data points: %d', n_r * n_n), 'FontSize', 10);
text(0.1, 0.36, 'Method:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.28, '  Monotonic splines (pchip)', 'FontSize', 9);
text(0.1, 0.20, '  Adaptive regional weighting', 'FontSize', 9);
text(0.1, 0.12, '  Asymptotic constraints for r >= 5.5', 'FontSize', 9);

print('-dpng', 'Fig1_Coefficient_Evolution.png', '-r300');
fprintf('  Saved: Fig1_Coefficient_Evolution.png\n');

%==========================================================================
% FIGURE 2: 2D COMPARISON - ORIGINAL VS FITTED DATA
%==========================================================================

fprintf('Figure 2: 2D comparison - original vs fitted data...\n');

figure(2);
set(gcf, 'Position', [120, 120, 1600, 900]);

% Original data
subplot(2,3,1);
imagesc(n_vals, r_vals, rbn_mean);
axis xy;
colorbar;
xlabel('n', 'FontSize', 11);
ylabel('r', 'FontSize', 11);
title('Original Data (rbn MC data)', 'FontSize', 12, 'FontWeight', 'bold');
colormap(jet);

% Fitted data
subplot(2,3,2);
imagesc(n_vals, r_vals, rbn_calc);
axis xy;
colorbar;
xlabel('n', 'FontSize', 11);
ylabel('r', 'FontSize', 11);
title('Fitted Data (Model)', 'FontSize', 12, 'FontWeight', 'bold');
colormap(jet);

% Error (absolute)
subplot(2,3,3);
imagesc(n_vals, r_vals, abs(err_all));
axis xy;
colorbar;
xlabel('n', 'FontSize', 11);
ylabel('r', 'FontSize', 11);
title('Absolute Error', 'FontSize', 12, 'FontWeight', 'bold');
colormap(hot);

% Relative error (percent)
rel_err = 100 * abs(err_all) ./ (abs(rbn_mean) + 1e-10);
subplot(2,3,4);
imagesc(n_vals, r_vals, rel_err);
axis xy;
colorbar;
xlabel('n', 'FontSize', 11);
ylabel('r', 'FontSize', 11);
title('Relative Error (%)', 'FontSize', 12, 'FontWeight', 'bold');
colormap(hot);
caxis([0 5]);

% Error distribution histogram
subplot(2,3,5);
hist(err_all(:), 50);
xlabel('Error', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
title('Error Distribution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Scatter plot: calculated vs original
subplot(2,3,6);
plot(rbn_mean(:), rbn_calc(:), 'b.', 'MarkerSize', 2); hold on;
plot([min(rbn_mean(:)) max(rbn_mean(:))], [min(rbn_mean(:)) max(rbn_mean(:))], ...
     'r-', 'LineWidth', 2);
xlabel('Original Data', 'FontSize', 11);
ylabel('Model Fit', 'FontSize', 11);
title('Calculated vs Original', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('Data', 'Perfect fit', 'Location', 'northwest');

print('-dpng', 'Fig2_Data_Comparison_2D.png', '-r300');
fprintf('  Saved: Fig2_Data_Comparison_2D.png\n');

%==========================================================================
% FIGURE 3: ERROR HEATMAP WITH PROBLEM REGIONS
%==========================================================================

fprintf('Figure 3: Error heatmap with problem regions...\n');

figure(3);
set(gcf, 'Position', [140, 140, 1400, 800]);

% Overall error heatmap
subplot(2,2,1);
imagesc(n_vals, r_vals, abs(err_all));
axis xy;
colorbar;
xlabel('n', 'FontSize', 12);
ylabel('r', 'FontSize', 12);
title('Absolute Error Heatmap', 'FontSize', 13, 'FontWeight', 'bold');
colormap(hot);
hold on;
% Mark problem regions
plot([min(n_vals) max(n_vals)], [2 2], 'c--', 'LineWidth', 2);
plot([min(n_vals) max(n_vals)], [4 4], 'g--', 'LineWidth', 2);
plot([min(n_vals) max(n_vals)], [5 5], 'g--', 'LineWidth', 2);
plot([min(n_vals) max(n_vals)], [5.5 5.5], 'm--', 'LineWidth', 2);

% Log-scale error
subplot(2,2,2);
imagesc(n_vals, r_vals, log10(abs(err_all) + 1e-10));
axis xy;
colorbar;
xlabel('n', 'FontSize', 12);
ylabel('r', 'FontSize', 12);
title('Log10 Absolute Error', 'FontSize', 13, 'FontWeight', 'bold');
colormap(hot);

% Relative error
subplot(2,2,3);
rel_err_map = 100 * abs(err_all) ./ (abs(rbn_mean) + 1e-10);
imagesc(n_vals, r_vals, rel_err_map);
axis xy;
colorbar;
xlabel('n', 'FontSize', 12);
ylabel('r', 'FontSize', 12);
title('Relative Error (%)', 'FontSize', 13, 'FontWeight', 'bold');
colormap(hot);
caxis([0 10]);

% Error by region bar chart
subplot(2,2,4);
regions_r = {'r<1', '1<=r<2', '2<=r<3', '3<=r<4', '4<=r<5', '5<=r<5.5', 'r>=5.5'};
idx_r1 = find(r_vals<1);
idx_r2 = find(r_vals>=1 & r_vals<2);
idx_r3 = find(r_vals>=2 & r_vals<3);
idx_r4 = find(r_vals>=3 & r_vals<4);
idx_r5 = find(r_vals>=4 & r_vals<5);
idx_r6 = find(r_vals>=5 & r_vals<5.5);
idx_r7 = find(r_vals>=5.5);

rmse_by_r = zeros(7, 1);
if ~isempty(idx_r1), rmse_by_r(1) = sqrt(mean(mean((err_all(idx_r1,:)).^2))); end
if ~isempty(idx_r2), rmse_by_r(2) = sqrt(mean(mean((err_all(idx_r2,:)).^2))); end
if ~isempty(idx_r3), rmse_by_r(3) = sqrt(mean(mean((err_all(idx_r3,:)).^2))); end
if ~isempty(idx_r4), rmse_by_r(4) = sqrt(mean(mean((err_all(idx_r4,:)).^2))); end
if ~isempty(idx_r5), rmse_by_r(5) = sqrt(mean(mean((err_all(idx_r5,:)).^2))); end
if ~isempty(idx_r6), rmse_by_r(6) = sqrt(mean(mean((err_all(idx_r6,:)).^2))); end
if ~isempty(idx_r7), rmse_by_r(7) = sqrt(mean(mean((err_all(idx_r7,:)).^2))); end

bar(rmse_by_r);
set(gca, 'XTickLabel', regions_r);
ylabel('RMSE', 'FontSize', 12);
title('RMSE by r Region', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

print('-dpng', 'Fig3_Error_Heatmap.png', '-r300');
fprintf('  Saved: Fig3_Error_Heatmap.png\n');

%==========================================================================
% FIGURE 4: REGIONAL ERROR ANALYSIS
%==========================================================================

fprintf('Figure 4: Regional error analysis...\n');

figure(4);
set(gcf, 'Position', [160, 160, 1600, 900]);

% Error vs n for different r ranges
subplot(2,3,1);
idx_r_crit = find(r_vals < 2);
if ~isempty(idx_r_crit)
    plot(n_vals, mean(abs(err_all(idx_r_crit, :)), 1), 'b-', 'LineWidth', 2);
    xlabel('n', 'FontSize', 11);
    ylabel('Mean Absolute Error', 'FontSize', 11);
    title('Critical r<2', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

subplot(2,3,2);
idx_r_trans = find(r_vals >= 2 & r_vals < 4);
if ~isempty(idx_r_trans)
    plot(n_vals, mean(abs(err_all(idx_r_trans, :)), 1), 'g-', 'LineWidth', 2);
    xlabel('n', 'FontSize', 11);
    ylabel('Mean Absolute Error', 'FontSize', 11);
    title('Transition 2<=r<4', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

subplot(2,3,3);
idx_r_prob = find(r_vals >= 4 & r_vals < 5);
if ~isempty(idx_r_prob)
    plot(n_vals, mean(abs(err_all(idx_r_prob, :)), 1), 'r-', 'LineWidth', 2);
    xlabel('n', 'FontSize', 11);
    ylabel('Mean Absolute Error', 'FontSize', 11);
    title('Problem 4<=r<5', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

% Error vs r for different n ranges
subplot(2,3,4);
idx_n_low = find(n_vals <= 4);
if ~isempty(idx_n_low)
    plot(r_vals, mean(abs(err_all(:, idx_n_low)), 2), 'b-', 'LineWidth', 2);
    xlabel('r', 'FontSize', 11);
    ylabel('Mean Absolute Error', 'FontSize', 11);
    title('Low n<=4', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

subplot(2,3,5);
idx_n_mid = find(n_vals > 4 & n_vals <= 40);
if ~isempty(idx_n_mid)
    plot(r_vals, mean(abs(err_all(:, idx_n_mid)), 2), 'g-', 'LineWidth', 2);
    xlabel('r', 'FontSize', 11);
    ylabel('Mean Absolute Error', 'FontSize', 11);
    title('Mid 4<n<=40', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

subplot(2,3,6);
idx_n_high = find(n_vals > 40);
if ~isempty(idx_n_high)
    plot(r_vals, mean(abs(err_all(:, idx_n_high)), 2), 'r-', 'LineWidth', 2);
    xlabel('r', 'FontSize', 11);
    ylabel('Mean Absolute Error', 'FontSize', 11);
    title('High n>40', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

print('-dpng', 'Fig4_Error_Analysis.png', '-r300');
fprintf('  Saved: Fig4_Error_Analysis.png\n');

%==========================================================================
% FIGURE 5: FIT COMPARISON AT SELECTED R VALUES
%==========================================================================

fprintf('Figure 5: Fit comparison at selected r values...\n');

figure(5);
set(gcf, 'Position', [180, 180, 1600, 1000]);

% Select representative r values
r_select = [0.5, 1.0, 2.0, 3.0, 4.5, 6.0, 8.0, 10.0];
n_plots = min(length(r_select), 8);

for k = 1:n_plots
    r_target = r_select(k);
    [~, r_idx] = min(abs(r_vals - r_target));
    
    subplot(2, 4, k);
    plot(n_vals, rbn_mean(r_idx, :), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
    plot(n_vals, rbn_calc(r_idx, :), 'r-', 'LineWidth', 2);
    xlabel('n', 'FontSize', 9);
    ylabel('rbn', 'FontSize', 9);
    
    if r_vals(r_idx) >= 4 && r_vals(r_idx) <= 5
        title_str = sprintf('r=%.2f (PROBLEM)', r_vals(r_idx));
    elseif r_vals(r_idx) < 2
        title_str = sprintf('r=%.2f (CRITICAL)', r_vals(r_idx));
    else
        title_str = sprintf('r=%.2f', r_vals(r_idx));
    end
    title(title_str, 'FontSize', 10);
    
    legend('Data', 'Model', 'Location', 'northeast', 'FontSize', 7);
    grid on;
    
    % Add RMSE
    err_slice = rbn_calc(r_idx, :) - rbn_mean(r_idx, :);
    rmse_slice = sqrt(mean(err_slice.^2));
    text(0.6, 0.1, sprintf('RMSE=%.5f', rmse_slice), 'Units', 'normalized', ...
         'FontSize', 7, 'BackgroundColor', 'w');
end

print('-dpng', 'Fig5_Fit_Comparison.png', '-r300');
fprintf('  Saved: Fig5_Fit_Comparison.png\n');

%==========================================================================
% FIGURE 6: MONOTONICITY AND ASYMPTOTIC BEHAVIOR
%==========================================================================

fprintf('Figure 6: Monotonicity and asymptotic behavior...\n');

figure(6);
set(gcf, 'Position', [200, 200, 1400, 800]);

% P and Q evolution (should decrease)
subplot(2,2,1);
plot(r_control_points, P_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b'); hold on;
plot(r_control_points, Q_vals, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r');
plot([min(r_control_points) max(r_control_points)], [0 0], 'k--', 'LineWidth', 1);
plot([4 4], [min([P_vals(:); Q_vals(:)]) max([P_vals(:); Q_vals(:)])], ...
     'g--', 'LineWidth', 2);
plot([5 5], [min([P_vals(:); Q_vals(:)]) max([P_vals(:); Q_vals(:)])], ...
     'g--', 'LineWidth', 2);
plot([5.5 5.5], [min([P_vals(:); Q_vals(:)]) max([P_vals(:); Q_vals(:)])], ...
     'm--', 'LineWidth', 2);
xlabel('r', 'FontSize', 12);
ylabel('Coefficient Value', 'FontSize', 12);
title('P and Q Evolution (should approach 0)', 'FontSize', 13, 'FontWeight', 'bold');
legend('P', 'Q', 'Zero', 'Problem 4-5', 'Problem 4-5', 'Saturation', 'Location', 'northeast');
grid on;

% A, B, C evolution with asymptotic targets
subplot(2,2,2);
plot(r_control_points, A_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b'); hold on;
plot(r_control_points, B_vals, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r');
plot(r_control_points, C_vals, 'g^-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'g');
plot([5.5 5.5], [min([A_vals(:); C_vals(:)]) max([A_vals(:); B_vals(:)])], ...
     'k--', 'LineWidth', 2);
% Mark measured asymptotic values
plot([5.5 max(r_control_points)], [A_inf A_inf], 'b:', 'LineWidth', 2);
plot([5.5 max(r_control_points)], [B_inf B_inf], 'r:', 'LineWidth', 2);
plot([5.5 max(r_control_points)], [C_inf C_inf], 'g:', 'LineWidth', 2);
if has_high_n
    plot([5.5 max(r_control_points)], [A_high_n A_high_n], 'b-.', 'LineWidth', 2);
    plot([5.5 max(r_control_points)], [C_high_n C_high_n], 'g-.', 'LineWidth', 2);
end
xlabel('r', 'FontSize', 12);
ylabel('Coefficient Value', 'FontSize', 12);
title('A, B, C Evolution (dotted: asymptotic, dash-dot: high-n)', ...
      'FontSize', 13, 'FontWeight', 'bold');
if has_high_n
    legend('A', 'B', 'C', 'r=5.5', 'A inf', 'B inf', 'C inf', 'A high-n', 'C high-n', ...
           'Location', 'northwest');
else
    legend('A', 'B', 'C', 'r=5.5', 'A inf', 'B inf', 'C inf', 'Location', 'northwest');
end
grid on;

% P and Q differences between consecutive points
subplot(2,2,3);
dP = diff(P_vals);
dQ = diff(Q_vals);
r_mid = (r_control_points(1:end-1) + r_control_points(2:end)) / 2;
bar_data = [dP(:), dQ(:)];
bar(r_mid, bar_data);
xlabel('r (midpoint)', 'FontSize', 12);
ylabel('Coefficient Change', 'FontSize', 12);
title('P and Q Changes Between Control Points', 'FontSize', 13, 'FontWeight', 'bold');
legend('dP', 'dQ', 'Location', 'northeast');
grid on;
hold on;
plot([min(r_mid) max(r_mid)], [0 0], 'k-', 'LineWidth', 2);

% Asymptotic convergence
subplot(2,2,4);
idx_sat = find(r_control_points >= 5.5);
if ~isempty(idx_sat)
    r_sat = r_control_points(idx_sat);
    plot(r_sat, P_vals(idx_sat), 'bo-', 'LineWidth', 2, 'MarkerSize', 10); hold on;
    plot(r_sat, Q_vals(idx_sat), 'rs-', 'LineWidth', 2, 'MarkerSize', 10);
    plot(r_sat, abs(A_vals(idx_sat) - A_inf), 'g^-', 'LineWidth', 2, 'MarkerSize', 10);
    plot(r_sat, abs(B_vals(idx_sat) - B_inf), 'md-', 'LineWidth', 2, 'MarkerSize', 10);
    plot(r_sat, abs(C_vals(idx_sat) - C_inf), 'cv-', 'LineWidth', 2, 'MarkerSize', 10);
    if has_high_n
        plot(r_sat, abs(A_vals(idx_sat) - A_high_n), 'b*-', 'LineWidth', 1.5, 'MarkerSize', 8);
        plot(r_sat, abs(C_vals(idx_sat) - C_high_n), 'g*-', 'LineWidth', 1.5, 'MarkerSize', 8);
    end
    plot([min(r_sat) max(r_sat)], [0.02 0.02], 'r--', 'LineWidth', 2);
    xlabel('r', 'FontSize', 12);
    ylabel('Value / Error', 'FontSize', 12);
    title('Asymptotic Convergence (r>=5.5)', 'FontSize', 13, 'FontWeight', 'bold');
    if has_high_n
        legend('P', 'Q', '|A-A inf|', '|B-B inf|', '|C-C inf|', ...
               '|A-A high-n|', '|C-C high-n|', 'Target 0.02', 'Location', 'northeast');
    else
        legend('P', 'Q', '|A-A inf|', '|B-B inf|', '|C-C inf|', 'Target 0.02', ...
               'Location', 'northeast');
    end
    grid on;
end

print('-dpng', 'Fig6_Monotonicity_Asymptotic.png', '-r300');
fprintf('  Saved: Fig6_Monotonicity_Asymptotic.png\n');

%==========================================================================
% FIGURE 7: PROBLEM REGIONS DETAILED ANALYSIS
%==========================================================================

fprintf('Figure 7: Problem regions detailed analysis...\n');

figure(7);
set(gcf, 'Position', [220, 220, 1400, 900]);

% Problem region 4<r<5 with low n
subplot(2,2,1);
idx_r_45 = find(r_vals > 4.0 & r_vals < 5.0);
idx_n_low = find(n_vals < 10);
if ~isempty(idx_r_45) && ~isempty(idx_n_low)
    imagesc(n_vals(idx_n_low), r_vals(idx_r_45), abs(err_all(idx_r_45, idx_n_low)));
    axis xy;
    colorbar;
    xlabel('n', 'FontSize', 12);
    ylabel('r', 'FontSize', 12);
    title('Problem: 4<r<5, n<10', 'FontSize', 13, 'FontWeight', 'bold');
    colormap(hot);
end

% Problem region 4<r<5 with high n
subplot(2,2,2);
idx_n_high_45 = find(n_vals > 60);
if ~isempty(idx_r_45) && ~isempty(idx_n_high_45)
    imagesc(n_vals(idx_n_high_45), r_vals(idx_r_45), abs(err_all(idx_r_45, idx_n_high_45)));
    axis xy;
    colorbar;
    xlabel('n', 'FontSize', 12);
    ylabel('r', 'FontSize', 12);
    title('Problem: 4<r<5, n>60', 'FontSize', 13, 'FontWeight', 'bold');
    colormap(hot);
end

% Extended critical region r<4, n<=4
subplot(2,2,3);
idx_r_crit = find(r_vals < 4.0);
idx_n_crit = find(n_vals <= 4);
if ~isempty(idx_r_crit) && ~isempty(idx_n_crit)
    imagesc(n_vals(idx_n_crit), r_vals(idx_r_crit), abs(err_all(idx_r_crit, idx_n_crit)));
    axis xy;
    colorbar;
    xlabel('n', 'FontSize', 12);
    ylabel('r', 'FontSize', 12);
    title('Extended Critical: r<4, n<=4', 'FontSize', 13, 'FontWeight', 'bold');
    colormap(hot);
end

% High-n, r~8 region
subplot(2,2,4);
idx_r_8 = find(r_vals >= 7.5 & r_vals <= 8.5);
idx_n_60 = find(n_vals > 60);
if ~isempty(idx_r_8) && ~isempty(idx_n_60)
    imagesc(n_vals(idx_n_60), r_vals(idx_r_8), abs(err_all(idx_r_8, idx_n_60)));
    axis xy;
    colorbar;
    xlabel('n', 'FontSize', 12);
    ylabel('r', 'FontSize', 12);
    title('High-n,r: n>60, r~8', 'FontSize', 13, 'FontWeight', 'bold');
    colormap(hot);
else
    axis off;
    text(0.5, 0.5, 'Region not available in dataset', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
end

print('-dpng', 'Fig7_Problem_Regions_Detail.png', '-r300');
fprintf('  Saved: Fig7_Problem_Regions_Detail.png\n');

%==========================================================================
% SUMMARY
%==========================================================================

fprintf('\n');
fprintf('================================================================\n');
fprintf('GRAPHICAL ANALYSIS COMPLETE\n');
fprintf('================================================================\n\n');

fprintf('Generated figures:\n');
fprintf('  1. Fig1_Coefficient_Evolution.png          - Coefficient transitions\n');
fprintf('  2. Fig2_Data_Comparison_2D.png             - 2D original vs fitted\n');
fprintf('  3. Fig3_Error_Heatmap.png                  - Error with problem regions\n');
fprintf('  4. Fig4_Error_Analysis.png                 - Regional error analysis\n');
fprintf('  5. Fig5_Fit_Comparison.png                 - Fit quality at r values\n');
fprintf('  6. Fig6_Monotonicity_Asymptotic.png        - Behavior verification\n');
fprintf('  7. Fig7_Problem_Regions_Detail.png         - Problem regions detail\n');
fprintf('\n');

fprintf('Summary:\n');
fprintf('  Overall RMSE: %.6f\n', rmse_all);
fprintf('  Critical RMSE: %.6f\n', rmse_crit);
fprintf('  R-squared: %.6f\n', R2);
fprintf('  Control points: %d\n', n_control);
fprintf('\n');
