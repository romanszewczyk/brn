clear all
more off;
page_screen_output(0);
page_output_immediately(1);

pkg load optim

% Fit the coefficient functions P(r), A(r), B(r), Q(r), C(r) of the
% second-degree rational correction b(n,r) (Eq. 15 of the paper) to the
% Monte Carlo surface in rbn_MC_data.mat, by constrained differential
% evolution. See MATHEMATICAL_BACKGROUND.md for the model and constraints.

log_file = fopen('optimization_log.txt', 'w');
fprintf(log_file, 'b(n,r) surface fit\n');
fprintf(log_file, 'Started: %s\n\n', datestr(now));

fprintf('b(n,r) surface fit\n\n');

%% LOAD DATA
fprintf('=== LOADING DATA ===\n');
fprintf(log_file, '--- DATA LOADING ---\n');

load('rbn_MC_data.mat');

if ~exist('n', 'var') || ~exist('r', 'var') || ~exist('rbn_MC_data', 'var')
    error('Required variables not found in rbn_MC_data.mat');
end

if size(n, 1) > size(n, 2)
    n = n';
end

r_values = r(:);
rbn_mean = rbn_MC_data;

fprintf('Dataset loaded successfully:\n');
fprintf('  r: %d points (%.2f to %.2f)\n', length(r_values), min(r_values), max(r_values));
fprintf('  n: %d points (%d to %d)\n', length(n), min(n), max(n));
fprintf('  Data: %d x %d = %d points\n\n', size(rbn_mean, 1), size(rbn_mean, 2), numel(rbn_mean));

fprintf(log_file, 'Data dimensions: %d x %d = %d points\n', ...
        size(rbn_mean, 1), size(rbn_mean, 2), numel(rbn_mean));
fprintf(log_file, 'r range: %.2f to %.2f (%d points)\n', min(r_values), max(r_values), length(r_values));
fprintf(log_file, 'n range: %d to %d (%d points)\n\n', min(n), max(n), length(n));

[nn, rr] = meshgrid(double(n), r_values);

%% CONTROL POINTS
r_control_points = [0.2, 0.4, 0.6, 0.8, 1.0, ...
                    1.2, 1.5, 1.8, ...
                    2.0, 2.3, 2.6, 2.9, ...
                    3.0, 3.3, 3.6, 4.0, ...
                    4.2, 4.4, 4.6, 4.8, 5.0, ...
                    5.5, 6.0, 7.0, 8.0, 10.0];
n_ctrl = length(r_control_points);

fprintf('=== CONTROL POINTS ===\n');
fprintf('Control points (%d): ', n_ctrl);
fprintf('%.1f ', r_control_points);
fprintf('\n\n');

fprintf(log_file, '--- CONTROL POINTS ---\n');
fprintf(log_file, 'Number: %d\n\n', n_ctrl);

%% MEASURE ASYMPTOTIC VALUES
fprintf('=== MEASURING ASYMPTOTIC VALUES ===\n');
fprintf(log_file, '--- ASYMPTOTIC MEASUREMENT ---\n');

idx_asym = find(r_values >= 5.5);
A_asym_vals = zeros(length(idx_asym), 1);
B_asym_vals = zeros(length(idx_asym), 1);
C_asym_vals = zeros(length(idx_asym), 1);

for k = 1:length(idx_asym)
    idx_r = idx_asym(k);
    bn_tmp = rbn_mean(idx_r, :);
    
    target_lin = @(a) sum(((a(1)*n + a(2)) ./ (a(3)*n + 1) - bn_tmp).^2);
    a0 = [-1.0, 1.5, -0.6];
    
    ctl_tmp.XVmin = a0 * 0.5;
    ctl_tmp.XVmax = a0 * 1.5;
    for j = 1:3
        if ctl_tmp.XVmin(j) > ctl_tmp.XVmax(j)
            tmp = ctl_tmp.XVmin(j);
            ctl_tmp.XVmin(j) = ctl_tmp.XVmax(j);
            ctl_tmp.XVmax(j) = tmp;
        end
    end
    ctl_tmp.constr = 1;
    ctl_tmp.NP = 200;
    ctl_tmp.refresh = 0;
    ctl_tmp.tol = 1e-8;
    ctl_tmp.maxiter = 2000;
    
    [a_tmp, ~, ~, ~] = de_min(target_lin, ctl_tmp);
    A_asym_vals(k) = a_tmp(1);
    B_asym_vals(k) = a_tmp(2);
    C_asym_vals(k) = a_tmp(3);
end

A_inf = mean(A_asym_vals);
B_inf = mean(B_asym_vals);
C_inf = mean(C_asym_vals);
A_inf_std = std(A_asym_vals);
B_inf_std = std(B_asym_vals);
C_inf_std = std(C_asym_vals);

fprintf('Measured from r >= 5.5:\n');
fprintf('  A_inf = %.6f +/- %.6f\n', A_inf, A_inf_std);
fprintf('  B_inf = %.6f +/- %.6f\n', B_inf, B_inf_std);
fprintf('  C_inf = %.6f +/- %.6f\n\n', C_inf, C_inf_std);

fprintf(log_file, 'Asymptotic values (from r>=5.5):\n');
fprintf(log_file, '  A_inf = %.6f (std=%.6f)\n', A_inf, A_inf_std);
fprintf(log_file, '  B_inf = %.6f (std=%.6f)\n', B_inf, B_inf_std);
fprintf(log_file, '  C_inf = %.6f (std=%.6f)\n\n', C_inf, C_inf_std);

%% MEASURE HIGH-N ASYMPTOTIC VALUES
fprintf('=== MEASURING HIGH-N ASYMPTOTIC VALUES ===\n');
fprintf(log_file, '--- HIGH-N ASYMPTOTIC MEASUREMENT ---\n');

idx_high_n = find(n > 50);
idx_high_r = find(r_values > 5.5);

if ~isempty(idx_high_n) && ~isempty(idx_high_r)
    A_high_n_vals = zeros(length(idx_high_r), 1);
    C_high_n_vals = zeros(length(idx_high_r), 1);
    
    for k = 1:length(idx_high_r)
        idx_r = idx_high_r(k);
        bn_tmp = rbn_mean(idx_r, idx_high_n);
        n_tmp = n(idx_high_n);
        
        target_lin = @(a) sum(((a(1)*n_tmp + a(2)) ./ (a(3)*n_tmp + 1) - bn_tmp).^2);
        a0 = [-1.5, 1.5, -1.0];
        
        ctl_tmp.XVmin = a0 * 0.7;
        ctl_tmp.XVmax = a0 * 1.3;
        for j = 1:3
            if ctl_tmp.XVmin(j) > ctl_tmp.XVmax(j)
                tmp = ctl_tmp.XVmin(j);
                ctl_tmp.XVmin(j) = ctl_tmp.XVmax(j);
                ctl_tmp.XVmax(j) = tmp;
            end
        end
        ctl_tmp.constr = 1;
        ctl_tmp.NP = 200;
        ctl_tmp.refresh = 0;
        ctl_tmp.tol = 1e-8;
        ctl_tmp.maxiter = 2000;
        
        [a_tmp, ~, ~, ~] = de_min(target_lin, ctl_tmp);
        A_high_n_vals(k) = a_tmp(1);
        C_high_n_vals(k) = a_tmp(3);
    end
    
    A_high_n = mean(A_high_n_vals);
    C_high_n = mean(C_high_n_vals);
    A_high_n_std = std(A_high_n_vals);
    C_high_n_std = std(C_high_n_vals);
    
    fprintf('Measured from n>50, r>5.5:\n');
    fprintf('  A_high_n = %.6f +/- %.6f\n', A_high_n, A_high_n_std);
    fprintf('  C_high_n = %.6f +/- %.6f\n\n', C_high_n, C_high_n_std);
    
    fprintf(log_file, 'High-n asymptotic values (from n>50, r>5.5):\n');
    fprintf(log_file, '  A_high_n = %.6f (std=%.6f)\n', A_high_n, A_high_n_std);
    fprintf(log_file, '  C_high_n = %.6f (std=%.6f)\n\n', C_high_n, C_high_n_std);
else
    A_high_n = -1.5;
    C_high_n = -1.0;
    fprintf('High-n data not available, using default targets:\n');
    fprintf('  A_high_n = %.6f\n', A_high_n);
    fprintf('  C_high_n = %.6f\n\n', C_high_n);
    fprintf(log_file, 'High-n targets (default):\n');
    fprintf(log_file, '  A_high_n = %.6f\n', A_high_n);
    fprintf(log_file, '  C_high_n = %.6f\n\n', C_high_n);
end

%% MODEL FUNCTION
function res = rbn_model_monotonic(nn, rr, params, r_points)
    n_points = length(r_points);
    P_vals = params(1:n_points);
    A_vals = params(n_points+1:2*n_points);
    B_vals = params(2*n_points+1:3*n_points);
    Q_vals = params(3*n_points+1:4*n_points);
    C_vals = params(4*n_points+1:5*n_points);
    
    rr_clamped = max(min(rr, r_points(end)), r_points(1));
    
    P = pchip(r_points, P_vals, rr_clamped);
    A = pchip(r_points, A_vals, rr_clamped);
    B = pchip(r_points, B_vals, rr_clamped);
    Q = pchip(r_points, Q_vals, rr_clamped);
    C = pchip(r_points, C_vals, rr_clamped);
    
    res = (P.*nn.^2 + A.*nn + B) ./ (Q.*nn.^2 + C.*nn + 1);
end

%% REGULARIZATION PENALTIES
function penalty = compute_penalties(params, r_points, A_inf, B_inf, C_inf, A_high_n, C_high_n)
    n_points = length(r_points);
    P_vals = params(1:n_points);
    A_vals = params(n_points+1:2*n_points);
    B_vals = params(2*n_points+1:3*n_points);
    Q_vals = params(3*n_points+1:4*n_points);
    C_vals = params(4*n_points+1:5*n_points);
    
    penalty = 0;

    % 1. Monotonicity of P(r) and Q(r)
    for i = 1:n_points-1
        dp = P_vals(i+1) - P_vals(i);
        dq = Q_vals(i+1) - Q_vals(i);
        
        if dp > 0
            penalty = penalty + 150 * dp^2;
        end
        if dq > 0
            penalty = penalty + 150 * dq^2;
        end
    end
    
    % 2. Asymptotic enforcement for r>=5.5: P,Q -> 0 and A,B,C -> their
    %    asymptotic values, with weight increasing toward large r
    for i = 1:n_points
        r_val = r_points(i);
        if r_val >= 5.5
            weight_asym = 100 * exp(0.25 * (r_val - 5.5));
            
            penalty = penalty + weight_asym * (P_vals(i)^2 + Q_vals(i)^2);
            penalty = penalty + weight_asym * 0.05 * ((A_vals(i) - A_inf)^2 + ...
                                                      (B_vals(i) - B_inf)^2 + ...
                                                      (C_vals(i) - C_inf)^2);
        end
    end
    
    % 3. High-n asymptotic guidance for A(r) and C(r) at r>5.5
    for i = 1:n_points
        r_val = r_points(i);
        if r_val > 5.5
            weight_high_n = 20 * exp(0.15 * (r_val - 5.5));
            
            penalty = penalty + weight_high_n * ((A_vals(i) - A_high_n)^2 + ...
                                                 (C_vals(i) - C_high_n)^2);
        end
    end
    
    % 4. Second-derivative smoothness of all coefficient functions
    for i = 2:n_points-1
        dr1 = r_points(i) - r_points(i-1);
        dr2 = r_points(i+1) - r_points(i);
        dr_avg = (dr1 + dr2) / 2;
        
        for coeff_idx = 0:4
            vals = params(coeff_idx*n_points+[i-1, i, i+1]);
            d2 = ((vals(3)-vals(2))/dr2 - (vals(2)-vals(1))/dr1) / dr_avg;
            penalty = penalty + 0.005 * d2^2;
        end
    end
    
    % 5. Low-r stability of the first three control points
    for coeff_idx = 0:4
        vals = params(coeff_idx*n_points+1:coeff_idx*n_points+3);
        penalty = penalty + 2.0 * sum((vals - mean(vals)).^2);
    end
    
    % 6. Extra smoothness across the 4 <= r <= 5.5 transition region
    idx_prob = find((r_points >= 4.0) & (r_points <= 5.5));
    if length(idx_prob) >= 3
        for i = 2:length(idx_prob)-1
            idx_i = idx_prob(i);
            for coeff_idx = 0:4
                vals = params(coeff_idx*n_points+[idx_i-1, idx_i, idx_i+1]);
                penalty = penalty + 0.02 * sum(diff(vals).^2);
            end
        end
    end
    
    % 7. Transition-zone smoothness (5.0 <= r <= 5.5)
    idx_trans = find((r_points >= 5.0) & (r_points <= 5.5));
    if length(idx_trans) >= 2
        for i = 1:length(idx_trans)-1
            idx_i = idx_trans(i);
            for coeff_idx = 0:4
                vals = params(coeff_idx*n_points+[idx_i, idx_i+1]);
                penalty = penalty + 0.05 * (vals(2) - vals(1))^2;
            end
        end
    end
end

%% ADAPTIVE REGIONAL WEIGHTS
fprintf('=== BUILDING ADAPTIVE REGIONAL WEIGHTS ===\n');
fprintf(log_file, '--- ADAPTIVE REGIONAL WEIGHTING ---\n');

weight_matrix = ones(size(nn));

idx_critical = (rr < 2.0) & (nn <= 4);
weight_matrix(idx_critical) = 5.0;
n_critical = sum(idx_critical(:));

idx_r_low = (rr < 2.0) & (nn > 4);
weight_matrix(idx_r_low) = 2.0;
n_r_low = sum(idx_r_low(:));

idx_n_low = (rr >= 2.0) & (rr < 4.0) & (nn <= 5);
weight_matrix(idx_n_low) = 3.0;
n_n_low = sum(idx_n_low(:));

idx_problem_45 = (rr > 4.0) & (rr < 5.0) & ((nn < 10) | (nn > 60));
weight_matrix(idx_problem_45) = 8.0;
n_problem_45 = sum(idx_problem_45(:));

idx_high_nr = (nn > 60) & (rr >= 7.5) & (rr <= 8.5);
weight_matrix(idx_high_nr) = 5.0;
n_high_nr = sum(idx_high_nr(:));

idx_extended_crit = (rr >= 2.0) & (rr < 4.0) & (nn <= 4);
weight_matrix(idx_extended_crit) = 6.0;
n_extended_crit = sum(idx_extended_crit(:));

idx_saturated = (rr >= 5.5);
weight_matrix(idx_saturated) = 0.3;
n_saturated = sum(idx_saturated(:));

fprintf('Weight distribution:\n');
fprintf('  Critical (r<2, n<=4):         %4d points, weight=5.0\n', n_critical);
fprintf('  Extended crit (2<=r<4, n<=4): %4d points, weight=6.0\n', n_extended_crit);
fprintf('  Problem 4-5 (extreme n):      %4d points, weight=8.0\n', n_problem_45);
fprintf('  High n,r (n>60, r~8):         %4d points, weight=5.0\n', n_high_nr);
fprintf('  Low-r (r<2, n>4):             %4d points, weight=2.0\n', n_r_low);
fprintf('  Low-n (2<=r<4, n<=5):         %4d points, weight=3.0\n', n_n_low);
fprintf('  Saturated (r>=5.5):           %4d points, weight=0.3\n', n_saturated);
fprintf('  Normal:                       %4d points, weight=1.0\n\n', ...
        numel(weight_matrix) - n_critical - n_r_low - n_n_low - n_saturated - ...
        n_problem_45 - n_high_nr - n_extended_crit);

fprintf(log_file, 'Weight distribution:\n');
fprintf(log_file, '  Critical region: %d points (5.0x)\n', n_critical);
fprintf(log_file, '  Extended critical: %d points (6.0x)\n', n_extended_crit);
fprintf(log_file, '  Problem 4-5 region: %d points (8.0x)\n', n_problem_45);
fprintf(log_file, '  High-n,r region: %d points (5.0x)\n', n_high_nr);
fprintf(log_file, '  Low-r region: %d points (2.0x)\n', n_r_low);
fprintf(log_file, '  Low-n region: %d points (3.0x)\n', n_n_low);
fprintf(log_file, '  Saturated region: %d points (0.3x)\n', n_saturated);
fprintf(log_file, '  Normal: %d points (1.0x)\n\n', ...
        numel(weight_matrix) - n_critical - n_r_low - n_n_low - n_saturated - ...
        n_problem_45 - n_high_nr - n_extended_crit);

%% INITIALIZATION
fprintf('=== INITIALIZATION AT CONTROL POINTS ===\n');
fprintf(log_file, '--- INITIALIZATION ---\n');

P_init = zeros(n_ctrl, 1);
A_init = zeros(n_ctrl, 1);
B_init = zeros(n_ctrl, 1);
Q_init = zeros(n_ctrl, 1);
C_init = zeros(n_ctrl, 1);

for k = 1:n_ctrl
    [~, idx_r] = min(abs(r_values - r_control_points(k)));
    bn_mean_tmp = rbn_mean(idx_r, :);
    
    if r_control_points(k) >= 5.5
        P_init(k) = 0.001;
        A_init(k) = A_inf;
        B_init(k) = B_inf;
        Q_init(k) = 0.001;
        C_init(k) = C_inf;
    elseif r_control_points(k) >= 5.0
        % Transition zone - blend between quadratic and asymptotic
        alpha = (r_control_points(k) - 5.0) / 0.5;
        target_lin = @(a) sum(((a(1)*n + a(2)) ./ (a(3)*n + 1) - bn_mean_tmp).^2);
        a0 = [A_high_n, 1.4, C_high_n];
        
        ctl_tmp.XVmin = a0 * 0.6;
        ctl_tmp.XVmax = a0 * 1.4;
        for j = 1:3
            if ctl_tmp.XVmin(j) > ctl_tmp.XVmax(j)
                tmp = ctl_tmp.XVmin(j);
                ctl_tmp.XVmin(j) = ctl_tmp.XVmax(j);
                ctl_tmp.XVmax(j) = tmp;
            end
        end
        ctl_tmp.constr = 1;
        ctl_tmp.NP = 300;
        ctl_tmp.refresh = 0;
        ctl_tmp.tol = 1e-8;
        ctl_tmp.maxiter = 3000;
        
        [a_tmp, ~, ~, ~] = de_min(target_lin, ctl_tmp);
        
        P_init(k) = 0.02 * (1 - alpha);
        A_init(k) = a_tmp(1);
        B_init(k) = a_tmp(2);
        Q_init(k) = 0.02 * (1 - alpha);
        C_init(k) = a_tmp(3);
    else
        target_quad = @(a) sum(((a(1)*n.^2 + a(2)*n + a(3)) ./ ...
                                 (a(4)*n.^2 + a(5)*n + 1) - bn_mean_tmp).^2);
        a0 = [0.8, -0.3, 1.3, 0.7, -0.5];
        
        ctl_tmp.XVmin = a0 * 0.3;
        ctl_tmp.XVmax = a0 * 1.7;
        for j = 1:5
            if ctl_tmp.XVmin(j) > ctl_tmp.XVmax(j)
                tmp = ctl_tmp.XVmin(j);
                ctl_tmp.XVmin(j) = ctl_tmp.XVmax(j);
                ctl_tmp.XVmax(j) = tmp;
            end
        end
        ctl_tmp.constr = 1;
        ctl_tmp.NP = 400;
        ctl_tmp.refresh = 0;
        ctl_tmp.tol = 1e-8;
        ctl_tmp.maxiter = 3000;
        
        [a_tmp, ~, ~, ~] = de_min(target_quad, ctl_tmp);
        
        P_init(k) = a_tmp(1);
        A_init(k) = a_tmp(2);
        B_init(k) = a_tmp(3);
        Q_init(k) = a_tmp(4);
        C_init(k) = a_tmp(5);
    end
    
    fprintf('  r=%4.1f: P=%7.4f A=%7.4f B=%7.4f Q=%7.4f C=%7.4f\n', ...
            r_control_points(k), P_init(k), A_init(k), B_init(k), Q_init(k), C_init(k));
end

params_init = [P_init(:); A_init(:); B_init(:); Q_init(:); C_init(:)];
fprintf('Initialization complete.\n\n');

fprintf(log_file, 'Sample initialization (first 5 control points):\n');
for k = 1:5
    fprintf(log_file, '  r=%.1f: P=%.4f A=%.4f B=%.4f Q=%.4f C=%.4f\n', ...
            r_control_points(k), P_init(k), A_init(k), B_init(k), Q_init(k), C_init(k));
end
fprintf(log_file, '\n');

%% MULTI-STAGE DIFFERENTIAL EVOLUTION
% Three stages of de_min with progressively tighter bounds around the
% incumbent: stage 1 fits the weighted data only; stages 2 and 3 add the
% regularization penalties and refine.
fprintf('================================================================\n');
fprintf('   MULTI-STAGE DIFFERENTIAL EVOLUTION\n');
fprintf('================================================================\n\n');
fprintf(log_file, '--- OPTIMIZATION STAGES ---\n\n');

% Stage 1: data fit only
fprintf('=== STAGE 1: data fit ===\n');
fprintf('NP=2500, MaxIter=25000\n');
fprintf(log_file, 'STAGE 1: data fit\n');
fflush(stdout);

target_s1 = @(params) sum(sum(((rbn_model_monotonic(nn, rr, params, r_control_points) ...
                                 - rbn_mean).^2) .* weight_matrix));

delta = 0.5;
params_min = params_init .* (1 - delta);
params_max = params_init .* (1 + delta);
for j = 1:numel(params_min)
    if params_min(j) > params_max(j)
        tmp = params_min(j);
        params_min(j) = params_max(j);
        params_max(j) = tmp;
    end
end

ctl.XVmin = params_min(:)';
ctl.XVmax = params_max(:)';
ctl.constr = 1;
ctl.NP = 2500;
ctl.maxiter = 25000;
ctl.tol = 1e-6;
ctl.refresh = 10;

tic;
[params_s1, obj1, nf1, ~] = de_min(target_s1, ctl);
time1 = toc;
fprintf('Stage 1: %.1f sec (%.1f min), obj=%.6e, nfeval=%d\n\n', time1, time1/60, obj1, nf1);
fprintf(log_file, '  Time: %.1f sec (%.1f min)\n', time1, time1/60);
fprintf(log_file, '  Objective: %.6e\n', obj1);
fprintf(log_file, '  Evaluations: %d\n\n', nf1);
fflush(stdout);

% Stage 2: add regularization penalties
fprintf('=== STAGE 2: data fit + penalties ===\n');
fprintf('NP=2000, MaxIter=20000\n');
fprintf(log_file, 'STAGE 2: data fit + penalties\n');
fflush(stdout);

target_s2 = @(params) sum(sum(((rbn_model_monotonic(nn, rr, params, r_control_points) ...
                                 - rbn_mean).^2) .* weight_matrix)) + ...
                          compute_penalties(params, r_control_points, A_inf, B_inf, C_inf, ...
                                          A_high_n, C_high_n);

delta = 0.3;
params_min = params_s1 .* (1 - delta);
params_max = params_s1 .* (1 + delta);
for j = 1:numel(params_min)
    if params_min(j) > params_max(j)
        tmp = params_min(j);
        params_min(j) = params_max(j);
        params_max(j) = tmp;
    end
end

ctl.XVmin = params_min(:)';
ctl.XVmax = params_max(:)';
ctl.NP = 2000;
ctl.maxiter = 20000;
ctl.tol = 1e-7;

tic;
[params_s2, obj2, nf2, ~] = de_min(target_s2, ctl);
time2 = toc;
fprintf('Stage 2: %.1f sec (%.1f min), obj=%.6e, nfeval=%d\n\n', time2, time2/60, obj2, nf2);
fprintf(log_file, '  Time: %.1f sec (%.1f min)\n', time2, time2/60);
fprintf(log_file, '  Objective: %.6e\n', obj2);
fprintf(log_file, '  Evaluations: %d\n\n', nf2);
fflush(stdout);

% Stage 3: final polish
fprintf('=== STAGE 3: final polish ===\n');
fprintf('NP=1500, MaxIter=15000\n');
fprintf(log_file, 'STAGE 3: final polish\n');
fflush(stdout);

delta = 0.15;
params_min = params_s2 .* (1 - delta);
params_max = params_s2 .* (1 + delta);
for j = 1:numel(params_min)
    if params_min(j) > params_max(j)
        tmp = params_min(j);
        params_min(j) = params_max(j);
        params_max(j) = tmp;
    end
end

ctl.XVmin = params_min(:)';
ctl.XVmax = params_max(:)';
ctl.NP = 1500;
ctl.maxiter = 15000;
ctl.tol = 1e-8;

tic;
[a_res, obj_final, nf3, ~] = de_min(target_s2, ctl);
time3 = toc;
fprintf('Stage 3: %.1f sec (%.1f min), obj=%.6e, nfeval=%d\n\n', time3, time3/60, obj_final, nf3);
fprintf(log_file, '  Time: %.1f sec (%.1f min)\n', time3, time3/60);
fprintf(log_file, '  Objective: %.6e\n', obj_final);
fprintf(log_file, '  Evaluations: %d\n\n', nf3);

total_time = time1 + time2 + time3;
fprintf('=== TOTAL TIME: %.1f sec (%.1f min / %.1f hours) ===\n\n', ...
        total_time, total_time/60, total_time/3600);
fprintf(log_file, 'Total time: %.1f sec (%.1f min / %.1f hours)\n\n', ...
        total_time, total_time/60, total_time/3600);
fflush(stdout);

rbn_calc = rbn_model_monotonic(nn, rr, a_res, r_control_points);

%% RESULTS ANALYSIS
fprintf('================================================================\n');
fprintf('   RESULTS ANALYSIS\n');
fprintf('================================================================\n\n');
fprintf(log_file, '--- RESULTS ANALYSIS ---\n\n');

err_all = rbn_calc - rbn_mean;
rmse_all = sqrt(mean(err_all(:).^2));
mae_all = mean(abs(err_all(:)));
max_err = max(abs(err_all(:)));
R2 = 1 - sum(err_all(:).^2) / sum((rbn_mean(:) - mean(rbn_mean(:))).^2);

fprintf('=== OVERALL FIT QUALITY ===\n');
fprintf('RMSE:      %.6f\n', rmse_all);
fprintf('MAE:       %.6f\n', mae_all);
fprintf('Max Error: %.6f\n', max_err);
fprintf('R-squared: %.6f\n\n', R2);

fprintf(log_file, 'OVERALL FIT QUALITY:\n');
fprintf(log_file, '  RMSE:      %.6f\n', rmse_all);
fprintf(log_file, '  MAE:       %.6f\n', mae_all);
fprintf(log_file, '  Max Error: %.6f\n', max_err);
fprintf(log_file, '  R-squared: %.6f\n\n', R2);

% Extract coefficients
P_vals = a_res(1:n_ctrl);
A_vals = a_res(n_ctrl+1:2*n_ctrl);
B_vals = a_res(2*n_ctrl+1:3*n_ctrl);
Q_vals = a_res(3*n_ctrl+1:4*n_ctrl);
C_vals = a_res(4*n_ctrl+1:5*n_ctrl);

P_vals = P_vals(:);
A_vals = A_vals(:);
B_vals = B_vals(:);
Q_vals = Q_vals(:);
C_vals = C_vals(:);

% Regional errors
fprintf('=== REGIONAL ERRORS ===\n');
fprintf(log_file, 'REGIONAL ERRORS:\n');

idx_crit = (rr < 2.0) & (nn <= 4);
rmse_crit = sqrt(mean(err_all(idx_crit).^2));
fprintf('Critical (r<2, n<=4):  RMSE=%.6f\n', rmse_crit);
fprintf(log_file, '  Critical (r<2, n<=4):  %.6f\n', rmse_crit);

idx_prob_r = find(r_values > 4.0 & r_values < 5.0);
idx_n_low = find(n <= 10);
idx_n_high = find(n > 60);
if ~isempty(idx_prob_r) && (~isempty(idx_n_low) || ~isempty(idx_n_high))
    idx_prob_n = [idx_n_low, idx_n_high];
    rmse_prob = sqrt(mean(mean((err_all(idx_prob_r, idx_prob_n)).^2)));
    fprintf('Problem 4-5 (extreme n): RMSE=%.6f\n', rmse_prob);
    fprintf(log_file, '  Problem 4-5 (extreme n): %.6f\n', rmse_prob);
end

idx_sat = (rr >= 5.5);
rmse_sat = sqrt(mean(err_all(idx_sat).^2));
fprintf('Saturated (r>=5.5):    RMSE=%.6f\n\n', rmse_sat);
fprintf(log_file, '  Saturated (r>=5.5):    %.6f\n\n', rmse_sat);

% Monotonicity
dP = diff(P_vals);
dQ = diff(Q_vals);
n_inc_P = sum(dP > 1e-6);
n_inc_Q = sum(dQ > 1e-6);

fprintf('=== MONOTONICITY ===\n');
fprintf('P increases: %d\n', n_inc_P);
fprintf('Q increases: %d\n\n', n_inc_Q);
fprintf(log_file, 'MONOTONICITY:\n');
fprintf(log_file, '  P increases: %d\n', n_inc_P);
fprintf(log_file, '  Q increases: %d\n\n', n_inc_Q);

fprintf('=== PARAMETERS AT KEY POINTS ===\n');
fprintf(log_file, 'PARAMETERS AT KEY POINTS:\n');
for i = 1:n_ctrl
    fprintf('r=%5.1f:  P=%9.6f  A=%9.6f  B=%9.6f  Q=%9.6f  C=%9.6f\n', ...
            r_control_points(i), P_vals(i), A_vals(i), B_vals(i), Q_vals(i), C_vals(i));
    fprintf(log_file, '  r=%5.1f:  P=%9.6f  A=%9.6f  B=%9.6f  Q=%9.6f  C=%9.6f\n', ...
            r_control_points(i), P_vals(i), A_vals(i), B_vals(i), Q_vals(i), C_vals(i));
end
fprintf('\n');
fprintf(log_file, '\n');

%% SAVE RESULTS
save('-v7', 'Results_surface_fit.mat', 'a_res', 'nn', 'rr', 'rbn_calc', ...
     'rbn_mean', 'r_control_points', 'weight_matrix', 'A_inf', 'B_inf', 'C_inf', ...
     'A_high_n', 'C_high_n', 'rmse_all', 'rmse_crit', 'R2', 'P_vals', 'A_vals', ...
     'B_vals', 'Q_vals', 'C_vals');

fprintf('================================================================\n');
fprintf('   DONE\n');
fprintf('================================================================\n');
fprintf('Results saved to: Results_surface_fit.mat\n');
fprintf('Log file: optimization_log.txt\n\n');

fprintf(log_file, '\nCompleted: %s\n', datestr(now));
fclose(log_file);

fprintf('Run analyze_results_PLOTS.m for the coefficient and error figures.\n\n');
