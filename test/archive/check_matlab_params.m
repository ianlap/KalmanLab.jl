% Check MATLAB initialization parameters 
addpath('../../../matlab/kf_functions');

% Load data exactly like Julia
data_matrix = readmatrix('../../../data/6krb25apr.txt');
phase_data = data_matrix(:, 2);

% Same initialization as Julia
n_fit = 100;
t = (0:n_fit-1)' * 1.0;
y = phase_data(1:n_fit);
A = [ones(n_fit, 1), t, t.^2/2];
coeffs = A \ y;
x0 = coeffs;
y_pred = A * coeffs;
residuals_fit = y - y_pred;
P0 = var(residuals_fit) * inv(A' * A);

fprintf('MATLAB Initialization Results:\n');
fprintf('x0 = [%.6f, %.6f, %.6f]\n', x0);
fprintf('P0 diagonal = [%.6e, %.6e, %.6e]\n', P0(1,1), P0(2,2), P0(3,3));
fprintf('P0 condition number = %.2e\n', cond(P0));
fprintf('P0 eigenvalues = [%.6e, %.6e, %.6e]\n', eig(P0));

% Test small run
data_trimmed = phase_data(n_fit+1:end);
test_data = data_trimmed(1:200);  % Small test

[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(test_data, 0.01, 1e-6, 100.0, 0.1, 0.01, 0.05, 3, 1.0, P0, x0, 0.0, 0.0, 86400);

fprintf('\nMATLAB Results (200 samples):\n');
fprintf('Phase RMS: %.3f ns\n', sqrt(mean(phase_est.^2)));
fprintf('Freq RMS: %.6f ns/s\n', sqrt(mean(freq_est.^2)));
fprintf('Final phase: %.6f ns\n', phase_est(end));
fprintf('Final freq: %.6f ns/s\n', freq_est(end));

% Check for NaN or Inf
if any(isnan(phase_est)) || any(isinf(phase_est))
    fprintf('WARNING: MATLAB has NaN or Inf values!\n');
else
    fprintf('MATLAB: No NaN/Inf detected\n');
end

% Save first 10 results for comparison
results_10 = [phase_est(1:10), freq_est(1:10)];
writematrix(results_10, 'matlab_debug_10.txt', 'Delimiter', '\t');