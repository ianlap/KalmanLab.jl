addpath('../../../../matlab/kf_functions');
addpath('../../../../matlab/AllanLab');
addpath('../../../../matlab/utilities');

% Load and initialize exactly like Julia
data_matrix = readmatrix('../../../../data/6krb25apr.txt');
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

% Use same subset as Julia
data_trimmed = phase_data(n_fit+1:end);
test_data = data_trimmed(1:1000);

% Run MATLAB filter on subset
[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(test_data, 0.01, 1e-6, 100.0, 0.1, 0.01, 0.05, 3, 1.0, P0, x0, 0.0, 0.0, 86400);

% Save results  
results_matrix = [phase_est(:), freq_est(:), drift_est(:), residuals(:), steers(:)];
writematrix(results_matrix, 'kalman_results_matlab_comparison.txt', 'Delimiter', '	');

fprintf('MATLAB completed - Phase RMS: %.3e, Residual RMS: %.3e\n', ...
        sqrt(mean(phase_est.^2)), sqrt(mean(residuals.^2)));
