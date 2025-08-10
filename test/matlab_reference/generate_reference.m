% Generate MATLAB Kalman filter reference for comparison with Julia
% Uses same data and parameters as Julia test

clear; clc;

% Add paths
addpath('../../../../matlab/kf_functions');
addpath('../../../../matlab/AllanLab');
addpath('../../../../matlab/utilities');

fprintf('=== MATLAB Kalman Filter Reference ===\n');

% Load same data file
data_file = '../../../../data/6krb25apr.txt';
fprintf('Loading data from: %s\n', data_file);

% Load phase data - phase error is in column 2
data_matrix = readmatrix(data_file);
phase_data = data_matrix(:, 2);  % Phase error data in column 2

fprintf('Loaded %d samples\n', length(phase_data));
fprintf('Data range: [%.2e, %.2e] ns\n', min(phase_data), max(phase_data));

% Use same parameters as Julia default_filter
g_p = 0.1;
g_i = 0.01; 
g_d = 0.05;
q_wpm = 100.0;    % Measurement noise (R)
q_wfm = 0.01;     % White frequency noise
q_rwfm = 1e-6;    % Random walk frequency noise
q_irwfm = 0.0;    % No integrated RWFM
q_diurnal = 0.0;  % No diurnal
nstates = 3;      % 3-state model
tau = 1.0;        % 1 second sampling
start_cov = 1e6;  % Large initial covariance

fprintf('\nMATLAB Filter configuration:\n');
fprintf('  PID gains: g_p=%.3f, g_i=%.3f, g_d=%.3f\n', g_p, g_i, g_d);
fprintf('  Noise params: q_wpm=%.1f, q_wfm=%.3f, q_rwfm=%.1e\n', q_wpm, q_wfm, q_rwfm);
fprintf('  States: %d, tau=%.1f\n', nstates, tau);

% Initialize using first 100 points (manual least squares like Julia)
n_fit = 100;
t = (0:n_fit-1)' * tau;
y = phase_data(1:n_fit);  % Column vector, not transposed

% Design matrix for 3-state model: [1 t t²/2]
A = [ones(n_fit, 1), t, t.^2/2];

% Least squares fit
coeffs = A \ y;
x0 = coeffs;

% Covariance from residuals (matching Julia)
y_pred = A * coeffs;
residuals_fit = y - y_pred;
P0 = var(residuals_fit) * inv(A' * A);

fprintf('\nAfter initialization (first %d points):\n', n_fit);
fprintf('  Initial state x0: [%.6f, %.6f, %.6f]\n', x0(1), x0(2), x0(3));
fprintf('  P0 diagonal: [%.3e, %.3e, %.3e]\n', P0(1,1), P0(2,2), P0(3,3));

% Use remaining data (trimmed like Julia)
data_trimmed = phase_data(n_fit+1:end);
fprintf('  Remaining data points: %d\n', length(data_trimmed));

% Run Kalman filter
fprintf('\nRunning MATLAB Kalman filter...\n');
tic;
[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(data_trimmed, q_wfm, q_rwfm, q_wpm, ...
               g_p, g_i, g_d, nstates, tau, P0, x0, q_irwfm, q_diurnal, 86400);
matlab_time = toc;

fprintf('MATLAB filter completed in %.3f seconds\n', matlab_time);

% Calculate statistics
phase_rms = sqrt(mean(phase_est.^2));
freq_rms = sqrt(mean(freq_est.^2));
drift_rms = sqrt(mean(drift_est.^2));
residual_rms = sqrt(mean(residuals.^2));

fprintf('\nMATLAB Estimate statistics:\n');
fprintf('  Phase RMS: %.3e ns\n', phase_rms);
fprintf('  Frequency RMS: %.3e ns/s\n', freq_rms);
fprintf('  Drift RMS: %.3e ns/s²\n', drift_rms);
fprintf('  Residual RMS: %.3e ns\n', residual_rms);

% Save results for Julia comparison
results_matrix = [phase_est(:), freq_est(:), drift_est(:), residuals(:), steers(:)];
writematrix(results_matrix, 'kalman_results_matlab.txt', 'Delimiter', '\t');
fprintf('\nResults saved to kalman_results_matlab.txt\n');

% Also save steering results
steering_matrix = [sumsteers(:), sumsumsteers(:)];
writematrix(steering_matrix, 'steering_results_matlab.txt', 'Delimiter', '\t');

fprintf('✅ MATLAB reference generated successfully!\n');