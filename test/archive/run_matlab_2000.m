% MATLAB comparison script - run identical test to Julia
clear; clc;

% Add paths
addpath('../../../../matlab/kf_functions');
addpath('../../../../matlab/AllanLab'); 
addpath('../../../../matlab/utilities');

% Load the exact same test data Julia used
test_data = readmatrix('test_data_2000.txt');
init_params = readmatrix('init_params_2000.txt');

% Extract initialization parameters (same as Julia)
x0 = init_params(1:3);
P0 = reshape(init_params(4:end), 3, 3);

fprintf('MATLAB running on %d samples\n', length(test_data));
fprintf('Initial state: [%.6f, %.6f, %.6f]\n', x0(1), x0(2), x0(3));

% Run MATLAB Kalman filter with identical parameters
[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(test_data, 0.01, 1e-6, 100.0, 0.1, 0.01, 0.05, 3, 1.0, P0, x0, 0.0, 0.0, 86400);

% Save results in same format as Julia
matlab_matrix = [phase_est(:), freq_est(:), drift_est(:), residuals(:), steers(:)];
writematrix(matlab_matrix, 'matlab_comparison_2000.txt', 'Delimiter', '\t');

fprintf('MATLAB results saved\n');
fprintf('Phase RMS: %.3e ns\n', sqrt(mean(phase_est.^2)));
fprintf('Residual RMS: %.3e ns\n', sqrt(mean(residuals.^2)));
