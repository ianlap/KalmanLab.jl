% MATLAB vs Julia Full Implementation Comparison
% Using the stable PID gains and noise parameters

addpath('../../../matlab/kf_functions');

fprintf('=== MATLAB vs Julia Full Implementation Comparison ===\n');

% Load data
data_matrix = readmatrix('../../../data/6krb25apr.txt');
phase_data = data_matrix(:, 2);
fprintf('Loaded %d samples\n', length(phase_data));

% Use the same stable parameters as Julia
T = 1.0;
tau0 = 1.0;

% Calculate PID gains using time constant (same formulas as Julia)
a = exp(-tau0 / T);
g_p = 1 - 3*a^2 + 2*a^3;
g_i = 1 - 3*a + 3*a^2 - a^3;
g_d = 1 - a^3;

fprintf('\nPID gains from time constant T=%.1f:\n', T);
fprintf('  a = exp(-tau/T) = %.6f\n', a);
fprintf('  g_p = %.6f\n', g_p);
fprintf('  g_i = %.6f\n', g_i);
fprintf('  g_d = %.6f\n', g_d);

% Stable noise parameters (matching Julia exactly)
q_wpm = 100.0;    % Measurement noise (R)
q_wfm = 6e-3;     % White frequency modulation
q_rwfm = 5e-9;    % Random walk frequency modulation
q_irwfm = 0.0;    % No integrated RWFM
q_diurnal = 0.0;  % No diurnal terms
period = 86400;   % Daily period (not used with q_diurnal=0)

fprintf('\nNoise parameters:\n');
fprintf('  q_wpm (R): %.1f\n', q_wpm);
fprintf('  q_wfm: %.3e\n', q_wfm);
fprintf('  q_rwfm: %.3e\n', q_rwfm);

% Conservative initialization (matching Julia)
nparams = 3;
init_state = [0.0; 0.0; 0.0];  % Zero initial state
start_cov = 1e6;               % Conservative initial covariance

fprintf('\nInitialization:\n');
fprintf('  Initial state: [%.1f, %.1f, %.1f]\n', init_state);
fprintf('  Initial covariance: %.0e * I\n', start_cov);

% Run MATLAB filter on full dataset
fprintf('\nRunning MATLAB Kalman filter on full dataset...\n');
fprintf('This may take a moment...\n');

tic;
[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(phase_data, q_wfm, q_rwfm, q_wpm, g_p, g_i, g_d, ...
               nparams, tau0, start_cov, init_state, q_irwfm, q_diurnal, period);
matlab_time = toc;

fprintf('MATLAB processing time: %.2f seconds\n', matlab_time);
fprintf('MATLAB rate: %.0f samples/second\n', length(phase_data)/matlab_time);

% Calculate MATLAB metrics
matlab_phase_rms = sqrt(mean(phase_est.^2));
matlab_freq_rms = sqrt(mean(freq_est.^2));
matlab_residual_rms = sqrt(mean(residuals.^2));

fprintf('\n=== MATLAB Results ===\n');
fprintf('Phase RMS: %.6f ns\n', matlab_phase_rms);
fprintf('Freq RMS: %.8f ns/s\n', matlab_freq_rms);
fprintf('Residual RMS: %.3f ns\n', matlab_residual_rms);
fprintf('Final phase: %.6f ns\n', phase_est(end));
fprintf('Final freq: %.8f ns/s\n', freq_est(end));

% Save MATLAB results
matlab_results = [phase_est, freq_est, drift_est, residuals, steers];

% Save sampled results (every 100th sample to match Julia)
sample_indices = 1:100:length(phase_est);
matlab_sampled = matlab_results(sample_indices, :);
writematrix(matlab_sampled, 'matlab_full_sampled.txt', 'Delimiter', '\t');

fprintf('Saved MATLAB sampled results to matlab_full_sampled.txt\n');

% Save summary
matlab_summary_file = 'matlab_full_summary.txt';
fid = fopen(matlab_summary_file, 'w');
fprintf(fid, 'MATLAB Full Dataset Results\n');
fprintf(fid, 'Samples: %d\n', length(phase_data));
fprintf(fid, 'Processing time: %.2f s\n', matlab_time);
fprintf(fid, 'Rate: %.0f samples/s\n', length(phase_data)/matlab_time);
fprintf(fid, 'Phase RMS: %.6f ns\n', matlab_phase_rms);
fprintf(fid, 'Freq RMS: %.8f ns/s\n', matlab_freq_rms);
fprintf(fid, 'Residual RMS: %.3f ns\n', matlab_residual_rms);
fprintf(fid, 'Final phase: %.6f ns\n', phase_est(end));
fprintf(fid, 'Final freq: %.8f ns/s\n', freq_est(end));
fclose(fid);

fprintf('Saved MATLAB summary to %s\n', matlab_summary_file);

% Load Julia results for comparison
if exist('julia_full_summary.txt', 'file')
    fprintf('\n=== Comparing with Julia Results ===\n');
    
    % Read Julia summary (this is a bit manual but works)
    julia_data = readmatrix('julia_full_sampled.txt');
    julia_phase = julia_data(:, 1);
    julia_freq = julia_data(:, 2);
    julia_residuals = julia_data(:, 4);
    
    julia_phase_rms = sqrt(mean(julia_phase.^2));
    julia_freq_rms = sqrt(mean(julia_freq.^2));
    julia_residual_rms = sqrt(mean(julia_residuals.^2));
    
    fprintf('Julia Phase RMS: %.6f ns\n', julia_phase_rms);
    fprintf('Julia Freq RMS: %.8f ns/s\n', julia_freq_rms);
    fprintf('Julia Residual RMS: %.3f ns\n', julia_residual_rms);
    
    % Compare sampled results
    n_compare = min(size(matlab_sampled, 1), size(julia_data, 1));
    
    phase_diff = matlab_sampled(1:n_compare, 1) - julia_phase(1:n_compare);
    freq_diff = matlab_sampled(1:n_compare, 2) - julia_freq(1:n_compare);
    residual_diff = matlab_sampled(1:n_compare, 4) - julia_residuals(1:n_compare);
    
    fprintf('\n=== Differences (MATLAB - Julia) ===\n');
    fprintf('Phase differences:\n');
    fprintf('  RMS: %.8f ns\n', sqrt(mean(phase_diff.^2)));
    fprintf('  Max: %.8f ns\n', max(abs(phase_diff)));
    fprintf('  Mean: %.8f ns\n', mean(phase_diff));
    
    fprintf('Frequency differences:\n');
    fprintf('  RMS: %.10f ns/s\n', sqrt(mean(freq_diff.^2)));
    fprintf('  Max: %.10f ns/s\n', max(abs(freq_diff)));
    fprintf('  Mean: %.10f ns/s\n', mean(freq_diff));
    
    % Save detailed comparison
    comparison_results = [matlab_sampled(1:n_compare, 1:2), julia_data(1:n_compare, 1:2), ...
                         phase_diff, freq_diff];
    comparison_header = {'MATLAB_Phase', 'MATLAB_Freq', 'Julia_Phase', 'Julia_Freq', ...
                        'Phase_Diff', 'Freq_Diff'};
    
    writematrix(comparison_results, 'matlab_julia_comparison.txt', 'Delimiter', '\t');
    fprintf('Saved detailed comparison to matlab_julia_comparison.txt\n');
    
    % Check if differences are acceptable
    if sqrt(mean(phase_diff.^2)) < 0.1 && sqrt(mean(freq_diff.^2)) < 0.001
        fprintf('\n✅ EXCELLENT AGREEMENT between MATLAB and Julia!\n');
    else
        fprintf('\n⚠️  Significant differences detected between implementations\n');
    end
    
else
    fprintf('\nJulia results not found. Run Julia test first.\n');
end

fprintf('\n🎉 MATLAB COMPARISON COMPLETE! 🎉\n');