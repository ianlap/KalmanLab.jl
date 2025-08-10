# Create comparison plots between MATLAB and Julia Kalman filter results

using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics
using Plots

println("=== Creating MATLAB vs Julia Comparison Plots ===")

# Load data
data_file = "/Users/ianlapinski/Desktop/masterclock-kflab/data/6krb25apr.txt"
data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]
println("Loaded $(length(phase_data)) samples")

# Create filter with same parameters
kf = KalmanFilter(0.1, 0.01, 0.05, 100.0, 0.01, 1e-6, 0.0, 0.0, 3, 1.0, zeros(3), 1e6)

# Initialize
data_trimmed = initialize!(kf, phase_data, 100)
println("Filter initialized")

# Use manageable subset for comparison
N_compare = 2000  # 2000 points for good visualization
test_data = data_trimmed[1:N_compare]
println("Using $(length(test_data)) samples for comparison")

# Run Julia filter
println("Running Julia Kalman filter...")
julia_results = KalmanLab.run(kf, test_data)

# Save Julia results for MATLAB comparison
julia_matrix = hcat(julia_results.phase_est, julia_results.freq_est, 
                   julia_results.drift_est, julia_results.residuals, julia_results.steers)
writedlm("test/julia_comparison_$(N_compare).txt", julia_matrix, '\t')

# Save the test data for MATLAB to use
writedlm("test/test_data_$(N_compare).txt", test_data, '\t')

# Save initialization parameters for MATLAB
init_params = [kf.x0; vec(kf.P0)]
writedlm("test/init_params_$(N_compare).txt", init_params, '\t')

println("Julia results and parameters saved for MATLAB comparison")

# Create MATLAB script for exact same test
matlab_script = """
% MATLAB comparison script - run identical test to Julia
clear; clc;

% Add paths
addpath('../../../../matlab/kf_functions');
addpath('../../../../matlab/AllanLab'); 
addpath('../../../../matlab/utilities');

% Load the exact same test data Julia used
test_data = readmatrix('test_data_$(N_compare).txt');
init_params = readmatrix('init_params_$(N_compare).txt');

% Extract initialization parameters (same as Julia)
x0 = init_params(1:3);
P0 = reshape(init_params(4:end), 3, 3);

fprintf('MATLAB running on %d samples\\n', length(test_data));
fprintf('Initial state: [%.6f, %.6f, %.6f]\\n', x0(1), x0(2), x0(3));

% Run MATLAB Kalman filter with identical parameters
[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(test_data, 0.01, 1e-6, 100.0, 0.1, 0.01, 0.05, 3, 1.0, P0, x0, 0.0, 0.0, 86400);

% Save results in same format as Julia
matlab_matrix = [phase_est(:), freq_est(:), drift_est(:), residuals(:), steers(:)];
writematrix(matlab_matrix, 'matlab_comparison_$(N_compare).txt', 'Delimiter', '\\t');

fprintf('MATLAB results saved\\n');
fprintf('Phase RMS: %.3e ns\\n', sqrt(mean(phase_est.^2)));
fprintf('Residual RMS: %.3e ns\\n', sqrt(mean(residuals.^2)));
"""

# Write MATLAB script
open("test/run_matlab_$(N_compare).m", "w") do f
    write(f, matlab_script)
end

println("MATLAB comparison script created: test/run_matlab_$(N_compare).m")
println("Run: cd test && matlab -nodisplay -r \"run_matlab_$(N_compare); exit\"")

# Check if MATLAB results already exist and compare
matlab_file = "test/matlab_comparison_$(N_compare).txt"
if isfile(matlab_file)
    println("\nMATLAB results found! Creating comparison plots...")
    
    matlab_matrix = readdlm(matlab_file, '\t', Float64)
    
    # Extract results
    j_phase = julia_results.phase_est
    j_freq = julia_results.freq_est
    j_drift = julia_results.drift_est
    j_residuals = julia_results.residuals
    
    m_phase = matlab_matrix[:, 1]
    m_freq = matlab_matrix[:, 2] 
    m_drift = matlab_matrix[:, 3]
    m_residuals = matlab_matrix[:, 4]
    
    # Calculate differences
    phase_diff = j_phase - m_phase
    freq_diff = j_freq - m_freq
    drift_diff = j_drift - m_drift
    residual_diff = j_residuals - m_residuals
    
    println("Comparison Statistics:")
    println("  Phase difference RMS: $(sqrt(mean(phase_diff.^2)):.3e) ns")
    println("  Max phase difference: $(maximum(abs.(phase_diff)):.3e) ns")
    println("  Freq difference RMS: $(sqrt(mean(freq_diff.^2)):.3e) ns/s")
    println("  Max freq difference: $(maximum(abs.(freq_diff)):.3e) ns/s")
    
    # Create comprehensive comparison plots
    sample_range = 1:min(500, length(j_phase))  # Show first 500 points clearly
    
    # Plot 1: State estimates comparison
    p1 = plot(sample_range, j_phase[sample_range], label="Julia Phase", color=:blue)
    plot!(p1, sample_range, m_phase[sample_range], label="MATLAB Phase", color=:red, linestyle=:dash)
    title!(p1, "Phase Estimates")
    ylabel!(p1, "Phase (ns)")
    
    p2 = plot(sample_range, j_freq[sample_range], label="Julia Freq", color=:blue)
    plot!(p2, sample_range, m_freq[sample_range], label="MATLAB Freq", color=:red, linestyle=:dash)
    title!(p2, "Frequency Estimates") 
    ylabel!(p2, "Frequency (ns/s)")
    
    # Plot 2: Differences
    p3 = plot(sample_range, phase_diff[sample_range], label="Phase Diff", color=:green)
    title!(p3, "Phase Difference (Julia - MATLAB)")
    ylabel!(p3, "Difference (ns)")
    
    p4 = plot(sample_range, freq_diff[sample_range], label="Freq Diff", color=:orange)
    title!(p4, "Frequency Difference (Julia - MATLAB)")
    ylabel!(p4, "Difference (ns/s)")
    xlabel!(p4, "Sample")
    
    # Combined plot
    comparison_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
    savefig(comparison_plot, "test/matlab_julia_comparison.png")
    
    # Plot 3: Full residuals comparison
    p5 = plot(1:length(j_residuals), j_residuals, label="Julia Residuals", alpha=0.7)
    plot!(p5, 1:length(m_residuals), m_residuals, label="MATLAB Residuals", alpha=0.7)
    title!(p5, "Measurement Residuals")
    xlabel!(p5, "Sample")
    ylabel!(p5, "Residuals (ns)")
    savefig(p5, "test/residuals_comparison.png")
    
    # Plot 4: Difference time series
    p6 = plot(1:length(phase_diff), phase_diff, title="Phase Difference Over Time")
    xlabel!(p6, "Sample")
    ylabel!(p6, "Julia - MATLAB (ns)")
    savefig(p6, "test/phase_difference_timeseries.png")
    
    println("Plots saved:")
    println("  - test/matlab_julia_comparison.png (main comparison)")
    println("  - test/residuals_comparison.png (residuals)")
    println("  - test/phase_difference_timeseries.png (differences)")
    
    # Summary statistics
    if maximum(abs.(phase_diff)) < 1e-10
        println("EXCELLENT: Phase estimates match to machine precision!")
    elseif maximum(abs.(phase_diff)) < 1e-6
        println("GOOD: Phase estimates match to microsecond precision")
    else
        println("WARNING: Significant differences detected")
    end
    
else
    println("\nMATLAB results not found.")
    println("To complete comparison:")
    println("1. cd test") 
    println("2. matlab -nodisplay -r \"run_matlab_$(N_compare); exit\"")
    println("3. Run this script again")
end