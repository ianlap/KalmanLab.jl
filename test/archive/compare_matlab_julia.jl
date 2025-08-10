# Compare MATLAB vs Julia Kalman filter results
# Run on same subset to avoid the NaN issue from MATLAB full dataset

using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics
using Plots  # For plotting differences

println("=== MATLAB vs Julia Kalman Filter Comparison ===")

# Load data
data_file = "/Users/ianlapinski/Desktop/masterclock-kflab/data/6krb25apr.txt"
data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]
println("Loaded $(length(phase_data)) samples")

# Use same parameters for both
kf = KalmanFilter(0.1, 0.01, 0.05, 100.0, 0.01, 1e-6, 0.0, 0.0, 3, 1.0, zeros(3), 1e6)

# Initialize
data_trimmed = initialize!(kf, phase_data, 100)
println("Filter initialized: x0 = $(round.(kf.x0, digits=6))")

# Run on smaller subset to match what works
test_size = 1000  # Use 1000 samples for comparison
test_data = data_trimmed[1:test_size]
println("Testing on $(length(test_data)) samples to avoid MATLAB NaN issue")

# Run Julia filter
println("\nRunning Julia filter...")
julia_results = KalmanLab.run(kf, test_data)

# Save Julia results
julia_matrix = hcat(julia_results.phase_est, julia_results.freq_est, 
                   julia_results.drift_est, julia_results.residuals, 
                   julia_results.steers)
writedlm("test/kalman_results_julia_comparison.txt", julia_matrix, '\t')

println("Julia results saved")
println("  Phase RMS: $(round(sqrt(mean(julia_results.phase_est.^2)), digits=3)) ns")
println("  Residual RMS: $(round(sqrt(mean(julia_results.residuals.^2)), digits=3)) ns")
println("  Final state: $(round.(julia_results.final_state, digits=6))")

# Generate MATLAB results for same subset
println("\nGenerating MATLAB results for comparison...")

# Create MATLAB script for this subset
matlab_script = """
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
coeffs = A \\ y;
x0 = coeffs;
y_pred = A * coeffs;
residuals_fit = y - y_pred;
P0 = var(residuals_fit) * inv(A' * A);

% Use same subset as Julia
data_trimmed = phase_data(n_fit+1:end);
test_data = data_trimmed(1:$test_size);

% Run MATLAB filter on subset
[phase_est, freq_est, drift_est, residuals, innovations, steers, ...
 rtP00, rtP11, rtP22, rtP01, rtP02, rtP12, sumsteers, sumsumsteers] = ...
 kalman_filter(test_data, 0.01, 1e-6, 100.0, 0.1, 0.01, 0.05, 3, 1.0, P0, x0, 0.0, 0.0, 86400);

% Save results  
results_matrix = [phase_est(:), freq_est(:), drift_est(:), residuals(:), steers(:)];
writematrix(results_matrix, 'kalman_results_matlab_comparison.txt', 'Delimiter', '\t');

fprintf('MATLAB completed - Phase RMS: %.3e, Residual RMS: %.3e\\n', ...
        sqrt(mean(phase_est.^2)), sqrt(mean(residuals.^2)));
"""

# Write and run MATLAB script
open("test/run_matlab_comparison.m", "w") do f
    write(f, matlab_script)
end

println("MATLAB script created. Run it manually to generate comparison data.")
println("\nTo complete the comparison:")
println("1. cd test")
println("2. matlab -nodisplay -r \"run_matlab_comparison; exit\"")
println("3. Then run this script again to compare results")

# Try to load MATLAB results if they exist
matlab_file = "test/kalman_results_matlab_comparison.txt"
if isfile(matlab_file)
    println("\nMATLAB results found - comparing...")
    matlab_results = readdlm(matlab_file, '\t', Float64)
    
    # Compare results
    phase_diff = julia_matrix[:, 1] - matlab_results[:, 1]
    freq_diff = julia_matrix[:, 2] - matlab_results[:, 2]
    drift_diff = julia_matrix[:, 3] - matlab_results[:, 3]
    residual_diff = julia_matrix[:, 4] - matlab_results[:, 4]
    
    println("Comparison (Julia - MATLAB):")
    println("  Phase difference RMS: $(sqrt(mean(phase_diff.^2)):.3e)")
    println("  Freq difference RMS: $(sqrt(mean(freq_diff.^2)):.3e)")
    println("  Drift difference RMS: $(sqrt(mean(drift_diff.^2)):.3e)")
    println("  Residual difference RMS: $(sqrt(mean(residual_diff.^2)):.3e)")
    
    # Plot differences if Plots is available
    try
        p1 = plot(phase_diff, title="Phase Estimate Difference (Julia - MATLAB)", xlabel="Sample", ylabel="Difference (ns)")
        p2 = plot(freq_diff, title="Frequency Estimate Difference", xlabel="Sample", ylabel="Difference (ns/s)")
        p3 = plot(drift_diff, title="Drift Estimate Difference", xlabel="Sample", ylabel="Difference (ns/s²)")
        p4 = plot(residual_diff, title="Residual Difference", xlabel="Sample", ylabel="Difference (ns)")
        
        plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))
        savefig("test/matlab_julia_comparison.png")
        println("Comparison plots saved to test/matlab_julia_comparison.png")
    catch e
        println("Could not create plots: $e")
    end
else
    println("\nMATLAB results not found. Run MATLAB comparison first.")
end