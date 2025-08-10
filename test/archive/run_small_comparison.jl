#!/usr/bin/env julia

# Run a smaller comparison test to avoid divergence issues
using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics

println("=== Small Julia vs MATLAB Comparison ===")

# Load data
data_file = "/Users/ianlapinski/Desktop/masterclock-kflab/data/6krb25apr.txt"
data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]
println("Loaded $(length(phase_data)) samples")

# Create filter with same parameters
kf = KalmanFilter(0.1, 0.01, 0.05, 100.0, 0.01, 1e-6, 0.0, 0.0, 3, 1.0, zeros(3), 1e6)

# Initialize
data_trimmed = initialize!(kf, phase_data, 100)
println("Filter initialized: x0 = $(round.(kf.x0, digits=6))")

# Run on small subset (100 samples to match our successful tests)
test_size = 100  
test_data = data_trimmed[1:test_size]
println("Testing on $(length(test_data)) samples")

# Run Julia filter
println("\nRunning Julia filter...")
julia_results = KalmanLab.run(kf, test_data)

# Save Julia results
julia_matrix = hcat(julia_results.phase_est, julia_results.freq_est, 
                   julia_results.drift_est, julia_results.residuals, 
                   julia_results.steers)
writedlm("test/julia_results_small.txt", julia_matrix, '\t')

println("Julia results saved")
println("  Phase RMS: $(round(sqrt(mean(julia_results.phase_est.^2)), digits=3)) ns")
println("  Residual RMS: $(round(sqrt(mean(julia_results.residuals.^2)), digits=3)) ns")
println("  Final state: $(round.(julia_results.final_state, digits=6))")

# Check if we have MATLAB results to compare
matlab_file = "test/m_phase_100.txt"
if isfile(matlab_file)
    println("\nComparing with existing MATLAB reference...")
    matlab_data = readdlm(matlab_file, Float64)
    matlab_phase = matlab_data[:, 2]  # Phase estimates are in column 2
    
    # Compare first 10 samples
    println("\nComparison (first 10 samples):")
    println("Index  MATLAB         Julia          |Difference|")
    println("-" ^ 50)
    
    for i in 1:min(10, length(matlab_phase), length(julia_results.phase_est))
        diff = abs(julia_results.phase_est[i] - matlab_phase[i])
        println("$i      $(rpad(round(matlab_phase[i], digits=6), 12)) $(rpad(round(julia_results.phase_est[i], digits=6), 12)) $(round(diff, digits=8))")
    end
    
    # Overall statistics  
    n_compare = min(length(matlab_phase), length(julia_results.phase_est))
    differences = abs.(julia_results.phase_est[1:n_compare] - matlab_phase[1:n_compare])
    
    println("\nSummary Statistics:")
    println("  Mean difference: $(round(mean(differences), digits=8)) ns")
    println("  Max difference:  $(round(maximum(differences), digits=8)) ns")
    println("  RMS difference:  $(round(sqrt(mean(differences.^2)), digits=8)) ns")
else
    println("No MATLAB reference found at $matlab_file")
end