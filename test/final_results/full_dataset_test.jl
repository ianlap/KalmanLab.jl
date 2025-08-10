#!/usr/bin/env julia

# Test the stable filter on the full ~400k dataset
using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics, LinearAlgebra

println("=== Full Dataset Test (~400k samples) ===")

# Load full dataset
data_file = "/Users/ianlapinski/Desktop/masterclock-kflab/data/6krb25apr.txt"
data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]

println("Loaded $(length(phase_data)) samples")

# Use the stable parameters from previous test
T = 1.0
tau0 = 1.0
a = exp(-tau0 / T)
g_p = 1 - 3*a^2 + 2*a^3
g_i = 1 - 3*a + 3*a^2 - a^3
g_d = 1 - a^3

# Known-good noise parameters
q_wpm = 100.0
q_wfm = 6e-3
q_rwfm = 5e-9

kf_full = KalmanFilter(
    g_p, g_i, g_d,                  # Stable PID gains
    q_wpm, q_wfm, q_rwfm,          # Known-good noise parameters
    0.0, 0.0,                      # No additional noise terms
    3, tau0,                       # 3-state, 1-second sampling
    [0.0, 0.0, 0.0],               # Zero initial state
    1e6                            # Conservative initial covariance
)

println("Running Kalman filter on full dataset...")
println("This may take a moment...")

# Time the operation
start_time = time()
results_full = KalmanLab.run(kf_full, phase_data)
end_time = time()

processing_time = end_time - start_time
samples_per_second = length(phase_data) / processing_time

println("\n=== Full Dataset Results ===")
println("Processing time: $(round(processing_time, digits=2)) seconds")
println("Rate: $(round(samples_per_second, digits=0)) samples/second")

# Calculate final metrics
phase_rms = sqrt(mean(results_full.phase_est.^2))
freq_rms = sqrt(mean(results_full.freq_est.^2))
residual_rms = sqrt(mean(results_full.residuals.^2))

# Check final covariance
P_final = results_full.final_P
eigenvals = real(eigvals(P_final))
cond_num = cond(P_final)

println("\nFinal Results:")
println("  Phase RMS: $(round(phase_rms, digits=6)) ns")
println("  Freq RMS: $(round(freq_rms, digits=8)) ns/s")
println("  Residual RMS: $(round(residual_rms, digits=3)) ns")
println("  Final phase: $(round(results_full.phase_est[end], digits=6)) ns")
println("  Final freq: $(round(results_full.freq_est[end], digits=8)) ns/s")
println("  P condition: $(round(cond_num, digits=2))")
println("  P eigenvalues: $(round.(eigenvals, digits=8))")

# Check for any issues
has_nan = any(isnan.(results_full.phase_est)) || any(isnan.(results_full.freq_est))
has_inf = any(isinf.(results_full.phase_est)) || any(isinf.(results_full.freq_est))
min_eigenval = minimum(eigenvals)

if has_nan
    println("  ⚠️  WARNING: NaN values detected!")
elseif has_inf
    println("  ⚠️  WARNING: Inf values detected!")
elseif phase_rms > 10 || freq_rms > 1
    println("  ⚠️  WARNING: Large RMS values - possible instability!")
else
    println("  ✅ SUCCESS: Filter is stable on full dataset!")
end

# Save results for comparison and plotting
println("\nSaving full dataset results...")

# Save complete results (may be large file)
julia_matrix = hcat(results_full.phase_est, results_full.freq_est, 
                   results_full.drift_est, results_full.residuals, 
                   results_full.steers)

# Save every 100th sample to make files manageable
sample_indices = 1:100:length(results_full.phase_est)
sampled_results = julia_matrix[sample_indices, :]
writedlm("test/julia_full_sampled.txt", sampled_results, '\t')

println("Saved sampled results (every 100th sample) to julia_full_sampled.txt")

# Save summary statistics
summary_file = "test/julia_full_summary.txt"
open(summary_file, "w") do f
    println(f, "Julia Full Dataset Results")
    println(f, "Samples: $(length(phase_data))")
    println(f, "Processing time: $(round(processing_time, digits=2)) s")
    println(f, "Rate: $(round(samples_per_second, digits=0)) samples/s")
    println(f, "Phase RMS: $(round(phase_rms, digits=6)) ns")
    println(f, "Freq RMS: $(round(freq_rms, digits=8)) ns/s")
    println(f, "Residual RMS: $(round(residual_rms, digits=3)) ns")
    println(f, "Final phase: $(round(results_full.phase_est[end], digits=6)) ns")
    println(f, "Final freq: $(round(results_full.freq_est[end], digits=8)) ns/s")
    println(f, "P condition: $(round(cond_num, digits=2))")
    println(f, "Min eigenval: $(round(min_eigenval, digits=8))")
end

println("Saved summary to julia_full_summary.txt")
println("\n🎉 FULL DATASET TEST COMPLETE! 🎉")