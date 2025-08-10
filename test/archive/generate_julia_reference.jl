# Generate Julia Kalman filter reference for comparison with MATLAB

using Pkg
Pkg.activate(".")

using KalmanLab
using DelimitedFiles
using Statistics

println("=== Julia Kalman Filter Reference ===")

# Load data - phase error is in column 2
data_file = "../../data/6krb25apr.txt"
println("Loading data from: $data_file")

data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]  # Phase error data in column 2

println("Loaded $(length(phase_data)) samples")
println("Data range: [$(minimum(phase_data):.2e), $(maximum(phase_data):.2e)] ns")

# Use same parameters as MATLAB
kf = KalmanFilter(
    0.1,        # g_p
    0.01,       # g_i  
    0.05,       # g_d
    100.0,      # q_wpm (measurement noise R)
    0.01,       # q_wfm
    1e-6,       # q_rwfm
    0.0,        # q_irwfm
    0.0,        # q_diurnal
    3,          # nstates
    1.0,        # tau
    zeros(3),   # x0 (will be updated by initialize!)
    1e6         # P0 (will be updated by initialize!)
)

println("\nJulia Filter configuration:")
println("  PID gains: g_p=$(kf.g_p), g_i=$(kf.g_i), g_d=$(kf.g_d)")
println("  Noise params: q_wpm=$(kf.q_wpm), q_wfm=$(kf.q_wfm), q_rwfm=$(kf.q_rwfm)")
println("  States: $(kf.nstates), tau=$(kf.tau)")

# Initialize using first 100 points (matching MATLAB exactly)
data_trimmed = initialize!(kf, phase_data, 100)
println("\nAfter initialization (first 100 points):")
println("  Initial state x0: $(round.(kf.x0, digits=6))")
println("  P0 diagonal: $(round.(diag(kf.P0), sigdigits=4))")
println("  Remaining data points: $(length(data_trimmed))")

# Run Julia Kalman filter on ALL remaining data (matching MATLAB)
println("\nRunning Julia Kalman filter...")
@time results = KalmanLab.run(kf, data_trimmed)
julia_time = @elapsed KalmanLab.run(kf, data_trimmed)  # More accurate timing

println("Julia filter completed in $julia_time seconds")

# Calculate statistics (same as MATLAB)
phase_rms = sqrt(mean(results.phase_est.^2))
freq_rms = sqrt(mean(results.freq_est.^2))
drift_rms = sqrt(mean(results.drift_est.^2))
residual_rms = sqrt(mean(results.residuals.^2))

println("\nJulia Estimate statistics:")
println("  Phase RMS: $(phase_rms:.3e) ns")
println("  Frequency RMS: $(freq_rms:.3e) ns/s")
println("  Drift RMS: $(drift_rms:.3e) ns/s²")
println("  Residual RMS: $(residual_rms:.3e) ns")

# Save results in same format as MATLAB
results_matrix = hcat(results.phase_est, results.freq_est, results.drift_est, 
                     results.residuals, results.steers)
writedlm("kalman_results_julia.txt", results_matrix, '\t')
println("\nResults saved to kalman_results_julia.txt")

# Also save steering results (matching MATLAB)
steering_matrix = hcat(results.sumsteers, results.sum2steers)
writedlm("steering_results_julia.txt", steering_matrix, '\t')

println("Julia reference generated successfully!")