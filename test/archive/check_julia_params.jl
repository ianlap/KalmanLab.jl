#!/usr/bin/env julia

# Check Julia initialization parameters to match MATLAB
using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics, LinearAlgebra

println("=== Julia Parameter Check ===")

# Load data exactly like MATLAB
data_file = "/Users/ianlapinski/Desktop/masterclock-kflab/data/6krb25apr.txt"
data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]

# Manual initialization to match MATLAB exactly
n_fit = 100
t = (0:n_fit-1) * 1.0
y = phase_data[1:n_fit]
A = [ones(n_fit) t t.^2/2]
coeffs = A \ y
x0 = coeffs
y_pred = A * coeffs
residuals_fit = y - y_pred
P0 = var(residuals_fit) * inv(A' * A)

println("Julia Initialization Results:")
println("x0 = $(round.(x0, digits=6))")
println("P0 diagonal = $(round.([P0[1,1], P0[2,2], P0[3,3]], digits=6))")
println("P0 condition number = $(round(cond(P0), digits=2))")
eigenvals = real(eigvals(P0))
println("P0 eigenvalues = $(round.(eigenvals, digits=6))")

# Create filter with these exact parameters
kf = KalmanFilter(0.1, 0.01, 0.05, 100.0, 0.01, 1e-6, 0.0, 0.0, 3, 1.0, x0, P0)

# Get trimmed data
data_trimmed = phase_data[n_fit+1:end]
test_data = data_trimmed[1:200]  # Same 200 samples as MATLAB

println("\nRunning Julia on 200 samples...")
results = KalmanLab.run(kf, test_data)

println("Julia Results (200 samples):")
println("  Phase RMS: $(round(sqrt(mean(results.phase_est.^2)), digits=3)) ns")
println("  Freq RMS: $(round(sqrt(mean(results.freq_est.^2)), digits=6)) ns/s")
println("  Final phase: $(round(results.phase_est[end], digits=6)) ns")
println("  Final freq: $(round(results.freq_est[end], digits=6)) ns/s")

# Check covariance condition
P_final = results.final_P
final_eigenvals = real(eigvals(P_final))
println("  Final P condition: $(round(cond(P_final), digits=2))")
println("  Final P eigenvals: $(round.(final_eigenvals, digits=6))")
println("  Min eigenval: $(round(minimum(final_eigenvals), digits=8))")

# Check for issues
has_nan = any(isnan.(results.phase_est)) || any(isnan.(results.freq_est))
has_inf = any(isinf.(results.phase_est)) || any(isinf.(results.freq_est))

if has_nan || has_inf
    println("  WARNING: Julia has NaN or Inf values!")
else
    println("  Julia: No NaN/Inf detected")
end

# Save first 10 for comparison
results_10 = hcat(results.phase_est[1:10], results.freq_est[1:10])
writedlm("julia_debug_10.txt", results_10, '\t')

# Compare with MATLAB
if isfile("matlab_debug_10.txt")
    matlab_10 = readdlm("matlab_debug_10.txt", Float64)
    println("\nFirst 5 sample comparison:")
    println("Idx  MATLAB Phase   Julia Phase    MATLAB Freq    Julia Freq")
    for i in 1:5
        println("$i    $(rpad(round(matlab_10[i,1], digits=6), 12)) $(rpad(round(results.phase_est[i], digits=6), 12)) $(rpad(round(matlab_10[i,2], digits=6), 12)) $(round(results.freq_est[i], digits=6))")
    end
end