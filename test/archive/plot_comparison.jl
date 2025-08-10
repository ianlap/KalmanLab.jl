#!/usr/bin/env julia

# Create comparison plots between Julia and MATLAB results
using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics, Plots

println("=== Creating Julia vs MATLAB Comparison Plots ===")

# Load the Julia results from our recent run
julia_file = "test/julia_results_small.txt"
julia_data = readdlm(julia_file, Float64)
julia_phase = julia_data[:, 1]
julia_freq = julia_data[:, 2]
julia_drift = julia_data[:, 3]
julia_residuals = julia_data[:, 4]

# Load MATLAB reference data
matlab_file = "test/m_phase_100.txt" 
matlab_data = readdlm(matlab_file, Float64)
matlab_indices = matlab_data[:, 1]
matlab_phase = matlab_data[:, 2]

# Load MATLAB frequency data (if available)
matlab_freq_file = "test/m_freq_100.txt"
if isfile(matlab_freq_file)
    matlab_freq_data = readdlm(matlab_freq_file, Float64) 
    matlab_freq = matlab_freq_data[:, 2]
else
    # Generate approximate frequency from phase differences (if no freq file)
    matlab_freq = diff([matlab_phase[1]; matlab_phase]) # Simple difference approximation
end

println("Loaded $(length(julia_phase)) Julia samples and $(length(matlab_phase)) MATLAB samples")

# Create time axis
n_samples = min(length(julia_phase), length(matlab_phase))
time_axis = 1:n_samples

# Create plots
theme(:default)

# Plot 1: Phase estimates comparison
p1 = plot(time_axis, matlab_phase[1:n_samples], 
          label="MATLAB Phase", linewidth=2, color=:red, linestyle=:solid,
          title="Phase Estimates Comparison",
          xlabel="Sample", ylabel="Phase (ns)")
plot!(p1, time_axis, julia_phase[1:n_samples], 
      label="Julia Phase", linewidth=2, color=:blue, linestyle=:dash)

# Plot 2: Phase difference 
differences = julia_phase[1:n_samples] - matlab_phase[1:n_samples]
p2 = plot(time_axis, differences, 
          label="Julia - MATLAB", linewidth=2, color=:green,
          title="Phase Estimate Differences", 
          xlabel="Sample", ylabel="Difference (ns)")
hline!(p2, [0], color=:black, linestyle=:dot, label="Zero")

# Plot 3: Frequency estimates comparison
p3 = plot(time_axis, matlab_freq[1:n_samples],
          label="MATLAB Frequency", linewidth=2, color=:red, linestyle=:solid,
          title="Frequency Estimates Comparison",
          xlabel="Sample", ylabel="Frequency (ns/s)")
plot!(p3, time_axis, julia_freq[1:n_samples],
      label="Julia Frequency", linewidth=2, color=:blue, linestyle=:dash)

# Plot 4: Residuals
p4 = plot(time_axis, julia_residuals[1:n_samples],
          label="Julia Residuals", linewidth=2, color=:purple,
          title="Residuals",
          xlabel="Sample", ylabel="Residual (ns)")

# Combine plots
combined_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))

# Save the plot
savefig(combined_plot, "test/julia_matlab_comparison.png")
println("Saved comparison plots to test/julia_matlab_comparison.png")

# Print summary statistics
println("\n=== Summary Statistics ===")
println("Phase Estimates:")
println("  MATLAB RMS: $(round(sqrt(mean(matlab_phase[1:n_samples].^2)), digits=3)) ns")
println("  Julia RMS:  $(round(sqrt(mean(julia_phase[1:n_samples].^2)), digits=3)) ns")
println("  Difference RMS: $(round(sqrt(mean(differences.^2)), digits=3)) ns")
println("  Max difference: $(round(maximum(abs.(differences)), digits=3)) ns")

println("\nFrequency Estimates:")
println("  MATLAB RMS: $(round(sqrt(mean(matlab_freq[1:n_samples].^2)), digits=6)) ns/s") 
println("  Julia RMS:  $(round(sqrt(mean(julia_freq[1:n_samples].^2)), digits=6)) ns/s")

freq_diff = julia_freq[1:n_samples] - matlab_freq[1:n_samples]
println("  Frequency difference RMS: $(round(sqrt(mean(freq_diff.^2)), digits=6)) ns/s")

println("\nResiduals:")
println("  Julia residual RMS: $(round(sqrt(mean(julia_residuals[1:n_samples].^2)), digits=3)) ns")

# Display the plot
display(combined_plot)