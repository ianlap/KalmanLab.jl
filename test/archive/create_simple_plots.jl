#!/usr/bin/env julia

# Create simple, focused comparison plots
using Pkg
Pkg.activate(".")

using DelimitedFiles, Plots

println("=== Creating Simple Comparison Plots ===")

# Load data
julia_data = readdlm("test/julia_results_small.txt", Float64)
matlab_data = readdlm("test/m_phase_100.txt", Float64)

julia_phase = julia_data[:, 1]
julia_freq = julia_data[:, 2]
matlab_phase = matlab_data[:, 2]

# Approximate MATLAB frequency from phase differences
matlab_freq = diff([matlab_phase[1]; matlab_phase])

n = length(julia_phase)
time = 1:n

# Set up plotting theme
theme(:default)

# Phase comparison plot
p1 = plot(time, matlab_phase, label="MATLAB", linewidth=2, color=:red, 
          title="Phase Estimates: Julia vs MATLAB",
          xlabel="Sample Number", ylabel="Phase Error (ns)")
plot!(p1, time, julia_phase, label="Julia", linewidth=2, color=:blue, linestyle=:dash)

# Frequency comparison plot  
p2 = plot(time, matlab_freq, label="MATLAB", linewidth=2, color=:red,
          title="Frequency Estimates: Julia vs MATLAB", 
          xlabel="Sample Number", ylabel="Frequency Error (ns/s)")
plot!(p2, time, julia_freq, label="Julia", linewidth=2, color=:blue, linestyle=:dash)

# Difference plots
phase_diff = julia_phase - matlab_phase
freq_diff = julia_freq - matlab_freq

p3 = plot(time, phase_diff, label="Julia - MATLAB", linewidth=2, color=:green,
          title="Phase Estimate Differences",
          xlabel="Sample Number", ylabel="Difference (ns)")
hline!(p3, [0], color=:black, linestyle=:dot, alpha=0.5, label="")

p4 = plot(time, freq_diff, label="Julia - MATLAB", linewidth=2, color=:orange,
          title="Frequency Estimate Differences", 
          xlabel="Sample Number", ylabel="Difference (ns/s)")
hline!(p4, [0], color=:black, linestyle=:dot, alpha=0.5, label="")

# Combined plot
combined = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))

# Save
savefig(combined, "test/simple_comparison_plots.png")
println("Saved to: test/simple_comparison_plots.png")

# Statistics
println("\nComparison Statistics:")
println("Phase differences:")
println("  RMS: $(round(sqrt(mean(phase_diff.^2)), digits=4)) ns")
println("  Max: $(round(maximum(abs.(phase_diff)), digits=4)) ns")

println("Frequency differences:")  
println("  RMS: $(round(sqrt(mean(freq_diff.^2)), digits=6)) ns/s")
println("  Max: $(round(maximum(abs.(freq_diff)), digits=6)) ns/s")

display(combined)