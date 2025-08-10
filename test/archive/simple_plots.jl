# Simple comparison plotting script

using Pkg; Pkg.activate(".")
using DelimitedFiles, Statistics, Plots

println("Loading comparison data...")

# Load Julia and MATLAB results  
julia_data = readdlm("test/julia_comparison_2000.txt", '\t', Float64)
matlab_data = readdlm("test/matlab_comparison_2000.txt", '\t', Float64)

println("Julia data shape: $(size(julia_data))")
println("MATLAB data shape: $(size(matlab_data))")

# Extract columns
j_phase = julia_data[:, 1]
j_freq = julia_data[:, 2] 
m_phase = matlab_data[:, 1]
m_freq = matlab_data[:, 2]

# Calculate differences
phase_diff = j_phase - m_phase
freq_diff = j_freq - m_freq

println("\nComparison Statistics:")
println("Julia Phase RMS: $(sqrt(mean(j_phase.^2)):.3e) ns")
println("MATLAB Phase RMS: $(sqrt(mean(m_phase.^2)):.3e) ns") 
println("Phase Difference RMS: $(sqrt(mean(phase_diff.^2)):.3e)")
println("Max Phase Difference: $(maximum(abs.(phase_diff)))") 

# Create plots
sample_range = 1:min(500, length(j_phase))

# Phase comparison
p1 = plot(sample_range, j_phase[sample_range], label="Julia", color=:blue, lw=2)
plot!(p1, sample_range, m_phase[sample_range], label="MATLAB", color=:red, lw=2, alpha=0.8)
title!(p1, "Phase Estimates Comparison")
ylabel!(p1, "Phase (ns)")

# Phase difference  
p2 = plot(sample_range, phase_diff[sample_range], color=:green, lw=2)
title!(p2, "Phase Difference (Julia - MATLAB)")
ylabel!(p2, "Difference (ns)")
xlabel!(p2, "Sample")

# Combined plot
final_plot = plot(p1, p2, layout=(2,1), size=(800,600))
savefig(final_plot, "test/comparison_plot.png")

println("\nComparison plot saved to test/comparison_plot.png")

# Check if they match
if maximum(abs.(phase_diff)) < 1e-6
    println("GOOD: Implementations match to reasonable precision")
else
    println("WARNING: Large differences detected - implementations may differ") 
    println("This could be due to:")
    println("  - Different matrix operations")
    println("  - Different numerical precision") 
    println("  - Different algorithm details")
end