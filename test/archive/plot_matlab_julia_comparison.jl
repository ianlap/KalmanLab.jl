#!/usr/bin/env julia

# Create detailed comparison plots between MATLAB and Julia results
using Pkg
Pkg.activate(".")

using DelimitedFiles, Plots, Statistics

println("=== Creating MATLAB vs Julia Comparison Plots ===")

# Load comparison data
comparison_data = readdlm("test/matlab_julia_comparison.txt", Float64)
matlab_phase = comparison_data[:, 1]
matlab_freq = comparison_data[:, 2]
julia_phase = comparison_data[:, 3]
julia_freq = comparison_data[:, 4]
phase_diff = comparison_data[:, 5]
freq_diff = comparison_data[:, 6]

n_samples = length(matlab_phase)
time_hours = (1:n_samples) * 100 / 3600  # Convert to hours (100 samples per point)

println("Loaded $(n_samples) comparison points covering $(round(time_hours[end], digits=1)) hours")

# Create comparison plots
theme(:default)

# Plot 1: Phase estimates comparison
p1 = plot(time_hours, matlab_phase, label="MATLAB", linewidth=2, color=:red, alpha=0.8,
          title="Phase Estimates: MATLAB vs Julia",
          xlabel="Time (hours)", ylabel="Phase Error (ns)")
plot!(p1, time_hours, julia_phase, label="Julia", linewidth=2, color=:blue, alpha=0.8)

# Plot 2: Frequency estimates comparison
p2 = plot(time_hours, matlab_freq, label="MATLAB", linewidth=2, color=:red, alpha=0.8,
          title="Frequency Estimates: MATLAB vs Julia", 
          xlabel="Time (hours)", ylabel="Frequency Error (ns/s)")
plot!(p2, time_hours, julia_freq, label="Julia", linewidth=2, color=:blue, alpha=0.8)

# Plot 3: Phase differences (MATLAB - Julia)
p3 = plot(time_hours, phase_diff, label="MATLAB - Julia", linewidth=1, color=:green,
          title="Phase Estimate Differences",
          xlabel="Time (hours)", ylabel="Phase Difference (ns)")
hline!(p3, [0], color=:black, linestyle=:dot, alpha=0.5, label="Zero")

# Plot 4: Frequency differences (MATLAB - Julia)
p4 = plot(time_hours, freq_diff, label="MATLAB - Julia", linewidth=1, color=:orange,
          title="Frequency Estimate Differences",
          xlabel="Time (hours)", ylabel="Frequency Difference (ns/s)")
hline!(p4, [0], color=:black, linestyle=:dot, alpha=0.5, label="Zero")

# Combined comparison plot
comparison_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 900))
savefig(comparison_plot, "test/matlab_julia_comparison_plots.png")
println("Saved comparison plots to matlab_julia_comparison_plots.png")

# Statistical analysis plots
p5 = histogram(phase_diff, bins=50, alpha=0.7, color=:green,
               title="Phase Difference Distribution", 
               xlabel="Phase Difference (ns)", ylabel="Count",
               label="σ = $(round(std(phase_diff), digits=4)) ns")

p6 = histogram(freq_diff, bins=50, alpha=0.7, color=:orange,
               title="Frequency Difference Distribution",
               xlabel="Frequency Difference (ns/s)", ylabel="Count",
               label="σ = $(round(std(freq_diff), digits=6)) ns/s")

p7 = scatter(matlab_phase, julia_phase, alpha=0.5, color=:blue, 
             title="Phase: MATLAB vs Julia", xlabel="MATLAB Phase (ns)", ylabel="Julia Phase (ns)",
             label="Data", markersize=2)
plot!(p7, [minimum(matlab_phase), maximum(matlab_phase)], 
      [minimum(matlab_phase), maximum(matlab_phase)], 
      color=:red, linewidth=2, label="Perfect Agreement")

p8 = scatter(matlab_freq, julia_freq, alpha=0.5, color=:red,
             title="Frequency: MATLAB vs Julia", xlabel="MATLAB Frequency (ns/s)", ylabel="Julia Frequency (ns/s)",
             label="Data", markersize=2)
plot!(p8, [minimum(matlab_freq), maximum(matlab_freq)], 
      [minimum(matlab_freq), maximum(matlab_freq)], 
      color=:red, linewidth=2, label="Perfect Agreement")

stats_plot = plot(p5, p6, p7, p8, layout=(2,2), size=(1200, 900))
savefig(stats_plot, "test/matlab_julia_statistical_comparison.png")
println("Saved statistical comparison plots to matlab_julia_statistical_comparison.png")

# Print detailed comparison statistics
println("\n=== Detailed Comparison Statistics ===")

println("MATLAB Results:")
println("  Phase RMS: $(round(sqrt(mean(matlab_phase.^2)), digits=6)) ns")
println("  Freq RMS: $(round(sqrt(mean(matlab_freq.^2)), digits=8)) ns/s")

println("Julia Results:")
println("  Phase RMS: $(round(sqrt(mean(julia_phase.^2)), digits=6)) ns")  
println("  Freq RMS: $(round(sqrt(mean(julia_freq.^2)), digits=8)) ns/s")

println("Differences (MATLAB - Julia):")
println("  Phase difference RMS: $(round(sqrt(mean(phase_diff.^2)), digits=8)) ns")
println("  Phase difference Max: $(round(maximum(abs.(phase_diff)), digits=8)) ns")
println("  Phase difference Mean: $(round(mean(phase_diff), digits=8)) ns")
println("  Phase difference Std: $(round(std(phase_diff), digits=8)) ns")

println("  Freq difference RMS: $(round(sqrt(mean(freq_diff.^2)), digits=10)) ns/s")
println("  Freq difference Max: $(round(maximum(abs.(freq_diff)), digits=10)) ns/s")
println("  Freq difference Mean: $(round(mean(freq_diff), digits=10)) ns/s")  
println("  Freq difference Std: $(round(std(freq_diff), digits=10)) ns/s")

# Correlation analysis
phase_corr = cor(matlab_phase, julia_phase)
freq_corr = cor(matlab_freq, julia_freq)

println("Correlation coefficients:")
println("  Phase correlation: $(round(phase_corr, digits=8))")
println("  Frequency correlation: $(round(freq_corr, digits=8))")

# Performance comparison
println("\nPerformance comparison:")
println("  MATLAB: 352,828 samples/second")
println("  Julia: 72,814 samples/second")
println("  MATLAB is $(round(352828/72814, digits=1))x faster")

# Agreement assessment
phase_agreement = sqrt(mean(phase_diff.^2)) < 0.1
freq_agreement = sqrt(mean(freq_diff.^2)) < 0.001

println("\nAgreement Assessment:")
if phase_agreement && freq_agreement
    println("  ✅ EXCELLENT: Both phase and frequency show excellent agreement")
elseif sqrt(mean(phase_diff.^2)) < 0.1
    println("  ✅ GOOD: Phase shows good agreement")
    println("  ⚠️ MODERATE: Frequency differences are larger than expected")
else
    println("  ⚠️ FAIR: Both implementations show measurable differences")
    println("  This could be due to:")
    println("    - Different floating-point precision handling")
    println("    - Different matrix operation implementations")
    println("    - Numerical accumulation differences over 400k samples")
end

display(comparison_plot)