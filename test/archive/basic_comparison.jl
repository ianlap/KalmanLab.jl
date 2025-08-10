# Basic numerical comparison without complex plots

using DelimitedFiles, Statistics

# Load results
julia_data = readdlm("test/julia_comparison_2000.txt", '\t', Float64)
matlab_data = readdlm("test/matlab_comparison_2000.txt", '\t', Float64)

println("=== MATLAB vs Julia Comparison Results ===")
println("Data shapes: Julia $(size(julia_data)), MATLAB $(size(matlab_data))")

# Column labels
labels = ["Phase", "Frequency", "Drift", "Residuals", "Steers"]

for col in 1:5
    j_data = julia_data[:, col]
    m_data = matlab_data[:, col]
    
    diff_data = j_data - m_data
    rel_diff = diff_data ./ (abs.(m_data) .+ 1e-15)  # Avoid division by zero
    
    println("\n$(labels[col]):")
    println("  Julia RMS:     $(sqrt(mean(j_data.^2)):.6e)")
    println("  MATLAB RMS:    $(sqrt(mean(m_data.^2)):.6e)")
    println("  Abs diff RMS:  $(sqrt(mean(diff_data.^2)):.6e)")
    println("  Max abs diff:  $(maximum(abs.(diff_data))):.6e")
    println("  Rel diff RMS:  $(sqrt(mean(rel_diff.^2))*100:.3f)%")
end

# Overall assessment
phase_diff = julia_data[:, 1] - matlab_data[:, 1]
max_phase_diff = maximum(abs.(phase_diff))

println("\n=== ASSESSMENT ===")
if max_phase_diff < 1e-12
    println("EXCELLENT: Results match to machine precision!")
elseif max_phase_diff < 1e-9  
    println("VERY GOOD: Results match to nanometer precision")
elseif max_phase_diff < 1e-6
    println("GOOD: Results match to micrometer precision")
elseif max_phase_diff < 1e-3
    println("ACCEPTABLE: Results match to millimeter precision") 
else
    println("WARNING: Significant differences detected")
end

println("Maximum phase difference: $(max_phase_diff:.2e) ns")