#!/usr/bin/env julia

# Compare the Simple Form Julia results with MATLAB and previous Joseph Form results

using DelimitedFiles, Statistics

println("=== Simple Form vs MATLAB vs Joseph Form Comparison ===")

# Simple form results (current run from test above)
simple_results = [
    1.527017, 1.360308, 1.034885, 0.58151, 0.046794,
    -0.511544, -1.030851, -1.45049, -1.71883, -1.799471
]

# MATLAB results
matlab_results = [
    1.52701670140534, 1.36355584208931, 1.04331495666785, 0.595436930220503, 0.06488614777913,
    -0.492019983768153, -1.01355270364178, -1.43937350373387, -1.71745894731139, -1.81036715154112
]

# Previous Julia Joseph form results  
joseph_results = [
    1.5270167014053444, 1.3603082373677884, 1.0348849544933418, 0.5815099256494289, 0.046794285865291496,
    -0.5115441965041603, -1.0308509475179612, -1.4504900683816457, -1.71882997569865, -1.7994714256730442
]

println("\nComparison for first 10 samples:")
println("Index  MATLAB           Joseph Form      Simple Form      |Simple-MATLAB|  |Simple-Joseph|")
println("-" ^ 95)

for i in 1:10
    diff_matlab = abs(simple_results[i] - matlab_results[i])
    diff_joseph = abs(simple_results[i] - joseph_results[i])
    
    println("$i      $(rpad(round(matlab_results[i], digits=6), 12)) $(rpad(round(joseph_results[i], digits=6), 12)) $(rpad(round(simple_results[i], digits=6), 12)) $(rpad(round(diff_matlab, digits=8), 12)) $(round(diff_joseph, digits=8))")
end

# Overall statistics
diff_matlab_all = abs.(simple_results - matlab_results)
diff_joseph_all = abs.(simple_results - joseph_results)

println("\nSummary Statistics:")
println("Simple vs MATLAB:")
println("  Mean difference: $(round(mean(diff_matlab_all), digits=8)) ns")
println("  Max difference:  $(round(maximum(diff_matlab_all), digits=8)) ns")
println("  RMS difference:  $(round(sqrt(mean(diff_matlab_all.^2)), digits=8)) ns")

println("\nSimple vs Joseph:")
println("  Mean difference: $(round(mean(diff_joseph_all), digits=8)) ns") 
println("  Max difference:  $(round(maximum(diff_joseph_all), digits=8)) ns")
println("  RMS difference:  $(round(sqrt(mean(diff_joseph_all.^2)), digits=8)) ns")