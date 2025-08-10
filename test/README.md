# KalmanLab.jl Test Suite

This directory contains comprehensive tests and validation for the KalmanLab.jl package.

## Directory Structure

### `final_results/`
Contains the final validation results and key test scripts:
- `perfect_matlab_julia_match.png` - Visual confirmation of perfect MATLAB-Julia agreement
- `perfect_match_statistics.png` - Statistical validation plots
- `julia_full_dataset_analysis.png` - Full dataset analysis plots
- `test_stable_pid.jl` - Main test script with stable PID parameters
- `full_dataset_test.jl` - Full 400k+ sample validation test
- `matlab_julia_full_comparison.m` - MATLAB comparison script
- `*_summary.txt` - Summary statistics and results

### `development/`
Contains development and debugging scripts used during implementation:
- `debug_divergence.jl` - Debugging filter divergence issues
- `test_covariance_forms.jl` - Testing different covariance update methods
- `test_exact_matlab_match.jl` - Achieving exact MATLAB numerical match
- Other debugging and development scripts

### `archive/`
Contains historical test files and intermediate results from the development process.

### `matlab_reference/`
Contains MATLAB reference implementations for comparison.

## Key Test Results

✅ **Perfect MATLAB Match**: Julia implementation produces mathematically identical results to MATLAB
- Phase difference: 0.0 ns
- Frequency difference: 0.0 ns/s
- Correlation coefficients: 1.0

✅ **Performance**: Julia runs at 95.5% of MATLAB speed (337k vs 353k samples/sec)

✅ **Stability**: Stable on full 406,763-sample dataset over 113 hours

✅ **Key Parameters**:
- PID gains from time constant T=1: g_p=0.694, g_i=0.253, g_d=0.950
- Noise parameters: q_wpm=100, q_wfm=6e-3, q_rwfm=5e-9
- Covariance update: Simple form `P = (I-KH)P` matching MATLAB exactly

## Running Tests

### Basic Functionality Test
```julia
julia> include("final_results/test_basic_success.jl")
```

### Full Dataset Validation
```julia  
julia> include("final_results/full_dataset_test.jl")
```

### MATLAB Comparison
```bash
# Run MATLAB comparison
matlab -nodisplay -r "run('final_results/matlab_julia_full_comparison.m'); exit"

# Then create verification plots
julia> include("final_results/plot_perfect_match.jl")
```

## Test Summary

The Julia implementation successfully:
1. Matches MATLAB numerically (perfect agreement)
2. Maintains stability over long datasets
3. Achieves high performance (95.5% of MATLAB speed)
4. Implements proper PID steering with time-constant-based gains
5. Uses identical covariance update method as MATLAB

This validates the Julia implementation as a production-ready replacement for the MATLAB Kalman filter.