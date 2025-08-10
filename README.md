# KalmanLab.jl

A Julia package for oscillator phase/frequency tracking using Kalman filters with PID steering control.

## Author

**Ian Lapinski** (ian.lapinski.99@gmail.com)

## Overview

This package provides a clean, performant implementation of Kalman filtering for time and frequency applications. It translates and improves upon a MATLAB implementation while leveraging Julia's type system and performance advantages.

### Key Features

- **Stateful Design**: Simple `KalmanFilter` type with intuitive API
- **Flexible Models**: Support for 2, 3, or 5 state clock models
- **PID Steering**: Integrated proportional-integral-derivative control
- **Holdover Analysis**: Easy transition from active steering to prediction
- **Parameter Optimization**: Built-in optimization of noise parameters
- **Performance**: Optimized for large datasets with minimal allocations

## Quick Start

```julia
using KalmanLab

# Create a 3-state filter with PID control
kf = KalmanFilter(
    g_p=0.1, g_i=0.01, g_d=0.05,           # PID gains
    q_wpm=100.0, q_wfm=0.01, q_rwfm=1e-6,  # Noise parameters
    q_irwfm=0.0, q_diurnal=0.0,            # Optional noise terms
    nstates=3, tau=1.0,                     # Model: 3 states, 1 sec sampling
    x0=zeros(3), P0=1e6                     # Initial conditions
)

# Run the filter on phase error data
results = run!(kf, phase_error_data)

# Access results
plot(results.phase_est)      # Estimated phase
plot(results.freq_est)       # Estimated frequency
plot(results.residuals)      # Measurement residuals

# Holdover prediction
kf.g_p = kf.g_i = kf.g_d = 0.0  # Turn off steering
predictions = holdoverpredict(kf, [10, 100, 1000])  # Predict 10, 100, 1000 steps
```

## Installation

This package is currently under development. To use:

```julia
] dev path/to/KalmanLab
```

## Documentation

- Mathematical background and theory
- API reference and examples
- Performance benchmarks vs MATLAB implementation
- Integration with stability analysis tools

## License

MIT License - see LICENSE file for details.