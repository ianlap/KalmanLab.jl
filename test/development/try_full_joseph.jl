#!/usr/bin/env julia

# Try the full Joseph form to see if that fixes the divergence
using Pkg
Pkg.activate(".")

using KalmanLab, DelimitedFiles, Statistics, LinearAlgebra

println("=== Testing Full Joseph Form ===")

# Temporarily modify the filter to use full Joseph form
function run_with_full_joseph(kf, data)
    # This is a copy of the run function but with full Joseph form
    
    N = length(data)
    nstates = kf.nstates
    τ = kf.tau
    g_p, g_i, g_d = kf.g_p, kf.g_i, kf.g_d
    
    # Process noise parameters
    q_wpm, q_wfm, q_rwfm = kf.q_wpm, kf.q_wfm, kf.q_rwfm
    q_irwfm, q_diurnal = kf.q_irwfm, kf.q_diurnal
    R = q_wpm  # Measurement noise
    
    # Initialize
    x = copy(kf.x0)
    P_init = isa(kf.P0, Matrix) ? copy(kf.P0) : Matrix{Float64}(kf.P0 * I, nstates, nstates)
    P = Symmetric(P_init)
    
    # System matrices
    Φ = Matrix{Float64}(I, nstates, nstates)
    Q = zeros(Float64, nstates, nstates)
    H = zeros(Float64, nstates)'
    H[1] = 1.0
    
    # Results
    phase_est = zeros(N)
    freq_est = zeros(N)
    residuals = zeros(N)
    sumsteers = zeros(N)
    sum2steers = zeros(N)
    
    # PID state
    pid_state = [0.0, 0.0]  # [sumx, last_steer]
    
    # Working copy of phase data
    phase = copy(data)
    
    for k in 1:N
        τₖ = τ
        absτ = abs(τₖ)
        
        # Update Φ
        Φ .= 0.0
        for i in 1:nstates
            Φ[i,i] = 1.0
        end
        if nstates >= 2
            Φ[1,2] = τₖ
        end
        if nstates >= 3
            Φ[1,3] = τₖ^2 / 2
            Φ[2,3] = τₖ
        end
        
        # Update Q
        Q .= 0.0
        τ² = τₖ^2
        τ³ = τₖ^3
        τ⁴ = τₖ^4
        τ⁵ = τₖ^5
        
        Q[1,1] = q_wfm*absτ + q_rwfm*τ³/3 + q_irwfm*τ⁵/20
        if nstates >= 2
            Q[1,2] = q_rwfm*τ²/2 + q_irwfm*τ⁴/8
            Q[2,1] = Q[1,2]
            Q[2,2] = q_rwfm*absτ + q_irwfm*τ³/3
        end
        if nstates >= 3
            Q[1,3] = q_irwfm*τ³/6
            Q[3,1] = Q[1,3]
            Q[2,3] = q_irwfm*τ²/2
            Q[3,2] = Q[2,3]
            Q[3,3] = q_irwfm*absτ
        end
        
        # Prediction step
        if k > 1
            x = Φ * x
            x[1] += pid_state[2] * τₖ
            if nstates >= 2
                x[2] += pid_state[2]
            end
            
            P = Φ * P * Φ' + Q
            P = Symmetric(P)
        end
        
        # Update measurement
        if k > 1
            phase[k] = sum2steers[k-1] + data[k] - data[k-1] + sumsteers[k-1]
        end
        
        # Innovation
        innovation = phase[k] - dot(H, x)
        
        # Update step
        S = (H * P * H')[1,1] + R
        K = (P * H') / S
        
        # State update
        x = x + K * innovation
        
        # FULL JOSEPH FORM COVARIANCE UPDATE
        I_KH = I - K * H
        P = I_KH * P * I_KH' + K * K' * R
        P = Symmetric(P)
        
        # Calculate residual
        residual = phase[k] - x[1]
        
        # PID steering
        pid_state[1] += x[1]
        
        if nstates >= 2
            steer = -g_p * x[1] - g_i * pid_state[1] - g_d * x[2]
        else
            steer = -g_p * x[1] - g_i * pid_state[1]
        end
        
        pid_state[2] = steer
        
        # Update cumulative steering
        if k == 1
            sumsteers[k] = pid_state[2]
            sum2steers[k] = sumsteers[k]
        else
            sumsteers[k] = sumsteers[k-1] + pid_state[2]
            sum2steers[k] = sum2steers[k-1] + sumsteers[k]
        end
        
        # Store results
        phase_est[k] = x[1]
        freq_est[k] = nstates >= 2 ? x[2] : 0.0
        residuals[k] = residual
    end
    
    return (phase_est=phase_est, freq_est=freq_est, residuals=residuals, 
            final_P=P, final_state=x)
end

# Test the full Joseph form
data_file = "/Users/ianlapinski/Desktop/masterclock-kflab/data/6krb25apr.txt"
data_matrix = readdlm(data_file, Float64)
phase_data = data_matrix[:, 2]

kf = KalmanFilter(0.1, 0.01, 0.05, 100.0, 0.01, 1e-6, 0.0, 0.0, 3, 1.0, zeros(3), 1e6)
data_trimmed = initialize!(kf, phase_data, 100)

println("Testing Full Joseph Form on 1000 samples:")
test_data = data_trimmed[1:1000]

results = run_with_full_joseph(kf, test_data)

# Check results
phase_rms = sqrt(mean(results.phase_est.^2))
freq_rms = sqrt(mean(results.freq_est.^2))
residual_rms = sqrt(mean(results.residuals.^2))

# Check final covariance
P_final = results.final_P
eigenvals = real(eigvals(P_final))
cond_num = cond(P_final)

println("Full Joseph Form Results:")
println("  Phase RMS: $(round(phase_rms, digits=3)) ns")
println("  Freq RMS: $(round(freq_rms, digits=6)) ns/s")
println("  Residual RMS: $(round(residual_rms, digits=3)) ns")
println("  Final phase: $(round(results.phase_est[end], digits=3)) ns")
println("  Final freq: $(round(results.freq_est[end], digits=6)) ns/s")
println("  P condition: $(round(cond_num, digits=2))")
println("  P eigenvalues: $(round.(eigenvals, digits=6))")
println("  Min eigenval: $(round(minimum(eigenvals), digits=8))")

if phase_rms < 1e6 && freq_rms < 1e3
    println("  Status: STABLE ✓")
else
    println("  Status: DIVERGED ✗")
end