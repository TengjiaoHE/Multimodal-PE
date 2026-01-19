# Multimodal Parabolic Equation (Multimodal-PE) MATLAB Code

## Citation

This code accompanies the following publication:

**He, T., & Guo, W. (2026).** Wide-angle multimodal parabolic equations: modeling of directional sound propagation in stratified range-dependent fluid waveguides. *Proceedings of the Royal Society A*,.

If you use this code in your research, please cite the above paper.

## Overview

This repository contains MATLAB code implementing the Multimodal-PE method for modeling underwater acoustic propagation in stratified, range-dependent environments. The method is particularly suited for:

- **Directional sources** with arbitrary prescribed directivity (which can be complex values to take into account phase information)
- **Multilayer seabeds** with density and sound speed discontinuities
- **Very-wide-angle propagation** (up to ±90°)
- **Range-dependent environments** with adaptive optimization

### Key Features

1. **Strong convergence**: O(1/M²) where M is the number of basis modes
2. **Analytical treatment** of interface discontinuities (no numerical instability)
3. **Adaptive optimization**: Fast modal propagation for range-invariant segments
4. **GPU acceleration** and precision control (single/double)
5. **Very-wide-angle capability**: Supports propagation angles up to ±90°

## File Descriptions

### Main Files

- `multimodal_pe.p` - Main Multimodal-PE solver (Eq. 2.15 in paper， p-code version for efficiency and stability considerations, the full version of M file is available upon reasonable request)
- `profile_upd_multilayer.m` - Updates K operator for sound speed profile (Eq. 2.23-2.26)
- `topo_upd_multilayer.m` - Updates C and L operators for topography (Eq. 2.22)
- `configPML_multilayer.m` - Configures PML for domain truncation (Eq. 2.18-2.20)

### Variable Naming Convention (Following Paper Notation)

| Code Variable | Paper Symbol | Description |
|--------------|--------------|-------------|
| `m`, `n` | m, n | Mode indices |
| `H` | H | Total depth (physical + PML) |
| `D` | D | Physical domain depth |
| `omega` | ω | Angular frequency (2πf) |
| `k0` | k₀ | Reference wavenumber (ω/c₀) |
| `c` | c(z) | Sound speed profile |
| `rho` | ρ(z) | Density profile |
| `atten` | α(z) | Attenuation profile |
| `h` | hⱼ | Layer interface depths |
| `Vp` | a | Projection coefficients |
| `gn` | φ(z) | Basis mode matrix |
| `psi` | ψ(r,z) | Pressure envelope |
| `C` | C | Second derivative operator |
| `K` | K | Profile operator |
| `L` | L | First derivative operator |

## Theoretical Background

### Multimodal Formulation

The method projects the acoustic pressure field onto basis modes:

**Eq. (2.6):**
```
ψ(r,z) = Σₘ aₘ(r)φₘ(z) = φᵀa
```

where φₘ(z) = √(2/H)sin(mπz/H) are sine basis modes satisfying pressure-release boundary conditions at z=0 and z=H.

### One-Way Equation

**Eq. (2.15):**
```
∂a/∂r = ik₀(√(A⁻¹B + I) - I)a
```

where:
- **A** (Eq. 2.9): Mass matrix
- **B** (Eq. 2.10): Stiffness matrix (B = C + K)
- **C** (Eq. 2.22): Second derivative operator
- **K** (Eq. 2.23-2.26): Profile operator

### Split-Step Padé Propagator

**Eq. (2.17):**
```
exp[ik₀Δr(√(A⁻¹B+I)-I)] ≈ ∏ₗ₌₁ᴸ (1+αₗA⁻¹B)/(1+βₗA⁻¹B)
```

### Directional Starter

**Eq. (2.35):**
```
a(r₀) = 2πW Λ(r₀) c(zₛ)
```

where:
- **W**: Projection matrix from normal modes to basis modes
- **Λ(r₀)**: Diagonal matrix with modal propagation factors
- **c(zₛ)**: Modal excitation coefficients (Eq. 2.30)

## Usage Example

```matlab
% Environment setup


% Source beam pattern (omnidirectional)
sbp.dir = @(phi, k) ones(size(phi));

% Run simulation
use_gpu = false;         % Set true for GPU acceleration
precision = 'double';    % 'single' or 'double'
grid = multimodal_pe(env, grid, sbp, use_gpu, precision);

% Results
TL = -20*log10(abs(grid.p));  % Transmission loss [dB]
```

### Directional Source Example

```matlab
% Directional source (e.g., double pistons, Eq. 3.2)
ra = 10;  % Piston radius [m]
phi_s = 15*pi/180;  % Steering angle [rad]

sbp.dir = @(phi, k) ...
    2*besselj(1, k*ra*sin(phi-phi_s))./(k*ra*sin(phi-phi_s)) + ...
    2*besselj(1, k*ra*sin(phi+phi_s))./(k*ra*sin(phi+phi_s));
```

## Algorithm Flow

1. **Initialization**
   - Set up grid and environmental parameters
   - Compute PML operators
   - Initialize basis modes φₘ(z)

2. **Starting Field** (Section 2(d))
   - Compute modal excitation coefficients (Eq. 2.30-2.31)
   - Project onto basis modes (Eq. 2.35)

3. **Environment Analysis**
   - Identify range-invariant segments
   - Select propagation method for each segment

4. **Range Marching**
   - **Invariant segments**: Modal propagation (eigendecomposition)
   - **Varying segments**: Padé method (Eq. 2.17)
   - Energy conservation via impedance matching (Section 2(e))

5. **Post-processing**
   - Reconstruct pressure: p = ψH₀⁽¹⁾(k₀r) (Eq. 2.2)
   - Compute transmission loss

## Performance Considerations

### Recommended Parameters

- **Number of modes**: M = 2Hf/c₀ (default)
- **Range step**: Δr ∈ [λ/2, λ] for Padé method
- **Padé order**: 8 (good balance of accuracy and efficiency)
- **Minimum segment length**: 20 steps for modal propagation

### Computational Complexity

- **Depth operator construction**: O(M³)
- **Padé propagation**: O(nₚM²) per range step
- **Modal propagation**: O(M³) for eigendecomposition (once per segment)

### Memory Usage

- **Double precision**: ~8M² bytes for operators
- **Single precision**: ~4M² bytes for operators
- **GPU**: Additional memory for large grids

## Validation Benchmarks

The code has been validated against:

1. **ASA benchmark wedge** (Section 3(a)(i))
2. **RAM range-dependent benchmark** (Section 3(a)(ii))
3. **Multilayer seabed** (Section 3(d))
4. **Directional propagation** (Section 3(c))

## Dependencies

### Required
- MATLAB R2019b or later

### Optional
- MATLAB Parallel Computing Toolbox (for GPU acceleration)

## Known Limitations

1. **2D only**: Current implementation is for range-depth (r-z) problems
2. **Fluid media**: Does not handle elastic seabeds (see elastic PE extensions)
3. **One-way**: No backscattering (inherent to PE methods)

## Troubleshooting

### Common Issues

**Numerical instability**
- Reduce range step Δr or depth step Δz
- Check PML parameters

**Slow performance**
- Use single precision
- Enable GPU acceleration
- Check segment identification (avoid too many small segments)

## Contact

For questions, bug reports, or contributions, please contact:

**Tengjiao He**  
College of Meteorology and Oceanography  
Shanghai Jiao Tong University  
Shanghai, China  
Email: hetengjiao@sjtu.edu.cn

**Wei Guo**  
College of Meteorology and Oceanography  
National University of Defense Technology  
Changsha, China  
Email: guowei23@nudt.edu.cn

## License

This code is provided as supplementary material to the publication. It is distributed "as-is" without warranty of any kind. Users are free to use and modify the code for research purposes, provided that proper citation is given to the original publication.


---

**Last updated**: January 2026
