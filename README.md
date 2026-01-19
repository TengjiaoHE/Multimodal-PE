## Cite

If you use this code, please cite:

- **Paper**: *Wide-angle multimodal parabolic equations: modeling of directional sound propagation in stratified range-dependent fluid waveguides*, Proceedings A, 2026. https://doi.org/XXXXXXXX

- **Code (paper-version release)**: Tengjiao He, *Multimodal-PE (v1.0.0)*, Zenodo. [https://doi.org/10.5281/zenodo.18298325](https://doi.org/10.5281/zenodo.18298325)

# Multimodal Parabolic Equation (Multimodal-PE) MATLAB Code

## Citation

This code accompanies the following publication:

**He, T., & Guo, W. (2026).** Wide-angle multimodal parabolic equations: modeling of directional sound propagation in stratified range-dependent fluid waveguides. *Proceedings of the Royal Society A*,.

If you use this code in your research, please cite the above paper.

---

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


---

## File Descriptions

### Main Files

- `multimodal_pe.p` - Main Multimodal-PE solver (Eq. 2.15 in paper, p-code version for efficiency and stability considerations, the full version of M file is available upon reasonable request)
- `profile_upd_multilayer.m` - Updates K operator for sound speed profile (Eq. 2.23-2.26)
- `topo_upd_multilayer.m` - Updates C and L operators for topography (Eq. 2.22)
- `configPML_multilayer.m` - Configures PML for domain truncation (Eq. 2.18-2.20)
- `read_env.m` - Loads environmental and computational parameters
- `fejer.m` - Weights of the Fejer2, Clenshaw-Curtis and Fejer1 quadratures by DFTs
- `ClenshawCurtis.m` - Clenshaw-Curtis scheme for numerical integration

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


---

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


---

## Usage Example

```matlab
% Environment setup
[grid, env, casename] = read_env('input.env'); % Load env file

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
### Explanation of env file
- the `*.env` file layout and conventions;
- the meaning (and units) of each field; and
- an example `input4.env` (downslope case).

#### 1. Comment syntax

- Comments use a **double backslash** delimiter: `\\`
- Anything after the first `\\` on a line is ignored by the parser.
- You can write **full-line comments** or **inline comments**.

Example:

```text
1600 \\ reference wave velocity [m/s]
```


#### 2. File layout (high level)

A valid `*.env` file is organized as:

1. **Case name** (string, one line)
2. **Header block** (global dimensions + physical/numerical parameters)
3. **Range nodes** `env.rg` (defines range dependence)
4. **Layer profile blocks** (one block per layer)
5. **Interface geometry / topography** (one block per interface)

> ⚠️ The parser reads values **strictly in order**. Do not insert extra tokens unless they are placed after `\\`.


#### 3. Case name

Line 1: free text.

```text
My test case name
```


#### 4. Header block (global settings)

##### 4.1 Dimensions

Immediately after the case name, provide:

1. `numLayer` (integer): number of layers  
2. `numRange` (integer): number of **range nodes** used to specify range-dependent profiles  
3. `nprofile` (array of `numLayer` integers): number of depth samples in each layer

Example (`numLayer=2`, `numRange=2`, each layer has 2 depth samples):

```text
2
2
2 2
```

##### 4.2 Physical & numerical parameters (in this exact order)

| In file | Variable | Meaning | Unit |
|---|---|---|---|
| `c0` | `env.c0` | reference wave velocity | m/s |
| `freq` | `env.freq` | frequency | Hz |
| `zs` | `grid.zs` | source depth | m |
| `zr` | `grid.zr` | receiver depth | m |
| `rmax` | `grid.rmax` | maximum range | m |
| `dr` | `grid.dr` | range step | m |
| `zmax` | `grid.zmax` | maximum computational depth | m |
| `dz` | `grid.dz` | depth step | m |
| `dpml` | `grid.dpml` | PML thickness | m |
| `npml` | `grid.npml` | number of PML grid points | count |
| `nP` | `grid.nP` | Padé order (wide-angle) | order |


#### 5. Range nodes: `env.rg`

After the header scalars, provide `numRange` numbers:

- `env.rg(1..numRange)` (meters)

These are the range locations at which layer profiles are specified.
The solver typically interpolates between them.

Example:

```text
0 4000
```


#### 6. Layer profile blocks

For each layer `i = 1..numLayer`, provide `nprofile(i)` **depth sample lines**.

##### 6.1 Per-line format

Each depth-sample line contains:

- `1 + numRange + 2` numeric fields, in this order:

1. `z` : depth coordinate for this sample  
2. `c(r1) ... c(rN)` : sound speed at each range node `env.rg`  
3. `rho` : density  
4. `alpha` : attenuation  

So the column layout is:

```text
z   c(r1)  c(r2) ... c(rN)   rho   alpha
```

##### 6.2 Depth coordinate convention (important)

The format itself does enforce `z` is:
- absolute depth from the sea surface.

##### 6.3 Units for `rho` and `alpha`

- `rho`: g/cm^3.
- `alpha`: dB/λ.  

#### 7. Interface geometry / topography

After all layer profile blocks, provide:

1. `toponr` : integer, number of control points per interface

Then for each interface `k = 1..(numLayer-1)` provide **2×toponr** numbers:

- first **all ranges**: `r1 r2 ... rM`
- then **all depths**: `z1 z2 ... zM`

Typically written as:

```text
r1 r2 ... rM
z1 z2 ... zM
```

Interfaces are generally treated as piecewise-linear between control points.


#### 8. Example: `input4.env` (downslope case)

This is an example with:

- 2 layers (water + bottom)
- 2 range nodes (`r=0` and `r=4000 m`)
- a sloping bottom from 100 m to 400 m water depth

```text
Example 4: Downslope double positon beam
2 \\ numLayer
2 \\ numRange
2 2 \\ nprofile for each layer

1600 \\ c0 [m/s]
1000 \\ frequency [Hz]
50 \\ source depth zs [m]
30 \\ receiver depth zr [m]
4000 \\ rmax [m]
2 \\ dr [m]
410 \\ zmax [m]
0.2 \\ dz [m]
10 \\ dpml [m]
32 \\ npml [grid points]
8 \\ Pade order nP

0 4000 \\ env.rg range nodes [m]

\\ Layer 1 (water): z, c(r=0), c(r=4000), rho, alpha
0   1500 1500  1.0 0.0
400 1500 1500  1.0 0.0

\\ Layer 2 (bottom): z, c(r=0), c(r=4000), rho, alpha
0   2200 2200  2.0 0.2
420 2200 2200  2.0 0.2

2 \\ toponr (interface control points)

\\ Interface 1: ranges then depths
0 4000
100 400
```

---

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


---

## Performance Considerations

### Recommended Parameters

- **Number of modes**: M = 3Hf/c₀ (default)
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


---

## Validation Benchmarks

The code has been validated against:

1. **ASA benchmark wedge** (running demo1.m, Fig.3 of the paper)
<img src='/images/demo1.png'>

2. **RAM range-dependent benchmark** (running demo2.m，Fig.6 of the paper)
<img src='/images/demo2.png'>

3. **Multilayer seabed** (running demo3.m, Fig.13 of the paper)
<img src='/images/demo3.png'>

4. **Directional propagation(Downslope)** (running demo4.m, Fig.9 of the paper)
<img src='/images/demo4.png'>

5. **Directional propagation(Dickins seamount)** (running demo5.m, Fig.12 of the paper)
<img src='/images/demo5.png'>

6. **Warm Core Eddy** (running demo6.m, highlighting high efficiency for range-independent segmentation scheme)
<img src='/images/demo6.png'>
   
---

## Dependencies

### Required
- MATLAB R2019b or later

### Optional
- MATLAB Parallel Computing Toolbox (for GPU acceleration)


---

## Known Limitations

1. **2D only**: Current implementation is for range-depth (r-z) problems
2. **Fluid media**: Does not handle elastic seabeds (see elastic PE extensions)
3. **One-way**: No backscattering (inherent to PE methods)


---

## Contact

For questions, bug reports, or contributions, please contact:

**Tengjiao He**  
School of Ocean and Civil Engineering  
Shanghai Jiao Tong University  
Shanghai, China  
Email: hetengjiao@sjtu.edu.cn

**Wei Guo**  
College of Meteorology and Oceanography  
National University of Defense Technology  
Changsha, China  
Email: guowei23@nudt.edu.cn


---

## License

This code is provided as supplementary material to the publication. It is distributed "as-is" without warranty of any kind. Users are free to use and modify the code for research purposes, provided that proper citation is given to the original publication.


---

**Last updated**: January 2026
