%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOPO_UPD_MULTILAYER - Update topography operators for multilayer medium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Citation:
% He, T., & Guo, W. (2026). Wide-angle multimodal parabolic equations: 
% modeling of directional sound propagation in stratified range-dependent 
% fluid waveguides. Proceedings of the Royal Society A
%
% Description:
% Computes the C (second derivative) and L (first derivative) operators
% that account for range-varying topography in the parabolic equation.
%
% Key Feature:
% Analytically handles density discontinuities at stratified interfaces
% through integration by parts (Eq. 2.11 in paper), eliminating the need
% for numerical treatment of boundary conditions.
%
% Paper References:
% - Eq. (2.9): A matrix (mass matrix)
% - Eq. (2.10): B = C + K matrix
% - Eq. (2.11): Integration by parts eliminating interface terms
% - Eq. (2.22): Analytical expressions for C
% - Eq. (2.37): Extension to J-layer seabed
%
% Last modified: January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, L] = topo_upd_multilayer(grid, rho, h, m, env, pml)
% TOPO_UPD_MULTILAYER Update topography operators for multilayer medium
%
% Computes C and L operators in projection space (Eq. 2.9-2.10, 2.22):
%   C = second derivative operator (diffraction)
%   L = first derivative operator (topography gradient)
%
% Inputs:
%   grid - Grid structure with depth discretization
%          .z    - Depth grid [0, H]
%          .H    - Total depth including PML
%          .D    - Physical domain depth
%   rho  - Density profile ρ(z) [nz x 1] [kg/m³]
%   h    - Layer interface depths [numLayer-1 x 1] [m]
%   m    - Mode indices [N x 1]
%   env  - Environment structure
%          .numLayer - Number of layers J
%   pml  - PML structure with precomputed operators
%          .C - PML contribution to C operator
%          .L - PML contribution to L operator
%
% Outputs:
%   C - Second derivative operator matrix [N x N]
%       Related to ∂²/∂z² in projection space
%   L - First derivative operator matrix [N x N]
%       Related to ∂/∂z in projection space
%
% Features:
%   - Precomputed coefficient matrices (only depend on m, n)
%   - Vectorized sinc computations
%   - Analytical treatment of interface discontinuities
%   - Efficient matrix operations

%% Extract grid parameters
z = grid.z.';      % Depth coordinate [nz x 1]
H = grid.H;        % Total depth H (including PML)
D = grid.D;        % Physical domain depth D
n = m.';           % Row vector of mode indices

% Precompute coefficient matrices (only depend on mode indices)
% These are constant for a given N and can be reused
% [Paper: Eq. 2.22, 2.25]
alp  = m + n;      % [N x N] sum of mode indices (m+n)
beta = m - n;      % [N x N] difference of mode indices (m-n)
E    = (pi^2/H^3) * (m .* n);  % [N x N] coefficient matrix π²mn/H³

% Extended layer boundaries: h₁=0, h₂, ..., h_J, h_{J+1}=D
% [Paper: Eq. 2.37]
h_ext = [0; h; D];

% Precompute constants
inv_H = 1/H;  % 1/H

% Initialize operators
C = 0;
L = 0;

%% Loop over layers
% [Paper: Eq. 2.37] Sum over all J layers
for idx_layer = 1:env.numLayer
    
    % Layer boundaries
    h_lower = h_ext(idx_layer);        % Lower interface z = hⱼ
    h_upper = h_ext(idx_layer + 1);    % Upper interface z = hⱼ₊₁
    
    % Layer mask and mean density
    % [Paper: Section 2(c)] Density ρⱼ assumed constant within each layer
    layer_mask = (z < h_upper) & (z >= h_lower);
    rho_layer = mean(rho(layer_mask));
    
    % Skip if no grid points in layer (numerical safety)
    if isnan(rho_layer) || rho_layer == 0
        continue;
    end
    
    % Precompute inverse density
    inv_rho = 1/rho_layer;
    
    %% Vectorized sinc computation
    % sinc(x) = sin(πx)/(πx) in MATLAB
    % Arguments are already scaled by H
    
    % Lower boundary sinc values [N x N]
    % [Paper: Eq. 2.22] sinc((m±n)hⱼ/H)
    sincA = sinc(alp * h_lower * inv_H);
    sincB = sinc(beta * h_lower * inv_H);
    
    % Upper boundary sinc values [N x N]
    % [Paper: Eq. 2.22] sinc((m±n)hⱼ₊₁/H)
    sincC = sinc(alp * h_upper * inv_H);
    sincD = sinc(beta * h_upper * inv_H);
    
    %% Compute F and G matrices
    % These represent integrated products of basis functions
    % [Paper: Eq. 2.22]
    % F and G arise from integration by parts
    F = h_upper*sincC - h_lower*sincA;  % [N x N]
    G = h_upper*sincD - h_lower*sincB;  % [N x N]
    
    %% Accumulate operators (vectorized)
    % Second derivative operator C
    % [Paper: Eq. 2.22, first equation]
    % C = -∑ⱼ (1/ρⱼ) * (π²mn/H³) * (F + G)
    C = C - inv_rho * E .* (F + G);
    
    % First derivative operator L
    % [Paper: Eq. 2.22, second equation]
    % L = -∑ⱼ (1/ρⱼ) * (1/H) * (F - G)
    L = L - inv_rho * inv_H * (F - G);
    
end

% Add PML contributions
% [Paper: Eq. 2.21, 2.24]
% PML integral from D to H with complex coordinate stretching
C = C + pml.C;
L = L + pml.L;

end