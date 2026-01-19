%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROFILE_UPD_MULTILAYER - Update profile operator K for multilayer medium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Citation:
% He, T., & Guo, W. (2026). Wide-angle multimodal parabolic equations: 
% modeling of directional sound propagation in stratified range-dependent 
% fluid waveguides. Proceedings of the Royal Society A,.
%
% Description:
% Computes the K operator (Eq. 2.23-2.26 in paper) which represents the 
% squared wavenumber operator ω²/c²(z) in the parabolic equation, 
% accounting for multilayer structure with density discontinuities.
%
% Paper References:
% - Eq. (2.23): K = C - E
% - Eq. (2.24): C contains depth derivative terms
% - Eq. (2.26): E contains refraction terms from ς²(z)
% - Eq. (2.38): Extension to J-layer seabed
%
% Last modified: January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = profile_upd_multilayer(grid, rho, c, atten, h, m, env, pml, omega)
% PROFILE_UPD_MULTILAYER Update profile operator for multilayer medium
%
% Computes K operator in projection space (Eq. 2.23-2.26 in paper):
%   K = C - E
% where C contains diffraction terms and E contains refraction terms.
%
% Inputs:
%   grid  - Grid structure with depth discretization
%           .z    - Depth grid [0, H]
%           .H    - Total depth including PML
%           .D    - Physical domain depth
%   rho   - Density profile ρ(z) [nz x 1] [kg/m³]
%   c     - Sound speed profile c(z) [nz x 1] [m/s]
%   atten - Attenuation profile α(z) [nz x 1] [dB/wavelength]
%   h     - Layer interface depths [numLayer-1 x 1] [m]
%   m     - Mode indices [N x 1]
%   env   - Environment structure
%           .numLayer - Number of layers J
%   pml   - PML structure with precomputed operators
%           .K - PML contribution to K operator
%   omega - Angular frequency ω = 2πf [rad/s]
%
% Output:
%   K - Profile operator matrix [N x N]
%       Hermitian matrix in projection space
%
% Features:
%   - Vectorized basis function computation
%   - Analytical integration for efficiency
%   - Handles density discontinuities at interfaces
%   - Complex sound speed for attenuation

%% Extract grid parameters
z = grid.z.';      % Depth coordinate [nz x 1]
H = grid.H;        % Total depth (including PML)
D = grid.D;        % Physical domain depth
N = numel(m);      % Number of basis modes M

% Precompute constants for efficiency
N_intz = N/2*5;  % Integration points: 5 per mode period
basis_coeff = sqrt(2/H);  % √(2/H) from Eq. (2.6)
pi_over_H = pi/H;         % π/H

% Extended layer boundaries: h₁=0, h₂, ..., h_J, h_{J+1}=D
% [Paper: Eq. 2.38]
h_ext = [0; h; D];

% Initialize operator
K = 0;

%% Loop over layers
% [Paper: Eq. 2.38] Sum over all J layers
for idx_layer = 1:env.numLayer
    
    % Layer boundaries
    h_lower = h_ext(idx_layer);        % Lower interface z = hⱼ
    h_upper = h_ext(idx_layer + 1);    % Upper interface z = hⱼ₊₁
    
    % Layer mask and mean density
    % [Paper: Section 2(c)] Density ρⱼ assumed constant within each layer
    layer_mask = (z <= h_upper) & (z >= h_lower);
    rho_layer = mean(rho(layer_mask));
    
    % Skip if no grid points in layer (numerical safety)
    if isnan(rho_layer) || rho_layer == 0
        continue;
    end
    
    % Number of Clenshaw-Curtis integration points for this layer
    layer_thickness = h_upper - h_lower;
    N_layer = floor(layer_thickness/H * N_intz) + 2;
    
    % Clenshaw-Curtis quadrature points and weights
    [z_wcc, wcc] = ClenshawCurtis(N_layer, [h_upper h_lower]);
    
    % Interpolate environmental parameters to quadrature points
    c_wcc = interp1(z, c, z_wcc, 'nearest');
    atten_wcc = interp1(z, atten, z_wcc, 'nearest');
    
    % Complex sound speed (include attenuation)
    % c_complex = c/(1 + iα/54.575)
    c_wcc = c_wcc ./ (1 + 1i*(atten_wcc)/54.575);
    
    % Squared wavenumber k²(z) = (ω/c(z))²
    % Index of refraction ς² = (c₀/c)²
    k_wcc = (omega./c_wcc).^2;
    
    %% Compute basis functions at quadrature points
    % [Paper: Eq. 2.6] φₘ(z) = √(2/H)sin(mπz/H)
    % Vectorized computation: psi is [N x N_layer]
    psi = basis_coeff * sin(m * pi_over_H .* z_wcc);
    
    %% Weighted integration (optimized)
    % Combine weights with wavenumber: [1 x N_layer]
    weighted_k = wcc.' .* k_wcc;
    
    % Matrix assembly: C += (1/ρ) * 0.5 * dh * ∫ k²φφᵀ dz
    % [Paper: Eq. 2.24, 2.38]
    % psi .* weighted_k broadcasts efficiently
    weighted_psi = psi .* weighted_k;  % [N x N_layer]
    
    % Add contribution to K operator
    % K += (1/ρⱼ) * 0.5 * (hⱼ₊₁-hⱼ) * (weighted_psi * psiᵀ)
    K = K + (0.5 * layer_thickness / rho_layer) * (weighted_psi * psi.');
    
end

% Add PML contribution
% [Paper: Eq. 2.24, Section 2(c)]
% PML integral from D to H with complex coordinate stretching
K = K + pml.K;

end