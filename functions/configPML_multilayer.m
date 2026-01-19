%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGPML_MULTILAYER - Configure Perfectly Matched Layer (PML)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Citation:
% He, T., & Guo, W. (2026). Wide-angle multimodal parabolic equations: 
% modeling of directional sound propagation in stratified range-dependent 
% fluid waveguides. Proceedings of the Royal Society A
%
% Description:
% Sets up the PML operators to absorb outgoing waves at the bottom boundary.
% Uses complex coordinate stretching to truncate the computational domain.
%
% Paper References:
% - Eq. (2.18): Complex coordinate transformation ẑ = z + i∫σ(ε)dε
% - Eq. (2.19): Derivative transformation
% - Eq. (2.20): Damping function s(z)
% - [36] Berenger (1994): Original PML concept
%
% Last modified: January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pml = configPML_multilayer(grid, m, c, rho, omega)
% CONFIGPML_MULTILAYER Set up PML operators for domain truncation
%
% Computes PML contributions to C, L, and K operators using complex
% coordinate stretching to absorb waves propagating toward infinity.
%
% [Paper: Eq. 2.18-2.20] PML formulation
% [Paper: Eq. 2.21, 2.24] PML contributions to operators
%
% Inputs:
%   grid  - Grid structure
%           .dpml - PML thickness d [m]
%           .H    - Total depth H = D + d [m]
%           .D    - Physical domain depth D [m]
%           .npml - Number of PML quadrature points
%   m     - Mode indices [N x 1]
%   c     - Sound speed in PML region [m/s]
%   rho   - Density in PML region [kg/m³]
%   omega - Angular frequency ω = 2πf [rad/s]
%
% Output:
%   pml - Structure with PML operator contributions:
%         .C - PML contribution to C operator [N x N]
%         .L - PML contribution to L operator [N x N]
%         .K - PML contribution to K operator [N x N]
%
% Note:
%   PML uses polynomial damping function (degree p=4, q=4)
%   Clenshaw-Curtis quadrature for efficient integration

%% Set-up of PMLs
d         = grid.dpml;  % PML thickness d
H         = grid.H;     % Total depth H
D         = grid.D;     % Physical domain depth D
N         = grid.npml;  % Number of quadrature points in PML

% Clenshaw-Curtis quadrature points and weights in PML [D, H]
[zcc, wcc]= ClenshawCurtis (N, [D H]);

% Sound speed in PML (assumed constant)
c_PML     = c*ones(1,N+1);

% Basis functions and derivatives at quadrature points
% [Paper: Eq. 2.6] φₘ(z) = √(2/H)sin(mπz/H)
% Derivative: φ'ₘ(z) = (mπ/H)√(2/H)cos(mπz/H)
psi_PML_p = m.*pi/H.*sqrt(2/H).*cos((m)*pi/H.*zcc);  % φ'ₘ at zcc
psi_PML   = sqrt(2/H)*sin((m)*pi/H.*zcc);            % φₘ at zcc

% Polynomial damping parameters
p         = 4;  % Degree for real part
q         = 4;  % Degree for imaginary part

% Normalized coordinate in PML: τ = (z-D)/d ∈ [0,1]
tao       = (zcc-D)./d;

% Complex damping function s(z)
% [Paper: Eq. 2.20] s(z) = 1 + iσ(z) for z ∈ [D, H]
% Here: s = exp(-pτ) - i(exp(-qτ) - 1)
s         = exp(-p*tao)-1i*(exp(-q*tao)-1);

% PML contributions to operators
% [Paper: Eq. 2.21, 2.24]

% C operator: second derivative contribution
% C_PML = (1/ρ) * 0.5 * d * ∫ (1/s) φ'ₘ φ'ₙ dz
% Matrix form: (wcc .* (-1./s) .* psi_PML_p) * psi_PML_p'
pml.C     = 1/rho*0.5*d*(wcc.'.*(-1./s).*psi_PML_p)*psi_PML_p.';

% L operator: first derivative contribution
% L_PML = (1/ρ) * 0.5 * d * ∫ s φₘ φₙ dz
pml.L     = 1/rho*0.5*d*(wcc.'.*s.*psi_PML)*psi_PML.';

% K operator: wavenumber contribution
% K_PML = (1/ρ) * 0.5 * d * ∫ s (ω²/c²) φₘ φₙ dz
pml.K     = 1/rho*0.5*d*(wcc.'.*s.*omega^2./c_PML.^2.*psi_PML)*psi_PML.';

end