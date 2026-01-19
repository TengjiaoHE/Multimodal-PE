%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTIMODAL PARABOLIC EQUATION (Multimodal-PE) SOLVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Citation:
% He, T., & Guo, W. (2026). Wide-angle multimodal parabolic equations: 
% modeling of directional sound propagation in stratified range-dependent 
% fluid waveguides. Proceedings of the Royal Society A,.
%
% Description:
% This code implements the Multimodal-PE method described in the paper.
% The method projects the acoustic pressure field onto a set of basis modes
% (sine functions satisfying boundary conditions at z=0 and z=H), and 
% advances the projection coefficients along range using a split-step Padé
% approximation.
%
% Key Features:
% - Handles stratified media with discontinuities in sound speed and density
% - Supports directional sources with arbitrary prescribed directivity
% - Adaptive optimization: modal propagation for range-invariant segments,
%   Padé method for range-varying segments
% - GPU acceleration and precision control (single/double)
% - Strong convergence: O(1/M^2) where M is number of basis modes
%
% Contact:
% Tengjiao He (hetengjiao@sjtu.edu.cn)
% Shanghai Jiao Tong Univerity, Shanghai, China
% Wei Guo (guowei23@nudt.edu.cn)
% National University of Defense Technology, Changsha, China
%
% License:
% This code is provided as supplementary material to the above publication.
% Use of this code should cite the paper. Provided "as-is" without warranty.
%
% Last modified: January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This demo is used to reproduce the sound field results presented in Fig.13 of the paper 

clear ; clc;

addpath(genpath('functions'));
addpath(genpath('envfiles'));
clear ; clc;

% Environment setup
[grid, env, casename] = read_env('input3.env'); % Load env file

% Source beam pattern (omnidirectional)
sbp.dir  = @(varphi, k) 1; % point source

% Run simulation
use_gpu = false;         % Set true for GPU acceleration
precision = 'double';    % 'single' or 'double'
grid = multimodal_pe(env, grid, sbp, use_gpu, precision);


% plotting results

tlmin = 30; tlmax = 80;

figure
pcolor(grid.r.'/1000,grid.z,-mag2db(abs(grid.p)));
view(2);
shading interp
colormap(flipud(jet))
set(gca,'ydir','reverse') ;
shading interp
clim([tlmin,tlmax])
hold on
plot(grid.r/1000,grid.h,'w','Linewidth',2)
xlabel('Range / km')
ylabel('Depth / m')
set(gca, 'Box', 'of', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XGrid', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1)
set(gca,'FontSize',18)
ylim([0 grid.D])
title(casename)
colorbar

