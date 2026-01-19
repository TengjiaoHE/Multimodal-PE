function [grid, env, casename] = read_env(env_file)

%READ_ENV  Read a Multimodal-PE environment input file (*.env).
%
%   [grid, env, casename] = read_env(env_file)
%
%   This parser reads a plain-text environment file describing a stratified,
%   range-dependent FLUID waveguide for the wide-angle Multimodal-PE solver.
%   The file supports:
%     - Multiple layers (env.numLayer)
%     - Range-dependent profiles specified at env.numRange range nodes (env.rg)
%     - Piecewise-linear interface topography between adjacent layers
%
%   Comment rule:
%     Any text following a double backslash \\ is treated as a comment and
%     removed (both full-line and inline comments).
%
%   Expected file layout (high level):
%     Line 1:  Case name (string)
%     Next:    numLayer, numRange, nprofile(1..numLayer)
%     Next:    c0, freq, zs, zr, rmax, dr, zmax, dz, dpml, npml, nP
%     Next:    env.rg (numRange values)
%     Next:    For each layer i=1..numLayer, provide nprofile(i) depth samples.
%              Each depth sample line has (1 + numRange + 2) numbers:
%                 z  c(r1) c(r2) ... c(rN)  rho  alpha
%              where N = numRange and r1..rN correspond to env.rg.
%     Next:    env.toponr (number of control points for interface topography)
%     Next:    For each interface i=1..numLayer-1, provide 2*env.toponr numbers
%              in the order: r1..rM then z1..zM (often written as two lines).
%
%   Outputs:
%     casename : case title string
%     grid     : struct of numerical/grid parameters (zs, zr, rmax, dr, zmax, dz, dpml, npml, nP)
%     env      : struct of environment profiles and geometry (rg, zg/ssp/rho/atten, topography)


    [fid, cleanup_tmp] = fopen_skip_comments(env_file); % open temp file with \\ comments stripped
    casename      = fgetl(fid); % case name / title (first line in .env)
    env.numLayer  = fscanf(fid, '%d', 1); % number of layers
    env.numRange  = fscanf(fid, '%d', 1); % number of range nodes for range-dependent profiles
    nprofile      = fscanf(fid, '%d', env.numLayer); 

    env.c0         = fscanf(fid, '%f', 1); % reference sound speed [m/s]
    env.freq       = fscanf(fid, '%f', 1); % frequency [Hz]
    grid.zs        = fscanf(fid, '%f', 1); % source depth [m]
    grid.zr        = fscanf(fid, '%f', 1); % receiver depth [m]
    grid.rmax      = fscanf(fid, '%f', 1); % maximum range [m]
    grid.dr        = fscanf(fid, '%f', 1); % range step [m]
    grid.zmax      = fscanf(fid, '%f', 1); % maximum depth in computation [m]
    grid.dz        = fscanf(fid, '%f', 1); % depth step [m]
    grid.dpml      = fscanf(fid, '%f', 1); % PML thickness [m]
    grid.npml      = fscanf(fid, '%f', 1); % number of PML grid points (discretisation)
    grid.nP        = fscanf(fid, '%f', 1); % Pade order (wide-angle)

    env.rg         = fscanf(fid, '%f', env.numRange ); 
    
    env.zg         = cell(env.numLayer,1);
    env.rho        = cell(env.numLayer,1);
    env.ssp        = cell(env.numLayer,1); 
    env.atten      = cell(env.numLayer,1);
    

    % ---------------------------------------------------------------------
    % Layer profiles:
    %   For each layer i, we read a matrix with size (3+numRange) x nprofile(i).
    %   For each depth sample, the expected columns are:
    %       z, c(r1..rN), rho, alpha
    %   where N = env.numRange and r1..rN correspond to env.rg.
    % ---------------------------------------------------------------------

    for i = 1 : env.numLayer

            profile     = fscanf(fid, '%f %f', [3+env.numRange , nprofile(i)]);
            env.zg(i)   = {profile(1, 1:nprofile(i))};
            env.ssp(i)  = {profile(2:2+env.numRange -1, 1:nprofile(i))};
            env.rho(i)  = {profile(2+env.numRange , 1:nprofile(i))};
            env.atten(i)= {profile(2+env.numRange +1, 1:nprofile(i))};

    end

    env.toponr     = fscanf(fid, '%f', 1); % number of control points for each interface topography


    % ---------------------------------------------------------------------
    % Topography / interface depths:
    %   The reader expects 2*env.toponr numbers per interface in this order:
    %     r1 ... rM  z1 ... zM
    %   (This is often provided as two lines: first ranges, then depths.)
    % ---------------------------------------------------------------------

    for i = 1 : env.numLayer-1

            topo        = fscanf(fid, '%f %f', [env.toponr , 2]);
            env.topor(i)= {topo(1:env.toponr, 1)};
            env.topoz(i)= {topo(1:env.toponr, 2)};

    end

    fclose(fid);

end


function [fid, cleanupObj] = fopen_skip_comments(srcFile)
    % Remove comments marked by the delimiter '\\' (double backslash): full-line or inline.
    tmp = [tempname, '.txt'];
    fin = fopen(srcFile, 'r');
    if fin < 0, error('Cannot open source file: %s', srcFile); end
    fout = fopen(tmp, 'w');
    if fout < 0, fclose(fin); error('Cannot create temp file.'); end

    while true
        ln = fgetl(fin);
        if ~ischar(ln), break; end
        idx = strfind(ln, '\\');          % match comment delimiter '\\'
        if ~isempty(idx), ln = ln(1:idx(1)-1); end
        fprintf(fout, '%s\n', ln);        % preserve line structure (including empty lines)
    end

    fclose(fin); fclose(fout);
    fid = fopen(tmp, 'r');
    if fid < 0, error('Internal error: cannot reopen temp file.'); end
    cleanupObj = onCleanup(@() delete_if_exists(tmp));
end

function delete_if_exists(p)
    if exist(p, 'file'), try, delete(p); end, end
end
