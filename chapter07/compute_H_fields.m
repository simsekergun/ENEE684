function modes = compute_H_fields(modes, varargin)
% COMPUTE_H_FIELDS  Reconstruct magnetic field components from a solved
%                   FEM mode struct and optionally plot Hx, Hy, Hz.
%
% ── Physics ────────────────
%   From Faraday's law in the time-harmonic convention e^(j*omega*t):
%
%     H = (1 / j*omega*mu0*mu_r) * curl(E)
%
%   With the waveguide ansatz  E(x,y,z) = Ẽ(x,y) * exp(j*beta*z),
%   all d/dz -> j*beta, giving:
%
%     Hx = prefac * ( j*beta*Ey  -  dEz/dy )
%     Hy = prefac * ( dEz/dx     -  j*beta*Ex )
%     Hz = prefac * ( dEy/dx     -  dEx/dy )
%
%   where  prefac = 1 / (j * omega * mu0 * mu_r).
%
%   For P1 (nodal Lagrange) elements Ez has a constant gradient over each
%   triangle, computed from Ez_dof nodal values and the element's
%   barycentric-coordinate gradients.
%
%   For Nedelec edge elements the curl components dEy/dx and dEx/dy are
%   also element-constant, computed analytically from the basis functions.
%
% ── Syntax ─────────────────
%   modes = compute_H_fields(modes)
%   modes = compute_H_fields(modes, Name, Value, ...)
%
% ── Optional parameters ────
%   'mu_r'        relative permeability  (default 1.0)
%   'plot_modes'  indices of modes to plot, e.g. [1 2 3]
%                 0 or []  -> no plots  (default [])
%                 'all'    -> plot every mode
%   'save_figs'   true/false — save each figure as PNG  (default false)
%
% ── Added fields to each mode struct ──────────────────────────────────────
%   .Hx_centroid  [Ne x 1] complex  — Hx at element centroids
%   .Hy_centroid  [Ne x 1] complex
%   .Hz_centroid  [Ne x 1] complex
%   .Ex_centroid  [Ne x 1] complex  — Ex (also stored for convenience)
%   .Ey_centroid  [Ne x 1] complex
%   .Ez_centroid  [Ne x 1] complex
%   .cx           [Ne x 1]          — centroid x coordinates  (µm)
%   .cy           [Ne x 1]          — centroid y coordinates  (µm)

%% ── Parse optional arguments ──────────────────────────────────────────────
p = inputParser;
addParameter(p, 'mu_r',       1.0,   @isnumeric);
addParameter(p, 'plot_modes', [],    @(x) isnumeric(x) || ischar(x));
addParameter(p, 'save_figs',  false, @islogical);
parse(p, varargin{:});
opt = p.Results;

% Resolve which mode indices to plot
Nmodes = numel(modes);
if ischar(opt.plot_modes) && strcmpi(opt.plot_modes, 'all')
    plot_idx = 1:Nmodes;
elseif isnumeric(opt.plot_modes) && ~isempty(opt.plot_modes) && any(opt.plot_modes > 0)
    plot_idx = opt.plot_modes(opt.plot_modes > 0);
else
    plot_idx = [];
end

%% ── Physical constants (µm unit system) ─────────────────────────────────
mu0 = 4*pi * 1e-13;   % H/µm

%% ── Local edge-node pair table ────────────────────────────────────────────
local_pairs = [2 3; 3 1; 1 2];  % local node indices for edge k

%% ── Loop over modes ────────
for m_idx = 1:Nmodes

    md  = modes(m_idx);
    Ne  = size(md.elems, 1);

    beta   = md.k;                                % complex β  (1/µm)
    prefac = 1 / (1j * md.omega * mu0 * opt.mu_r);

    Hx_c = zeros(Ne,1,'like',1j);
    Hy_c = zeros(Ne,1,'like',1j);
    Hz_c = zeros(Ne,1,'like',1j);
    Ex_c = zeros(Ne,1,'like',1j);
    Ey_c = zeros(Ne,1,'like',1j);
    Ez_c = zeros(Ne,1,'like',1j);
    Sz_c = zeros(Ne,1,'like',1j);
    
    for ie = 1:Ne
        nds = md.elems(ie,:);           % [1x3] global node indices
        xy  = md.nodes(nds,:);          % [3x2] node coordinates
        edg = md.elem2edge(ie,:);       % [1x3] global edge DOF indices
        sgn = md.edge_sign(ie,:);       % [1x3] ±1

        % Barycentric gradients (constant over element)
        [~, detJ, grads] = element_geometry(xy);   % grads [3x2]
        elem_area = abs(detJ)/2;
        
        % Edge lengths
        edge_len = zeros(1,3);
        for k = 1:3
            ni = local_pairs(k,1);  nj = local_pairs(k,2);
            edge_len(k) = norm(xy(nj,:) - xy(ni,:));
        end

        % Centroid barycentric coordinates
        lam_c = [1/3, 1/3, 1/3];

        % ── Accumulate E-field and its relevant derivatives ───────────────
        Ex_q   = 0;   Ey_q   = 0;
        dEx_dy = 0;   dEy_dx = 0;   % needed for Hz

        for k = 1:3
            ni  = local_pairs(k,1);
            nj  = local_pairs(k,2);

            % Nedelec vector basis at centroid: W = s*l*(lam_i*grad_j - lam_j*grad_i)
            Wk = sgn(k) * edge_len(k) * ...
                 (lam_c(ni)*grads(nj,:) - lam_c(nj)*grads(ni,:));  % [1x2]

            coeff  = md.Et_dof(edg(k));
            Ex_q   = Ex_q + coeff * Wk(1);
            Ey_q   = Ey_q + coeff * Wk(2);

            % Curl of Nedelec basis (element-constant):
            gi = grads(ni,:);
            gj = grads(nj,:);
            dWkx_dy = sgn(k)*edge_len(k)*(gi(2)*gj(1) - gj(2)*gi(1));
            dWky_dx = sgn(k)*edge_len(k)*(gi(1)*gj(2) - gj(1)*gi(2));

            dEx_dy = dEx_dy + coeff * dWkx_dy;
            dEy_dx = dEy_dx + coeff * dWky_dx;
        end

        % Ez at centroid (average of P1 nodal values)
        Ez_q = mean(md.Ez_dof(nds));

        % Gradient of Ez (constant for P1 elements)
        %   Ez = sum_k  Ez_k * phi_k  =>  grad(Ez) = sum_k  Ez_k * grads(k,:)
        dEz_dx = md.Ez_dof(nds(1))*grads(1,1) + ...
                 md.Ez_dof(nds(2))*grads(2,1) + ...
                 md.Ez_dof(nds(3))*grads(3,1);
        dEz_dy = md.Ez_dof(nds(1))*grads(1,2) + ...
                 md.Ez_dof(nds(2))*grads(2,2) + ...
                 md.Ez_dof(nds(3))*grads(3,2);

        % ── Faraday's law ──
        %   curl(E)_x =  j*beta*Ey  -  dEz/dy
        %   curl(E)_y =  dEz/dx     -  j*beta*Ex
        %   curl(E)_z =  dEy/dx     -  dEx/dy
        Hx_c(ie) = prefac * (1j*beta*Ey_q  - dEz_dy);
        Hy_c(ie) = prefac * (dEz_dx        - 1j*beta*Ex_q);
        Hz_c(ie) = prefac * (dEy_dx        - dEx_dy);

        Ex_c(ie) = Ex_q;
        Ey_c(ie) = Ey_q;
        Ez_c(ie) = Ez_q;
        
        Sz_c(ie) = 0.5 * real(Ex_c(ie) * conj(Hy_c(ie))- Ey_c(ie) * conj(Hx_c(ie)) )*elem_area; % multiplied with area gives the power
    end % element loop

    % ── Normalise by peak |E| of this mode ────────────────────────────────
    Normalizer = sum(Sz_c);
    if Normalizer > 0
        scale = 1 / Normalizer;
        Ex_c = Ex_c * scale;  Ey_c = Ey_c * scale;  Ez_c = Ez_c * scale;
        Hx_c = Hx_c * scale;  Hy_c = Hy_c * scale;  Hz_c = Hz_c * scale;
    end

    % ── Centroid coordinates ──────────────────────────────────────────────
    cx = (md.nodes(md.elems(:,1),1) + md.nodes(md.elems(:,2),1) + ...
          md.nodes(md.elems(:,3),1)) / 3;
    cy = (md.nodes(md.elems(:,1),2) + md.nodes(md.elems(:,2),2) + ...
          md.nodes(md.elems(:,3),2)) / 3;

    % ── Store in struct ────
    modes(m_idx).Ex_centroid = Ex_c;
    modes(m_idx).Ey_centroid = Ey_c;
    modes(m_idx).Ez_centroid = Ez_c;
    modes(m_idx).Hx_centroid = Hx_c;
    modes(m_idx).Hy_centroid = Hy_c;
    modes(m_idx).Hz_centroid = Hz_c;
    modes(m_idx).cx          = cx;
    modes(m_idx).cy          = cy;

end % mode loop

fprintf('[compute_H_fields]  H-fields computed for %d mode(s).\n', Nmodes);

%% ── Plotting ───────────────
for m_idx = plot_idx
    if m_idx > Nmodes
        warning('compute_H_fields: mode index %d out of range (only %d modes).', ...
            m_idx, Nmodes);
        continue
    end
    plot_EH_mode(modes(m_idx), m_idx, opt.save_figs);
end

end % ── compute_H_fields ─────


