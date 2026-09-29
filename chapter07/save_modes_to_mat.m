%% ============================================================
%  FUNCTION: save_modes_to_mat
%  Saves all mode fields + mesh info to waveguide_modes.mat
%
%  The .mat file contains:
%    mesh          - struct with nodes, elems, epsilon_r
%    mode_fields   - struct array, one entry per mode, with:
%                      n_eff, beta, te_fraction
%                      Ex, Ey, Ez, Hx, Hy, Hz  (at element centroids)
%                      cx, cy                  (centroid coordinates µm)
%                      Et_dof, Ez_dof          (raw DOF vectors)
%    sim_params    - struct with wavelength, k0, geometry
%% ============================================================
function save_modes_to_mat(modes, wavelength, w_core, h_core, filename)

    % ---- Mesh (same for all modes) ----
    mesh.nodes     = modes(1).nodes;       % [Nn x 2]  µm
    mesh.elems     = modes(1).elems;       % [Ne x 3]  1-based
    mesh.epsilon_r = modes(1).epsilon_r;   % [Ne x 1]
    mesh.edges     = modes(1).edges;       % [Nedge x 2]

    % ---- Simulation parameters ----
    sim_params.wavelength_um = wavelength;
    sim_params.k0_per_um     = modes(1).k0;
    sim_params.omega_rad_s   = modes(1).omega;
    sim_params.w_core_um     = w_core;
    sim_params.h_core_um     = h_core;
    sim_params.num_modes     = length(modes);
    sim_params.units         = 'All lengths in micrometres (µm)';
    sim_params.field_note    = ['Fields evaluated at element centroids. ', ...
        'Normalised so that max(|E|)=1 over all centroids of mode 1.'];

    % ---- Mode fields ----
    mode_fields = struct();
    for m = 1:length(modes)
        md = modes(m);
        mode_fields(m).mode_index   = m;
        mode_fields(m).n_eff        = md.n_eff;
        mode_fields(m).beta_per_um  = md.k;          % complex β (1/µm)
        mode_fields(m).te_fraction  = md.te_fraction;
        mode_fields(m).tm_fraction  = md.tm_fraction;

        % Electric field components at centroids [Ne x 1] complex
        mode_fields(m).Ex = md.Ex_centroid;
        mode_fields(m).Ey = md.Ey_centroid;
        mode_fields(m).Ez = md.Ez_centroid;

        % Magnetic field components at centroids [Ne x 1] complex
        mode_fields(m).Hx = md.Hx_centroid;
        mode_fields(m).Hy = md.Hy_centroid;
        mode_fields(m).Hz = md.Hz_centroid;

        % Element centroid coordinates [Ne x 1]
        mode_fields(m).cx = md.cx;   % µm
        mode_fields(m).cy = md.cy;   % µm

        % Raw DOF vectors (useful for post-processing overlap integrals)
        mode_fields(m).Et_dof  = md.Et_dof;    % [Nedge x 1]
        mode_fields(m).Ez_dof  = md.Ez_dof;    % [Nn x 1]
    end

    % ---- Save ----
    fname = [filename '.mat'];
    save(fname, 'mesh', 'mode_fields', 'sim_params', '-v7.3');
    fprintf('Saved: %s\n', fname);
    fprintf('  Contents:\n');
    fprintf('    mesh         – nodes (%dx%d), elems (%dx%d), epsilon_r (%dx1)\n', ...
        size(mesh.nodes,1), size(mesh.nodes,2), ...
        size(mesh.elems,1), size(mesh.elems,2), size(mesh.epsilon_r,1));
    fprintf('    mode_fields  – %d modes, each with Ex/Ey/Ez/Hx/Hy/Hz at %d centroids\n', ...
        length(mode_fields), size(mesh.elems,1));
    fprintf('    sim_params   – wavelength=%.3f µm, k0=%.4f µm^-1\n', ...
        wavelength, modes(1).k0);
end