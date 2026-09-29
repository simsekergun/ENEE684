%% test 1
% FEM solver for EM modes in dielectric waveguides
% Original formulation is based on: https://doi.org/10.1080/02726340290084012
%
% Usage:
%   Run the script directly for a demo 
%   at 1550 nm, or call compute_modes() from your own script.
%
% Requirements: MATLAB R2020b+ (for complex sparse eigs)

%%
clear; 
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

wavelength = 1.55; % µm
% ---- Geometry (all lengths in micrometres) ----
w_core   = 1.6;   % core width
h_core   = 0.7;   % core height
h_sim = 4.0;     % height of  thesimulation width
w_sim    = 6.00;   % total simulation width
% --- materials
n_core   = get_refractive_index('Si3N4', wavelength); % refractive index of the core
n_clad   = get_refractive_index('SiO2', wavelength);  % refractive index of the cladding

% create a string for automation
filename = ['W' int2str(round(w_core*1000)) 'H' int2str(round(h_core*1000)) 'L' int2str(round(wavelength*1000))];

num_modes  = 2;         % number of modes to search for
mesh_res   = round(w_sim/wavelength*20);   % approximate triangles along the longest side

% compute & plot settings
compute_H_field_profiles = 1;  % set it to 1, if you want to calculate the H fields
plot_mesh = 1;
save_data = 0;
% end of inputs

%% Call mesh builder
build_mesh;

%%
tic
fprintf('Assembling & solving eigenvalue problem ...\n');
modes = compute_modes(nodes, elems, epsilon_r, wavelength, w_core, h_core, save_data, ...
     compute_H_field_profiles, filename, 'num_modes', num_modes, 'mu_r', 1.0);

toc
%% ---- Print results ----
fprintf('\n--- Guided modes ---\n');
for m = 1:length(modes)
    fprintf('Mode %d:  n_eff = %.6f + %.2ei,  TE-frac = %.3f\n', ...
        m, real(modes(m).n_eff), imag(modes(m).n_eff), modes(m).te_fraction);
end
