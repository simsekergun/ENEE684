%=========================================================
%  phase_vs_group_velocity.m
%
%  Animates a Gaussian wave packet propagating in a dispersive
%  medium where v_phase ≠ v_group, replicating the classic
%  figure showing E(z,t) and its envelope F(z,t).
%
%  Physics:
%    E(z,t) = F(z,t) * cos(k0*z - omega0*t)   [full field]
%    F(z,t) = exp( -(z - v_g*t)^2 / (2*sigma^2) ) [envelope]
%    v_p = omega0 / k0   (phase velocity — carrier fringes)
%    v_g = d(omega)/dk   (group velocity — envelope peak)
%% =========================================================
clear;%  clc; close all;
 
%% ---- Mac / renderer fix ----------------------------------
set(0, 'DefaultFigureRenderer', 'painters');
 
%% =========================================================
%  USER PARAMETERS  — edit here
%% =========================================================
k0      = 8;        % carrier wavenumber [rad/m]
omega0  = 12;       % carrier angular frequency [rad/s]
v_g     = 0.8;      % group velocity [m/s]   (envelope speed)
% NOTE: v_p = omega0/k0 = 1.5 m/s  — set by k0 and omega0 above.
%       Change k0 or omega0 to adjust v_p.
%       v_p > v_g  → fringes overtake the envelope (shown below)
%       v_p < v_g  → envelope overtakes the fringes (swap values)
 
sigma   = 2.5;      % Gaussian envelope half-width [m]
t_end   = 140;       % total animation time [s]
dt      = 0.01;     % time step [s]  (smaller = smoother, slower)
 
SAVE_GIF      = false;           % true → write 'phase_vs_group.gif'
GIF_FILENAME  = 'phase_vs_group_velocity.gif';
GIF_DELAY     = dt;              % frame delay in GIF [s]
 
%% =========================================================
%  Derived quantities
%% =========================================================
v_p     = omega0 / k0;           % phase velocity
t_vec   = 0 : dt : t_end;
Nt      = numel(t_vec);
z_half  = 15;                    % half-width of displayed window
Nz      = 6000;
 
%% =========================================================
%  Colours
%% =========================================================
C_field  = [0.10 0.10 0.10];   % near-black — full field E
C_env    = [0.00 0.45 0.74];   % blue       — envelope F
C_vg     = [0.00 0.45 0.74];   % blue       — group velocity marker
C_vp     = [0.85 0.15 0.10];   % red        — phase velocity marker
C_trail  = [0.70 0.70 0.70];   % grey       — position trails
 
%% =========================================================
%  Figure layout
%% =========================================================
fig = figure('Name','Phase vs Group Velocity', ...
             'NumberTitle','off', ...
             'Position',[60 60 1050 600], ...
             'Color','w');
 
%% ---- Top panel: wave packet (moving window) --------------
ax1 = subplot(3,1,[1 2]);
hold(ax1,'on');  grid(ax1,'on');  box(ax1,'on');
ax1.FontSize = 11;
ax1.XTickLabel = [];           % shared x-axis with bottom panel
 
% Initial data
t0    = t_vec(1);
z_ctr = v_g * t0;             % envelope centre
z     = linspace(z_ctr - z_half, z_ctr + z_half, Nz);
F0    = gaussian_env(z, t0, v_g, sigma);
E0    = F0 .* cos(k0*z - omega0*t0);
 
h_E      = plot(ax1, z, E0,  '-',  'Color', C_field, 'LineWidth', 1.0);
h_Fp     = plot(ax1, z,  F0, '--', 'Color', C_env,   'LineWidth', 2.2);
h_Fn     = plot(ax1, z, -F0, '--', 'Color', C_env,   'LineWidth', 2.2);
 
% Envelope peak marker (blue square)
h_vg_mk  = plot(ax1, z_ctr, 1.0,  's', 'Color', C_vg, ...
                'MarkerSize', 11, 'MarkerFaceColor', C_vg);
 
% Phase crest marker (red triangle)
z_pc0    = nearest_phase_crest(k0, omega0, t0, z_ctr);
h_vp_mk  = plot(ax1, z_pc0, cos(k0*z_pc0 - omega0*t0), 'v', ...
                'Color', C_vp, 'MarkerSize', 11, 'MarkerFaceColor', C_vp);
 
% Velocity arrows using annotation (positioned relative to axes)
ylabel(ax1, 'Amplitude', 'FontSize', 12);
ax1.YLim = [-1.35 1.55];
 
% Text labels that move with the markers
h_vg_lbl = text(ax1, z_ctr + 0.4, 1.22, ...
    sprintf('\\bf v_g = %.2f m/s', v_g), ...
    'Color', C_vg, 'FontSize', 10.5, 'FontWeight','bold');
h_vp_lbl = text(ax1, z_pc0 + 0.4, cos(k0*z_pc0-omega0*t0) - 0.22, ...
    sprintf('\\bf v_p = %.2f m/s', v_p), ...
    'Color', C_vp, 'FontSize', 10.5, 'FontWeight','bold');
 
% Annotations: labels for E and F
h_E_lbl  = text(ax1, z_ctr - z_half*0.72, 0.85, 'E(z,t)', ...
    'FontSize', 13, 'FontAngle','italic', 'Color', C_field);
h_F_lbl  = text(ax1, z_ctr + z_half*0.55, 0.78, 'F(z,t)', ...
    'FontSize', 13, 'FontAngle','italic', 'Color', C_env);
 
h_title  = title(ax1, ...
    sprintf('Gaussian Wave Packet   |   v_p / v_g = %.2f   |   t = %.2f s', ...
            v_p/v_g, t0), 'FontSize', 12, 'FontWeight','bold');
 
legend(ax1, [h_E, h_Fp, h_vg_mk, h_vp_mk], ...
    {'E(z,t) — full field', 'F(z,t) — envelope (\pmF)', ...
     sprintf('Envelope peak  (v_g = %.2f m/s)', v_g), ...
     sprintf('Phase crest    (v_p = %.2f m/s)', v_p)}, ...
    'Location','northeast','FontSize',9.5,'NumColumns',2);
 
%% ---- Bottom panel: position trails vs time ---------------
ax2 = subplot(3,1,3);
hold(ax2,'on');  grid(ax2,'on');  box(ax2,'on');
ax2.FontSize = 11;
 
% Pre-allocate trail history
z_vg_hist = NaN(1, Nt);
z_vp_hist = NaN(1, Nt);
t_hist     = NaN(1, Nt);
 
h_vg_trail = plot(ax2, NaN, NaN, '-', 'Color', C_vg, 'LineWidth', 2);
h_vp_trail = plot(ax2, NaN, NaN, '-', 'Color', C_vp, 'LineWidth', 2);
 
% Ideal slope lines (drawn once)
t_line = [0, t_end];
plot(ax2, t_line, v_g * t_line, '--', 'Color', C_vg, 'LineWidth', 0.8);
plot(ax2, t_line, v_p * t_line, '--', 'Color', C_vp, 'LineWidth', 0.8);
 
xlabel(ax2, 'Time  t  [s]',     'FontSize', 12);
ylabel(ax2, 'Position  z  [m]', 'FontSize', 12);
title(ax2, 'Position Trails of Envelope Peak (v_g) and Phase Crest (v_p)', ...
      'FontSize', 11);
ax2.XLim = [0, t_end];
ax2.YLim = [-2, v_p * t_end * 1.05];
legend(ax2, [h_vg_trail, h_vp_trail], ...
    {sprintf('Envelope peak  v_g = %.2f m/s', v_g), ...
     sprintf('Phase crest    v_p = %.2f m/s', v_p)}, ...
    'Location','northwest','FontSize',9.5);
 
%% =========================================================
%  Animation loop
%% =========================================================
for n = 1:Nt
    t = t_vec(n);
    if ~ishandle(fig), break; end   % exit if window closed
 
    % Moving window centred on envelope peak
    z_ctr = v_g * t;
    z     = linspace(z_ctr - z_half, z_ctr + z_half, Nz);
 
    % Field and envelope
    F     = gaussian_env(z, t, v_g, sigma);
    E     = F .* cos(k0*z - omega0*t);
 
    % Phase crest nearest to envelope centre
    z_pc  = nearest_phase_crest(k0, omega0, t, z_ctr);
    % Clamp to window
    z_pc  = max(z(1)+0.5, min(z(end)-0.5, z_pc));
 
    %% Update top panel
    set(h_E,  'XData', z, 'YData', E);
    set(h_Fp, 'XData', z, 'YData', F);
    set(h_Fn, 'XData', z, 'YData', -F);
 
    set(h_vg_mk, 'XData', z_ctr, 'YData', 1.0);
    set(h_vp_mk, 'XData', z_pc,  'YData', cos(k0*z_pc - omega0*t));
 
    ax1.XLim = [z(1), z(end)];
 
    set(h_vg_lbl, 'Position', [z_ctr + 0.4, 1.22, 0]);
    set(h_vp_lbl, 'Position', [z_pc  + 0.4, cos(k0*z_pc-omega0*t) - 0.22, 0]);
    set(h_E_lbl,  'Position', [z_ctr - z_half*0.72, 0.85, 0]);
    set(h_F_lbl,  'Position', [z_ctr + z_half*0.55, 0.78, 0]);
 
    set(h_title,  'String', ...
        sprintf('Gaussian Wave Packet   |   v_p / v_g = %.2f   |   t = %.2f s', ...
                v_p/v_g, t));
 
    %% Update bottom panel (position trails)
    z_vg_hist(n) = z_ctr;
    z_vp_hist(n) = z_pc;
    t_hist(n)    = t;
 
    set(h_vg_trail, 'XData', t_hist(1:n),     'YData', z_vg_hist(1:n));
    set(h_vp_trail, 'XData', t_hist(1:n),     'YData', z_vp_hist(1:n));
 
    drawnow limitrate;
 
    %% Optional: save GIF
    if SAVE_GIF
        frame = getframe(fig);
        im    = frame2im(frame);
        [A, map] = rgb2ind(im, 256);
        if n == 1
            imwrite(A, map, GIF_FILENAME, 'gif', ...
                    'LoopCount', Inf, 'DelayTime', GIF_DELAY);
        else
            imwrite(A, map, GIF_FILENAME, 'gif', ...
                    'WriteMode','append', 'DelayTime', GIF_DELAY);
        end
    end
end
 
%% ========================================================= 
function F = gaussian_env(z, t, v_g, sigma)
% GAUSSIAN_ENV  Gaussian envelope centred at v_g*t.
    F = exp(-(z - v_g*t).^2 / (2*sigma^2));
end
 
function z_crest = nearest_phase_crest(k0, omega0, t, z_ref)
% NEAREST_PHASE_CREST  Returns z of the cosine maximum nearest to z_ref.
%   Crests occur where  k0*z - omega0*t = 2*pi*n  =>  z = (2*pi*n + omega0*t)/k0
    n_approx = round((k0*z_ref - omega0*t) / (2*pi));
    z_crest  = (2*pi*n_approx + omega0*t) / k0;
end
 