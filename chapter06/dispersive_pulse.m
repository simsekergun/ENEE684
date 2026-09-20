% Dispersive pulse propagation visualizer
%  Solves A(Z,T) analytically for a Gaussian pulse in a dispersive medium.
%
%  A(Z,T) = A0/sqrt(1 + j*Z/LD) * exp( -T^2 / (2*T0^2*(1 + j*Z/LD)) )
%  LD = T0^2 / |beta2|   (dispersion length)
 
clear; clc; close all;
 
%% --- Parameters (change these to explore) ---
T0      = 0.1e-12;        % Initial pulse half-width at 1/e^2 intensity [s] (1 ps)
beta2   = 2.17e-26;        % GVD parameter [s^2/m]  (normal dispersion, beta2 > 0)
A0      = 1;            % Peak amplitude [sqrt(W)]
Nz      = 200;          % Number of propagation steps
NT      = 2000;         % Number of time samples
 
LD      = T0^2 / abs(beta2);          % Dispersion length [m]
Z_max   = 4 * LD;                     % Propagate to 4 dispersion lengths
 
%% --- Grids ---
T  = linspace(-6*T0, 6*T0, NT);       % Time grid [s]
Zv = linspace(0, Z_max, Nz);          % Propagation distance grid [m]
 
%% --- Analytical solution ---
% Preallocate intensity matrix  (rows = T, cols = Z)
intensity = zeros(NT, Nz);
fwhm_vec  = zeros(1, Nz);
chirp_mat = zeros(NT, Nz);            % Instantaneous frequency offset
 
for k = 1:Nz
    Z    = Zv(k);
    zeta = 1 + 1j * Z / LD;          % Complex denominator
 
    A = (A0 / sqrt(zeta)) .* exp( -T.^2 ./ (2 * T0^2 * zeta) );
 
    intensity(:, k) = abs(A).^2;
 
    % Instantaneous frequency: delta_omega = -d(phase)/dT
    phase            = angle(A);
    d_phase          = gradient(phase, T);
    chirp_mat(:, k)  = -d_phase;      % [rad/s]
 
    % FWHM of intensity profile (find half-maximum crossings)
    I     = intensity(:, k);
    I_max = max(I);
    idx   = find(I >= I_max/2);
    if numel(idx) >= 2
        fwhm_vec(k) = T(idx(end)) - T(idx(1));
    else
        fwhm_vec(k) = 0;
    end
end
 
% Analytical FWHM: FWHM = 2*sqrt(2*ln2)*T0*sqrt(1+(Z/LD)^2)
fwhm_analytic = 2*sqrt(2*log(2)) * T0 * sqrt(1 + (Zv/LD).^2);
 
%% --- Normalise axes for cleaner plots ---
T_ps   = T  / 1e-12;          % [ps]
Z_norm = Zv / LD;             % [units of LD]
fwhm_ps       = fwhm_vec      / 1e-12;
fwhm_an_ps    = fwhm_analytic / 1e-12;
chirp_THz     = chirp_mat     / (2*pi*1e12);   % [THz]
 
%% =========================================================
%  Figure 1 — 2-D false-colour map of intensity vs (T, Z)
%% =========================================================
fig = figure(1); clf;
fig.Position = [400 400 1000 500];
imagesc(Z_norm, T_ps, intensity);
set(gca,'YDir','normal');
colormap(parula);
cb = colorbar;
cb.Label.String = 'Intensity |A|^2  (W)';
xlabel('Propagation distance  Z / L_D');
ylabel('Retarded time  T  (ps)');
title(sprintf('Gaussian pulse broadening  (\\beta_2 = %.2f ps^2/km,  T_0 = %.1f fs)', ...
              beta2*1e27, T0/1e-15));
 
% Overlay FWHM edges on the map
hold on;
plot(Z_norm, +fwhm_an_ps/2, 'w--', 'LineWidth', 1.5, 'DisplayName','FWHM / 2');
plot(Z_norm, -fwhm_an_ps/2, 'w--', 'LineWidth', 1.5, 'HandleVisibility','off');
% legend('Location','northeast','TextColor','w');
hold off;
pause(1)
print -dpng ../../figures/gaussian_pulse_broadening1
pause(1)
% 
% %% =========================================================
% %  Figure 2 — stacked pulse profiles at selected distances
% %% =========================================================
% fig = figure(2); clf;
% fig.Position = [400 800 600 400];
% 
% Z_plot_norm = [0, 0.5, 1, 2, 3, 4];           % in units of LD
% colors = parula(numel(Z_plot_norm));
% 
% ax = axes;
% hold on;
% for m = 1:numel(Z_plot_norm)
%     Ztarget = Z_plot_norm(m) * LD;
%     [~, ki] = min(abs(Zv - Ztarget));
%     I = intensity(:, ki);
%     % Offset successive curves vertically for clarity
%     offset = (m-1) * 0.25;
%     plot(T_ps, I + offset, 'Color', colors(m,:), 'LineWidth', 1.8, ...
%          'DisplayName', sprintf('Z = %.1f L_D', Z_plot_norm(m)));
% end
% hold off;
% xlabel('Retarded time  T  (ps)');
% ylabel('|A|^2 + offset  (a.u.)');
% title('Pulse profiles at increasing propagation distances');
% legend('Location','northeast');
% grid on; box on;
%  print -dpng ../../figures/gaussian_pulse_broadening3
% 
% %% =========================================================
% %  Figure 3 — FWHM growth vs distance
% %% =========================================================
% fig = figure(2); clf;
% fig.Position = [400 800 1000 300];
% 
% plot(Z_norm, fwhm_ps, 'b.', 'MarkerSize', 6, 'DisplayName','Numerical (half-max)');
% hold on;
% plot(Z_norm, fwhm_an_ps, 'r-', 'LineWidth', 2, ...
%      'DisplayName','Analytic: FWHM_0 \cdot \surd(1+(Z/L_D)^2)');
% xlabel('Propagation distance  Z / L_D');
% ylabel('Intensity FWHM  (ps)');
% title('Full-width at half-maximum vs propagation distance');
% legend('Location','northwest');
% grid on; box on;
% 
% % Mark Z = LD (FWHM = sqrt(2) * FWHM_0)
% FWHM0 = 2*sqrt(2*log(2)) * T0 / 1e-12;
% yline(sqrt(2)*FWHM0, 'k--', '\surd2 \times FWHM_0  at  Z = L_D', ...
%       'LabelHorizontalAlignment','right');
% xline(1, 'k:', 'Z = L_D', 'LabelVerticalAlignment','bottom');
% print -dpng ../../figures/gaussian_pulse_broadening2 
% % %% =========================================================
% % %  Figure 4 — chirp (instantaneous frequency) at select Z
% % %% =========================================================
% % figure(4);
% % hold on;
% % for m = 1:numel(Z_plot_norm)
% %     Ztarget = Z_plot_norm(m) * LD;
% %     [~, ki] = min(abs(Zv - Ztarget));
% %     plot(T_ps, chirp_THz(:, ki), 'Color', colors(m,:), 'LineWidth', 1.8, ...
% %          'DisplayName', sprintf('Z = %.1f L_D', Z_plot_norm(m)));
% % end
% % hold off;
% % xlabel('Retarded time  T  (ps)');
% % ylabel('\delta\omega / 2\pi  (THz)');
% % title('Instantaneous frequency chirp  (normal dispersion, \beta_2 > 0)');
% % legend('Location','northeast');
% % yline(0,'k--');
% % grid on; box on;
% 
% %% --- Print summary to command window ---
% fprintf('\n=== Dispersion parameters ===\n');
% fprintf('  T0      = %.2f ps\n',  T0/1e-12);
% fprintf('  beta2   = %.2e s^2/m\n', beta2);
% fprintf('  LD      = %.4f m\n',   LD);
% fprintf('  FWHM_0  = %.4f ps\n',  FWHM0);
% fprintf('  FWHM at Z=LD : %.4f ps  (= sqrt(2) x FWHM_0 = %.4f ps)\n', ...
%         sqrt(2)*FWHM0, sqrt(2)*FWHM0);
% fprintf('  FWHM at Z=4LD: %.4f ps\n', fwhm_an_ps(end));
% 