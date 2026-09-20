% =========================================================
%  legend.m
%
%  Computes from the Sellmeier equation:
%    (1) Phase refractive index  n(λ)
%    (2) Group index              ng(λ) = n - λ·dn/dλ
%    (3) Group-velocity dispersion parameter  D(λ)  [ps/(nm·km)]
%    (4) GVD coefficient         β₂(λ)  [ps²/km]
%
%  Reference wavelength: 1550 nm (telecom C-band)
%  Zero-dispersion wavelength (ZDW) of fused silica ≈ 1270 nm
%% =========================================================
clear; close all;
 
%% ---- Mac / renderer fix ----------------------------------
set(0, 'DefaultFigureRenderer', 'painters');
 
%% =========================================================
%  Sellmeier coefficients for fused silica (SiO₂)
%  λ must be in micrometres [μm]
%% =========================================================
B1 = 0.6962;   C1 = 0.0684^2;   % [μm²]
B2 = 0.4079;   C2 = 0.1162^2;
B3 = 0.8975;   C3 = 9.896^2;
 
%% =========================================================
%  Wavelength axis and refractive index
%% =========================================================
lam     = linspace(0.35, 2.2, 50000);   % λ [μm], fine grid for accurate derivatives
lam_nm  = lam * 1e3;                     % λ [nm]  (for plots)
 
% Sellmeier equation
n2 = 1 + B1*lam.^2./(lam.^2 - C1) ...
       + B2*lam.^2./(lam.^2 - C2) ...
       + B3*lam.^2./(lam.^2 - C3);
n  = sqrt(n2);                            % phase refractive index
 
%% =========================================================
%  Numerical derivatives  (central differences via gradient)
%  dlam  in μm  →  derivatives in μm⁻¹, μm⁻²
%% =========================================================
dlam    = lam(2) - lam(1);
dndl    = gradient(n,    dlam);           % dn/dλ   [μm⁻¹]
d2ndl2  = gradient(dndl, dlam);          % d²n/dλ² [μm⁻²]
 
%% =========================================================
%  Derived quantities
%% =========================================================
% Group index
ng      = n - lam .* dndl;
 
% Dispersion parameter D [ps/(nm·km)]
%   D = -(λ/c) · d²n/dλ²
%   Unit conversion: λ[μm]×10⁻⁶ → [m],  d²n/dλ²[μm⁻²]×10¹² → [m⁻²]
%   D[s/m²] × 10⁶ → [ps/(nm·km)]
%   Combined factor: ×10¹² × 10⁻⁶ × 10⁶ = ×10¹²
c_ms    = 2.998e8;                        % speed of light [m/s]
D       = -(lam * 1e12 / c_ms) .* d2ndl2;  % [ps/(nm·km)]
 
% GVD coefficient β₂ [ps²/km]
%   β₂ = -D·λ²/(2πc),  with λ[nm], c[nm/ps] = 2.998×10⁵ nm/ps
c_nmps  = 2.998e5;                        % c in [nm/ps]
beta2   = -D .* (lam_nm).^2 ./ (2*pi*c_nmps);   % [ps²/km]
 
%% =========================================================
%  Key values at 1550 nm
%% =========================================================
lam0_nm = 1550;                           % reference wavelength [nm]
lam0    = lam0_nm * 1e-3;                % [μm]
[~, idx1550] = min(abs(lam_nm - lam0_nm));
 
n_1550    = n(idx1550);
ng_1550   = ng(idx1550);
D_1550    = D(idx1550);
beta2_1550 = beta2(idx1550);
 
% Zero-dispersion wavelength (ZDW)
[~, zdw_idx] = min(abs(D));              % index where D ≈ 0
lam_ZDW_nm   = lam_nm(zdw_idx);
 
%% =========================================================
%  Colours & line styles
%% =========================================================
C_n     = [0.00 0.45 0.74];   % blue
C_ng    = [0.85 0.33 0.10];   % red-orange
C_D     = [0.47 0.67 0.19];   % green
C_b2    = [0.49 0.18 0.56];   % purple
C_ref   = [0.93 0.69 0.13];   % gold  (1550 nm marker)
C_zdw   = [0.40 0.40 0.40];   % grey  (ZDW marker)
 
%% =========================================================
%  FIGURE 1 — Four-panel overview (wide wavelength range)
%% =========================================================
fig1 = figure('Name','Fused Silica Dispersion','NumberTitle','off', ...
              'Position',[40 40 0.9*1100 0.9*850],'Color','w');
 
% ---- Panel 1: Refractive index n(λ) ----------------------
ax1 = subplot(4,1,1);
plot(lam_nm, n, '-', 'Color', C_n, 'LineWidth', 2);
hold on;
xline(lam0_nm,  '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':', 'Color', C_zdw, 'LineWidth', 1.2);
plot(lam0_nm, n_1550, 'o', 'Color', C_ref, 'MarkerSize', 8, ...
     'MarkerFaceColor', C_ref);
text(lam0_nm+20, n_1550+0.003, sprintf('n = %.5f', n_1550), ...
     'Color', C_ref, 'FontSize',11, 'FontWeight','bold');
ylabel('n(\lambda)', 'FontSize', 12);
title('Phase Refractive Index — Fused Silica (SiO_2)', 'FontSize', 11);
xlim([350 2200]);
ylim([min(n)-0.01, max(n)+0.02]);
grid on; box on;
ax1.XTickLabel = [];
legend('n(\lambda)', '\lambda = 1550 nm', 'ZDW', ...
       'Location','northeast','FontSize',11);
 
% ---- Panel 2: Group index ng(λ) --------------------------
ax2 = subplot(4,1,2);
plot(lam_nm, ng, '-', 'Color', C_ng, 'LineWidth', 2);
hold on;
xline(lam0_nm,    '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':',  'Color', C_zdw, 'LineWidth', 1.2);
plot(lam0_nm, ng_1550, 'o', 'Color', C_ref, 'MarkerSize', 8, ...
     'MarkerFaceColor', C_ref);
text(lam0_nm+20, ng_1550+0.003, sprintf('n_g = %.5f', ng_1550), ...
     'Color', C_ref, 'FontSize',11, 'FontWeight','bold');
ylabel('n_g(\lambda)', 'FontSize', 12);
title('Group Index', 'FontSize', 11);
xlim([350 2200]);
grid on; box on;
ax2.XTickLabel = [];
legend('n_g(\lambda)', '\lambda = 1550 nm', 'ZDW', ...
       'Location','northeast','FontSize',11);
 
% ---- Panel 3: Dispersion parameter D(λ) -----------------
ax3 = subplot(4,1,3);
plot(lam_nm, D, '-', 'Color', C_D, 'LineWidth', 2);
hold on;
yline(0,    '-k',  'LineWidth', 1.0);
xline(lam0_nm,    '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':',  'Color', C_zdw, 'LineWidth', 1.2);
plot(lam0_nm, D_1550, 'o', 'Color', C_ref, 'MarkerSize', 8, ...
     'MarkerFaceColor', C_ref);
plot(lam_ZDW_nm, 0, 'p', 'Color', C_zdw, 'MarkerSize', 11, ...
     'MarkerFaceColor', C_zdw);
text(lam0_nm+20,    D_1550+2, sprintf('D = %.1f ps/(nm·km)', D_1550), ...
     'Color', C_ref, 'FontSize',11, 'FontWeight','bold');
text(lam_ZDW_nm+20, 6,        sprintf('ZDW = %.0f nm', lam_ZDW_nm), ...
     'Color', C_zdw, 'FontSize',11, 'FontWeight','bold');
ylabel('D  [ps/(nm{\cdot}km)]', 'FontSize', 12);
title('Dispersion Parameter  D(\lambda) = -(\lambda/c)\cdot d^2n/d\lambda^2', 'FontSize', 11);
xlim([350 2200]);
ylim([-250 80]);
grid on; box on;
ax3.XTickLabel = [];
legend('D(\lambda)', 'D = 0', '\lambda = 1550 nm', 'ZDW', ...
       'Location','southeast','FontSize',11);
 
% ---- Panel 4: GVD coefficient β₂(λ) ---------------------
ax4 = subplot(4,1,4);
plot(lam_nm, beta2, '-', 'Color', C_b2, 'LineWidth', 2);
hold on;
yline(0,    '-k',  'LineWidth', 1.0);
xline(lam0_nm,    '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':',  'Color', C_zdw, 'LineWidth', 1.2);
plot(lam0_nm, beta2_1550, 'o', 'Color', C_ref, 'MarkerSize', 8, ...
     'MarkerFaceColor', C_ref);
text(lam0_nm+20, beta2_1550-5, ...
     sprintf('\\beta_2 = %.1f ps^2/km', beta2_1550), ...
     'Color', C_ref, 'FontSize',11, 'FontWeight','bold');
xlabel('Wavelength  \lambda  [nm]', 'FontSize', 12);
ylabel('\beta_2  [ps^2/km]', 'FontSize', 12);
title('GVD Coefficient  \beta_2(\lambda) = -D\lambda^2/(2\pic)       (anomalous: \beta_2 < 0)', ...
      'FontSize', 11);
xlim([350 2200]);
ylim([-250 150]);
grid on; box on;
legend('\beta_2(\lambda)', '\beta_2 = 0', '\lambda = 1550 nm', 'ZDW', ...
       'Location','southeast','FontSize',11);
 
sgtitle('Fused Silica (SiO_2) — Sellmeier Dispersion Model', ...
        'FontSize', 13, 'FontWeight','bold');
print -dpng ../../figures/Fused_Silica_Dispersion 
pause(1)
%% =========================================================
%  FIGURE 2 — Zoom around 1550 nm (telecom C-band window)
%% =========================================================
lam_win = [1200 1900];                   % display window [nm]
mask = lam_nm >= lam_win(1) & lam_nm <= lam_win(2);
 
fig2 = figure('Name','Fused Silica — Near 1550 nm','NumberTitle','off', ...
              'Position',[80 80 1000 750],'Color','w');
 
% ---- n and ng together ------------------------------------
ax_a = subplot(3,1,1);
yyaxis left
plot(lam_nm(mask), n(mask),  '-', 'Color', C_n,  'LineWidth', 2.2);
ylabel('Phase index  n(\lambda)', 'FontSize', 11, 'Color', C_n);
ax_a.YColor = C_n;
yyaxis right
plot(lam_nm(mask), ng(mask), '-', 'Color', C_ng, 'LineWidth', 2.2);
ylabel('Group index  n_g(\lambda)', 'FontSize', 11, 'Color', C_ng);
ax_a.YColor = C_ng;
hold on;
yyaxis left
xline(lam0_nm,    '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':',  'Color', C_zdw, 'LineWidth', 1.2);
xlim(lam_win);  grid on;  box on;
ax_a.XTickLabel = [];
title('Refractive Index and Group Index near 1550 nm', 'FontSize', 11);
legend('n(\lambda)', 'n_g(\lambda)', '\lambda_{ref}=1550 nm', 'ZDW', ...
       'Location','East','FontSize',11);
 
% ---- D(λ) zoom -------------------------------------------
ax_b = subplot(3,1,2);
plot(lam_nm(mask), D(mask), '-', 'Color', C_D, 'LineWidth', 2.2);
hold on;
yline(0, '-k', 'LineWidth', 1.0);
xline(lam0_nm,    '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':',  'Color', C_zdw, 'LineWidth', 1.2);
plot(lam0_nm, D_1550, 'o', 'MarkerSize', 9, ...
     'Color', C_ref, 'MarkerFaceColor', C_ref);
plot(lam_ZDW_nm, 0, 'p', 'MarkerSize', 12, ...
     'Color', C_zdw, 'MarkerFaceColor', C_zdw);
 
% Shade anomalous dispersion region (D > 0 → β₂ < 0)
idx_anom = lam_nm(mask) >= lam_ZDW_nm;
lam_anom = lam_nm(mask);   D_anom = D(mask);
fill([lam_anom(idx_anom), fliplr(lam_anom(idx_anom))], ...
     [D_anom(idx_anom),   zeros(1,sum(idx_anom))], ...
     C_D, 'FaceAlpha', 0.10, 'EdgeColor','none');
 
% Labels
text(1420, D_1550+2.5, sprintf('D(1550) = %.1f ps/(nm·km)', D_1550), ...
     'FontSize', 12, 'Color', C_ref, 'FontWeight','bold');
text(lam_ZDW_nm-180, -15, sprintf('ZDW = %.0f nm', lam_ZDW_nm), ...
     'FontSize', 12, 'Color', C_zdw, 'FontWeight','bold');
 
% Region labels
text(1260, -17, 'Normal  (D<0)', 'FontSize', 12, 'Color',[0.5 0.5 0.5], ...
     'HorizontalAlignment','center');
text(1700, 15, 'Anomalous  (D>0)', 'FontSize', 12, 'Color', C_D*0.7, ...
     'HorizontalAlignment','center');
 
ylabel('D  [ps/(nm{\cdot}km)]', 'FontSize', 11);
xlim(lam_win);  ylim([-60 40]);
grid on;  box on;
ax_b.XTickLabel = [];
title('Dispersion Parameter', 'FontSize', 11);
 
% ---- β₂(λ) zoom ------------------------------------------
ax_c = subplot(3,1,3);
plot(lam_nm(mask), beta2(mask), '-', 'Color', C_b2, 'LineWidth', 2.2);
hold on;
yline(0, '-k', 'LineWidth', 1.0);
xline(lam0_nm,    '--', 'Color', C_ref, 'LineWidth', 1.5);
xline(lam_ZDW_nm, ':',  'Color', C_zdw, 'LineWidth', 1.2);
plot(lam0_nm, beta2_1550, 'o', 'MarkerSize', 9, ...
     'Color', C_ref, 'MarkerFaceColor', C_ref);
 
% Shade anomalous (β₂ < 0)
b2_anom = beta2(mask);
fill([lam_anom(idx_anom), fliplr(lam_anom(idx_anom))], ...
     [b2_anom(idx_anom),  zeros(1,sum(idx_anom))], ...
     C_b2, 'FaceAlpha', 0.10, 'EdgeColor','none');
 
text(1560, beta2_1550+6, sprintf('\\beta_2(1550) = %.1f ps^2/km', beta2_1550), ...
     'FontSize', 10, 'Color', C_ref, 'FontWeight','bold');
 
xlabel('Wavelength  \lambda  [nm]', 'FontSize', 11);
ylabel('\beta_2  [ps^2/km]', 'FontSize', 11);
xlim(lam_win);  ylim([-80 10]);
grid on;  box on;
title('GVD Coefficient  \beta_2  (anomalous: shaded)', 'FontSize', 11);
 
sgtitle('Fused Silica (SiO_2) near 1550 nm — Telecom C-band', ...
        'FontSize', 13, 'FontWeight','bold');
print -dpng ../../figures/Fused_Silica_Dispersion2 
%% =========================================================
%  Print summary table to Command Window
%% =========================================================
fprintf('\n========================================================\n');
fprintf('  Fused Silica (SiO2) — Key Values at %.0f nm\n', lam0_nm);
fprintf('========================================================\n');
fprintf('  Phase index          n    = %.6f\n',    n_1550);
fprintf('  Group index          ng   = %.6f\n',    ng_1550);
fprintf('  Dispersion param.    D    = %+.3f  ps/(nm·km)\n', D_1550);
fprintf('  GVD coefficient      β₂   = %+.3f  ps²/km\n',    beta2_1550);
fprintf('  Zero-disp. wavelen.  ZDW  = %.1f  nm\n',          lam_ZDW_nm);
fprintf('  Dispersion regime         = %s\n', ...
        ternary_str(D_1550 > 0, 'Anomalous (D>0, β₂<0)', 'Normal (D<0, β₂>0)'));
fprintf('========================================================\n\n');
 
%% =========================================================
function s = ternary_str(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end
 