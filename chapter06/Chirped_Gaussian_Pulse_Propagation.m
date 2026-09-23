%% Chirped Gaussian Pulse Propagation — Cbeta2 < 0 (compression regime)
% Visualizes T(z)/T0 = sqrt( (1 + C*z/LD)^2 + (z/LD)^2 )
% for several negative-product (C*beta2 < 0) chirp values.

clear; clc; 

%% ── Parameters ────────────────────────────────────────────────────────────
% Normalised propagation axis  ξ = z / L_D
xi   = linspace(0, 1.2, 1000);   % 0 … 2 L_D

% Chirp values with C*beta2 < 0  (here beta2 < 0, so C > 0)
% C > 0 paired with anomalous dispersion (beta2 < 0)
C_vals   = [0, 1, 2, 3];    % unchirped + four compression cases
colors   = ["#555555","#0072BD","#D95319","#77AC30","#7E2F8E"];
lstyle   = {'-','--','--','--','--'};

%% ── Broadening factor ─────────────────────────────────────────────────────
% T(z)/T0  for C*beta2 < 0  → sign in cross-term is negative → use –|C|
bf = @(C, xi) sqrt( (1 - C.*xi).^2 + xi.^2 );   % beta2 < 0, C > 0

%% ── z_min (normalised) for each C ────────────────────────────────────────
z_min_norm = @(C) abs(C) ./ (1 + C.^2);          % |C|/(1+C²) in units of L_D
T_min_norm = @(C) 1 ./ sqrt(1 + C.^2);           % T_min / T0

%% ── Plot 1 : broadening factor vs ξ ─────────────────────────────────────
figure('Name','Pulse Compression (Cβ₂ < 0)', ...
       'Position',[100 100 1280 520],'Color','w');
hold on; grid on; box on;

for k = 1:numel(C_vals)
    C = C_vals(k);
    if C == 0
        y = sqrt(1 + xi.^2);           % unchirped reference
        lbl = 'C = 0 (unchirped)';
    else
        y   = bf(C, xi);
        lbl = sprintf('C = +%d  (Cβ₂ < 0)', C);
    end
    plot(xi, y, lstyle{k}, 'Color', colors(k), ...
         'LineWidth', 2.0, 'DisplayName', lbl);

    % Mark compression minimum for chirped cases
    if C ~= 0
        xm = z_min_norm(C);
        Tm = T_min_norm(C);
        plot(xm, Tm, 'o', 'MarkerSize', 7, 'MarkerFaceColor', colors(k), ...
             'MarkerEdgeColor','k','HandleVisibility','off');
        text(xm + 0.02, Tm + 0.04, ...
             sprintf('z_{min}/L_D = %.2f\nT_{min}/T_0 = %.2f', xm, Tm), ...
             'Color', colors(k), 'FontSize', 10, 'FontWeight','bold');
    end
end

xlabel('Propagation distance  z / L_D',  'FontSize', 14);
ylabel('Pulse-width ratio  T(z) / T_0',  'FontSize', 14);
title('Chirped Gaussian in Dispersive Medium  (C\beta_2 < 0)', ...
      'FontSize', 14, 'FontWeight','bold');
legend('Location','northwest','FontSize',14);
ylim([0 2]);
xlim([0 1.2]);
yline(1,'k:','LineWidth',1.2);          % T/T0 = 1 reference line
print -dpng Chirped_Gaussian_Pulse_Propagation

%% ── Plot 2 : 2-D map  T(z)/T0  over (C, ξ) ────────────────────────────
C_grid  = linspace(0.5, 6, 200);
xi_grid = linspace(0, 2, 300);
[XI, CC] = meshgrid(xi_grid, C_grid);
BF = bf(CC, XI);                        % broadening factor grid

figure('Name','2-D Compression Map','Position',[900 100 700 480],'Color','w');
imagesc(xi_grid, C_grid, BF);
set(gca,'YDir','normal');
colormap(flipud(hot));
cb = colorbar;
cb.Label.String = 'T(z) / T_0';
cb.Label.FontSize = 11;
caxis([0 2.5]);
hold on;

% Locus of compression minima:  z_min = |C|/(1+C²)
xi_min_locus = z_min_norm(C_grid);
plot(xi_min_locus, C_grid, 'c--', 'LineWidth', 2.0, ...
     'DisplayName', 'z_{min}(C) locus');

xlabel('z / L_D', 'FontSize', 13);
ylabel('Chirp parameter  C  (with Cβ₂ < 0)', 'FontSize', 13);
title('Pulse-Width Ratio Map  (Cβ₂ < 0)', 'FontSize', 14, 'FontWeight','bold');
legend('Location','northeast','FontSize',10,'TextColor','w');

%% ── Console summary ──────────────────────────────────────────────────────
fprintf('\n%-6s  %-14s  %-14s\n','C','z_min / L_D','T_min / T_0');
fprintf('%s\n', repmat('-',1,38));
for C = C_vals(2:end)
    fprintf('%-6d  %-14.4f  %-14.4f\n', C, z_min_norm(C), T_min_norm(C));
end