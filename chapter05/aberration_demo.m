%% aberration_demo.m
% =========================================================================
%  LENS ABERRATION DEMONSTRATION
%  Shows chromatic and monochromatic aberrations for three lens types:
%    1. Singlet      — BK7 biconvex,  f ≈ 100 mm
%    2. Achromat     — BK7+F2 doublet, f ≈ 100 mm (cemented)
%    3. Cooke Triplet— SK16-F2-SK16,  f ≈ 100 mm, f/4
%
%  Wavelengths: Blue 400 nm | Green 532 nm | Red 633 nm
%
%  Required files (same directory):
%    sellmeier.m, raytrace_sequential.m,
%    lens_singlet.m, lens_achromat.m, lens_triplet.m
%
%  Usage:  >> aberration_demo
% =========================================================================

clear; close all; clc;

%% ── Wavelength & colour setup ──────────────────────────────────────────
LAMBDA  = [400,  532,  633];            % nm
COLORS  = {[0.15 0.15 1], ...          % blue
           [0.0  0.72 0.0], ...        % green
           [0.9  0.0  0.0]};           % red
LABELS  = {'Blue (400 nm)', 'Green (532 nm)', 'Red (633 nm)'};

% Rays-per-analysis
N_MERIDIONAL = 9;     % meridional rays for cross-section plot
N_SPOT_RINGS = 6;     % rings in spot diagram
N_SPOT_PTS   = 24;    % points per ring in spot diagram
N_SA_RAYS    = 60;    % rays for spherical-aberration curve

%% ── Define lenses ───────────────────────────────────────────────────────
lenses = {lens_singlet(), lens_achromat(), lens_triplet()};
L_names = {'Singlet', 'Achromat', 'Cooke Triplet'};

%% ── Main loop over lens types ───────────────────────────────────────────
all_results = cell(3,1);

for L = 1:numel(lenses)
    lens = lenses{L};
    fprintf('\n════════════════════════════════════════\n');
    fprintf(' %s\n', lens.name);
    fprintf(' %s\n', lens.description);
    fprintf('════════════════════════════════════════\n');

    %-- Find paraxial focal length (each wavelength) ----------------------
    f_paraxial = zeros(1,3);
    for w = 1:3
        f_paraxial(w) = find_paraxial_focus(lens, LAMBDA(w));
        fprintf('  f(%d nm) = %.3f mm\n', LAMBDA(w), f_paraxial(w));
    end

    %-- Spherical aberration curve (green) --------------------------------
    [h_arr, dz_sa] = spherical_aberration(lens, LAMBDA(2), N_SA_RAYS);
    fprintf('  LSA (green): marginal - paraxial = %.3f mm\n', dz_sa(end)-dz_sa(1));

    %-- Spot diagrams at green paraxial focus ----------------------------
    z_img = f_paraxial(2);
    spots = struct();
    for w = 1:3
        [sx, sy] = compute_spot(lens, LAMBDA(w), z_img, N_SPOT_RINGS, N_SPOT_PTS);
        spots(w).x = sx;  spots(w).y = sy;
    end

    %-- Ray-fan (transverse aberration vs normalised pupil coord) --------
    [p_fan, ta_fan] = ray_fan(lens, LAMBDA, z_img, 40);

    %-- Store results for comparison ------------------------------------
    all_results{L}.lens        = lens;
    all_results{L}.f_paraxial  = f_paraxial;
    all_results{L}.h_sa        = h_arr;
    all_results{L}.dz_sa       = dz_sa;
    all_results{L}.spots       = spots;
    all_results{L}.p_fan       = p_fan;
    all_results{L}.ta_fan      = ta_fan;

    %-- Plot this lens ---------------------------------------------------
    figure('Name', sprintf('Lens %d: %s', L, lens.name), ...
           'NumberTitle','off', 'Position', [80*(L-1)+20, 50, 1200, 1200]);

    plot_lens_analysis(lens, LAMBDA, COLORS, LABELS, ...
        h_arr, dz_sa, spots, p_fan, ta_fan, f_paraxial, N_MERIDIONAL);

    sgtitle(lens.name, 'FontSize', 14, 'FontWeight','bold');
end

%% ── Comparison summary figure ──────────────────────────────────────────
figure('Name','Comparison: Chromatic & Spherical Aberration', ...
       'NumberTitle','off','Position',[150 100 1100 500]);


subplot(1,2,1);
hold on; grid on;
for L = 1:3
    f_p = all_results{L}.f_paraxial;
    for w = 1:3
        scatter(L + 0.2*(w-2), f_p(w), 80, COLORS{w}, 'filled');
    end
    % Range bar
    plot([L L], [min(f_p) max(f_p)], 'k-', 'LineWidth',2);
end
set(gca,'XTick',1:3,'XTickLabel', L_names, 'FontSize',11);
ylabel('Paraxial focal length (mm)','FontSize',11);
title('Chromatic Aberration: Focal Length vs Wavelength','FontSize',12);
legend(LABELS{:},'Location','northwest','FontSize',9);
xlim([0.5 3.5]); box on;

subplot(1,2,2);
hold on; grid on; box on;
for L = 1:3
    h   = all_results{L}.h_sa;
    dz  = all_results{L}.dz_sa;
    plot(dz - dz(1), h / (lenses{L}.D/2), '-', ...
        'LineWidth', 2, 'DisplayName', L_names{L});
end
xlabel('Longitudinal Spherical Aberration (mm)','FontSize',11);
ylabel('Normalised Pupil Height  y/y_{max}','FontSize',11);
title('Spherical Aberration Curve — 532 nm','FontSize',12);
legend('Location','northwest','FontSize',10);

sgtitle('Lens Performance Comparison','FontSize',13,'FontWeight','bold');

fprintf('\nDone. Figures generated.\n');

%% =========================================================================
%%  LOCAL FUNCTIONS
%% =========================================================================

% ──────────────────────────────────────────────────────────────────────────
function f = find_paraxial_focus(lens, lambda_nm)
% Trace a paraxial ray (1% of aperture height) and find z-intercept.
    h_par = 0.01 * lens.D / 2;
    ray.pos = [0; h_par; -1e4];
    ray.dir = [0; 0; 1];
    [~, ex] = raytrace_sequential(lens.surfaces, ray, lambda_nm);
    if ex.vignetted
        f = NaN; return;
    end
    % y-intercept: y0 + t*dy = 0  =>  t = -y0/dy
    t = -ex.pos(2) / ex.dir(2);
    f = ex.pos(3) + t * ex.dir(3);
end

% ──────────────────────────────────────────────────────────────────────────
function [h_arr, dz_arr] = spherical_aberration(lens, lambda_nm, N)
% Longitudinal spherical aberration: focus position vs ray height.
    h_max  = 0.98 * lens.D / 2;
    h_arr  = linspace(0.01 * lens.D/2, h_max, N);
    dz_arr = zeros(size(h_arr));
    for k = 1:N
        ray.pos = [0; h_arr(k); -1e4];
        ray.dir = [0; 0; 1];
        [~, ex] = raytrace_sequential(lens.surfaces, ray, lambda_nm);
        if ex.vignetted || abs(ex.dir(2)) < 1e-10
            dz_arr(k) = NaN;
        else
            t = -ex.pos(2) / ex.dir(2);
            dz_arr(k) = ex.pos(3) + t * ex.dir(3);
        end
    end
end

% ──────────────────────────────────────────────────────────────────────────
function [sx, sy] = compute_spot(lens, lambda_nm, z_img, N_rings, N_pts)
% Compute spot diagram at z_img plane by tracing ring bundles.
    sx = []; sy = [];
    h_max = 0.97 * lens.D / 2;
    % Also add a central ray
    ray0.pos = [0;0;-1e4]; ray0.dir = [0;0;1];
    [~,ex0] = raytrace_sequential(lens.surfaces, ray0, lambda_nm);
    if ~ex0.vignetted
        t0 = (z_img - ex0.pos(3)) / ex0.dir(3);
        sx(end+1) = ex0.pos(1) + t0*ex0.dir(1);
        sy(end+1) = ex0.pos(2) + t0*ex0.dir(2);
    end
    for r = 1:N_rings
        rho = h_max * r / N_rings;
        for j = 0:N_pts-1
            phi  = 2*pi*j/N_pts;
            yin  = rho * sin(phi);
            xin  = rho * cos(phi);
            ray.pos = [xin; yin; -1e4];
            ray.dir = [0; 0; 1];
            [~, ex] = raytrace_sequential(lens.surfaces, ray, lambda_nm);
            if ex.vignetted, continue; end
            t = (z_img - ex.pos(3)) / ex.dir(3);
            sx(end+1) = ex.pos(1) + t*ex.dir(1); %#ok<AGROW>
            sy(end+1) = ex.pos(2) + t*ex.dir(2); %#ok<AGROW>
        end
    end
end

% ──────────────────────────────────────────────────────────────────────────
function [p_out, ta_out] = ray_fan(lens, lambdas, z_img, N)
% Transverse ray fan (meridional): TA vs normalised pupil position.
    h_max   = 0.97 * lens.D / 2;
    p_arr   = linspace(-1, 1, N);
    ta_out  = zeros(numel(lambdas), N);
    for w = 1:numel(lambdas)
        for k = 1:N
            yin = p_arr(k) * h_max;
            ray.pos = [0; yin; -1e4];
            ray.dir = [0; 0; 1];
            [~, ex] = raytrace_sequential(lens.surfaces, ray, lambdas(w));
            if ex.vignetted || abs(ex.dir(3)) < 1e-12
                ta_out(w,k) = NaN; continue;
            end
            t = (z_img - ex.pos(3)) / ex.dir(3);
            ta_out(w,k) = ex.pos(2) + t*ex.dir(2);  % y at image plane
        end
    end
    % Centre (reference image height from paraxial green ray)
    ref = ta_out(2, round(N/2) + 1);   % near-axis green
    for w = 1:numel(lambdas)
        ta_out(w,:) = ta_out(w,:) - ref;  % transverse aberration
    end
    p_out = p_arr;
end

% ──────────────────────────────────────────────────────────────────────────
function plot_lens_analysis(lens, LAMBDA, COLORS, LABELS, ...
    h_arr, dz_arr, spots, p_fan, ta_fan, f_paraxial, N_MER)
% Arrange 6 subplots for one lens.

    %──── Subplot 1: Cross-section + ray traces ────────────────────────
    subplot(3,1,1);
    draw_cross_section(lens, LAMBDA, COLORS, N_MER);

    %──── Subplot 2: Spot diagram ──────────────────────────────────────
    subplot(3,2,3);
    hold on; axis equal; grid on; box on;
    for w = 1:3
        scatter(spots(w).x*1e3, spots(w).y*1e3, 10, COLORS{w}, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end
    xlabel('x (\mum)','FontSize',10);
    ylabel('y (\mum)','FontSize',10);
    title('Spot Diagram at Green Paraxial Focus','FontSize',10);
    legend(LABELS{:},'FontSize',8,'Location','best');
    ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';

    %──── Subplot 3: Longitudinal Spherical Aberration ────────────────
    subplot(3,2,4);
    valid = ~isnan(dz_arr);
    dz_ref = dz_arr(1);    % paraxial focus (smallest ray height)
    plot(dz_arr(valid) - dz_ref, h_arr(valid), '-', ...
        'Color', COLORS{2}, 'LineWidth', 2.5);
    hold on; xline(0,'k--','LineWidth',1);
    xlabel('Focus shift \DeltaZ (mm)','FontSize',10);
    ylabel('Ray height h (mm)','FontSize',10);
    title('Longitudinal Spherical Aberration  [532 nm]','FontSize',10);
    grid on; box on;

    %──── Subplot 4: Chromatic focal shift ────────────────────────────
    subplot(3,2,5);
    hold on; grid on; box on;
    f_ref = f_paraxial(2);    % green reference
    for w = 1:3
        scatter(LAMBDA(w), f_paraxial(w) - f_ref, 80, COLORS{w}, ...
            'filled', 'DisplayName', LABELS{w});
        text(LAMBDA(w)+8, f_paraxial(w)-f_ref, ...
            sprintf('%.3f mm', f_paraxial(w)-f_ref), 'Color', COLORS{w}, 'FontSize',9);
    end
    yline(0,'k--');
    xlim([360 680]); xlabel('Wavelength (nm)','FontSize',10);
    ylabel('\DeltaF  (mm)','FontSize',10);
    title('Chromatic Focal Shift  (relative to 532 nm)','FontSize',10);
    legend('Location','best','FontSize',8);

    %──── Subplot 5: Ray fan (transverse aberration) ──────────────────
    subplot(3,2,6);
    hold on; grid on; box on;
    for w = 1:3
        ta_mm = ta_fan(w,:) * 1e3;   % convert to µm... actually to µm:
        % ta_fan is in mm, multiply by 1000 for µm
        plot(p_fan, ta_fan(w,:)*1e3, '-', 'Color', COLORS{w}, ...
            'LineWidth', 2.0, 'DisplayName', LABELS{w});
    end
    yline(0,'k--','LineWidth',1);
    xlabel('Normalised pupil coordinate','FontSize',10);
    ylabel('Transverse aberration (\mum)','FontSize',10);
    title('Ray Fan Plot  (on-axis)','FontSize',10);
    legend('Location','best','FontSize',8);


    figure('Position',[300 300 1200 300]);   
    draw_cross_section(lens, LAMBDA, COLORS, N_MER);    

end

% ──────────────────────────────────────────────────────────────────────────
function draw_cross_section(lens, LAMBDA, COLORS, N_MER)
% Draw the 2-D meridional cross-section and trace coloured rays.

    hold on; axis equal; grid on; box on;
    surfaces = lens.surfaces;
    D_ap     = lens.D;

    %-- Background shading for glass elements --
    glass_fill = containers.Map({'BK7','N-BK7','F2','N-F2','SK16','SF11','BAF10','LAFN7'}, ...
        {[0.75 0.88 1.0], [0.75 0.88 1.0], [1.0 0.82 0.75], [1.0 0.82 0.75], ...
         [0.78 0.92 0.78],[0.95 0.80 0.95],[0.95 0.95 0.75],[0.85 0.75 0.95]});

    % Draw filled glass elements (connect consecutive surfaces of same glass)
    for k = 1:numel(surfaces)-1
        g_k  = surfaces{k}.glass_after;
        g_k1 = surfaces{k+1}.glass_after;
        if ~(strcmpi(g_k,'air') || strcmpi(g_k,'vacuum'))
            % same element: fill between surface k and k+1
            Dmin = min(surfaces{k}.D, surfaces{k+1}.D) * 0.99;
            y_pts = linspace(-Dmin/2, Dmin/2, 120);
            z_left  = surf_z_profile(surfaces{k},   y_pts);
            z_right = surf_z_profile(surfaces{k+1}, y_pts);
            col = [0.80 0.90 1.00];   % default light blue
            if isKey(glass_fill, g_k)
                col = glass_fill(g_k);
            end
            patch([z_left, fliplr(z_right)], [y_pts, fliplr(y_pts)], col, ...
                'EdgeColor','none','FaceAlpha',0.6);
        end
    end

    %-- Draw surface outlines --
    for k = 1:numel(surfaces)
        s    = surfaces{k};
        y_s  = linspace(-s.D/2, s.D/2, 200);
        z_s  = surf_z_profile(s, y_s);
        plot(z_s, y_s, 'k-', 'LineWidth', 1.5);
        % Rim lines
        for sgn = [-1, 1]
            z_r = surf_z_profile(s, sgn * s.D/2);
            if k < numel(surfaces)
                s2  = surfaces{k+1};
                g_k = s.glass_after;
                if strcmpi(g_k,'air') || strcmpi(g_k,'vacuum')
                    continue;
                end
                z_r2 = surf_z_profile(s2, sgn * min(s.D,s2.D)/2);
                plot([z_r z_r2], [sgn*s.D/2, sgn*min(s.D,s2.D)/2], 'k-','LineWidth',1.5);
            end
        end
    end

    %-- Optical axis --
    z_start = surfaces{1}.z   - 25;
    z_end   = surfaces{end}.z + 25;
    plot([z_start z_end], [0 0], 'k--', 'LineWidth', 0.8);

    %-- Trace meridional rays for each wavelength -------------------------
    y_rays = linspace(-0.90*D_ap/2, 0.90*D_ap/2, N_MER);
    z_obj  = surfaces{1}.z - 20;    % object plane (to the left of lens)

    for w = 1:3
        for yi = 1:numel(y_rays)
            ray.pos = [0; y_rays(yi); z_obj];
            ray.dir = [0; 0; 1];
            [hist, ex] = raytrace_sequential(surfaces, ray, LAMBDA(w));
            if ex.vignetted, continue; end

            % Collect z,y at each surface
            z_pts = [z_obj;  arrayfun(@(h) h.pos(3), hist(:))];
            y_pts = [y_rays(yi); arrayfun(@(h) h.pos(2), hist(:))];

            % Propagate to paraxial focus plane
            f_z = find_paraxial_focus(lens, LAMBDA(2));
            t_ext = (f_z - ex.pos(3)) / ex.dir(3);
            z_end_r = ex.pos(3) + t_ext * ex.dir(3);
            y_end_r = ex.pos(2) + t_ext * ex.dir(2);
            z_pts = [z_pts; z_end_r];
            y_pts = [y_pts; y_end_r];

            alpha_val = 0.55 + 0.1*(yi == round(numel(y_rays)/2));
            plot(z_pts, y_pts, '-', 'Color', [COLORS{w}, alpha_val], 'LineWidth', 1.0);
        end
    end

    % Focal length annotation (green)
    f_g = find_paraxial_focus(lens, LAMBDA(2));
    xline(f_g, '--', 'Color', [0 0.6 0], 'LineWidth', 1.2);
    text(f_g+1, D_ap/2*0.85, sprintf('f_{532}=%.1f mm',f_g), ...
        'Color',[0 0.5 0],'FontSize',8);

    xlabel('z (mm)','FontSize',10);
    ylabel('y (mm)','FontSize',10);
    title('Meridional Ray Trace (on-axis)','FontSize',10);
    legend([{'\color{blue}400 nm','\color[rgb]{0 0.72 0}532 nm','\color{red}633 nm'}], ...
        'FontSize', 8, 'Location','northwest');

    % Tidy axis
    z_all = cellfun(@(s) s.z, surfaces);
    xlim([min(z_all)-30, f_g+30]);
    ylim([-D_ap/2*1.3, D_ap/2*1.3]);    
end

% ──────────────────────────────────────────────────────────────────────────
function z = surf_z_profile(s, y)
% Compute z-coordinate on surface s at height(s) y.
%   z = z_v + R - sign(R)*sqrt(R^2 - y^2)   for spherical
%   z = z_v                                   for flat

    if isinf(s.R)
        z = s.z * ones(size(y));
    else
        R = s.R;
        z = s.z + R - sign(R) * sqrt(max(0, R^2 - y.^2));
    end
end
