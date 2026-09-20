function [history, exit_ray] = raytrace_sequential(surfaces, ray_in, lambda_nm)
% RAYTRACE_SEQUENTIAL  Trace a 3-D ray through a sequential surface list
%
%   [history, exit_ray] = RAYTRACE_SEQUENTIAL(surfaces, ray_in, lambda_nm)
%
%   Inputs:
%     surfaces  - 1×N cell array of surface structs, each with fields:
%                   .R          radius of curvature (mm); +Inf for flat
%                   .z          vertex z-position (mm)
%                   .D          clear-aperture diameter (mm)
%                   .glass_after  glass name string after this surface
%     ray_in    - struct with:
%                   .pos  [3×1]  start position [x;y;z] (mm)
%                   .dir  [3×1]  unit direction vector [dx;dy;dz]
%     lambda_nm - wavelength in nm (scalar)
%
%   Outputs:
%     history   - struct array (1 entry per surface hit) with .pos, .dir
%     exit_ray  - struct with .pos, .dir, .vignetted (logical)
%
%   Sign convention (standard optical):
%     - Optical axis along +z
%     - R > 0 : center of curvature to right of vertex
%     - R < 0 : center of curvature to left  of vertex
%
%   Calls: sellmeier.m

    pos = ray_in.pos(:);
    dir = ray_in.dir(:);
    dir = dir / norm(dir);

    history   = struct('pos', {}, 'dir', {});
    vignetted = false;

    glass_before = 'air';

    for k = 1:numel(surfaces)
        s = surfaces{k};

        % Refractive indices on both sides
        n1 = get_n(glass_before, lambda_nm);
        n2 = get_n(s.glass_after, lambda_nm);

        % Intersect ray with surface
        [t, Q, N_out] = surf_intersect(pos, dir, s);

        if isnan(t) || t < 1e-9
            vignetted = true;
            break;
        end

        % Aperture check
        r_hit = sqrt(Q(1)^2 + Q(2)^2);
        if r_hit > s.D / 2
            vignetted = true;
            break;
        end

        % Refract
        [dir, tir] = snell_refract(dir, N_out, n1, n2);
        if tir
            vignetted = true;
            break;
        end

        pos = Q;
        glass_before = s.glass_after;
        history(end+1).pos = pos;   %#ok<AGROW>
        history(end).dir   = dir;
    end

    exit_ray.pos       = pos;
    exit_ray.dir       = dir;
    exit_ray.vignetted = vignetted;
end

% -----------------------------------------------------------------------
function n = get_n(glass, lambda_nm)
    if strcmpi(glass,'air') || strcmpi(glass,'vacuum') || isempty(glass)
        n = 1.0;
    else
        n = sellmeier(glass, lambda_nm);
    end
end

% -----------------------------------------------------------------------
function [t, Q, N_out] = surf_intersect(pos, dir, s)
% Find the first valid intersection of the ray with surface s.
% N_out is the outward normal (pointing toward the incoming medium).

    z_v = s.z;
    R   = s.R;

    if isinf(R)
        %--- Flat surface ---
        if abs(dir(3)) < 1e-14
            t = NaN; Q = []; N_out = [];
            return;
        end
        t = (z_v - pos(3)) / dir(3);
        Q = pos + t * dir;
        N_out = [0; 0; -sign(dir(3))];   % points toward incoming ray

    else
        %--- Spherical surface ---
        C    = [0; 0; z_v + R];           % centre of curvature
        Prel = pos - C;

        b    = dot(Prel, dir);
        c    = dot(Prel, Prel) - R^2;
        disc = b^2 - c;

        if disc < 0
            t = NaN; Q = []; N_out = [];
            return;
        end

        sq = sqrt(disc);
        t1 = -b - sq;
        t2 = -b + sq;

        % Pick smallest positive t (first hit from the left)
        if t1 > 1e-9
            t = t1;
        elseif t2 > 1e-9
            t = t2;
        else
            t = NaN; Q = []; N_out = [];
            return;
        end

        Q = pos + t * dir;
        % Outward normal: from centre → Q, normalised, pointing toward incoming ray
        N_out = (Q - C) / abs(R);
        if dot(N_out, dir) > 0
            N_out = -N_out;    % flip so it faces the incoming ray
        end
    end
end

% -----------------------------------------------------------------------
function [dir_out, tir] = snell_refract(dir_in, N_out, n1, n2)
% Vector form of Snell's law.
% N_out points into the medium of n1 (toward incoming ray).

    tir = false;
    mu  = n1 / n2;

    cos_i = -dot(dir_in, N_out);
    if cos_i < 0
        % Ray hits back face — flip normal
        N_out = -N_out;
        cos_i = -dot(dir_in, N_out);
    end

    sin2_r = mu^2 * (1 - cos_i^2);
    if sin2_r > 1
        tir = true;
        dir_out = dir_in;
        return;
    end

    cos_r   = sqrt(1 - sin2_r);
    dir_out = mu * dir_in + (mu * cos_i - cos_r) * N_out;
    dir_out = dir_out / norm(dir_out);
end
