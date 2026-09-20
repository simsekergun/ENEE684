function lens = lens_singlet()
% LENS_SINGLET  Biconvex BK7 singlet lens, f ≈ 100 mm at 532 nm
%
%   lens = LENS_SINGLET()
%
%   Design parameters:
%     Material   : Schott N-BK7 (borosilicate crown)
%     Shape      : Symmetric biconvex
%     Diameter   : 25.4 mm  (1 inch)
%     Focal length: ~100 mm at 532 nm (green)
%     Center thickness: 5.8 mm
%
%   Radius derivation (thin-lens, symmetric biconvex):
%     1/f = (n-1) * 2/R  =>  R = 2*(n-1)*f
%     n_BK7(532nm) ≈ 1.5191  =>  R ≈ 103.8 mm
%
%   Returns a struct with fields:
%     .name        descriptive string
%     .surfaces    cell array of surface structs for raytrace_sequential
%     .D           clear aperture diameter (mm)
%     .f_nominal   nominal focal length (mm)
%     .description short description string

    n_design = sellmeier('BK7', 532);   % ≈ 1.5191
    f_target = 100;                      % mm
    R = 2 * (n_design - 1) * f_target;  % symmetric biconvex radius ≈ 103.8 mm

    D_ap = 25.4;     % clear aperture diameter (mm)
    d_ct = 5.8;      % centre thickness (mm)

    % Surface 1 — front convex (R > 0 : centre to the right)
    s1.R           = R;
    s1.z           = 0;
    s1.D           = D_ap;
    s1.glass_after = 'BK7';

    % Surface 2 — rear convex (R < 0 : centre to the left)
    s2.R           = -R;
    s2.z           = d_ct;
    s2.D           = D_ap;
    s2.glass_after = 'air';

    lens.name             = 'Singlet — BK7 Biconvex';
    lens.surfaces         = {s1, s2};
    lens.D                = D_ap;
    lens.f_nominal        = f_target;
    lens.center_thickness = d_ct;
    lens.glasses          = {'BK7'};
    lens.description      = sprintf(['BK7 symmetric biconvex  |  ' ...
        'R1 = %.1f mm  |  R2 = %.1f mm  |  ' ...
        'd = %.1f mm  |  D = %.1f mm  |  f ≈ %d mm'], R, -R, d_ct, D_ap, f_target);
end
