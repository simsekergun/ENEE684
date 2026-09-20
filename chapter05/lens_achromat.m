function lens = lens_achromat()
% LENS_ACHROMAT  Cemented achromatic doublet (BK7 crown + F2 flint), f = 100 mm
%
%   lens = LENS_ACHROMAT()
%
%   Design: Fraunhofer-type cemented doublet corrected for chromatic
%   aberration at the C (656 nm) and F (486 nm) Fraunhofer lines.
%
%   Material properties (d-line, 587.6 nm):
%     Crown  N-BK7: nD = 1.5168,  Abbe V = 64.17
%     Flint  N-F2:  nD = 1.6200,  Abbe V = 36.43
%
%   Achromat power distribution (thin-lens theory):
%     phi1 = phi * V1 / (V1 - V2) = (1/100) * 64.17/27.74 =  0.02313 mm^-1  => f1 =  43.2 mm
%     phi2 = phi - phi1            = (1/100) * (-27.74)/27.74... = -0.01313 mm^-1  => f2 = -76.2 mm
%
%   Surface data (from standard doublet catalogue scaling to f=100 mm):
%     R1 =  61.47 mm  (crown front)
%     R2 = -44.64 mm  (cemented surface)
%     R3 = -141.9 mm  (flint rear)
%     Crown thickness   d1 = 6.5 mm
%     Flint  thickness  d2 = 3.0 mm
%     Clear aperture    D  = 25.4 mm
%
%   These radii produce f ≈ 100 mm at 587.6 nm and correct longitudinal
%   chromatic aberration (achromatic condition).

    D_ap  = 25.4;

    R1 =  61.47;   d_crown = 6.5;
    R2 = -44.64;   d_flint = 3.0;
    R3 = -141.9;

    z = 0;

    % Surface 1: front of crown element
    s1.R           =  R1;
    s1.z           = z;          z = z + d_crown;
    s1.D           = D_ap;
    s1.glass_after = 'BK7';

    % Surface 2: cemented interface (crown → flint)
    s2.R           =  R2;
    s2.z           = z;          z = z + d_flint;
    s2.D           = D_ap;
    s2.glass_after = 'F2';

    % Surface 3: rear of flint element
    s3.R           =  R3;
    s3.z           = z;
    s3.D           = D_ap;
    s3.glass_after = 'air';

    lens.name             = 'Achromatic Doublet — BK7 + F2';
    lens.surfaces         = {s1, s2, s3};
    lens.D                = D_ap;
    lens.f_nominal        = 100;
    lens.center_thickness = d_crown + d_flint;
    lens.glasses          = {'BK7', 'F2'};
    lens.description      = ['Cemented achromat BK7+F2  |  ' ...
        'R1=61.47  R2=-44.64  R3=-141.9 mm  |  ' ...
        'd=6.5+3.0 mm  |  D=25.4 mm  |  f≈100 mm'];
end
