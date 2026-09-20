function lens = lens_triplet()
% LENS_TRIPLET  Cooke triplet: SK16 – F2 – SK16, f ≈ 100 mm, f/4
%
%   lens = LENS_TRIPLET()
%
%   The Cooke triplet (H.D. Taylor, 1893) uses three separated elements
%   to correct spherical aberration, coma, astigmatism, and field curvature
%   while also achieving partial chromatic correction.
%
%   Configuration:  Biconvex crown | Biconcave flint | Biconvex crown
%
%   Materials (Schott):
%     Element 1 & 3: SK16  (dense crown)   nD = 1.6203, V = 60.41
%     Element 2    : F2    (dense flint)    nD = 1.6200, V = 36.43
%
%   Design data (scaled to f = 100 mm at 532 nm):
%     R1  =  58.52 mm   E1 front
%     R2  = -56.78 mm   E1 rear
%     d1  =   5.5  mm   E1 thickness
%     gap1=   5.0  mm   air gap after E1
%     R3  = -35.30 mm   E2 front (concave)
%     R4  =  35.30 mm   E2 rear  (concave)
%     d2  =   2.5  mm   E2 thickness
%     gap2=   5.5  mm   air gap after E2
%     R5  =  96.54 mm   E3 front
%     R6  = -59.12 mm   E3 rear
%     d3  =   5.0  mm   E3 thickness
%     Clear apertures: E1,E3 = 25 mm; E2 (stop) = 18 mm

    D_outer = 25.0;   % E1 and E3 aperture (mm)
    D_stop  = 18.0;   % E2 / stop aperture (mm)

    R1 =  58.52;  d1   = 5.5;
    R2 = -56.78;  gap1 = 5.0;
    R3 = -35.30;  d2   = 2.5;
    R4 =  35.30;  gap2 = 5.5;
    R5 =  96.54;  d3   = 5.0;
    R6 = -59.12;

    z = 0;

    % Element 1: SK16 biconvex (positive crown)
    s1.R = R1;  s1.z = z;  s1.D = D_outer;  s1.glass_after = 'SK16';  z = z + d1;
    s2.R = R2;  s2.z = z;  s2.D = D_outer;  s2.glass_after = 'air';   z = z + gap1;

    % Element 2: F2 biconcave (negative flint, acts as aperture stop)
    s3.R = R3;  s3.z = z;  s3.D = D_stop;   s3.glass_after = 'F2';    z = z + d2;
    s4.R = R4;  s4.z = z;  s4.D = D_stop;   s4.glass_after = 'air';   z = z + gap2;

    % Element 3: SK16 biconvex (positive crown)
    s5.R = R5;  s5.z = z;  s5.D = D_outer;  s5.glass_after = 'SK16';  z = z + d3;
    s6.R = R6;  s6.z = z;  s6.D = D_outer;  s6.glass_after = 'air';

    lens.name             = 'Cooke Triplet — SK16 · F2 · SK16';
    lens.surfaces         = {s1, s2, s3, s4, s5, s6};
    lens.D                = D_outer;
    lens.f_nominal        = 100;
    lens.center_thickness = d1 + gap1 + d2 + gap2 + d3;
    lens.glasses          = {'SK16', 'F2', 'SK16'};
    lens.description      = ['Cooke triplet SK16-F2-SK16  |  f/4  |  ' ...
        'd = 5.5|5.0|2.5|5.5|5.0 mm  |  f ≈ 100 mm'];
end
