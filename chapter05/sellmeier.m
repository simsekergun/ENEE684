function n = sellmeier(glass, lambda_nm)
% SELLMEIER  Refractive index via Sellmeier dispersion equation
%
%   n = SELLMEIER(glass, lambda_nm)
%
%   Inputs:
%     glass     - Glass type string: 'BK7','N-BK7','F2','N-F2','SF11',
%                 'SK16','BAF10','LAFN7','air'
%     lambda_nm - Wavelength(s) in nanometers (scalar or array)
%
%   Output:
%     n         - Refractive index (same size as lambda_nm)
%
%   Sellmeier equation:
%     n^2 = 1 + B1*L2/(L2-C1) + B2*L2/(L2-C2) + B3*L2/(L2-C3)
%   where L2 = (lambda in micrometers)^2, C in um^2
%
%   Glass data from Schott catalog (2022)
%
%   Examples:
%     n_BK7_green = sellmeier('BK7', 532)     % -> 1.5191
%     n_F2_blue   = sellmeier('F2', 400)      % -> 1.6600
%     n_arr       = sellmeier('SK16', [400 532 633])

lambda_um2 = (lambda_nm / 1000).^2;   % wavelength^2 in um^2

switch upper(strtrim(glass))

    case {'BK7', 'N-BK7'}
        % Borosilicate crown — most common optical glass, nD=1.5168, V=64.17
        B = [1.03961212,  0.231792344, 1.01046945];
        C = [6.00069867e-3, 2.00179144e-2, 1.03560653e2];

    case {'F2', 'N-F2'}
        % Dense flint — nD=1.6200, V=36.43; used in achromats
        B = [1.34533359,  0.209073176, 0.937357162];
        C = [9.97743871e-3, 4.70450767e-2, 1.11886764e2];

    case {'SF11', 'N-SF11'}
        % Extra-dense flint — nD=1.7847, V=25.68; high dispersion
        B = [1.73759695,  0.313747346, 1.89878101];
        C = [1.31887070e-2, 6.23068142e-2, 1.55236290e2];

    case 'SK16'
        % Dense crown — nD=1.6203, V=60.41; Cooke triplet crown
        B = [1.34317774,  0.241144399, 0.994317969];
        C = [7.04687339e-3, 2.29005e-2,   9.27508526e1];

    case 'BAF10'
        % Barium flint — nD=1.6700, V=47.11
        B = [1.5851495,   0.143559385, 1.08521269];
        C = [9.26681282e-3, 4.24489805e-2, 1.05613573e2];

    case 'LAFN7'
        % Lanthanum flint — nD=1.7495, V=34.95
        B = [1.66842615,  0.298512803, 1.07743760];
        C = [1.03159999e-2, 4.69216348e-2, 8.27592000e1];

    case {'AIR', 'VACUUM', ''}
        n = ones(size(lambda_nm));
        return;

    otherwise
        error('sellmeier: Unknown glass type ''%s''. Supported: BK7, F2, SF11, SK16, BAF10, LAFN7, air.', glass);
end

% Compute refractive index
n2 = 1 + B(1)*lambda_um2./(lambda_um2 - C(1)) ...
       + B(2)*lambda_um2./(lambda_um2 - C(2)) ...
       + B(3)*lambda_um2./(lambda_um2 - C(3));
n = sqrt(n2);
end
