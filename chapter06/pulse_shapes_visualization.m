% =========================================================
%  pulse_shapes_visualization.m
%
%  Visualizes five common pulse shapes and their analytic
%  Fourier transforms (intensity and spectral density).
%
%  Pulses:  Gaussian | Sech | Rect | Lorentzian | Double-Exp.
%% =========================================================
clear; 
%% ---- Parameters ------------------------------------------
tau   = 1;                              % pulse width parameter [normalised]
 
% Time axis
Nt    = 10000;
t_max = 7 * tau;
t     = linspace(-t_max, t_max, Nt);   % [s / normalised units]
 
% Angular-frequency axis
Nw    = 1000;
w_max = 30 / tau;
omega = linspace(-w_max, w_max, Nw);   % [rad/s / normalised units]
 
%% ---- Compute pulses and their Fourier transforms ---------
[Ft_G, Fw_G] = pulse_gaussian   (t, omega, tau);
[Ft_S, Fw_S] = pulse_sech       (t, omega, tau);
[Ft_R, Fw_R] = pulse_rect       (t, omega, tau);
[Ft_L, Fw_L] = pulse_lorentzian (t, omega, tau);
[Ft_D, Fw_D] = pulse_doubleexp  (t, omega, tau);
 
% Collect into cell arrays for looped plotting
labels  = {'Gaussian',   'Sech',        'Rect', ...
           'Lorentzian', 'Double-Exp.'};
Ft_all  = {Ft_G, Ft_S, Ft_R, Ft_L, Ft_D};
Fw_all  = {Fw_G, Fw_S, Fw_R, Fw_L, Fw_D};
 
% FWHM pulse widths (in units of tau) and time-bandwidth products
Delta_t_coeff = [2*sqrt(log(2)), 1.7627, 1.0, 1.2872, log(2)]; % * tau
TBP           = [0.441,          0.315,  0.886, 0.142,  0.142 ];
 
% Colour palette (matches MATLAB default but in a fixed order)
clr = [0.00 0.45 0.74;   % blue   – Gaussian
       0.85 0.33 0.10;   % red    – Sech
       0.47 0.67 0.19;   % green  – Rect
       0.49 0.18 0.56;   % purple – Lorentzian
       0.93 0.69 0.13];  % gold   – Double-Exp.
 
%% =========================================================
%  FIGURE 1 — Individual subplots (5 rows × 2 columns)
%% =========================================================
fig1 = figure('Name','Pulse Shapes & Spectra (Individual)', ...
              'NumberTitle','off','Position',[30 30 1300 980]);
 
for k = 1:5
    % Normalised intensity (peak = 1)
    It = abs(Ft_all{k}).^2;   It = It / max(It);
    Iw = abs(Fw_all{k}).^2;   Iw = Iw / max(Iw);
 
    % ---- Time domain ----------------------------------------
    subplot(5, 2, 2*k-1);
    plot(t/tau, It, 'Color', clr(k,:), 'LineWidth', 2);
    hold on;
 
    % Shade FWHM region
    mask_t = It >= 0.5;
    t_fill  = t(mask_t)/tau;
    It_fill = It(mask_t);
    fill([t_fill, fliplr(t_fill)], ...
         [It_fill, zeros(1,sum(mask_t))], ...
         clr(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
 
    yline(0.5, '--k', 'FWHM', 'LineWidth', 0.9, ...
          'LabelHorizontalAlignment','right','FontSize',8);
 
    % Annotate FWHM span
    if any(mask_t)
        t_lo = t(find(mask_t,1,'first')) / tau;
        t_hi = t(find(mask_t,1,'last'))  / tau;
        dt_fwhm = t_hi - t_lo;
        text(0, 0.57, sprintf('\\Deltat/\\tau = %.4f', dt_fwhm), ...
             'HorizontalAlignment','center','FontSize',8,'Color',clr(k,:));
    end
 
    title(sprintf('%s   |   \\Deltat = %.4f\\tau   |   TBP = %.3f', ...
                  labels{k}, Delta_t_coeff(k), TBP(k)), 'FontSize', 11);
    xlabel('t / \tau', 'FontSize', 11);
    ylabel('|F(t)|^2  (norm.)', 'FontSize', 11);
    xlim([-5 5]);  ylim([0 1.15]);
    grid on;  box on;
    if k == 1
        text(-4.8, 1.10, 'TIME DOMAIN', 'FontWeight','bold','FontSize',10, ...
             'Color',[0.3 0.3 0.3]);
    end
 
    % ---- Frequency domain -----------------------------------
    subplot(5, 2, 2*k);
 
    f_norm = omega / (2*pi);          % convert rad/s → cycles/s, still normalised by tau
    plot(f_norm*tau, Iw, 'Color', clr(k,:), 'LineWidth', 2);
    hold on;
 
    mask_w  = Iw >= 0.5;
    fw_fill = f_norm(mask_w)*tau;
    Iw_fill = Iw(mask_w);
    fill([fw_fill, fliplr(fw_fill)], ...
         [Iw_fill, zeros(1,sum(mask_w))], ...
         clr(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
 
    yline(0.5, '--k', 'FWHM', 'LineWidth', 0.9, ...
          'LabelHorizontalAlignment','right','FontSize',8);
 
    % Annotate FWHM span
    if any(mask_w)
        fw_lo   = f_norm(find(mask_w,1,'first'))*tau;
        fw_hi   = f_norm(find(mask_w,1,'last' ))*tau;
        df_fwhm = fw_hi - fw_lo;
        text(0, 0.57, sprintf('\\Deltaf\\cdot\\tau = %.4f', df_fwhm), ...
             'HorizontalAlignment','center','FontSize',8,'Color',clr(k,:));
    end
 
    title(sprintf('Spectrum: %s', labels{k}), 'FontSize', 11);
    xlabel('f\tau  ( = \omega\tau / 2\pi )', 'FontSize', 11);
    ylabel('|F(\omega)|^2  (norm.)', 'FontSize', 11);
    xlim([-3 3]);  ylim([0 1.15]);
    grid on;  box on;
    if k == 1
        text(-2.88, 1.10, 'FREQUENCY DOMAIN', 'FontWeight','bold','FontSize',10, ...
             'Color',[0.3 0.3 0.3]);
    end
end
 
sgtitle('Pulse Shapes and Their Fourier Transforms', ...
        'FontSize', 14, 'FontWeight', 'bold');
pause(1)
% print -dpng ../../figures/pulse_shapes_and_spectra
 % keyboard
%% =========================================================
%  FIGURE 2 — Overlay comparison (linear scale)
%% =========================================================
fig2 = figure('Name','Pulse Shapes — Overlay (Linear)', ...
              'NumberTitle','off','Position',[60 60 1050 440]);
 
subplot(1,2,1);
for k = 1:5
    It = abs(Ft_all{k}).^2;
    plot(t/tau, It/max(It), 'Color', clr(k,:), 'LineWidth', 2);
    hold on;
end
yline(0.5,'--k','LineWidth',0.8);
xlabel('t / \tau','FontSize',11);
ylabel('|F(t)|^2  (norm.)','FontSize',11);
title('Pulse Shapes — Time Domain','FontSize',12);
legend(labels,'Location','north','NumColumns',2,'FontSize',9);
xlim([-5 5]);  ylim([0 1.12]);
grid on;  box on;
 
subplot(1,2,2);
for k = 1:5
    Iw = abs(Fw_all{k}).^2;
    plot(omega/(2*pi)*tau, Iw/max(Iw), 'Color', clr(k,:), 'LineWidth', 2);
    hold on;
end
yline(0.5,'--k','LineWidth',0.8);
xlabel('f\tau  ( = \omega\tau / 2\pi )','FontSize',11);
ylabel('|F(\omega)|^2  (norm.)','FontSize',11);
title('Spectra — Frequency Domain','FontSize',12);
legend(labels,'Location','north','NumColumns',2,'FontSize',9);
xlim([-3 3]);  ylim([0 1.12]);
grid on;  box on;
 
sgtitle('Overlay Comparison of Pulse Shapes','FontSize',13,'FontWeight','bold');
 
%% =========================================================
%  FIGURE 3 — Overlay comparison (log scale for spectrum)
%% =========================================================
fig3 = figure('Name','Pulse Shapes — Overlay (Log Spectrum)', ...
              'NumberTitle','off','Position',[90 90 1050 440]);
 
subplot(1,2,1);
for k = 1:5
    It = abs(Ft_all{k}).^2;
    plot(t/tau, It/max(It), 'Color', clr(k,:), 'LineWidth', 2);
    hold on;
end
yline(0.5,'--k','LineWidth',0.8);
xlabel('t / \tau','FontSize',11);
ylabel('|F(t)|^2  (norm.)','FontSize',11);
title('Pulse Shapes — Time Domain','FontSize',12);
legend(labels,'Location','north','NumColumns',2,'FontSize',9);
xlim([-5 5]);  ylim([0 1.12]);
grid on;  box on;
 
subplot(1,2,2);
for k = 1:5
    Iw = abs(Fw_all{k}).^2;
    semilogy(omega/(2*pi)*tau, Iw/max(Iw), 'Color', clr(k,:), 'LineWidth', 2);
    hold on;
end
yline(0.5,'--k','LineWidth',0.8);
xlabel('f\tau  ( = \omega\tau / 2\pi )','FontSize',11);
ylabel('|F(\omega)|^2  (norm., log scale)','FontSize',11);
title('Spectra — Log Scale','FontSize',12);
legend(labels,'Location','southwest','NumColumns',1,'FontSize',9);
xlim([-5 5]);  ylim([1e-4 1.5]);
grid on;  box on;
 
sgtitle('Overlay Comparison — Log Spectrum','FontSize',13,'FontWeight','bold');
 
%% =========================================================
%  FIGURE 4 — Time-Bandwidth Product Summary
%% =========================================================
fig4 = figure('Name','Time-Bandwidth Products','NumberTitle','off', ...
              'Position',[120 120 520 380]);
 
bar_h = bar(TBP, 0.55, 'FaceColor','flat');
for k = 1:5
    bar_h.CData(k,:) = clr(k,:);
end
set(gca, 'XTickLabel', labels, 'FontSize', 11);
ylabel('\Deltat \cdot \Deltaf  (TBP)','FontSize',12);
title('Time-Bandwidth Products (FWHM)','FontSize',12,'FontWeight','bold');
ylim([0 1.1]);
grid on;  box on;
text(1:5, TBP+0.03, arrayfun(@(x)sprintf('%.3f',x), TBP, 'UniformOutput',false), ...
     'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
pause(1)
 print -dpng ../../figures/pulse_dt_df_products
%% =========================================================
 
function [Ft, Fw] = pulse_gaussian(t, omega, tau)
% PULSE_GAUSSIAN  Gaussian pulse envelope and its Fourier transform.
%
%   F(t) = exp( -t^2 / (2*tau^2) )
%   F(w) = sqrt(2*pi)*tau * exp( -tau^2*omega^2/2 )
%
%   FWHM:  Delta_t = 2*sqrt(ln2)*tau  ≈ 1.6651*tau
%   TBP:   Delta_t * Delta_f = 0.4413
 
    Ft = exp(-0.5 * (t ./ tau).^2);
    Fw = sqrt(2*pi) * tau * exp(-0.5 * (tau .* omega).^2);
end
 
function [Ft, Fw] = pulse_sech(t, omega, tau)
% PULSE_SECH  Hyperbolic-secant pulse and its Fourier transform.
%
%   F(t) = sech( t/tau )
%   F(w) = (pi*tau/2) * sech( pi*tau*omega/2 )
%
%   FWHM:  Delta_t = 2*acosh(sqrt(2))*tau/pi ... = 1.7627*tau
%   TBP:   Delta_t * Delta_f = 0.3148
 
    Ft = sech(t ./ tau);
    Fw = (pi * tau / 2) .* sech(pi .* tau .* omega ./ 2);
end
 
function [Ft, Fw] = pulse_rect(t, omega, tau)
% PULSE_RECT  Rectangular (top-hat) pulse and its Fourier transform.
%
%   F(t) = Pi(t/tau) = 1 for |t| <= tau/2,  0 otherwise
%   F(w) = tau * sin(tau*omega/2) / (tau*omega/2)
%        = tau * sinc( tau*omega / (2*pi) )   [MATLAB sinc convention]
%
%   FWHM:  Delta_t = tau
%   TBP:   Delta_t * Delta_f = 0.8859
%
%   Note: MATLAB's sinc(x) = sin(pi*x)/(pi*x), so the argument is omega/(2*pi).
 
    Ft = double(abs(t) <= tau/2);
    Fw = tau .* sinc(tau .* omega ./ (2*pi));   % sinc(x)=sin(πx)/(πx) in MATLAB
end
 
function [Ft, Fw] = pulse_lorentzian(t, omega, tau)
% PULSE_LORENTZIAN  Lorentzian pulse and its Fourier transform.
%
%   F(t) = 1 / ( 1 + (t/tau)^2 )
%   F(w) = pi*tau * exp( -|tau*omega| )
%
%   FWHM:  Delta_t = 2*sqrt(sqrt(2)-1)*tau ≈ 1.2872*tau
%   TBP:   Delta_t * Delta_f = 0.1421
%
%   Note: some texts quote the prefactor as 2*pi*tau (different convention).
%   The normalised plots are identical either way.
 
    Ft = 1 ./ (1 + (t ./ tau).^2);
    Fw = pi .* tau .* exp(-abs(tau .* omega));
end
 
function [Ft, Fw] = pulse_doubleexp(t, omega, tau)
% PULSE_DOUBLEEXP  Double-sided exponential pulse and its Fourier transform.
%
%   F(t) = exp( -|t|/tau )
%   F(w) = 2*tau / ( 1 + (omega*tau)^2 )
%
%   FWHM:  Delta_t = ln(2)*tau ≈ 0.6931*tau
%   TBP:   Delta_t * Delta_f = 0.1421
 
    Ft = exp(-abs(t) ./ tau);
    Fw = 2 .* tau ./ (1 + (omega .* tau).^2);
end
 