% LLE time evolution solver using split-step Fourier method
% ----------------------------+
% Raonaqul Islam, UMBC        |
% Date started: March 2, 2025 |
% Last updated: March 9, 2025 |
% ----------------------------+

clear
clc
close

% Parameters
% -----------
alpha   = 3.2;                            % Detuning
beta    = -0.01;                          % 2nd-order dispersion parameter
F       = 2;                              % Pump amplitude
gamma   = 1;                              % Normalized nonlinear coefficient
tph     = 8.23e-10;                       % Photon lifetime
kappa   = 1/tph;                          % Cavity loss

% Spatial domain discretization
% ------------------------------
N       = 512;                            % No. of points on axis
naxis   = (-N/2:N/2-1).';                 % General axis
dtheta  = 2*pi/N;                         % Spatial step-size
theta   = naxis*dtheta;                   % Spatial domain
dmu     = 1;                              % Mode number domain step-size
mu      = fftshift(dmu*naxis);            % Mode number domain

% Avoided-mode crossing calculations
% ----------------------------------
a = 0.1;                                  % Mode split strength
b = 0;                                    % Mode number at which avoided crossing occurs
Dint = (1i.*(beta/2).*(mu).^2);           % 2nd-order Integrated dispersion
Dint = Dint + a;                          % Avoided-mode crossing
Dint(mu==b) = 0;                          % Set avoided crossing point

% Power and dispersion tuning
% ----------------------------
F_tilde = fft(F*ones(1,N)).';             % Fourier transform of power

% Time domain discretization
% ---------------------------
dt      = 1e-4;                           % Temporal step size
N_tau_p = 100;                            % No. of photon lifetimes
N_steps = N_tau_p/dt;                     % No. of steps

% Generate figure to plot on
figure('Units','inches','Position',[6 4 10 4],'Color','white'); 

% Time evolution
% ---------------
psi_0   = 2*sech(5*theta);                % Initial guess

for ix = 1:N_steps

    % Nonlinear part (half-step)        
    psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

    % Linear part (full-step)
    psi_0_tilde = fft(psi_nl);    
    A_tilde     = -(1+1i*alpha)+Dint;      
    psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
    psi_l       = ifft(psi_l_tilde);
    
    % Nonlinear part (half-step)
    psi_0       = psi_l;
    psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
    psi_0       = psi_nl2; 

    % Output after one-step
    psi_out     = psi_nl2;
    
    % Plot solution in azimuthal domain at every photon lifetime
    if rem(ix,1/dt) == 0
        plot(theta,abs(psi_out).^2,'LineWidth',2.5,'Color','#4D8DBA')
        xlabel('\theta')
        ylabel('|\psi|^2')
        title(['No. of \tau_{p} = ' num2str(ix*dt)])
        set(gca,'FontName','Arial','FontSize',16,'LineWidth',1.2,'TickDir','Out') 
        xlim([-pi pi]);
        xticklabels(strrep(xticklabels,'-','−'))               
        box off
        drawnow
    end
    print('./figures/LLE_AMX_soliton','-dpng','-r600'); % Save figure 
end

% POWER CONVERSION EFFICIENCY
% ---------------------------
w_fft = abs(fftshift(fft(psi_out))).^2;
u     = fftshift(mu);
eta   = sum(w_fft(u~=0))/(F^2*sqrt(kappa))*100;

% Plot soliton in mode number domain
% -----------------------------------
% Prepare the solution
psi_fft  = abs(fftshift(fft(psi_out))).^2;  % Obtain fourier transform of the solution
psi_norm = psi_fft/max(psi_fft);            % Normalize the output
psi_db   = 10*log10(psi_norm);              % Calculate normalized output power in dB

% Plot
figure('Units','inches','Position',[6 4 10 4],'Color','White');
soliton = stem(fftshift(mu),psi_db,'LineWidth',2,'Marker','None','Color','#89A482');
soliton.BaseValue = -400; % Turn the plot upside down
xlabel('Relative mode number, \mu')
ylabel('Power (dB)')
set(gca,'FontName','Arial','FontSize',16,'LineWidth',1.2,'TickDir','Out')  
xlim([-100 100])
ylim([-80 0]);
xticklabels(strrep(xticklabels,'-','−'))
yticklabels(strrep(yticklabels,'-','−'))
box off
print('./figures/LLE_AMX_comb','-dpng','-r600'); % Save figure 

% Display conversion efficiency
fprintf('For alpha = %.2f, F = %.2f, power conversion efficiency = %.2f%% \n',alpha,F,eta);