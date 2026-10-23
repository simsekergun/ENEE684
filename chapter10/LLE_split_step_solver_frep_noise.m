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
sigma   = 1;                              % Noise variance

% Spatial domain discretization
% ------------------------------
N       = 512;                            % No. of points on axis
naxis   = (-N/2:N/2-1).';                 % General axis
dtheta  = 2*pi/N;                         % Spatial step-size
theta   = naxis*dtheta;                   % Spatial domain
dmu     = 1;                              % Mode number domain step-size
mu      = fftshift(dmu*naxis);            % Mode number domain


% Power and dispersion tuning
% ----------------------------
F_tilde = fft(F*ones(1,N)).';             % Fourier transform of power

% Time domain discretization
% ---------------------------
dt      = 1e-4;                           % Temporal step size
N_tau_p = 100;                            % No. of photon lifetimes
N_steps = N_tau_p/dt;                     % No. of steps
theta0  = zeros(N_steps,1);               % Initialize pulse-center-position array

% Generate figure to plot on
figure('Units','inches','Position',[6 4 10 8],'Color','white'); 

% Time evolution
% ---------------
psi_0   = 2*sech(5*theta);                % Initial guess

for ix = 1:N_steps

    % Nonlinear part (half-step)        
    psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

    % Linear part (full-step)
    psi_0_tilde = fft(psi_nl);    
    A_tilde     = -(1+1i*alpha)+1i*sigma*randn*mu+1i.*beta.*(mu).^2;      
    psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
    psi_l       = ifft(psi_l_tilde);
    
    % Nonlinear part (half-step)
    psi_0       = psi_l;
    psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
    psi_0       = psi_nl2; 

    % Output after one-step
    psi_out     = psi_nl2;

    % Calculate and save pulse-centroid
    I           = abs(psi_out).^2;               % Intensity
    theta0(ix)  = sum(I.*theta)/sum(I);          % Centroid
    
    % Plot solution in azimuthal domain at every photon lifetime
    if rem(ix,1/dt) == 0
        % Plot soliton pulse
        subplot(211)
        plot(theta,abs(psi_out).^2,'LineWidth',2.5,'Color','#4D8DBA')
        xlabel('\theta')
        ylabel('|\psi|^2')
        title(['No. of \tau_{p} = ' num2str(ix*dt)])
        set(gca,'FontName','Arial','FontSize',16,'LineWidth',1.2,'TickDir','Out') 
        xlim([-pi pi]);
        xticklabels(strrep(xticklabels,'-','−'))               
        box off
        xline(0,'--')
        drawnow
        % Plot position fluctuations
        subplot(212)
        plot(1:ix,theta0(1:ix),'LineWidth',2.5,'Color','#2391FF')
        xlabel('Time Step')
        ylabel('Azimuthal Position')
        set(gca,'FontName','Arial','FontSize',16,'LineWidth',1.2,'TickDir','Out')  
        box off
        yline(0,'--')
        xlim([1 N_steps])
        drawnow
    end
end
print('./figures/LLE_frep_noise','-dpng','-r600'); % Save figure 