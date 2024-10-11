clc;clear all;close all

%% Parameters
visc_nu = 1.56*10^(-5); 
rho = 1.225;
visc_mu = visc_nu*rho;
H = 1;  

%% Constants
A_Plus = 26;
kappa = 0.41;

 
%% DNS
Re_tau_values = [180 395 590]; % Pressure grad
U_tau = Re_tau_values.*visc_nu/H;

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);


%% Solver
dy = 0.01; 

height = input('Enter the height condition (full/half): ', 's');
neumann_scheme = input('Enter the Neumann scheme (1/2): ');

if strcmp(height, 'full')
y = 0:dy:2*H;  
elseif strcmp(height, 'half')
y = 0:dy:H;   
end
            
N = length(y); 

%% Forward Euler

U_Turbulent_undamped = zeros(N, length(Re_tau_values));
U_Turbulent_damped = zeros(N, length(Re_tau_values));
U_Plus_undamped = zeros(N, length(Re_tau_values));
U_Plus_damped = zeros(N, length(Re_tau_values));


%% figures
for j = 1:length(Re_tau_values)
    % Solve for the current Re_tau
    [U_Turbulent_undamped(:, j), U_Turbulent_damped(:, j)] = ...
        solveTurbulence(Re_tau_values(j), height, neumann_scheme, visc_nu, H, kappa, A_Plus, N, dy, y);


    y_Plus_Turbulent(:, j) = y * U_tau(j) / visc_nu;
    U_Plus_undamped(:, j) = U_Turbulent_undamped(:, j) / U_tau(j);
    U_Plus_damped(:, j) = U_Turbulent_damped(:, j) / U_tau(j);


    figure;
    semilogx(y_Plus_Turbulent(:, j), U_Plus_undamped(:, j), 'DisplayName', ['Turbulent Undamped (Re_{\tau} = ', num2str(Re_tau_values(j)), ')']);
    hold on;
    semilogx(y_Plus_Turbulent(:, j), U_Plus_damped(:, j), 'DisplayName', ['Turbulent Damped (Re_{\tau} = ', num2str(Re_tau_values(j)), ')']);
    semilogx(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), DNSData.(['U_mean_', num2str(Re_tau_values(j))]), 'k--', 'DisplayName', ['DNS (Re_{\tau} = ', num2str(Re_tau_values(j)), ')']);

    title(['Poiseuille Flow Laminar vs Turbulent (Re_{\tau} = ', num2str(Re_tau_values(j)), ')']);
    grid on;
    xlabel('y^+');
    ylabel('U^+');
    legend('show');
    hold off;

    height_max = max(y);

    figure;
    plot(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max, DNSData.(['R_uu_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['u_{rms} = ', num2str(Re_tau_values(j)), ')']);
    hold on;
    plot(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max, DNSData.(['R_vv_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['v_{rms} = ', num2str(Re_tau_values(j)), ')']);
    plot(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max, DNSData.(['R_ww_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['w_{rms} = ', num2str(Re_tau_values(j)), ')']);
    plot(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max, -DNSData.(['R_uv_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['-uv_{rms} = ', num2str(Re_tau_values(j)), ')']);
    grid on;
    xlabel('y/H');
    %ylabel('$\sqrt{\overline{-u^{\prime} v^{\prime}}}/U_\tau$', 'Interpreter', 'latex');
    legend('show');
    hold off;

    figure;
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), DNSData.(['R_uu_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['u_{rms} = ', num2str(Re_tau_values(j)), ')']);
    hold on;
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), DNSData.(['R_vv_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['v_{rms} = ', num2str(Re_tau_values(j)), ')']);
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), DNSData.(['R_ww_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['w_{rms} = ', num2str(Re_tau_values(j)), ')']);
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), -DNSData.(['R_uv_', num2str(Re_tau_values(j))]),'k--', 'DisplayName', ['-uv_{rms} = ', num2str(Re_tau_values(j)), ')']);
    grid on;
    xlabel('y^+', 'Interpreter', 'latex');
    %ylabel('$\sqrt{\overline{-u^{\prime} v^{\prime}}}/U_\tau$', 'Interpreter', 'latex');
    legend('show');
    hold off;
    

 %% Normalized REynolds stress
    mu_du_dy = DNSData.(['DU_dy_', num2str(Re_tau_values(j))]) * visc_mu;  % Viscous stress
    rho_uv = -DNSData.(['R_uv_', num2str(Re_tau_values(j))]) * rho;        % Reynolds stress
    tau = mu_du_dy + rho_uv;

    mu_du_dy_normalized = mu_du_dy ./ tau;
    rho_uv_normalized = rho_uv ./ tau;


    figure;
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), mu_du_dy_normalized, 'k--', ...
    'DisplayName', ['$\frac{\mu \frac{d \overline{u}}{dy}}{\tau} = ', num2str(Re_tau_values(j)), '$']);
    hold on;
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), rho_uv_normalized, 'b--', ...
    'DisplayName', ['$\frac{\rho \overline{u^{\prime}v^{\prime}}}{\tau} = ', num2str(Re_tau_values(j)), '$']);
    grid on;

    xlabel('$y^+$', 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex');
    xlim([0 12]);
    hold off;


   


end

    figure;
    hold on;

for j = 1:length(Re_tau_values)
    % Solve for the current Re_tau
    [U_Turbulent_undamped(:, j), U_Turbulent_damped(:, j)] = ...
        solveTurbulence(Re_tau_values(j), height, neumann_scheme, visc_nu, H, kappa, A_Plus, N, dy, y);


    y_Plus_Turbulent(:, j) = y * U_tau(j) / visc_nu;
    U_Plus_undamped(:, j) = U_Turbulent_undamped(:, j) / U_tau(j);
    U_Plus_damped(:, j) = U_Turbulent_damped(:, j) / U_tau(j);

  %% Kinetic energy
     k = 0.5*( DNSData.(['R_vv_', num2str(Re_tau_values(j))]).^2 + DNSData.(['R_uu_', num2str(Re_tau_values(j))]).^2+DNSData.(['R_ww_', num2str(Re_tau_values(j))]).^2);
    
    plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), k, '--', ...
        'DisplayName', ['$k = ', num2str(Re_tau_values(j)), '$']);
end

grid on;
xlabel('$y^+$', 'Interpreter', 'latex');
ylabel('$k$', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex');
hold off;


