function plotReynoldsStress(DNSData, Re_tau, visc_mu, U_Turbulent_damped, U_Turbulent_undamped, dy, kappa, A_Plus, y_Plus_Turbulent,U_tau,height)
    height_max = max(DNSData.(['y_', num2str(Re_tau)]));

    H=1;
if strcmp(height, 'full')
y = 0:dy:H;  
elseif strcmp(height, 'half')
y = 0:dy:H;   
end
    
    du_dy_exp_damped = gradient(U_Turbulent_damped, y);
    l_t_damped = -(kappa * y' .* (1 - exp(-y_Plus_Turbulent / A_Plus)));
    Reynolds_stress_damped = -l_t_damped.^2 .* abs(du_dy_exp_damped).^2;
    uv_rms_damped = sqrt(-Reynolds_stress_damped)/U_tau;

    du_dy_exp_undamped = gradient(U_Turbulent_undamped, y);
    l_t_undamped = kappa * y';
    Reynolds_stress_undamped = -l_t_undamped.^2 .* abs(du_dy_exp_undamped).^2;
    uv_rms_undamped = sqrt(-Reynolds_stress_undamped)/U_tau;
    
y_plot = -H:2*dy:H;

    for j = 1:length(Re_tau)
        % Reynolds stress plot normalized by y/H
figure;
% Plotting the DNS and model data as per the original code
plot(DNSData.(['y_', num2str(Re_tau)])/height_max-1, -DNSData.(['R_uv_', num2str(Re_tau)]), 'k--', 'DisplayName', '-uv_{rms}');
hold on;
plot(y/H-1, uv_rms_undamped, 'r', 'DisplayName', 'Prandlt model');
plot(y/H-1, uv_rms_damped, 'b', 'DisplayName', 'van Driest damping');

% Plotting flipped data for the last three plots
plot(-flipud(fliplr(DNSData.(['y_', num2str(Re_tau)])/height_max))+1, -flipud(fliplr(-DNSData.(['R_uv_', num2str(Re_tau)]))), 'k--');
plot(y/H, -flipud(fliplr(uv_rms_undamped)), 'r');  % Flipping horizontally
plot(y/H, -flipud(fliplr(uv_rms_damped)), 'b'); 

        grid on;
        ylabel('-uv_{rms}');
        xlabel('y/H');
        legend('show');
        ylim([-1 1]);
        hold off;
    end

        for j = 1:length(Re_tau)
        % Reynolds stress plot normalized by y/H
        figure;
        %plot(DNSData.(['y_plus_', num2str(Re_tau)])/height_max, DNSData.(['R_uu_', num2str(Re_tau)]), 'k--', 'DisplayName', 'u_{rms}');
        hold on;
        %plot(DNSData.(['y_plus_', num2str(Re_tau)])/height_max, DNSData.(['R_vv_', num2str(Re_tau)]), 'k--', 'DisplayName', 'v_{rms}');
        %plot(DNSData.(['y_plus_', num2str(Re_tau)])/height_max, DNSData.(['R_ww_', num2str(Re_tau)]), 'k--', 'DisplayName', 'w_{rms}');
        plot(DNSData.(['y_plus_', num2str(Re_tau)])/height_max, -DNSData.(['R_uv_', num2str(Re_tau)]), 'k--', 'DisplayName', 'DNS');
        hold on;
        plot(y_Plus_Turbulent, uv_rms_undamped, 'DisplayName', 'Prandlt model');
        plot(y_Plus_Turbulent, uv_rms_damped, 'DisplayName', 'van Driest damping');
        grid on;
        ylabel('-uv_{rms}');
        xlabel('y^+');
        xlim([0 100]);
        legend('show','Location','southeast');
        hold off;
    end
end

