function plotKineticEnergy(Re_tau_values, DNSData)
    figure;
    hold on;
    for j = 1:length(Re_tau_values)
        k = 0.5 * (DNSData.(['R_vv_', num2str(Re_tau_values(j))]).^2 + ...
            DNSData.(['R_uu_', num2str(Re_tau_values(j))]).^2 + ...
            DNSData.(['R_ww_', num2str(Re_tau_values(j))]).^2);
        plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), k, '--', 'DisplayName', ['$k = ', num2str(Re_tau_values(j)), '$']);
    end
    grid on;
    xlabel('$y^+$', 'Interpreter', 'latex');
    ylabel('$k$', 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex');
    hold off;
end
