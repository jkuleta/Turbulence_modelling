function plotReynoldsStressNearWall(DNSData, Re_tau, visc_mu, U_Turbulent_damped, U_Turbulent_undamped, y, kappa, A_Plus, y_Plus_Turbulent)
    % Function to plot normalized Reynolds stress and viscous stress
    
    %% Extract required data
    mu_du_dy = DNSData.(['DU_dy_', num2str(Re_tau)]) * visc_mu;  % Viscous stress
    rho_uv = -DNSData.(['R_uv_', num2str(Re_tau)]) * 1.225;        % Reynolds stress
    tau = mu_du_dy + rho_uv;

    mu_du_dy_normalized = mu_du_dy ./ tau;
    rho_uv_normalized = rho_uv ./ tau;

    %% Turbulence model calculations
    du_dy_exp_damped = gradient(U_Turbulent_damped, y);
    l_t_damped = -(kappa * y' .* (1 - exp(-y_Plus_Turbulent / A_Plus)));
    Reynolds_stress_damped = -l_t_damped.^2 .* abs(du_dy_exp_damped).^2;

    du_dy_exp_undamped = gradient(U_Turbulent_undamped, y);
    l_t_undamped = kappa * y';
    Reynolds_stress_undamped = -l_t_undamped.^2 .* abs(du_dy_exp_undamped).^2;

    % Calculate normalized stresses for damped and undamped
    tau_damped = visc_mu * du_dy_exp_damped - 1.225 * Reynolds_stress_damped;
    rho_uv_normalized_damped = -1.225 * Reynolds_stress_damped ./ tau_damped;
    mu_du_dy_normalized_damped = visc_mu * du_dy_exp_damped ./ tau_damped;

    tau_undamped = visc_mu * du_dy_exp_undamped - 1.225 * Reynolds_stress_undamped;
    rho_uv_normalized_undamped = -1.225 * Reynolds_stress_undamped ./ tau_undamped;
    mu_du_dy_normalized_undamped = visc_mu * du_dy_exp_undamped ./ tau_undamped;

    % Set line width
    lineWidth = 2.5;  % Increased line width for better visibility

    % Plotting
    figure;
    hold on;

    % Plot DNS data
    plot(DNSData.(['y_plus_', num2str(Re_tau)]), mu_du_dy_normalized, 'k--', ...
        'DisplayName', 'DNS', 'LineWidth', lineWidth);  % DNS data
    plot(DNSData.(['y_plus_', num2str(Re_tau)]), rho_uv_normalized, 'k', ...
        'DisplayName', '', 'LineWidth', lineWidth);  % Hide display name

    % Plot normalized stresses for damped and undamped
    plot(y_Plus_Turbulent, mu_du_dy_normalized_damped, 'b--', ...
        'DisplayName', 'van Driest damping', 'LineWidth', lineWidth);  % Damped
    plot(y_Plus_Turbulent, rho_uv_normalized_damped, 'b', ...
        'DisplayName', '', 'LineWidth', lineWidth);  % Hide display name
    plot(y_Plus_Turbulent, mu_du_dy_normalized_undamped, 'r--', ...
        'DisplayName', 'Prandtl model', 'LineWidth', lineWidth);  % Undamped
    plot(y_Plus_Turbulent, rho_uv_normalized_undamped, 'r', ...
        'DisplayName', '', 'LineWidth', lineWidth);  % Hide display name

    % Add a text box to describe the line styles

    dim = [0.75, 0.4, 0.2, 0.1]; % Position of the textbox [x, y, width, height]
str = {'Dashed line: $\frac{\mu}{\tau(y)} \frac{d \overline{u}}{dy}$', ...
       'Solid line: $\frac{\rho \overline{u^{\prime}v^{\prime}}}{\tau(y)}$'};
annotation('textbox', dim, 'String', str, 'Interpreter', 'latex', ...
    'FontSize', 12, 'FitBoxToText', 'on', 'LineStyle', 'none');


    % Enhance grid and axes appearance
    ax = gca; % Get current axes
    ax.XAxis.LineWidth = 1.5; % Thicker x-axis
    ax.YAxis.LineWidth = 1.5; % Thicker y-axis
    ax.Box = 'on'; % Add a box around the plot
    
    grid on;
    grid minor; % Enable minor grid lines for better visibility
    ax.GridLineStyle = '-'; % Solid grid lines
    ax.MinorGridLineStyle = ':'; % Dotted minor grid lines

    % Set font size for better visibility
    ax.FontSize = 14;  % Increase font size for axes
    xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 16);  % Increase label font size
    %ylabel('$\frac{\mu}{\tau(y)} \frac{d \overline{u}}{dy}$ and $\frac{\rho \overline{u^{\prime}v^{\prime}}}{\tau(y)}$', ...
    %    'Interpreter', 'latex', 'FontSize', 16);  % Increase label font size

    % Display legend without vertical lines
           lgd = legend('show', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 14);  % Increase legend font size
        
        % Adjust legend position to move it outside the axes
        lgd.Position(1) = lgd.Position(1) + 0.15; % Move to the right
        lgd.Position(2) = lgd.Position(2) + 0.1; % Move up

    hold off;

    % Set x-axis limits
    xlim([0 100]);

    % Set y-axis limits
    %ylim([-1, 1]); % Set appropriate y limits based on your data

    % Set figure size
    fig = gcf; % Get current figure
    fig.Position(3) = 1200; % Set width of the figure (in pixels)
    fig.Position(4) = 600;  % Set height of the figure (in pixels)
end
