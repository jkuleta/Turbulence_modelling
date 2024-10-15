function plotVelocityComparison(Re_tau_values, DNSData, U_Prandalt, U_kepsilon_damped, H, visc_nu)
    lineWidth = 2.5;
    numReTau = length(Re_tau_values);
    
    for j = 1:numReTau
        U_tau = Re_tau_values(j) * visc_nu / H;
        
        %% Load Data
        % Mesh
        dy = 0.001; 
        y = linspace(0, H, 1000)';     
        U_Prandalt_damped = U_Prandalt(j).U_damped;    
        U_kepsilon_damped_current = U_kepsilon_damped(j).U; 

        %% Normalize Prandalt and K-epsilon
        y_Plus(:, j) = y * U_tau / visc_nu;
        U_Prandalt_Plus_damped(:, j) = U_Prandalt_damped / U_tau;     
        U_kepsilon_Plus_damped(:, j) = U_kepsilon_damped_current / U_tau;     

        %% Create a new figure for each Re_tau value
        figure;
        hold on;
        
        % Plot Prandalt (blue)
        semilogx(y_Plus(:, j), U_Prandalt_Plus_damped(:, j), 'b-', ...
            'DisplayName', 'Prandtl van Driest', 'LineWidth', lineWidth);  

        % Plot k-epsilon (red)
        semilogx(y_Plus(:, j), U_kepsilon_Plus_damped(:, j), 'r-', ...
            'DisplayName', 'k-\epsilon van Driest', 'LineWidth', lineWidth);

        % Plot DNS data (black dashed line)
        semilogx(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), ...
            DNSData.(['U_mean_', num2str(Re_tau_values(j))]), 'k--', ...
            'DisplayName', 'DNS', 'LineWidth', lineWidth);

        %% Add vertical lines for layer boundaries without adding to the legend
        y_viscous = 5;
        y_buffer = 30;
        y_log = 0.3 * H * U_tau / visc_nu;  % Adjust log layer according to y/H

        % Plot vertical lines (without adding to legend)
        xline(y_viscous, 'k-', 'y^+ = 5', 'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'FontSize', 16);
        xline(y_buffer, 'k-', 'y^+ = 30', 'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'FontSize', 16);
        xline(y_log, 'k-', 'y^+ = 0.3y/H', 'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'FontSize', 16);

        %% Enhance grid and axes appearance
        ax = gca;  % Get current axes
        ax.XAxis.LineWidth = 1.5;  % Thicker x-axis
        ax.YAxis.LineWidth = 1.5;  % Thicker y-axis
        ax.Box = 'on';  % Add a box around the plot
        
        grid on;
        grid minor;  % Enable minor grid lines for better visibility
        ax.GridLineStyle = '-';  % Solid grid lines
        ax.MinorGridLineStyle = ':';  % Dotted minor grid lines

        % Set font size for better visibility
        ax.FontSize = 14;  % Increase font size for axes
        xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 16);  % Increase label font size
        ylabel('$U^+$', 'Interpreter', 'latex', 'FontSize', 16);  % Increase label font size

        % Display legend without vertical lines
        legend('show', 'Location', 'northwest', 'FontSize', 14);  % Increase legend font size
        hold off;
        
        % Set x-axis limits
        xlim([10^(-1) 10^3]);

        % Set semilog scale on the x-axis
        set(gca, 'XScale', 'log');

        fig = gcf;  % Get current figure
        fig.Position(3) = 1200;  % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)
        hold off;
    end
end
