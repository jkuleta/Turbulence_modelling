function plotVelocityProfiles(Re_tau_values, DNSData, U_Prandalt,U_kepsilon_damped,U_kepsilon_undamped, H,visc_nu)
    lineWidth = 2.5;  

    numReTau = length(Re_tau_values);
    
    for j = 1:numReTau
        U_tau = Re_tau_values(j) * visc_nu / H;
        %% Load Data
        % Mesh
        dy = 0.001; 
        y = linspace(0, H, 1000)';     
        N = length(y);
        U_Prandalt_damped = U_Prandalt(j).U_damped;    
        U_Prandalt_undamped = U_Prandalt(j).U_undamped; 

        %% Normalize Prandalt
        y_Prandalt_Plus(:, j) = y * U_tau / visc_nu;
        U_Prandalt_Plus_undamped(:, j) = U_Prandalt_undamped / U_tau; 
        U_Prandalt_Plus_damped(:, j) = U_Prandalt_damped / U_tau;     

        %% Plot prandalt
        figure;
        hold on;
        semilogx(y_Prandalt_Plus(:, j), U_Prandalt_Plus_undamped(:, j), ...
            'DisplayName', 'Prandtl', ...
            'LineWidth', lineWidth);  
        semilogx(y_Prandalt_Plus(:, j), U_Prandalt_Plus_damped(:, j), ...
            'DisplayName', 'Prandalt van Driest', ...
            'LineWidth', lineWidth);  
        semilogx(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), ...
            DNSData.(['U_mean_', num2str(Re_tau_values(j))]), ...
            'k--', 'DisplayName', 'DNS', ...
            'LineWidth', lineWidth);

        %% Figure - Format from ChatGPT
        % Add vertical lines for layer boundaries without adding to the legend
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
        ylabel('$U^+$', 'Interpreter', 'latex', 'FontSize', 16);  % Increase label font size

        % Display legend without vertical lines
        legend('show', 'Location', 'northwest', 'FontSize', 14);  % Increase legend font size
        hold off;
        
        % Set x-axis limits
        xlim([10^(-1) 10^3]);

        % Set semilog scale on the x-axis
        set(gca, 'XScale', 'log');

        fig = gcf; % Get current figure
        fig.Position(3) = 1200; % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)
        hold off;

        %% Load Data
        % Mesh
        U_kepsilon_damped_current = U_kepsilon_damped(j).U;    
        U_kepsilon_undamped_current = U_kepsilon_undamped(j).U; 

        %% Normalize Kepsilon
        y_kepsilon_Plus_undamped(:, j) = U_kepsilon_undamped(j).y_Plus;
        y_kepsilon_Plus_damped(:, j) = U_kepsilon_damped(j).y_Plus;
        U_kepsilon_Plus_undamped(:, j) = U_kepsilon_undamped_current / U_tau; 
        U_kepsilon_Plus_damped(:, j) = U_kepsilon_damped_current / U_tau;     

        %% Plot prandalt
        figure;
        hold on;
        semilogx(y_kepsilon_Plus_undamped(:, j), U_kepsilon_Plus_undamped(:, j), ...
            'DisplayName', 'k -\epsilon', ...
            'LineWidth', lineWidth);  
        semilogx(y_kepsilon_Plus_damped(:, j), U_kepsilon_Plus_damped(:, j), ...
            'DisplayName', 'k -\epsilon van Driest', ...
            'LineWidth', lineWidth);  
        semilogx(DNSData.(['y_plus_', num2str(Re_tau_values(j))]), ...
            DNSData.(['U_mean_', num2str(Re_tau_values(j))]), ...
            'k--', 'DisplayName', 'DNS', ...
            'LineWidth', lineWidth);

        %% Figure - Format from ChatGPT
        % Add vertical lines for layer boundaries without adding to the legend
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
        ylabel('$U^+$', 'Interpreter', 'latex', 'FontSize', 16);  % Increase label font size

        % Display legend without vertical lines
        legend('show', 'Location', 'northwest', 'FontSize', 14);  % Increase legend font size
        hold off;
        
        % Set x-axis limits
        xlim([10^(-1) 10^3]);

        % Set semilog scale on the x-axis
        set(gca, 'XScale', 'log');

        fig = gcf; % Get current figure
        fig.Position(3) = 1200; % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)
        hold off;
    end
end
