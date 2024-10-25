function plotReynoldsStress(Re_tau_values,DNSData, visc_nu, U_Prandalt, U_kepsilon_damped, U_kepsilon_undamped,H);

    %% Constants
    H = 1;
    A_Plus = 26;
    kappa = 0.41;
    
    lineWidth = 2.5;  
    numReTau = length(Re_tau_values);
    
    %% Reynolds Stress calc
    for j = 1:numReTau
        U_tau = Re_tau_values(j) * visc_nu / H;
        height_max = max(DNSData.(['y_', num2str(Re_tau_values(j))]));

        %% Prandalt
        U_Prandalt_damped = U_Prandalt(j).U_damped;    
        U_Prandalt_undamped = U_Prandalt(j).U_undamped; 

        %% Mesh
        y = linspace(0, H, 1000)';     
        y_Prandalt_Plus = y * U_tau / visc_nu; % Normalized y+
        
        % Normalize velocities
        U_Prandalt_Plus_damped = U_Prandalt_damped / U_tau; 
        U_Prandalt_Plus_undamped = U_Prandalt_undamped / U_tau;     

        %% Prandtl Stress Calculation  
        % For damped model
        du_dy_exp_damped = gradient(U_Prandalt_damped, y);
        l_t_damped = -(kappa * y .* (1 - exp(-y_Prandalt_Plus / A_Plus)));
        Prandalt_Reynolds_stress_damped = -l_t_damped.^2 .* abs(du_dy_exp_damped).^2;
        Prandalt_uv_rms_damped = sqrt(-Prandalt_Reynolds_stress_damped) / U_tau;
        
        % For undamped model
        du_dy_exp_undamped = gradient(U_Prandalt_undamped, y);
        l_t_undamped = kappa * y; % Without damping
        Prandalt_Reynolds_stress_undamped = -l_t_undamped.^2 .* abs(du_dy_exp_undamped).^2;
        Prandalt_uv_rms_undamped = sqrt(-Prandalt_Reynolds_stress_undamped) / U_tau;

        %% full channel depth y/2H    
        figure;
        plot(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max-1, -DNSData.(['R_uv_', num2str(Re_tau_values(j))]), 'k--', 'DisplayName', 'DNS', 'LineWidth', lineWidth);
        hold on;
        h1 = plot(y/H-1, Prandalt_uv_rms_undamped, 'DisplayName', 'Prandtl model', 'LineWidth', lineWidth);
        color_undamped = get(h1, 'Color');
        h2 = plot(y/H-1, Prandalt_uv_rms_damped, 'DisplayName', 'van Driest damping', 'LineWidth', lineWidth);
        color_damped = get(h2, 'Color');
        % flip side mirror
        plot(-flipud(fliplr(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max))+1, -flipud(fliplr(-DNSData.(['R_uv_', num2str(Re_tau_values(j))]))), 'k--', 'LineWidth', lineWidth,'HandleVisibility', 'off');
        plot(y/H, -flipud(fliplr(Prandalt_uv_rms_undamped)), 'LineWidth', lineWidth,'HandleVisibility', 'off', 'Color', color_undamped);  
        plot(y/H, -flipud(fliplr(Prandalt_uv_rms_damped)), 'LineWidth', lineWidth,'HandleVisibility', 'off', 'Color', color_damped);
        plot([-1 1], [1 -1], 'b', 'DisplayName', '$\frac{\tau}{\tau_w}$', 'LineWidth', lineWidth*0.75); 
        yline(0, 'k--', 'LineWidth', 1,'HandleVisibility', 'off');
        grid on;
        ylabel('-uv_{rms}');
        xlabel('y/H');
        legend('show', 'Interpreter', 'latex');
        ylim([-1 1]);
        hold off;
        
        %% Enhance grid and axes appearance
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
        xlabel('y/H', 'FontSize', 16);  % Increase label font size
        ylabel('-uv_{rms}', 'FontSize', 16);  % Increase label font size

        % Display legend with LaTeX interpretation
        legend('show', 'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex');  % Increase legend font size
        hold off;
        
        % Set figure size
        fig = gcf; % Get current figure
        fig.Position(3) = 1200; % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)

        %% Figure near wall region
       
        figure;
        plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))])/height_max, -DNSData.(['R_uv_', num2str(Re_tau_values(j))]), 'k--', 'DisplayName', 'DNS', 'LineWidth', lineWidth);
        hold on;
        plot(y_Prandalt_Plus, Prandalt_uv_rms_undamped, 'DisplayName', 'Prandtl model', 'LineWidth', lineWidth);
        plot(y_Prandalt_Plus, Prandalt_uv_rms_damped, 'DisplayName', 'van Driest damping', 'LineWidth', lineWidth);
        grid on;
        ylabel('-uv_{rms}');
        xlabel('y^+');
        xlim([0 100]);
        ylim([0 1]);
        legend('show', 'Location', 'southeast');
        hold off;
        
        %% Enhance grid and axes appearance
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
        xlabel('y^+', 'FontSize', 16);  % Increase label font size
        ylabel('-uv_{rms}', 'FontSize', 16);  % Increase label font size

        % Display legend with LaTeX interpretation
        legend('show', 'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'latex');  % Increase legend font size
        hold off;
        
        % Set figure size
        fig = gcf; % Get current figure
        fig.Position(3) = 1200; % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)

        
        %% Kepsilon
        U_kepsilon_damped_current = U_kepsilon_damped(j).U;    
        U_kepsilon_undamped_current = U_kepsilon_undamped(j).U; 
        k_kepsilon_damped_current = U_kepsilon_damped(j).k;    
        k_kepsilon_undamped_current = U_kepsilon_undamped(j).k; 
        epsilon_kepsilon_damped_current = U_kepsilon_damped(j).epsilon;    
        epsilon_kepsilon_undamped_current = U_kepsilon_undamped(j).epsilon; 
        C_mu = 0.09;
         

        %% Mesh
        y_kepsilon_Plus_undamped = U_kepsilon_undamped(j).y_Plus;
        y_kepsilon_Plus_damped = U_kepsilon_damped(j).y_Plus; 
        y_kepsilon_undamped = y_kepsilon_Plus_undamped * visc_nu / U_tau;
        y_kepsilon_damped = y_kepsilon_Plus_damped * visc_nu / U_tau; 

        %% Prandtl Stress Calculation  
        % For damped model
        du_dy_exp_damped = gradient(U_kepsilon_damped_current, y_kepsilon_damped);
        f_u_damped = (1 - exp(-y_kepsilon_Plus_damped / A_Plus)).^2;
        kepsilon_Reynolds_stress_damped = -(f_u_damped.*C_mu.*k_kepsilon_damped_current.^2./epsilon_kepsilon_damped_current) .* abs(du_dy_exp_damped);
        kepsilon_uv_rms_damped = sqrt(-kepsilon_Reynolds_stress_damped) / U_tau;
        
        % For undamped model
        du_dy_exp_undamped = gradient(U_kepsilon_undamped_current,y_kepsilon_undamped);
        f_u_undamped = (1 - y_kepsilon_Plus_undamped *0);
        kepsilon_Reynolds_stress_undamped = -(f_u_undamped.*C_mu.*k_kepsilon_undamped_current.^2./epsilon_kepsilon_undamped_current) .* abs(du_dy_exp_undamped);
        kepsilon_uv_rms_undamped = sqrt(-kepsilon_Reynolds_stress_undamped) / U_tau;

                %% full channel depth y/2H    
        figure;
        plot(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max-1, -DNSData.(['R_uv_', num2str(Re_tau_values(j))]), 'k--', 'DisplayName', 'DNS', 'LineWidth', lineWidth);
        hold on;
        h1 = plot(y_kepsilon_undamped/H-1, kepsilon_uv_rms_undamped, 'DisplayName', '$k-\epsilon$', 'LineWidth', lineWidth);
        color_undamped = get(h1, 'Color');
        h2 = plot(y_kepsilon_damped/H-1, kepsilon_uv_rms_damped, 'DisplayName', '$k-\epsilon$ van Driest', 'LineWidth', lineWidth);
        color_damped = get(h2, 'Color');
        % flip side mirror
        plot(-flipud(fliplr(DNSData.(['y_', num2str(Re_tau_values(j))])/height_max))+1, -flipud(fliplr(-DNSData.(['R_uv_', num2str(Re_tau_values(j))]))), 'k--', 'LineWidth', lineWidth,'HandleVisibility', 'off');
        plot(y_kepsilon_undamped/H, -flipud(fliplr(kepsilon_uv_rms_undamped)), 'LineWidth', lineWidth,'HandleVisibility', 'off', 'Color', color_undamped);  
        plot(y_kepsilon_damped/H, -flipud(fliplr(kepsilon_uv_rms_damped)), 'LineWidth', lineWidth,'HandleVisibility', 'off', 'Color', color_damped);
        plot([-1 1], [1 -1], 'b', 'DisplayName', '$\frac{\tau}{\tau_w}$', 'LineWidth', lineWidth*0.75); 
        yline(0, 'k--', 'LineWidth', 1,'HandleVisibility', 'off');
        grid on;
        ylabel('-uv_{rms}');
        xlabel('y/H');
        legend('show', 'Interpreter', 'latex');
        ylim([-1 1]);
        hold off;
        
        %% Enhance grid and axes appearance
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
        xlabel('y/H', 'FontSize', 16);  % Increase label font size
        ylabel('-uv_{rms}', 'FontSize', 16);  % Increase label font size

        % Display legend with LaTeX interpretation
        legend('show', 'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'latex');  % Increase legend font size
        hold off;
        
        % Set figure size
        fig = gcf; % Get current figure
        fig.Position(3) = 1200; % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)

        %% Figure near wall region
       
        figure;
        plot(DNSData.(['y_plus_', num2str(Re_tau_values(j))])/height_max, -DNSData.(['R_uv_', num2str(Re_tau_values(j))]), 'k--', 'DisplayName', 'DNS', 'LineWidth', lineWidth);
        hold on;
        plot(y_kepsilon_Plus_undamped, kepsilon_uv_rms_undamped, 'DisplayName', 'k-\epsilon', 'LineWidth', lineWidth);
        plot(y_kepsilon_Plus_damped, kepsilon_uv_rms_damped, 'DisplayName', 'k-\epsilon van Driest', 'LineWidth', lineWidth);
        grid on;
        ylabel('-uv_{rms}');
        xlabel('y^+');
        xlim([0 100]);
        ylim([0 1]);
        legend('show', 'Location', 'southeast');
        hold off;
        
        %% Enhance grid and axes appearance
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
        xlabel('y^+', 'FontSize', 16);  % Increase label font size
        ylabel('-uv_{rms}', 'FontSize', 16);  % Increase label font size

        % Display legend with LaTeX interpretation
        legend('show', 'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'latex');  % Increase legend font size
        hold off;
        
        % Set figure size
        fig = gcf; % Get current figure
        fig.Position(3) = 1200; % Set width of the figure (in pixels)
        fig.Position(4) = 600;  % Set height of the figure (in pixels)

    end

end
