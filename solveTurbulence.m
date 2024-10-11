function [U_Turbulent_undamped U_Turbulent_damped] = solveTurbulence(Re_tau_values, height, neumann_scheme, visc_nu, H, kappa, A_Plus, N, dy,y)

    U_Turbulent_all = cell(length(Re_tau_values), 2); % Column for damped and undamped

    for j = 1:length(Re_tau_values)
        % Current Re_tau and U_tau calculation
        Re_tau = Re_tau_values(j);
        U_tau = Re_tau * visc_nu / H;
        DP_dx = -Re_tau^2*visc_nu^2/H^3;
        
        % Define y and y_plus arrays
      
        y_Plus_Turbulent = y * U_tau / visc_nu;

        %% Damped Model %%
        % Initialize A_Turbulent matrix and B_Turbulent vector for damped model
        A_Turbulent = zeros(N, N);
        l_t = kappa * y .* (1 - exp(-y_Plus_Turbulent / A_Plus));  % Damped model

        % Calculate B_Turbulent for the damped model
        inside_sqrt = visc_nu^2 - 4 * DP_dx * (H - y) .* (l_t).^2;
        B_Turbulent = [(-visc_nu + sqrt(inside_sqrt)) ./ (2 * (l_t).^2) .* dy]';

        % Fill the A_Turbulent matrix for interior points
        for i = 2:N-1
            A_Turbulent(i, i-1) = -1;    % Coefficient U(y-dy)
            A_Turbulent(i, i) = 1;       % Coefficient U(y)
            A_Turbulent(i, i+1) = 0;     % Coefficient U(y+dy)
        end

        % Boundary conditions for B_Turbulent
        B_Turbulent(1) = 0;             % Initial condition
        B_Turbulent(end) = 0;           % Boundary condition

        % Boundary conditions for A_Turbulent based on height and Neumann scheme
        if strcmp(height, 'full')
            A_Turbulent(1, 1) = 1;
            A_Turbulent(1, 2) = 0;
            A_Turbulent(end, end-1) = 0;
            A_Turbulent(end, end) = 1;
        elseif strcmp(height, 'half')
            A_Turbulent(1, 1) = 1;
            A_Turbulent(1, 2) = 0;
            if neumann_scheme == 1
                A_Turbulent(end, end-1) = -1;
                A_Turbulent(end, end) = 1;
            elseif neumann_scheme == 2
                A_Turbulent(end, end-2) = -3/2;
                A_Turbulent(end, end-1) = 2;
                A_Turbulent(end, end) = -0.5;
            else
                error('Invalid Neumann scheme');
            end
        else
            error('Invalid height condition. Choose either "full" or "half".');
        end

        % Solve the system for U_Turbulent (damped model)

        U_Turbulent_damped = inv(A_Turbulent) * B_Turbulent;

        %% Undamped Model %%
        % Initialize A_Turbulent matrix and B_Turbulent vector for undamped model
        A_Turbulent = zeros(N, N);
        l_t = kappa * y;  % Undamped model

        % Calculate B_Turbulent for the undamped model
        inside_sqrt = visc_nu^2 - 4 * DP_dx * (H - y) .* (l_t).^2;
        B_Turbulent = [(-visc_nu + sqrt(inside_sqrt)) ./ (2 * (l_t).^2) .* dy]';

        % Fill the A_Turbulent matrix for interior points
        for i = 2:N-1
            A_Turbulent(i, i-1) = -1;    % Coefficient U(y-dy)
            A_Turbulent(i, i) = 1;       % Coefficient U(y)
            A_Turbulent(i, i+1) = 0;     % Coefficient U(y+dy)
        end

        % Boundary conditions for B_Turbulent
        B_Turbulent(1) = 0;             % Initial condition
        B_Turbulent(end) = 0;           % Boundary condition

        % Boundary conditions for A_Turbulent based on height and Neumann scheme
        if strcmp(height, 'full')
            A_Turbulent(1, 1) = 1;
            A_Turbulent(1, 2) = 0;
            A_Turbulent(end, end-1) = 0;
            A_Turbulent(end, end) = 1;
        elseif strcmp(height, 'half')
            A_Turbulent(1, 1) = 1;
            A_Turbulent(1, 2) = 0;
            if neumann_scheme == 1
                A_Turbulent(end, end-1) = -1;
                A_Turbulent(end, end) = 1;
            elseif neumann_scheme == 2
                A_Turbulent(end, end-2) = -3/2;
                A_Turbulent(end, end-1) = 2;
                A_Turbulent(end, end) = -0.5;
            else
                error('Invalid Neumann scheme');
            end
        else
            error('Invalid height condition. Choose either "full" or "half".');
        end

        % Solve the system for U_Turbulent (undamped model)
        U_Turbulent_undamped = inv(A_Turbulent) * B_Turbulent;

end
