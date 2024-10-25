function U_Prandalt = PrandaltSolver(Re_tau_values, visc_nu, H)

%% Constants
A_Plus = 26;
kappa = 0.41;

% sturct
U_Prandalt = struct();

%% Mesh
dy = 0.001; 
y = linspace(0, H, 1000)';     
N = length(y);

for j = 1:length(Re_tau_values)

    %% Re_tau, U_tau parameters
    Re_tau = Re_tau_values(j);
    U_tau = Re_tau * visc_nu / H;
    DP_dx = -Re_tau^2 * visc_nu^2 / H^3;
    
    % Y plus
    y_Plus = y * U_tau / visc_nu;
    
    %% damping loop
    model_types = {'damped', 'undamped'};
    U_Output = zeros(N, 2); 

    for model_index = 1:length(model_types)
        model_type = model_types{model_index};
        

        A = zeros(N, N);
        

        if strcmp(model_type, 'damped')
            l_t = kappa * y .* (1 - exp(-y_Plus / A_Plus));  % damped 
        else
            l_t = kappa * y;  % undamped
        end
        
        %% Find Matrix B
        inside_sqrt = visc_nu^2 - 4 * DP_dx * (H - y) .* (l_t).^2;
        B = (-visc_nu + sqrt(inside_sqrt)) ./ (2 * (l_t).^2) .* dy;

        %% Forward Scheme
        for i = 2:N-1
            A(i, i-1) = -1;    % coeff U(y-dy)
            A(i, i) = 1;       % coeff U(y)
            A(i, i+1) = 0;     % coeff U(y+dy)
        end

        %% BC
        % B - Matrix
        B(1) = 0;             
        B(end) = 0;          

        % A - Matrix
        A(1, 1) = 1;
        A(1, 2) = 0;
        A(end, end-1) = -1;
        A(end, end) = 1;

        %% Solve
        U_Output(:, model_index) = inv(A) * B;
    end

    %% struct
    U_Prandalt(j).Re_tau = Re_tau; 
    U_Prandalt(j).U_damped = U_Output(:, 1); 
    U_Prandalt(j).U_undamped = U_Output(:, 2); 
end

end
