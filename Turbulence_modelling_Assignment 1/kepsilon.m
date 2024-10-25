function [U_kepsilon_damped, U_kepsilon_undamped] = kepsilon(Re_tau, visc_nu, rho, H)
    N = 1000;
    dy = H / (N - 1); 
    MESH = meshGenerator(H, N);
    U_tau = Re_tau * visc_nu / H;
    y_Plus = MESH.y * U_tau / visc_nu;

    %% Damping function
    f_u_func_damped = @(y_Plus, A_Plus) (1 - exp(-y_Plus / A_Plus)).^2;
    f2_func_damped = @(Re_tau, y_Plus) (1 - 2/9 * exp(-(Re_tau / 6).^2)) .* (1 - exp(-y_Plus / 5)).^2;
 
    f_u_func_undamped = @(y_Plus, A_Plus) (y_Plus * A_Plus * 0 + 1);
    f2_func_undamped = @(Re_tau, y_Plus) (y_Plus * Re_tau * 0 + 1);   

    %% Solver
    % initial guess damped

    visc_nu_eddy = zeros(N, 1);
    U = 180*ones(N, 1); 
    k = 180*ones(N, 1); 
    epsilon = 180*ones(N, 1); 
    alpha_U = 1;  
    alpha_k = 1;  
    alpha_epsilon = 1;  
    
    U_damped_struct = KEpsilonSolver(Re_tau, visc_nu, rho, H, f_u_func_damped, f2_func_damped,visc_nu_eddy,U,k,epsilon,alpha_U,alpha_k,alpha_epsilon,MESH,dy,N);

    %% Struct damped
    U_kepsilon_damped.Re_tau = Re_tau;
    U_kepsilon_damped.y_Plus = y_Plus;
    U_kepsilon_damped.U = U_damped_struct.U;
    U_kepsilon_damped.k = U_damped_struct.k;
    U_kepsilon_damped.epsilon = U_damped_struct.epsilon;
    U_kepsilon_damped.visc_nu_eddy = U_damped_struct.visc_nu_eddy;

    % initial guess undamped
    N = 2000;
    dy = H / (N - 1); 
    MESH = meshGenerator(H, N);
    y_Plus = MESH.y * U_tau / visc_nu;
    visc_nu_eddy = zeros(N, 1);
    U = 200*ones(N, 1); 
    k = 200*ones(N, 1); 
    epsilon = 200*ones(N, 1);
    alpha_U = 1;  
    alpha_k = 1;  
    alpha_epsilon = 1; 
    U_undamped_struct = KEpsilonSolver(Re_tau, visc_nu, rho, H, f_u_func_undamped, f2_func_undamped,visc_nu_eddy,U,k,epsilon,alpha_U,alpha_k,alpha_epsilon,MESH,dy,N);

    % Re_tau 180 Converges for Setting: N = 1500; initial values 200; Relaxation factor 1
    % Re_tau 395 Sor of converges for: N = 1700; Initial values 300: Relaxation factor 1
    % Re_tau 590 Converges for Setting: N = 1000; initial values 200; Relaxation factor 1

    %% Struct undamped
    U_kepsilon_undamped.Re_tau = Re_tau;
    U_kepsilon_undamped.y_Plus = y_Plus;
    U_kepsilon_undamped.U = U_undamped_struct.U;
    U_kepsilon_undamped.k = U_undamped_struct.k;
    U_kepsilon_undamped.epsilon = U_undamped_struct.epsilon;
    U_kepsilon_undamped.visc_nu_eddy = U_undamped_struct.visc_nu_eddy;

end
