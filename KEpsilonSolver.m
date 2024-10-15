function U_output = KEpsilonSolver(Re_tau, visc_nu, rho, H, f_u_func, f2_func,visc_nu_eddy,U,k,epsilon,alpha_U,alpha_k,alpha_epsilon,MESH,dy,N)

%% Mesh
y = MESH.y;
ddy = MESH.ddy;
d2dy2 = MESH.d2dy2;

%% Constants
A_Plus = 26;
kappa = 0.41;

%k-epslion
C_mu = 0.09;
C_1 = 1.44;
C_2 = 1.92;
sigma_k = 1;
sigma_epsilon = 1.3;


U_tau = Re_tau .* visc_nu / H; 
DP_dx = -rho * Re_tau^2 * visc_nu^2 / H^3;

%% Boundary conditions
visc_nu_eddy(1) = 0;
visc_nu_eddy(end) = visc_nu_eddy(end - 1);
k(1) = 0;
k(end) = k(end - 1);
epsilon(1) = 0;
epsilon(end) = epsilon(end - 1);
U(1) = 0;
U(end) = U(end - 1);

%% Solver 
maxiter = 500;
tolerance = 1e-6;
converged = false;

U_kepsilon = struct();
for i = 1:maxiter
    U_old = U; 
    k_old = k;
    epsilon_old = epsilon;

    dUdy = ddy * U;
    y_Plus = y*U_tau/visc_nu;

    f_u = f_u_func(y_Plus, A_Plus); 
    f2 = f2_func(Re_tau, y_Plus); 

    %% Epsilon loop
    % initialize
    visc_nu_eddy = f_u.* C_mu .* (k.^2) ./ epsilon;

    % Pk
    Pk = visc_nu_eddy .* dUdy.^2;
    visc_term = visc_nu + visc_nu_eddy / sigma_epsilon;

    % A epsilon
    A_epsilon = bsxfun(@times, visc_term, d2dy2) + bsxfun(@times, ddy * visc_term, ddy);
    A_epsilon(1, :) = 0;
    A_epsilon(1, 1) = 1;
    A_epsilon(end, :) = 0;
    A_epsilon(end, end - 1) = -1;
    A_epsilon(end, end) = 1;

    for index = 2:N-1
        A_epsilon(index, index) = A_epsilon(index, index) - C_2* f2(index) * epsilon(index) / max(k(index), 1e-10); 
    end
    
    % B epsilon
    B_epsilon = -C_1 * epsilon ./ max(k, 1e-10) .* Pk; 
    B_epsilon(1) = visc_nu * (k(1) - 2 * k(2) + k(3)) / dy^2;
    B_epsilon(end) = 0;

    % Solve epsilon
    epsilon_new = linsolve(A_epsilon, B_epsilon);
    epsilon_new = max(epsilon_new, 1e-8); 
    epsilon = alpha_epsilon * epsilon_new + (1 - alpha_epsilon) * epsilon_old;

    %% k loop 
    visc_nu_eddy = f_u .* C_mu .* k.^2 ./ epsilon;
    visc_term = visc_nu + visc_nu_eddy / sigma_k;

    % Matrix k
    A_k = bsxfun(@times, visc_term, d2dy2) + bsxfun(@times, ddy * visc_term, ddy);
    A_k(1, :) = 0;
    A_k(1, 1) = 1;
    A_k(end, :) = 0;
    A_k(end, end - 1) = -1;
    A_k(end, end) = 1;

    for j = 2:N-1
        A_k(j, j) = A_k(j, j) - epsilon(j) / k(j); 
    end
    
    % B matrix
    B_k = -Pk;
    B_k(1) = 0;
    B_k(end) = 0;

    % Solve k
    k_new = linsolve(A_k, B_k);
    k_new = max(k_new, 1e-8); 
    k = alpha_k * k_new + (1 - alpha_k) * k_old;
    
    %% U loop
    visc_nu_eddy = f_u.* C_mu .* k.^2 ./ epsilon; 
    visc_term = visc_nu + visc_nu_eddy;

    % A U
    A_U = bsxfun(@times, visc_term, d2dy2) + bsxfun(@times, ddy * visc_term, ddy);
    A_U(1, :) = 0;
    A_U(1, 1) = 1;
    A_U(end, :) = 0;
    A_U(end, end - 1) = -1;
    A_U(end, end) = 1;

    % B U
    B_U = ones(N, 1) * 1 / rho * DP_dx;
    B_U(1) = 0;
    B_U(end) = 0;

    % Solve U
    U_new = linsolve(A_U, B_U);
    U = alpha_U * U_new + (1 - alpha_U) * U_old;

    %% Convergance
    residual = norm(U - U_old);
    residual_history(mod(i-1, 10) + 1) = residual;  % Store the latest residual
    disp(['Iteration: ', num2str(i), ' Residual: ', num2str(residual)]);

    if residual < tolerance
        converged = true;
        disp('Converged');
        break; 
    end

    if i > 10 && max(residual_history) - min(residual_history) < 1e-8
        converged = true;
        disp('Converged');
        break;
    end
end
if ~converged
    warning('Solution did not converge');
end 
    U_output = struct();
    U_output.U = U;                     
    U_output.k = k;                     
    U_output.epsilon = epsilon;        
    U_output.visc_nu_eddy = visc_nu_eddy;  
end
