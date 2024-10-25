clc; clear; close all;

%% Parameters
visc_nu = 1.56 * 10^(-5); 
rho = 1.225;
visc_mu = visc_nu * rho;
H = 1;  

%% Constants
A_Plus = 26;
kappa = 0.41;

%% DNS
Re_tau_values = [180 395 590];  
U_tau = Re_tau_values .* visc_nu / H;

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);

%% Solver
N = 1000;
dy = H / (N - 1); 

MESH = mesh_chat(H, N);
y = MESH.y;
ddy = MESH.ddy;
d2dy2 = MESH.d2dy2;

%% Troubleshoot
visc_nu_eddy = zeros(N, 1);
Re_tau = Re_tau_values(1);
DP_dx = -rho * Re_tau^2 * visc_nu^2 / H^3;

%% k-epsilon model
% constants
C_mu = 0.09;
C_1 = 1.44;
C_2 = 1.92;
sigma_k = 1;
sigma_epsilon = 1.3;

%% Equation2
% model initiation
U = 180*ones(N, 1);
k = 180*ones(N, 1);
epsilon = ones(N, 1);

%% Boundary conditions
visc_nu_eddy(1) = 0;
visc_nu_eddy(end) = visc_nu_eddy(end - 1);
k(1) = 0;
k(end) = k(end - 1);
epsilon(1) = 0;
epsilon(end) = epsilon(end - 1);
U(1) = 0;
U(end) = U(end - 1);

% Loop parameters
maxiter = 1000; % Maximum iterations
tolerance = 1e-9; % Convergence tolerance
converged = false;

for i = 1:maxiter
    U_old = U; % Save the previous velocity for convergence check

    dUdy = ddy * U;

    y_Plus = y*U_tau(1)/visc_nu;

    f2 = (1-2/9*exp(-(Re_tau/6).^2)).*(1-exp(-y_Plus/5)).^2;
    f_u = (1 - exp(-y_Plus/ A_Plus)).^2;

    visc_nu_eddy = f_u.* C_mu .* (k.^2) ./ epsilon;

    % Compute Pk
    Pk = visc_nu_eddy .* dUdy.^2;

    visc_term = visc_nu + visc_nu_eddy / sigma_epsilon;

    % Matrix for epsilon
    A_epsilon = bsxfun(@times, visc_term, d2dy2) + bsxfun(@times, ddy * visc_term, ddy);
    A_epsilon(1, :) = 0;
    A_epsilon(1, 1) = 1;
    A_epsilon(end, :) = 0;
    A_epsilon(end, end - 1) = -1;
    A_epsilon(end, end) = 1;

    for index = 2:N-1
        A_epsilon(index, index) = A_epsilon(index, index) - C_2* f2(index) * epsilon(index) / max(k(index), 1e-10); % Avoid division by zero
    end

    B_epsilon = -C_1 * epsilon ./ max(k, 1e-10) .* Pk; % Avoid division by zero
    B_epsilon(1) = visc_nu * (k(1) - 2 * k(2) + k(3)) / dy^2;
    B_epsilon(end) = 0;

    % Solve for epsilon
    epsilon = linsolve(A_epsilon, B_epsilon);

    % Apply bounds to epsilon
    epsilon = max(epsilon, 1e-8); % Ensure epsilon is non-negative

    %disp(['Iteration: ', num2str(i), ' Max Epsilon: ', num2str(max(epsilon))]);

    % Calculate k 
    visc_nu_eddy = f_u .* C_mu .* k.^2 ./ epsilon;
    visc_term = visc_nu + visc_nu_eddy / sigma_k;

    % Matrix for k
    A_k = bsxfun(@times, visc_term, d2dy2) + bsxfun(@times, ddy * visc_term, ddy);
    A_k(1, :) = 0;
    A_k(1, 1) = 1;
    A_k(end, :) = 0;
    A_k(end, end - 1) = -1;
    A_k(end, end) = 1;

    for j = 2:N-1
        A_k(j, j) = A_k(j, j) - epsilon(j) / k(j); % Avoid division by zero
    end

    B_k = -Pk;
    B_k(1) = 0;
    B_k(end) = 0;

    % Solve for k
    k = linsolve(A_k, B_k);

    % Apply bounds to k
    k = max(k, 1e-8); % Ensure k is non-negative

    visc_nu_eddy = f_u.* C_mu .* k.^2 ./ epsilon; % Avoid division by zero
    visc_term = visc_nu + visc_nu_eddy;

    % Matrix for U
    A_U = bsxfun(@times, visc_term, d2dy2) + bsxfun(@times, ddy * visc_term, ddy);
    A_U(1, :) = 0;
    A_U(1, 1) = 1;
    A_U(end, :) = 0;
    A_U(end, end - 1) = -1;
    A_U(end, end) = 1;

    B_U = ones(N, 1) * 1 / rho * DP_dx;
    B_U(1) = 0;
    B_U(end) = 0;

    U = linsolve(A_U, B_U);
    %U = max(30, 1e-8);


    % Check convergence
    residual = norm(U - U_old);
    disp(['Iteration: ', num2str(i), ' Residual: ', num2str(residual)]);

    if residual < tolerance
        converged = true;
        disp('Converged successfully.');
        break; % Exit the loop if converged
    end
end

if ~converged
    disp('Did not converge within maximum iterations.');
end

figure;
semilogx(y*U_tau(1)/visc_nu, U/U_tau(1), 'b-', 'LineWidth', 2); % Plot U in semi-log scale
hold on;
semilogx(DNSData.(['y_plus_', num2str(Re_tau_values(1))]), ...
            DNSData.(['U_mean_', num2str(Re_tau_values(1))]), ...
            'k--', 'DisplayName', ['DNS'], ...
            'LineWidth', 1.5);
xlabel('y+');
ylabel('Velocity (U/U_\tau)');
title('Velocity Profile in Semi-log Format');
grid on;

