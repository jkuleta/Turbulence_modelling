               % Number of grid points
H = 1;                 % Total height of the channel
L_half = H / 2;        % Half-channel height (used as characteristic length)
dy = 0.0001;        % Grid spacing
nu = 1e-6;             % Kinematic viscosity
rho = 1.225;           % Density
Cmu = 0.09;            % Model constant for turbulent viscosity
sigma_k = 1.0;         % Model constant for k
sigma_eps = 1.3;       % Model constant for epsilon
Ce1 = 1.45;            % Dissipation constant C_e1
Ce2 = 1.9;             % Dissipation constant C_e2

y = 0: dy :H;
N = length(y)

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);


% Target Reynolds number (Re_tau)
Re_tau = 180;  % Example: Re_tau = 180

% Calculate friction velocity u_tau from Re_tau
u_tau_target = Re_tau * nu / L_half;

% Calculate wall shear stress (tau_w)
%tau_w = rho * u_tau_target^2;

% Calculate the pressure gradient (dP/dx) required to achieve target Re_tau
%dPdx = 2 * tau_w / H;
dPdx = -Re_tau^2*nu^2/L_half^3;

% Lower bounds for stability
k_min = 1e-8;          % Minimum value for k
eps_min = 1e-8;        % Minimum value for epsilon
nu_T_max = 1;          % Maximum allowed turbulent viscosity for stability

% Initialize variables
U = zeros(N,1);        % Velocity profile
k = ones(N,1) * 1e-3;  % Small initial value for k (Turbulent kinetic energy)
eps = ones(N,1) * 1e-3; % Small initial value for epsilon
nu_T = zeros(N,1);     % Turbulent viscosity
Pi = zeros(N,1);       % Production term for turbulence

% Boundary conditions
U(1) = 0;              % No-slip condition at the wall (velocity U = 0)
k(1) = 0;              % k = 0 at the wall
eps(1) = (k(1) -2*k(2)+k(3))/dy^2; % Wall dissipation rate condition based on viscosity

% Boundary conditions at the center
U(N) = U(N-1);         % Symmetry at the center
k(N) = k(N-1);         % Symmetry for k at the center
eps(N) = eps(N-1);     % Symmetry for epsilon at the center

% Iteration parameters
maxIter = 10000;         % Maximum number of iterations for convergence
tol = 1e-6;            % Tolerance for convergence

% Main iterative solver
for iter = 1:maxIter
    U_old = U;
    k_old = k;
    eps_old = eps;
    
    % Update Turbulent viscosity
    for i = 2:N-1
        nu_T(i) = Cmu * (k(i)^2) / eps(i);  % nu_T = Cmu * (k^2 / epsilon)
        nu_T(i) = min(nu_T(i), nu_T_max);   % Limit turbulent viscosity
    end
    
    % Update Production term Pi (Pi = nu_T * (dU/dy)^2)
    for i = 2:N-1
        Pi(i) = nu_T(i) * ((U(i+1) - U(i-1)) / (2*dy))^2;
    end
    
    % Solve Momentum equation (Finite Difference)
    for i = 2:N-1
        U(i) = ( -(2*nu + nu_T(i+1) + nu_T(i)) * U(i+1) - ...
                 (2*nu + nu_T(i-1) + nu_T(i)) * U(i-1) + ...
                 (2 * dy^2 * (1/rho) * dPdx)) / ...
                 (-4*nu - nu_T(i+1) - 2*nu_T(i) - nu_T(i-1));
    end
    
    % Solve k equation (Finite Difference)
    for i = 2:N-1
        k(i) = ( -(2*nu + nu_T(i+1) + nu_T(i)) * k(i+1) - ...
                 (2*nu + nu_T(i-1) + nu_T(i)) * k(i-1) + ...
                 (2 * dy^2 * sigma_k * (eps(i) - Pi(i)))) / ...
                 (-4*nu - nu_T(i+1) - 2*nu_T(i) - nu_T(i-1));
        k(i) = max(k(i), k_min);  % Ensure k is not negative
    end
    
    % Solve epsilon equation (Finite Difference)
    for i = 2:N-1
        eps(i) = ( -(2*nu + nu_T(i+1) + nu_T(i)) * eps(i+1) - ...
                     (2*nu + nu_T(i-1) + nu_T(i)) * eps(i-1) + ...
                     (2 * dy^2 * sigma_eps * (Ce2 * (eps(i)^2 / k(i)) - Ce1 * (eps(i) / k(i)) * Pi(i)))) / ...
                     (-4*nu - nu_T(i+1) - 2*nu_T(i) - nu_T(i-1));
        eps(i) = max(eps(i), eps_min);  % Ensure eps is not negative
    end
    
    % Enforce boundary conditions at the wall
    U(1) = 0;            % No-slip condition at the wall (U = 0)
    k(1) = 0;            % k = 0 at the wall
    eps(1) = (k(1) -2*k(2)+k(3))/dy^2;  % Wall dissipation rate condition based on viscosity
    
    % Enforce boundary conditions at the center
    U(N) = U(N-1);       % Symmetry at the center
    k(N) = k(N-1);       % Symmetry at the center for k
    eps(N) = eps(N-1);   % Symmetry at the center for epsilon
    
    % Convergence check
    if max(abs(U - U_old)) < tol && max(abs(k - k_old)) < tol && max(abs(eps - eps_old)) < tol
        fprintf('Convergence achieved after %d iterations\n', iter);
        break;
    end
end

% Display the target Re_tau
fprintf('Target Re_tau: %f\n', Re_tau);

% Plot results
y = linspace(0, H, N);  % Discretized y-axis

U_plus = U /u_tau_target;

y_plus = y*u_tau_target/nu;

figure(1);
semilogx(y_plus, U_plus);
grid on;
hold on;
%semilogx(DNSData.(['y_plus_', num2str(Re_tau)]), ...
            %DNSData.(['U_mean_', num2str(Re_tau)]))
% Plot the velocity profile U
figure;
subplot(3,1,1);
plot(U, y, 'r');
xlabel('U');
ylabel('y');
title('Velocity Profile (U)');

% Plot the turbulent kinetic energy (k)
subplot(3,1,2);
plot(k, y, 'g');
xlabel('k');
ylabel('y');
title('Turbulent Kinetic Energy (k)');

% Plot the dissipation rate (epsilon)
subplot(3,1,3);
plot(eps, y, 'b');
xlabel('\epsilon');
ylabel('y');
title('Dissipation Rate (\epsilon)');