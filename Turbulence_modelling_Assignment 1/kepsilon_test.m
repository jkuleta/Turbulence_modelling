N = 500;               % Number of grid points
H = 1;                 % Total height of half channel
dy = H / (N-1);        % Grid spacing
nu = 1.56*10^(-5);     % Kinematic viscosity
rho = 1.225;           % Density
Cmu = 0.09;            % Model constant for turbulent viscosity
sigma_k = 1.0;         % Model constant for k
sigma_eps = 1.3;       % Model constant for epsilon
Ce1 = 1.44;            % Dissipation constant C_e1
Ce2 = 1.92;             % Dissipation constant C_e2
A_plus = 100;

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);

% Reynolds number (Re_tau)
Re_tau = 180;  

% Calculate friction velocity u_tau from Re_tau
u_tau = Re_tau * nu / H;

% Calculate the pressure gradient (dP/dx) required to achieve target Re_tau
dPdx = -Re_tau^2*nu^2/H^3;

% Lower bounds for stability
k_min = 1e-12;          % Minimum value for k
eps_min = 1e-12;        % Minimum value for epsilon
nu_T_max = 1;          % Maximum allowed turbulent viscosity for stability

% Initialize variables
U = linspace(0, u_tau, N)';        % Velocity profile
k = ones(N,1) * (0.03 * u_tau^2);  % Small initial value for k (Turbulent kinetic energy)
eps = ones(N,1) * (0.03 * u_tau^3 / H); % Small initial value for epsilon
nu_T = zeros(N,1);     % Turbulent viscosity
Pi = zeros(N,1);       % Production term for turbulence

% Iteration parameters
maxIter = 10000000;  % Maximum number of iterations
tol = 1e-12;  % Tolerance for convergence

% Calculate y_plus
y = linspace(0, H, N);
y_vert = y';
y_plus = y_vert * u_tau / nu;

% Main iterative solver with relaxation
for iter = 1:maxIter
    U_old = U;   % Store the old velocity for convergence check

    % Apply Van Driest correction
    %f_mu = 1;
    f_mu = (1 - exp(-y_plus / A_plus)).^2;

    % Update turbulent viscosity nu_T
    nu_T = Cmu * f_mu .* (k.^2) ./ eps;

        % Update Production term Pi
    for i = 2:N-1
        Pi(i) = nu_T(i) * ((U(i+1) - U(i-1)) / (2*dy))^2;
        
        % Update velocity U(i) 
        U(i) = ( -(2*nu + nu_T(i+1) + nu_T(i)) * U(i+1) - ...
                   (2*nu + nu_T(i-1) + nu_T(i)) * U(i-1) + ...
                   (2 * dy^2 * (1/rho) * dPdx)) / ...
                   (-4*nu - nu_T(i+1) - 2*nu_T(i) - nu_T(i-1));
       
        % Update k(i) 
        k(i) = ( -(2*nu + nu_T(i+1) + nu_T(i)) * k(i+1) - ...
                   (2*nu + nu_T(i-1) + nu_T(i)) * k(i-1) + ...
                   (2 * dy^2 * sigma_k * (eps(i) - Pi(i)))) / ...
                   (-4*nu - nu_T(i+1) - 2*nu_T(i) - nu_T(i-1));
        
        % Ensure k does not go below minimum
        k(i) = max(k(i), k_min);

        % Update eps(i) 
        eps(i) = ( -(2*nu + nu_T(i+1) + nu_T(i)) * eps(i+1) - ...
                     (2*nu + nu_T(i-1) + nu_T(i)) * eps(i-1) + ...
                     (2 * dy^2 * sigma_eps * (Ce2 * (eps(i)^2 / k(i)) - Ce1 * (eps(i) / k(i)) * Pi(i)))) / ...
                     (-4*nu - nu_T(i+1) - 2*nu_T(i) - nu_T(i-1));

        % Ensure epsilon does not go below minimum
        eps(i) = max(eps(i), eps_min);  
    end

    % Enforce boundary conditions at the wall and center
    U(1) = 0;   % No-slip condition at the wall
    k(1) = 0;   % k = 0 at the wall
    eps(1) = 2 * nu * (k(1) -2*k(2)+k(3))/dy^2;
    %eps(1) = 1e-6;
    nu_T(1) = 0;

    U(N) = U(N-1);   % Symmetry at the center for U
    k(N) = k(N-1);   % Symmetry at the center for k
    eps(N) = eps(N-1);  % Symmetry at the center for epsilon
    nu_T(N) = nu_T(N-1);

    % Convergence check
    if max(abs(U - U_old)) < tol
        fprintf('Convergence achieved after %d iterations\n', iter);
        break;
    end
end


% Display the target Re_tau
fprintf('Re_tau: %f\n', Re_tau);

U_plus = U /u_tau;

figure(1);
semilogx(y_plus, U_plus);
grid on;
hold on;
semilogx(DNSData.(['y_plus_', num2str(Re_tau)]), ...
            DNSData.(['U_mean_', num2str(Re_tau)]))

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


% Calculate k_plus and epsilon_plus
k_plus = k * H / (u_tau^2);
epsilon_plus = eps * H * nu / (u_tau^4);

% Create a new figure for k+ and epsilon+ plots
figure;

% Plot k+
subplot(2, 1, 1);
plot(y_plus, k_plus, 'm');
xlabel('k^+');
ylabel('y^+');
title('Turbulent Kinetic Energy (k^+) vs y^+');
grid on;

% Plot epsilon+
subplot(2, 1, 2);
plot(y_plus, epsilon_plus, 'c');
xlabel('y^+');
ylabel('\epsilon^+');
title('Dissipation Rate (\epsilon^+) vs y^+');
grid on;

% Adjust figure layout for better visibility
sgtitle('Turbulent Kinetic Energy and Dissipation Rate in Wall Units');
