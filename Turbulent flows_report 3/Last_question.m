clc; clear; close all;

% Constants
nu_visc = 1.51e-5; % Kinematic viscosity (m^2/s)
D = 0.03; % Example characteristic length (m), adjust as needed
num_tests = 12; % Number of tests (adjust if different)

% Initialize arrays for storing results
Re = zeros(1, num_tests);
Lambda_f = zeros(1, num_tests);
lambda_f = zeros(1, num_tests);
eta_K = zeros(1, num_tests);
L = zeros(1, num_tests);
var_u = zeros(1, num_tests);

% Load data and compute quantities for each test
load("Exercise3.mat"); % Ensure correct file path
for i = 1:num_tests
    data = Jet(i);
    
    % Extract fluctuations
    fluc_u = data.u - mean(data.u);
    u_bar = mean(data.u); % Mean velocity
    var_u(i) = var(data.u);

    % Time parameters
    dt = data.t(2) - data.t(1);
    
    % Turbulent integral scale
    n = 1; 
    R_E(n) = mean(fluc_u.^2) / var_u(i);
    while R_E(n) > 0
        n = n + 1;
        R_E(n) = mean(fluc_u(1:end-n+1) .* fluc_u(n:end)) / var_u(i);
    end
    
    % Ensure R_E and data.t(1:n) have the same length
    R_E = R_E(:);
    T_E = trapz(data.t(1:length(R_E)), R_E); % Integral timescale
    
    % Taylor micro length scale
    tau_E = sqrt(2 * var(fluc_u) / mean((diff(fluc_u) / dt).^2));
    lambda_f(i) = u_bar * tau_E;
    
    % Taylor macro length scale
    Lambda_f(i) = u_bar * T_E;
    
    % Kolmogorov scales
    epsilon = 30 * nu_visc * var(fluc_u) / (lambda_f(i)^2); % Energy dissipation rate
    eta_K(i) = (nu_visc^3 / epsilon)^(1/4); % Kolmogorov length scale
    
    % Reynolds number
    Re(i) = u_bar * D / nu_visc; % Reynolds number

    L(i) = var_u(i)^(3/2) / epsilon;
end

% Compute scale ratios
ratio_Lambda_eta = Lambda_f ./ eta_K; % Ratio of Lambda_f to eta_K
ratio_Lambda_lambda = Lambda_f ./ lambda_f; % Ratio of Lambda_f to lambda_f

% Log-log plot of ratios vs Reynolds number
figure;
loglog(Re, ratio_Lambda_eta, 'o-', 'LineWidth', 1.5, 'DisplayName', '\Lambda_f / \eta_K');
hold on;
loglog(Re, ratio_Lambda_lambda, 'o-', 'LineWidth', 1.5, 'DisplayName', '\Lambda_f / \lambda_f');

Re_t = L .* sqrt(var_u) / nu_visc;

% Add trendlines
%Re_trend = logspace(log10(min(Re)), log10(max(Re)), 12); % Generate Re values for trendline
t_eta = (Lambda_f ./ L) .* Re_t.^(3/4); % Trendline for \Lambda_f / \eta_K ~ Re^(3/4)
t_lambda = (Lambda_f ./ L) .* Re_t.^(1/2) ./ sqrt(30); % Trendline for \Lambda_f / \lambda_f ~ Re^(1/2)

loglog(Re, t_eta, '--', 'LineWidth', 1.5, 'DisplayName', 'Re^{3/4}');
loglog(Re, t_lambda, '--', 'LineWidth', 1.5, 'DisplayName', 'Re^{1/2}');

% Add trendlines
Re_trend = logspace(log10(min(Re)), log10(max(Re)), 100); % Generate Re values for trendline
trend_eta = Re_trend.^(3/4); % Trendline for \Lambda_f / \eta_K ~ Re^(3/4)
trend_lambda = Re_trend.^(1/2) ./ sqrt(30); % Trendline for \Lambda_f / \lambda_f ~ Re^(1/2)

loglog(Re_trend, 0.5*trend_eta, '--', 'LineWidth', 1.5, 'DisplayName', '~Re^{3/4}');
loglog(Re_trend, trend_lambda, '--', 'LineWidth', 1.5, 'DisplayName', '~Re^{1/2}');

% Annotate the plot
grid on;
xlabel('Reynolds number (Re)');
ylabel('Length scale ratios');
legend('show', 'Location', 'northwest');
title('Length Scale Ratios vs Reynolds Number');
hold off;
