clc;clear all;close all

% Parameters
visc = 1.56*10^(-5); 
rho = 1.225;
H = 1;  
kappa = 0.41;
Re_tau = 180; % Pressure grad
U_tau = Re_tau*visc/H;
dy = 0.001;    % step

filename_180 = 'DNS_180.txt'; 
Data_DNS_180 = readmatrix(filename_180); 

height = input('Enter the height condition (full/half): ', 's');
neumann_scheme = input('Enter the Neumann scheme (1/2): ');

if strcmp(height, 'full')
y = 0:dy:2*H;  
elseif strcmp(height, 'half')
y = 0:dy:H;   
end
            
N = length(y);   

DP_dx = -Re_tau^2*visc^2/H^3;

%% Turbulent 
A_Turbulent = zeros(N, N);   
%B_Turbulent = ones(N, 1) * (-visc + sqrt(visc^2 - 4 * 1/rho * DP_dx * (H - y) .* (kappa * y).^2)) / (2 * (kappa * y).^2) .* dy;


inside_sqrt = visc^2 - 4 * DP_dx * (H - y) .* (kappa * y).^2;
B_Turbulent = [(-visc + sqrt(inside_sqrt)) ./ (2 * (kappa * y).^2) .* dy]';

%{
sqrt_term = sqrt((4 * DP_dx * kappa^2 * (H - y) .* y.^2)  + visc^2);

% Calculate the numerator
numerator =  visc * sqrt_term + DP_dx * kappa^2 * y.^3 - 2 * H * DP_dx * kappa^2 * y.^2 - visc^2;

% Calculate the denominator
denominator = kappa^2  * y.^3 .* sqrt_term;

% Final result
result = numerator / denominator;

B_Turbulent = ones(N, 1) * result * dy^2;
%}
for i = 2:N-1
    %A_Turbulent(i,i-1) = 0;    % coeff U(y-dy)
    A_Turbulent(i,i) = -1;     % Cohaleff U(y)
    A_Turbulent(i,i+1) = 1;    % Coeff U(y+dy)
end

% Boundary conditions
B_Turbulent(1) = 0;            % Initial condition
B_Turbulent(end) = 0;          % Boundary condition


if strcmp(height, 'full')
    A_Turbulent = diag(-2 * ones(N, 1)) + diag(ones(N-1, 1), 1) + diag(ones(N-1, 1), -1);
    A_Turbulent(1, 1) = 1;
    A_Turbulent(1, 2) = 0;
    A_Turbulent(end, end-1) = 0;
    A_Turbulent(end, end) = 1;
elseif strcmp(height, 'half')
    A_Turbulent = diag(-2 * ones(N, 1)) + diag(ones(N-1, 1), 1) + diag(ones(N-1, 1), -1);
    A_Turbulent(1, 1) = 1;
    A_Turbulent(1, 2) = 0;
    if neumann_scheme == 1
        A_Turbulent(end, end-1) = -1;
        A_Turbulent(end, end) = 1;
    elseif neumann_scheme == 2
        A_Turbulent(end, end-2) = -0.5;
        A_Turbulent(end, end-1) = 2;
        A_Turbulent(end, end) = -1.5;
    else
        error('Invalid Neumann scheme');
    end
else
    error('Invalid height condition. Choose either "full" or "half".');
end

% Solve
U_Turbulent = inv(A_Turbulent)*B_Turbulent;


y_Plus_Turbulent = y * U_tau/visc;
U_Plus_Turbulent = U_Turbulent/U_tau;

y_Plus_DNS_180 = Data_DNS_180(:,2);
U_Plus_DNS_180 = Data_DNS_180(:,3);

figure(2);
semilogx(y_Plus_Turbulent, U_Plus_Turbulent);
hold on;
semilogx(y_Plus_DNS_180, U_Plus_DNS_180);
title('Poiseuille Flow Laminar vs Turbulent');
grid on;
xlabel('y^+');
ylabel('U^+');
legend('Turbulent');
hold off;


