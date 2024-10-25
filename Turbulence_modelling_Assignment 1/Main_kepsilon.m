clc;clear all;close all

%% Parameters
visc_nu = 1.56*10^(-5); 
rho = 1.225;
visc_mu = visc_nu*rho;
H = 1;  
N = 1001;
dy = H/(N-1);

%% Constants
A_Plus = 26;
kappa = 0.41;

 
%% DNS
Re_tau_values = [180 395 590]; % Pressure grad
U_tau = Re_tau_values.*visc_nu/H;

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);


%% Solver


MESH = mesh(H, N);

y = MESH.y;
ddy = MESH.ddy;
d2dy2 = MESH.d2dy2;

visc_nu_eddy = 0 * ones(N,1);

Re_tau = Re_tau_values(1);
DP_dx = -Re_tau^2*visc_nu^2/H^3;

%Get laminar velocity
U = ones(N,1);
A = zeros(N,N);
B = ones(N,1);
    
B = ((1./((visc_nu + visc_nu_eddy) .* rho)) .* DP_dx .* dy^2) .* B;
%Boundary conditions
B(1) = 0;
B(end) = 0;
 
A = d2dy2;
A(1, :) = 0;
A(1,1) = 1;
A(end, :)=0;
A(end, end-1) = -1;
A(end, end) = 1;
    
U_laminar = linsolve(A,B);

%% k-epsilon model
%constants
C_mu = 0.09;
C_1 = 1.44;
C_2 = 1.92;
sigma_k = 1;
sigma_epsilon = 1.3;

%model initiation
U_new = zeros(N,1);
U = U_laminar;

k = 0.001*ones(N,1);
epsilon = 0.001*ones(N,1);

tolerance = 10^-6;
error = norm(U - U_new, 2) / length(U);
maxit = 20000;
it =0;

% Avoid division by zero or very small values
k(k < 1e-6) = 1e-6;
epsilon(epsilon < 1e-6) = 1e-6;
 
% Estimate eddy viscosity (C_nu * k^2 / epsilon)
visc_nu_eddy = C_mu .* (k.^2) ./ epsilon;
visc_nu_eddy(1) =0;
visc_nu_eddy(end) = visc_nu_eddy(end-1);


U = U_new;

%Model loop
while error > tolerance && it < 5
    k(k < 1e-6) = 1e-6;
    epsilon(epsilon < 1e-6) = 1e-6;

    dUdy = ddy * U_new;


    % dUdy = zeros(N,1);
    % dUdy(1) = (U_new(2) - U_new(1))/dy;
    % dUdy(end) = 0;
    % for i=2:N-1
    %     dUdy(i) = (U_new(i+1) - U_new(i-1))/(2*dy);
    % end
    % A_epsilon = ones(N);
    % for i = 2:N-1
    %     A_epsilon(i, i-1) = 4*sigma_epsilon*visc_nu + 4*visc_nu_eddy(i) + visc_nu_eddy (i-1) + visc_nu_eddy (i+1);  % Coefficient U(y-dy)
    %     A_epsilon (i, i) = -8*visc_nu*sigma_epsilon - 8*visc_nu_eddy(i);       % Coefficient U(y)
    %     A_epsilon(i, i+1) =  4*sigma_epsilon*visc_nu + 4*visc_nu_eddy(i) + visc_nu_eddy (i-1) + visc_nu_eddy (i+1);    % Coefficient U(y+dy)
    % end
    % 
    % A_epsilon(1, 1) = 1;
    % A_epsilon(end, end-1) = -1;
    % A_epsilon(end, end) = 1;


    % Calculate Pk 
    Pk = visc_nu_eddy .* dUdy^2;

    % Compute epsilon
    f_epsilon = sigma_epsilon * (C_2.* epsilon.^2 ./ k - C_1 * epsilon .* Pk ./k);
    
    B_epsilon = 4*dy^2 .* f_epsilon;
    B_epsilon(1) = visc_nu *k(1) -2*k(2)+k(3)/dy^2;
    B_epsilon(end) = 0;
 
    
    epsilon_old = epsilon;
    epsilon = linsolve(A_epsilon, B_epsilon);

    % Compute k
    A_k = ones(N);
    for i = 2:N-1
        A_k(i, i-1) = 4*sigma_k*visc_nu + 4*visc_nu_eddy(i) + visc_nu_eddy (i-1) + visc_nu_eddy (i+1);  % Coefficient U(y-dy)
        A_k(i, i) = -8*visc_nu*sigma_k - 8*visc_nu_eddy(i);       % Coefficient U(y)
        A_k(i, i+1) =  4*sigma_k*visc_nu + 4*visc_nu_eddy(i) + visc_nu_eddy (i-1) + visc_nu_eddy (i+1);    % Coefficient U(y+dy)
    end
    
    A_k(1, 1) = 1;
    A_k(end, end-1) = -1;
    A_k(end, end) = 1;

    f_k = sigma_k .* (epsilon - Pk);
    B_k = 4*dy^2 .*f_k;
    B_k(1) = 0;
    B_k(end) = 0;
    
    k_old = k;
    k = linsolve(A_k, B_k);

    %Calculate eddy viscosity 

    visc_nu_eddy = C_mu * k.^2 ./epsilon;
    visc_nu_eddy(1) =0;
    visc_nu_eddy(end) = visc_nu_eddy(end-1);

    B = 1/rho .*DP_dx*ones(N,1);
    B(1) = 0;
    B(end) = 0;
       
    for i = 2:N-1
            A(i, i-1) = 1;    % Coefficient U(y-dy)
            A(i, i) = -2;       % Coefficient U(y)
            A(i, i+1) = 1;     % Coefficient U(y+dy)
    end
        
    A(1, 1) = 1;
    A(end, end-1) = -1;
    A(end, end) = 1;
    U_new = linsolve(A, B);

    error = norm(U - U_new, 2) / length(U);
    U = U_new;
    it = it+1;
end

figure(1);
plot(U, y);
grid on;



U_plus = U / U_tau(1);

y_plus = y*U_tau(1)/visc_nu;

figure(2);
semilogx(y_plus, U_plus);
grid on;
hold on;
semilogx(DNSData.(['y_plus_', num2str(Re_tau_values(1))]), ...
            DNSData.(['U_mean_', num2str(Re_tau_values(1))]))

