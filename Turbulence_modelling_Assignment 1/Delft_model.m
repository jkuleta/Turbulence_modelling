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


%% k-epsilon model
%constants
C_mu = 0.09;
C_1 = 1.44;
C_2 = 1.92;
sigma_k = 1;
sigma_epsilon = 1.3;

%model initiation
U= zeros(N,1);


k = 0.001*ones(N,1);
epsilon = 0.001*ones(N,1);

tolerance = 10^-6;
error = 1000;
maxit = 20000;
it =0;

% Avoid division by zero or very small values
k(k < 1e-6) = 1e-6;
epsilon(epsilon < 1e-6) = 1e-6;
 
% Estimate eddy viscosity (C_nu * k^2 / epsilon)
visc_nu_eddy = C_mu .* (k.^2) ./ epsilon;
visc_nu_eddy(1) =0;
visc_nu_eddy(end) = visc_nu_eddy(end-1);


%Model loop
for i=1:1000
    fprintf('iteration number: %d.\n',i);
    epsilon_old = epsilon;
    k_old = k;
    U_old = U;
    visc_nu_eddy_oldb = visc_nu_eddy;

    k(k < 1e-6) = 1e-6;
    epsilon(epsilon < 1e-6) = 1e-6;

    dUdy = ddy * U;


    % Calculate Pk 
    Pk = visc_nu_eddy .* dUdy.^2;

    %A-matrix for epsilon

    visc_mu_tot = rho * (visc_nu + visc_nu_eddy/sigma_epsilon);

    A_epsilon =   bsxfun(@times, visc_mu_tot, d2dy2) ... 
        + bsxfun(@times, (ddy*visc_mu_tot), ddy);

    A_epsilon(1,:) = 0;
    A_epsilon(end, :) =0;
    A_epsilon(1,1) = 1;
    A_epsilon(end, end) = 1;
    A_epsilon(end, end-1) = -1;

    for index=2:N-1
        A_epsilon(index,index) = A_epsilon(index,index) - C_2*rho*epsilon(index)/k(index);
    end


    % Compute epsilon
    B_epsilon = - C_1 * epsilon .* Pk ./k;
    
    B_epsilon(1) = visc_nu *k(1) -2*k(2)+k(3)/dy^2;
    B_epsilon(end) = 0;
 
    
    
    epsilon = linsolve(A_epsilon, B_epsilon);

    % Compute k
    %A-matrix for k

    visc_mu_tot = rho*(visc_nu + visc_nu_eddy/sigma_k);

    A_k = bsxfun(@times, visc_mu_tot, d2dy2) ... 
        + bsxfun(@times, (ddy*visc_mu_tot), ddy);
    
    A_k(1,:) = 0;
    A_k(end, :) =0;
    A_k(1,1) = 1;
    A_k(end, end) = 1;
    A_k(end, end-1) = -1;

    for index=2:N-1
        A_k(index,index) = A_k(index,index) - rho.*epsilon(index)./k(index);
    end

    B_k = - Pk;
    B_k(1) = 0;
    B_k(end) = 0;

    k = linsolve(A_k, B_k);
    

    %Calculate eddy viscosity 

    visc_nu_eddy = C_mu * k.^2 ./epsilon;
    visc_nu_eddy(1) =0;
    visc_nu_eddy(end) = visc_nu_eddy(end-1);
    B = 1/rho.*(1./(visc_nu+visc_nu_eddy)) .*DP_dx.*ones(N,1);
    B(1) = 0;
    B(end) = 0;
    
    A = zeros(N);
    for i = 2:N-1
            A(i, i-1) = 1;    % Coefficient U(y-dy)
            A(i, i) = -2;       % Coefficient U(y)
            A(i, i+1) = 1;     % Coefficient U(y+dy)
    end
        
    A(1, 1) = 1;
    A(end, end-1) = -1;
    A(end, end) = 1;
    U = linsolve(A, B);

    error = norm(U_old - U, 2) / length(U);
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

