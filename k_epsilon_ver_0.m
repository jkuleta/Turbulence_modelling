clc;clear all;close all

%% Parameters
visc_nu = 1.56*10^(-5); 
rho = 1.225;
visc_mu = visc_nu*rho;
H = 1;  
N = 101;
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


MESH = mesh_Delft(H, N);

y = MESH.y;
ddy = MESH.ddy;
d2dy2 = MESH.d2dy2;


visc_nu_eddy = zeros(N,1);
Re_tau = Re_tau_values(1);
DP_dx = -rho*Re_tau^2*visc_nu^2/H^3;

%% k-epsilon model
%constants
C_mu = 0.09;
C_1 = 1.44;
C_2 = 1.92;
sigma_k = 1;
sigma_epsilon = 1.3;

%model initiation
U= zeros(N,1);


k = 100*ones(N,1);
epsilon = 100*ones(N,1);
epsilon_min = 1e-12;
k_min = 1e-12;

visc_nu_eddy = C_mu .* (k.^2) ./ epsilon;


%loop parameters
maxiter = 1000;



%% K-epsilon model

for i=1:maxiter

    dUdy = ddy*U;

    visc_nu_eddy = C_mu .* (k.^2) ./ epsilon;

    %compute Pk
    Pk = visc_nu_eddy .* dUdy.^2;

    %solve for epsilon
    visc = visc_nu + visc_nu_eddy / sigma_epsilon;
    A_epsilon = bsxfun(@times, visc, d2dy2) + bsxfun(@times,ddy*visc, ddy);

    A_epsilon(1, :) = 0;
    A_epsilon(1,1) = 1;
    A_epsilon(end, :) =0;
    A_epsilon(end, end-1) = -1;
    A_epsilon(end, end) = 1;

    for index=2:N-1
        A_epsilon(index,index) = A_epsilon(index,index) - C_2*epsilon(index)/k(index);
    end

    B_epsilon = - C_1* epsilon ./ k .* Pk;
    B_epsilon(1) = visc_nu *(k(1) -2*k(2)+k(3))/dy^2;
    B_epsilon(end)=0;
    
    epsilon = linsolve(A_epsilon, B_epsilon);

    epsilon(epsilon < epsilon_min) = epsilon_min;
    
    %calculate k 
    
    visc_nu_eddy = C_mu .* (k.^2) ./ epsilon;

    visc = visc_nu + visc_nu_eddy / sigma_k;

    A_k = bsxfun(@times, visc, d2dy2) + bsxfun(@times, ddy*visc, ddy);
    A_k(1, :) = 0;
    A_k(1,1) = 1;
    A_k(end, :) =0;
    A_k(end, end-1) = -1;
    A_k(end, end) = 1;

    for j=2:N-1
        A_k(j,j) = A_k(j,j) - epsilon(j) / k(j);
    end
    B_k = -Pk;
    B_k(1) = 0;
    B_k(end) = 0;

    k = linsolve(A_k, B_k);
    k(k < k_min) = k_min;
    
    %update eddy viscosity
    visc_nu_eddy = C_mu .* (k.^2) ./ epsilon;
    
    %solve rans
    viscosity = visc_nu + visc_nu_eddy;
    A_U = bsxfun(@times, viscosity, d2dy2) + bsxfun(@times, ddy * viscosity, ddy);

    A_U(1, :) = 0;
    A_U(1,1) = 1;
    A_U(end, :) =0;
    A_U(end, end-1) = -1;
    A_U(end, end) = 1;

    B_U = 1/rho * DP_dx * ones(N, 1);
    B_U(1) = 0;
    B_U(end) = 0;

    U_old = U;

    U = linsolve(A_U, B_U);


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
    


