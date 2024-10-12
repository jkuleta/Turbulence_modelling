clc;clear all;close all

%% Parameters
visc_nu = 1.56*10^(-5); 
rho = 1.225;
visc_mu = visc_nu*rho;
H = 1;  

%% Constants
A_Plus = 26;
kappa = 0.41;

 
%% DNS
Re_tau_values = [180 395 590]; % Pressure grad
U_tau = Re_tau_values.*visc_nu/H;

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);


%% Solver
dy = 0.001; 

height = input('Enter the height condition (full/half): ', 's');
neumann_scheme = input('Enter the Neumann scheme (1/2): ');

if strcmp(height, 'full')
y = 0:dy:H;  
elseif strcmp(height, 'half')
y = 0:dy:H;   
end
            
N = length(y); 

%% Forward Euler

U_Turbulent_undamped = zeros(N, length(Re_tau_values));
U_Turbulent_damped = zeros(N, length(Re_tau_values));
U_Plus_undamped = zeros(N, length(Re_tau_values));
U_Plus_damped = zeros(N, length(Re_tau_values));


%% figures
for j = 1:length(Re_tau_values)
    %% Solve for the current Re_tau
    [U_Turbulent_undamped(:, j), U_Turbulent_damped(:, j)] = ...
        solveTurbulence(Re_tau_values(j), height, neumann_scheme, visc_nu, H, kappa, A_Plus, N, dy, y);


    y_Plus_Turbulent(:, j) = y * U_tau(j) / visc_nu;
    U_Plus_undamped(:, j) = U_Turbulent_undamped(:, j) / U_tau(j);
    U_Plus_damped(:, j) = U_Turbulent_damped(:, j) / U_tau(j);

    %% Near wall plot
%plotReynoldsStressNearWall(DNSData, Re_tau_values(j), visc_mu, U_Turbulent_damped(:, j), U_Turbulent_undamped(:, j), y, kappa, A_Plus, y_Plus_Turbulent(:, j));


    %% Reynolds stress plot
    %plotReynoldsStress(DNSData, Re_tau_values(j), visc_mu, U_Turbulent_damped(:, j), U_Turbulent_undamped(:, j), dy, kappa, A_Plus, y_Plus_Turbulent(:, j),U_tau(j),height);

end

%% Velocity plot
plotVelocityProfiles(Re_tau_values, DNSData, U_Plus_undamped, U_Plus_damped, y_Plus_Turbulent, H)


%% Kinetic energy plot
%plotKineticEnergy(Re_tau_values, DNSData);




