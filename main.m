clc;clear all;close all

%% Parameters
visc_nu = 1.56*10^(-5); 
rho = 1.225;
visc_mu = visc_nu*rho;
H = 1;  

%% DNS
Re_tau_values = [180 395 590]; 
U_tau = Re_tau_values.*visc_nu/H;

filenames = {'DNS_180.txt', 'DNS_395.txt', 'DNS_590.txt'};
DNSData = LoadDNSData(filenames);


%% Loop
for j = 1:length(Re_tau_values)

    %% Prandalt
    U_Prandalt = PrandaltSolver(Re_tau_values, visc_nu, H);

    %% K-Epsilon
    Re_tau = Re_tau_values(j);
    [U_kepsilon_damped(j), U_kepsilon_undamped(j)] = kepsilon(Re_tau, visc_nu, rho, H);

end

%% Velocity plot
%plotVelocityProfiles(Re_tau_values, DNSData, U_Prandalt,U_kepsilon_damped,U_kepsilon_undamped, H,visc_nu);
%% Stress Plot
%plotReynoldsStress(Re_tau_values,DNSData,visc_mu, U_Prandalt, U_kepsilon_damped, U_kepsilon_undamped,H);
%% Model comp
%plotVelocityComparison(Re_tau_values, DNSData, U_Prandalt, U_kepsilon_damped, H, visc_nu);



