clc; clear; close all;

%% Load case info
load('CaseInfo.mat')
T = CaseInfo(1).T;
nu = 1.14e-6;
wt = [0 45 90 135]; 

%% Load RANS model
load('out_MatRANS_case1.mat') % Make sure its correct

%% Laminar theory solution 

t = 0:0.1:1; % idk about this check
y = 0:0.1:1;


figure;
hold on;
for i = 1:4
    % input velocity
    U_0 = CaseInfo(i).U0m *sin(2*pi/T*t); % eqn 10.20

    % Boundary condtions:
    % at y= inf u = U_0
    % at y = 0 u = 0;

    omega = 2*pi/T;
    a = CaseInfo(i).U0m/omega; % eqn 5.16
    delta = a * 3*pi/4*(2/CaseInfo(i).Re)^(0.5); % eqn 5.15
    delta_1 = 4/(3*pi)*delta; %  eqn 5.14

    u_theory = zeros(size(y)); 
    for j = 1:length(y)
        u_theory(j) = CaseInfo(i).U0m * sin(deg2rad(wt(i))) - CaseInfo(i).U0m * exp(-y(j) / delta_1) * sin(deg2rad(wt(i)) - y(j) / delta_1); % eqn 5.12
    end
    
    plot( u_theory/CaseInfo(i).U0m,y/a, 'DisplayName', ['\omega t = ' num2str(wt(i)) '°']);
end

ylabel('y/a');
xlabel('u/U0m');
%title('Velocity Profile for Laminar Flow at Different Phases');
legend show;
grid on;
hold off;