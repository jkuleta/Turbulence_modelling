clc; clear; close all;

%% case info
Cases = getAllWaveConditions();
nu = 1.14e-6;

%Cases
% 1 - Laminar 
% 2 - Transition 
% 3 - Turbulent 
% 4 - Turbulent (rough)

%% Question 1
for i = 1:4
    %% Friction priori 

    % Rough turbulent Case 4
    if Cases(i).ks > 0
        f_w = exp(5.5 * (Cases(i).a_over_ks^(-0.16)) - 6.7); % eqn 5.69

    % Laminar flow
    % Re <= 1.5 ^10^5 
    elseif Cases(i).Re <= 1.5e5
        f_w = 2 / sqrt(Cases(i).Re); % eqn 5.59
        
    % Turbulent flow
    % Re >= 5 ^10^5
    elseif Cases(i).Re >= 5e5
        f_w = 0.035 / (Cases(i).Re^(0.16)); % eqn 5.60
        
    % Transitional flow
    % 1.5 ^10^5 <= Re <= 5 ^10^5 
    elseif Cases(i).Re > 1.5e5 && Cases(i).Re < 5e5 % eqn 5.61
        f_w = 0.004 - 0.005; 
    end

    %% maximum velocity
    U_fm = sqrt(f_w / 2) * Cases(i).U0m; % eqn 5.57

   % update struct
    Cases(i).f_w = f_w;
    Cases(i).U_fm = U_fm;
end

%% Question 2
for i = 1:4

    % Rough turbulent Case 4
    if Cases(i).caseNumber == 4
        % delta y
        delta_y = Cases(i).ks*0.02; % eqn 9.54
        % k_plus
        ks = Cases(i).ks;
        ks_plus = ks*Cases(i).U_fm/nu; %eqn 9.24 
        % Ks_plus is more than 70 so it is hyraulicaly rough
    else
        delta_y = nu/Cases(i).U_fm; % eqn 9.49
        ks = 0.1* nu/Cases(i).U_fm; %eqn 9.24
        ks_plus = ks*Cases(i).U_fm/nu; %eqn 9.24 % should be 0.1
    end

   % update struct
    Cases(i).delta_y = delta_y;
    Cases(i).ks = ks;
    Cases(i).ks_plus = ks_plus;
end
 
CaseInfo = Cases;
save('CaseInfo.mat', 'CaseInfo');


