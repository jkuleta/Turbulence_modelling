function waveConditions = getAllWaveConditions()
    % Define the data from the table
    cases = [1, 2, 3, 4];
    testNumber = [3, 7, 10, 13];            % Test no. from Jensen et al. (1989)
    T = [9.72, 9.72, 9.72, 9.72];           % Period (s)
    U0m = [0.23, 0.68, 2.0, 2.0];           % Maximum velocity (m/s)
    ks = [0, 0, 0, 0.00084];                % Roughness height (m)
    Re = [7.2e4, 6.3e5, 5.4e6, 5.4e6];      % Reynolds number
    a_over_ks = [NaN, NaN, NaN, 3700];      % a/k_s parameter (NaN for the dash)

    % Preallocate structure array
    waveConditions(4) = struct();  % assuming 4 cases
    
    % Loop through each case and assign data
    for i = 1:4
        waveConditions(i).caseNumber = cases(i);
        waveConditions(i).testNumber = testNumber(i);
        waveConditions(i).T = T(i);
        waveConditions(i).U0m = U0m(i);
        waveConditions(i).ks = ks(i);
        waveConditions(i).Re = Re(i);
        waveConditions(i).a_over_ks = a_over_ks(i);
    end
end
