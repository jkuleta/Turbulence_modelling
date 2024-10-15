function DNSData = LoadDNSData(filenames)
    % Initialize a structure to hold the output data
    DNSData = struct();
    
    % Define Re_tau values corresponding to the input files
    Re_tau_values = [180, 395, 590];
    
    % Loop through each file and read the data
    for i = 1:length(filenames)
        data = readmatrix(filenames{i});

        % Create field names dynamically for each Re_tau value
        DNSData.(['y_', num2str(Re_tau_values(i))]) = data(:, 1);
        DNSData.(['y_plus_', num2str(Re_tau_values(i))]) = data(:, 2);   % Column 1 = y^+
        DNSData.(['U_mean_', num2str(Re_tau_values(i))]) = data(:, 3);  
        DNSData.(['DU_dy_', num2str(Re_tau_values(i))]) = data(:, 4); 
        DNSData.(['R_uu_', num2str(Re_tau_values(i))]) = data(:, 5);     % Column 3 = R_uu
        DNSData.(['R_vv_', num2str(Re_tau_values(i))]) = data(:, 6);     % Column 4 = R_vv
        DNSData.(['R_ww_', num2str(Re_tau_values(i))]) = data(:, 7);     % Column 5 = R_ww
        DNSData.(['R_uv_', num2str(Re_tau_values(i))]) = data(:, 8);     % Column 6 = R_uv
        DNSData.(['R_uw_', num2str(Re_tau_values(i))]) = data(:, 9);     % Column 7 = R_uw
        DNSData.(['R_vw_', num2str(Re_tau_values(i))]) = data(:, 10);     % Column 8 = R_vw
    end
end
