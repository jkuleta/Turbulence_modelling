clc; clear; close all;

nu_visc = 1.51 * 1e-5;
f_sample = 50000;
load("Exercise3.mat");

data = Jet(12);

fluc_u = data.u - mean(data.u);

figure(1);
plot(data.t, fluc_u);
grid on;
xlabel("Time");
ylabel("Fluctuation");

%% Compute statistics
u_bar = mean(data.u);
sigma_u = std(data.u);
var_u = var(data.u);
S_u = skewness(data.u);
K_u = kurtosis(data.u);
turb_intesity = sigma_u /u_bar;

%% Compare with normal distribution
[pdf_u,x] = ksdensity(fluc_u);
pn = 1/sigma_u/sqrt(2*pi)*exp(-fluc_u.^2/(2*var_u));

figure(2);
plot(x, pdf_u, 'LineWidth', 2);
hold on;
grid on;
plot(fluc_u, pn, '.');

%Time correlation
n=1; 
R_E(n) = mean(fluc_u.^2)/var_u;
while R_E(n)>0
    n=n+1;
    R_E(n) = mean(fluc_u(1:end-n+1).*fluc_u(n:end))/var_u;
end

dt = data.t(2) - data.t(1);


T_E = trapz(data.t(1:n),R_E);
%These two calculations are equivalent according to Furhmann.
R_tt = 2*(R_E(2) - R_E(1))/dt^2;
tau_E2 = 1/sqrt(-0.5*R_tt);
tau_E = sqrt(2*var_u/mean((diff(fluc_u)/dt).^2));

figure(3);
plot(data.t(1:n),R_E);
ylim([0, inf]);
hold on;
grid on;
plot(data.t(1:n), 1-data.t(1:n).^2/tau_E^2)
ylim([0, inf]);


%% Taylor macro length scale
Lambda_f = u_bar * T_E;
lambda_f = u_bar * tau_E;

%% Kolmogorov scale
epsilon = 30 * nu_visc * var_u / (lambda_f^2);
tau_K = sqrt(nu_visc / epsilon);
eta_K = (nu_visc^3/epsilon)^(1/4);

%%% Parameters (ensure these are defined in your script)
N = length(fluc_u);         % Length of the fluctuating signal
f_sample = 1/dt;            % Sampling frequency
f_N = f_sample / 2;         % Nyquist frequency
df = f_sample / N;          % Frequency resolution

%% FFT Analysis
fft_values = fft(fluc_u);    % Compute FFT of the signal
Pxx_fft = (1/(f_sample * N)) * abs(fft_values).^2; % Normalize FFT to match PSD
Pxx_fft(2:end-1) = 2 * Pxx_fft(2:end-1);  % Double values except for DC and Nyquist
fft_freq = (0:(N/2)) * df;  % Frequency vector for FFT (0 to Nyquist)

% Welch's Method
% Use pwelch to calculate PSD
[welch_values, welch_freq] = pwelch(fluc_u, [], [], [], f_sample);

integral_signal = trapz(fft_freq(1:(N/2)),Pxx_fft(1:(N/2)));

integral_welch = trapz(welch_freq, welch_values);

% Plot the spectra
figure(4);
loglog(fft_freq, Pxx_fft(1:N/2+1), 'LineWidth',1.25); % Plot FFT power spectrum
hold on;
grid on;
loglog(welch_freq, welch_values, 'LineWidth',1.25);   % Plot Welch's PSD
xlabel('Frequency f (Hz)');
ylabel('Power Spectral Density S(f)'); %I am not sure whether this should power or energy we can ask about it
legend('FFT', 'Welch method');
hold off;


%% Wavenumber spectrum 
S_f = Pxx_fft(1:(N/2));
F_k = u_bar /(4*pi) .*S_f;

k = 2*pi*fft_freq(1:(N/2)) / u_bar;
integral_f = 2*trapz(k, F_k); % Multiplying by 2 = variance

%Welch
S_f_Welch = welch_values;
F_k_Welch = u_bar /(4*pi) .*S_f_Welch;

k_Welch = 2*pi*welch_freq / u_bar;
integral_f_Welch = 2*trapz(k_Welch, F_k_Welch);


%von Karman
F_k_karman = (Lambda_f * var_u / pi) ./ (1 + 70.78 .* (k .* Lambda_f / (2 * pi)).^2).^(5/6);

% Plot wave number spectrum
figure(5);
loglog(k, F_k, '.');
hold on;
grid on;
loglog(k_Welch, F_k_Welch, '.');
loglog(k, F_k_karman, 'LineWidth', 1.5);


% Find index where von Karman wavenumber is 10^4
[~, end_idx] = min(abs(k - 1e4)); % Find closest index to k = 10^4
k_end = k(end_idx); % Actual wavenumber at the closest index
F_k_end = F_k_karman(end_idx); % Corresponding F(k) value at this point

% Define the inertial subrange: go backward from k_end to a lower bound (e.g., k_min = 1e1)
k_min = 1e1;
inertial_idx = k >= k_min & k <= k_end; % Indices between k_min and k_end
k_inertial = k(inertial_idx); % Subset of wavenumbers for inertial subrange

% Compute the -5/3 slope line
F_k_ref = F_k_end * (k_inertial / k_end).^(-5/3); % Scale using -5/3 slope
loglog(k_inertial, F_k_ref, 'k-.', 'LineWidth', 1.5); % Overlay the line

% Place vertical reference lines with labels at the bottom
xline(1e2, 'r--', 'Inertial', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(4500, 'r--', 'Dissipation', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

% Characteristic wave numbers - question 10
k_macro = 1 / Lambda_f;
k_micro = 1 / lambda_f;
k_kolmogorov = 1 / eta_K;

xline(k_macro, 'm--', 'k_{macro}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(k_micro, 'm--', 'k_{micro}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(k_kolmogorov, 'm--', 'k_{kolmogorov}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

% Add label for the -5/3 slope
%text(k_inertial(floor(end/2)), F_k_ref(floor(end/2)), '-5/3 slope', 'Rotation', -39, 'FontSize', 10);

% Add legend and labels
legend('Experimental', 'Welch', 'von Karman', '-5/3 slope','Location','southwest');
ylabel('Wave number spectrum F(k)');
xlabel('Wavenumber k');
hold off;


%%





