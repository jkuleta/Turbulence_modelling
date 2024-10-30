clear all; close all;

load("out_MatRANS_case4.mat");
load("Exercise4.mat");
load("CaseInfo.mat");

n = [1, 4, 7, 10];
c = 4; % case number from the data

U0m = MatRANS.U0m;
T = CaseInfo.T;

omega = 2 * pi / T;

a = U0m / omega;

phase_angle = [0, 45, 90, 135];
phase_angle_rad = phase_angle .* pi / 180 + 8 * pi; % using data from the fifth period
time_list = phase_angle_rad ./ omega;

tolerance = 1e-5;
index_list = zeros(size(time_list));

for i = 1:length(time_list)
    index = find(abs(MatRANS.t - time_list(i)) < tolerance);
     if ~isempty(index)
        index_list(i) = index(1);
    else
        index_list(i) = 0;
    end
end

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);

for i = 1:length(phase_angle_rad)
    index = index_list(i);
    u_plot = MatRANS.u(index, :) ./ U0m;
    y_plot = MatRANS.y ./ a;

    u_comparison = WBL(c).u(:, n(i)) ./ WBL(c).U0m;
    y_comparison = WBL(c).y_u ./ a;

    subplot(1, length(phase_angle_rad), i);
    plot(u_plot, y_plot, "LineWidth", 1);
    grid on;
    hold on;
    plot(u_comparison, y_comparison, "o", "MarkerFaceColor", [1, 0.5, 0], "MarkerSize", 3);

    if i == 4
        legend("k-$\omega$ model", "Experimental", 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
    end

    ylabel("$y / a $", 'Interpreter', 'latex', 'FontSize', 12);
    xlabel("$\overline{U} / U_{0m}$", 'Interpreter', 'latex', 'FontSize', 12);

    title(sprintf('$\\omega t = %d^\\circ$', phase_angle(i)), 'FontSize', 12, 'Interpreter', 'latex');
end

saveas(gcf, sprintf('case%d_u', c), 'png');


% Create a new figure for k/U0m^2 vs. y/a plots

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);

for i = 1:length(phase_angle_rad)
    index = index_list(i);

    k_comparison = 0.65 * (WBL(c).uu(:, n(i)) + WBL(c).vv(:, n(i)));
    k_comparison = k_comparison ./ (WBL(c).U0m^2);

    y_comparison = WBL(c).y_uuvv ./ a;
    k_plot = MatRANS.k(index, :) ./ (U0m^2);

    subplot(1, length(phase_angle_rad), i);
    plot(k_plot, y_plot, "LineWidth", 1);
    grid on;
    hold on;
    plot(k_comparison, y_comparison, "o", "MarkerFaceColor", [1, 0.5, 0], "MarkerSize", 3);

    if i == 4
        legend("k-$\omega$ model", "Experimental", 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
    end

    ylabel("$y / a $", 'Interpreter', 'latex', 'FontSize', 12);
    xlabel("$k / U_{0m}^2$", 'Interpreter', 'latex', 'FontSize', 12);

    title(sprintf('$\\omega t = %d^\\circ$', phase_angle(i)), 'FontSize', 12, 'Interpreter', 'latex');
end

saveas(gcf, sprintf('case%d_k', c), 'png');

figure;
set(gcf, 'Position', [100, 100, 1200, 400]);

for i = 1:length(phase_angle)
    index = index_list(i);

    uv = MatRANS.nu_t(index, :) .* gradient(MatRANS.u(index, :), MatRANS.y);
    uv_plot = uv ./ (U0m^2);
    y_plot = MatRANS.y ./ a;

    uv_comparison = -WBL(c).uv(:, n(i)) ./ (WBL(c).U0m^2);
    y_comparison = WBL(c).y_uv ./ a;

    subplot(1, length(phase_angle), i);
    plot(uv_plot, y_plot, "LineWidth", 1);
    grid on;
    hold on;
    plot(uv_comparison, y_comparison, "o", "MarkerFaceColor", [1, 0.5, 0], "MarkerSize", 3);

    if i == 4
        legend("k-$\omega$ model", "Experimental", 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
    end

    ylabel("$y / a$", 'Interpreter', 'latex', 'FontSize', 12);
    xlabel("$-\overline{u'v'} / U_{0m}^2$", 'Interpreter', 'latex', 'FontSize', 12);

    title(sprintf('$\\omega t = %d^\\circ$', phase_angle(i)), 'FontSize', 12, 'Interpreter', 'latex');
end

saveas(gcf, sprintf('case%d_uv', c), 'png');


time = MatRANS.t(index_list(1):end);
omegat = (omega * time - 8 * pi) * 180 / pi;
tau0 = MatRANS.tau0(index_list(1):end);
tau0_plot = tau0 ./ (MatRANS.rho * U0m^2);

tau0_comparison = WBL(c).tau0 / (MatRANS.rho * U0m^2);
omegat_comparison = WBL(c).omegat_tau0;

figure;
plot(omegat, tau0_plot, "LineWidth", 1.75);
grid on;
hold on;
plot(omegat_comparison, tau0_comparison, "o", "MarkerFaceColor", [1, 0.5, 0], "MarkerSize", 3);
legend("k-$\omega$ model", "Experimental", 'Interpreter', 'latex', 'FontSize', 12);
ylabel("$\tau_0 / \rho U_{0m}^2$", 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$\omega t\ [^\circ]$', 'Interpreter', 'latex', 'FontSize', 12);

saveas(gcf, sprintf('case%d_tau0',c) , 'png');

