clear all; close all;

load("out_MatRANS3.mat");
load("Exercise4.mat");
load("CaseInfo.mat");


U0m = MatRANS.U0m;
T = CaseInfo.T;

omega = 2 * pi / T;

a = U0m / omega;

phase_angle = 90;
phase_angle_rad = phase_angle .* pi / 180 + 8 * pi; % using data from the fifth period
time = phase_angle_rad ./ omega;

index = 307;
u_case3 = MatRANS.u(index, :) ./ U0m;
k_case3 = MatRANS.k(index, :) ./ U0m^2;
y_case3 = MatRANS.y ./a;

load("out_MatRANS4.mat");
u_case4 = MatRANS.u(index,:) ./ U0m;
k_case4 = MatRANS.k(index, :) ./ U0m^2;
y_case4 = MatRANS.y ./a;

subplot(1,2,1);
plot(u_case3, y_case3,"LineWidth", 1);
grid on;
hold on;
plot(u_case4, y_case4, "LineWidth", 1);
legend("Case 3: smooth turbulent", "Case 4: rough turbulent");
ylabel("$y / a $", 'Interpreter', 'latex', 'FontSize', 12);
xlabel("$\overline{U} / U_{0m}$", 'Interpreter', 'latex', 'FontSize', 12);

subplot(1,2,2);
plot(k_case3, y_case3, "LineWidth",1);
grid on;
hold on;
plot(k_case4, y_case4, "LineWidth", 1);
legend("Case 3: smooth turbulent", "Case 4: rough turbulent");
ylabel("$y / a $", 'Interpreter', 'latex', 'FontSize', 12);
xlabel("$k / U_{0m}^2$", 'Interpreter', 'latex', 'FontSize', 12);


