clear all;

load("out_MatRANS_case3.mat");
load("Exercise4.mat");

n = [1, 4, 7, 10];
c = 3; %case number from the data


figure;
plot(WBL(c).omegat_tau0, WBL(c).tau0, ".");
grid on;
