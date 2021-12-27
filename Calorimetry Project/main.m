%purpose: The purpose of this script is to take in a matrix of data,
%separate the data, fit lines using linear least squares regression,
%combine least squares lines using weighted average, propagate error
%through sigma y, and determine certain temperatures necessary to calculate
%the specific heat of an unknown sample.
%inputs: A matrix of data with time, water temperature, and termocouple
%data: time (s), water temp (degrees C), termocouple data (degrees C)
%outputs: multiple plots of data along with a value for specific heat of
%the sample and corresponding propagated error,
%assumptions: - well insulated
%             - Adiabatic
%             - isolated system
%Author names: Matt Bechtel, Pax Armon
%ID numbers: 109802403, 109799973
%Date created - 10/4/21
%Date modified - 10/14/21

%% Housekeeping
close all; clear; clc;

%% read in data
data = readmatrix('Sample_D.txt');

% splits data into desired vectors
time = data(:, 2);
water_temp = data(:,3);
tc_1 = data(:, 4);
tc_2 = data(:, 5);

%% Fitting Line 1
t_f_1 = time(1:840); % time vector before sample was added for thermocouple 1 and 2 
tc_1_f1 = tc_1(1:840); %tc 1 vector for fit 1
tc_2_f1 = tc_2(1:840); %tc 2 vector for fit 1

[tc_1_line, m1_f1_new, b1_f1_new, sigma_m1_f1, sigma_b1_f1, sigma_y1_f1, Q1_f1] = fitline(t_f_1, tc_1_f1, t_f_1); %fits line to first thermocouple data
[tc_2_line, m2_f1_new, b2_f1_new, sigma_m2_f1, sigma_b2_f1, sigma_y2_f1, Q2_f1] = fitline(t_f_1, tc_2_f1, t_f_1); %fits line to second thermocouple data

[fit1_m_new, sigma_m_new_f1] = weighted_avg(m1_f1_new, m2_f1_new, sigma_m1_f1, sigma_m2_f1); %weighted average of slopes of both lines
[fit1_b_new, sigma_b_new_f1] = weighted_avg(b1_f1_new, b2_f1_new, sigma_b1_f1, sigma_b2_f1); %weighted average of intercepts of both lines

y_f1_new = fit1_m_new * time + fit1_b_new; %final weighted average linear best fit line 1

% time associatd with T0
t_new = 871.881; %[s]

% final T0
T0 = fit1_m_new * t_new + fit1_b_new; %[degrees C]  

% combine sigma y from each line in quadrature to produce a sigma y (book
% version)
sigma_y1_1 = sigma_y1_f1 .* ones(size(time));
sigma_y2_1 = sigma_y2_f1 .* ones(size(time));
sigma_y1_final = sqrt((sigma_y1_1).^2 + (sigma_y2_1).^2);

sigma_T0_initial = err(Q1_f1, Q2_f1, t_new, sigma_y1_final(1));

% final error for T0
sigma_T0 = sqrt((2)^2 + (sigma_T0_initial)^2);

fprintf("T0 = %f with error %f \n", T0, sigma_T0);
%% Fitting Line 2

% max temp from each thermocouple data set
[max_temp_tc1, index1] = max(tc_1);
[max_temp_tc2, index2] = max(tc_2);

% time for t_max fit
t_f_2 = time(974:end); 

tc_1_f2 = tc_1(974:end); %tc 1 data in time range
tc_2_f2 = tc_2(974:end); %tc 2 data in time range

% function calls to fit lines using linear least squares regression
[tc_1_line2, m1_f2_new, b1_f2_new, sigma_m1_f2, sigma_b1_f2, sigma_y1_f2, Q1_f2] = fitline(t_f_2, tc_1_f2, t_f_2); 
[tc_2_line2, m2_f2_new, b2_f2_new, sigma_m2_f2, sigma_b2_f2, sigma_y2_f2, Q2_f2] = fitline(t_f_2, tc_2_f2, t_f_2);

% weighted average 
[fit2_m_new, sigma_m_new_f2] = weighted_avg(m1_f2_new, m2_f2_new, sigma_m1_f2, sigma_m2_f2);
[fit2_b_new, sigma_b_new_f2] = weighted_avg(b1_f2_new, b2_f2_new, sigma_b1_f2, sigma_b2_f2);

y_f2_new = fit2_m_new * time + fit2_b_new;

sigma_y1_2 = sigma_y1_f2 .* ones(size(time));
sigma_y2_2 = sigma_y2_f2 .* ones(size(time));

sigma_y2_final = sqrt((sigma_y1_2).^2 + (sigma_y2_2).^2);

%% Average temp
T_h = fit2_m_new * time(974) + fit2_b_new;
T_l = T0;

T_avg = (T_h + T_l) / 2;

%% Fit third line to T_avg
%time avg temp occurs
t_f_3 = time(840:859); %eye ball using rounded digits to four sigfigs and picked closest index corresponding to temperature

time_2 = time(800:900);

tc_1_f3 = tc_1(840:859);
tc_2_f3 = tc_2(840:859);

[tc_1_line3, m1_f3_new, b1_f3_new, sigma_m1_f3, sigma_b1_f3, sigma_y1_f3] = fitline(t_f_3, tc_1_f3, t_f_3); 
[tc_2_line3, m2_f3_new, b2_f3_new, sigma_m2_f3, sigma_b2_f3, sigma_y2_f3] = fitline(t_f_3, tc_2_f3, t_f_3);

[fit3_m_new, sigma_m_new_f3] = weighted_avg(m1_f3_new, m2_f3_new, sigma_m1_f3, sigma_m2_f3);
[fit3_b_new, sigma_b_new_f3] = weighted_avg(b1_f3_new, b2_f3_new, sigma_b1_f3, sigma_b2_f3);

y_f3_new = fit3_m_new * time_2 + fit3_b_new;

T_avg_vec = ones(1,length(time)) .* T_avg;

%solve for t_avg -- (y - b) / m = t
t_avg = (T_avg - fit3_b_new) / fit3_m_new;

%% Finding T_2
T2 = fit2_m_new * t_avg + fit2_b_new;

sigma_T2_1 = err(Q1_f2, Q2_f2, t_avg, sigma_y2_final(1));
sigma_T2 = sqrt((sigma_T2_1)^2 + (2)^2);

fprintf("T2 = %f with error %f \n", T2, sigma_T2);

%% Finding T_1
time_T1 = time(1:868);
water_temp_T1 = water_temp(1:868);
T1 = mean(water_temp_T1);
sigma_T1_initial = std(water_temp_T1);
sigma_T1 = sqrt((2)^2 + (sigma_T1_initial)^2);

fprintf("T1 = %f with error %f \n", T1, sigma_T1);

%% Solving for Specific Heat
m_s = 88.9; %[g]
sigma_ms = .001; %[g]

m_c = 510; %[g]
sigma_mc = .05; %[g]

C_c = .895; %[J/gC]
num = m_c * C_c * (T2 - T0);
denom = m_s * (T1 - T2);

C_s = num / denom;

fprintf("The specific heat of this sample is: %f J/gC \n", C_s)

%% Error Associated with specific heat

partial_m_c = (C_c * (T2 - T0)) / (m_s * (T1-T2));
partial_T2 = (m_c * C_c * (T1 - T0)) / (m_s * (T1-T2)^2);
partial_T0 = (-m_c * C_c) / (m_s * (T1 - T2));
partial_m_s = (-m_c * C_c * (T2 - T0)) / (m_s^2 * (T1 - T2));
partial_T1 = (-m_c * C_c * (T2 - T0)) / (m_s * (T1 - T2)^2);

sigma_cs = sqrt((partial_m_c * sigma_mc)^2 + (partial_T2 * sigma_T2)^2 + (partial_T0 * sigma_T0)^2 + (partial_m_s * sigma_ms)^2 + (partial_T1 * sigma_T1)^2);

fprintf("The error associated with the specific heat is: %f J/gC", sigma_cs)

%% plotting

figure(1)
%Raw data
subplot(2,2,1)
plot(time, tc_1, 'r.', 'MarkerSize', 6)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 6)
title('Raw Data')
xlabel('Time [s]')
ylabel('Temperature [\circC]')
legend('Thermocouple 1', 'Thermocouple 2', 'FontSize', 12, 'Location', 'NW')
hold off

%T0 fit
subplot(2,2,2)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
plot(time, y_f1_new, 'k', 'LineWidth', 2)
errorbar(time, y_f1_new, sigma_y1_final, 'Color', '#EDB120');
plot(t_new, T0, 'r.', 'MarkerSize', 25);
title('Raw Data with Line T_0 Fit Line', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 15)
ylabel('Temperature [\circC]', 'FontSize', 15)
legend('TC 1', 'TC 2', 'LLS line', 'Error', 'T0', 'FontSize', 12)
set(gca,'FontSize',13)
hold off

%T_max fit
subplot(2,2,3)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
plot(time, y_f2_new, 'k', 'LineWidth', 2)
errorbar(time, y_f2_new, sigma_y2_final, 'Color', '#EDB120');
plot(time(974), T_h, 'm.', 'MarkerSize', 25);
legend('TC 1', 'TC 2', 'T_{max} line fit', 'Error', 'T_{max}', 'FontSize', 12)
title('Raw Data with Max Temp Fit Line', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 15)
ylabel('Temperature [\circC]', 'FontSize', 15)
set(gca,'FontSize',13)
hold off

%T_avg fit
subplot(2,2,4)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
plot(time_2, y_f3_new, 'k', 'LineWidth', 2)
plot(time, T_avg_vec, 'm', 'LineWidth', 2)
title('Raw Data with T_{avg} Fit Line', 'FontSize', 16)
legend('TC 1', 'TC 2', 'LLS Line', 'T_{avg}', 'FontSize', 12)
xlabel('Time [s]', 'FontSize', 15)
ylabel('Temperature [\circC]', 'FontSize', 15)
set(gca,'FontSize',13)
hold off

figure(2)
subplot(2,2,1)
plot(time, water_temp, 'b.', 'MarkerSize', 2)
title('Water Temperature', 'FontSize', 16)
legend('Water Temp', 'FontSize', 12)
xlabel('Time [s]', 'FontSize', 15)
ylabel('Temperature [\circC]', 'FontSize', 15)
set(gca,'FontSize',13)
 
subplot(2,2,2)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
plot(time, y_f1_new, 'g-', 'LineWidth', 2)
plot(time, y_f2_new, 'k-', 'LineWidth', 2)
plot(time_2, y_f3_new, 'k--', 'LineWidth', 2)
title('Raw Data with Final Fits', 'FontSize', 16)
legend('TC 1', 'TC 2', 'T_0 fit', 'T_{max} fit', 'T_{avg} fit', 'FontSize', 12)
xlabel('Time [s]', 'FontSize', 15)
ylabel('Temperature [\circC]', 'FontSize', 15)
set(gca,'FontSize',13)

subplot(2,2,3)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
plot(t_new,T0, 'c.', 'MarkerSize', 15)
plot(t_avg, T2,'k.', 'MarkerSize', 15) 
title('Raw Data with Desired Values')
legend('TC1', 'TC2', 'T0', 'T2')
xlabel('Time [s]')
ylabel('Temperature [\circC]')
hold off

subplot(2,2,4)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
plot(time, m1_f1_new * time + b1_f1_new, 'k-', 'LineWidth', 1)
plot(time, m2_f1_new * time + b2_f1_new, 'k-', 'LineWidth', 1)
plot(time, m1_f2_new * time + b1_f2_new, 'k-', 'LineWidth', 1)
plot(time, m2_f2_new * time + b2_f2_new, 'k-', 'LineWidth', 1)
plot(time_2, m1_f3_new * time_2 + b1_f3_new, 'k-', 'LineWidth', 1)
plot(time_2, m2_f3_new * time_2 + b2_f3_new, 'k-', 'LineWidth', 1)
title('Raw Data with Each Thermocouple Fit')
legend('Thermocouple 1', 'Thermocouple 2', 'Best Fit Lines')
xlabel('Time [s]')
ylabel('Temperature [\circC]')
hold off

figure(3)
plot(time, tc_1, 'r.', 'MarkerSize', 4)
hold on
plot(time, tc_2, 'b.', 'MarkerSize', 4)
errorbar(time, y_f2_new, sigma_y2_final, 'Color', '#EDB120');
plot(time, y_f2_new, 'k', 'LineWidth', 2)
plot(t_avg, T2, 'm.', 'MarkerSize', 15)
title('Raw Data with Best Fit Line Extrapolated to T_2', 'FontSize', 16)
legend('Thermocouple 1', 'Thermocouple 2', 'Error', 'Best Fit Line', 'T2')
xlabel('Time [s]')
ylabel('Temperature [\circC]')
hold off