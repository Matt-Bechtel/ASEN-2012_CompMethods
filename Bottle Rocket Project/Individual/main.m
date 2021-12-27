% Purpose: The purpose of this script is to calculate and plot the
%          trajectory of a bottle rocket using ODE45  
% 
% Inputs: None
%
% Outputs: None
%
% Assumptions: NA
%
% Author: Matt Bechtel
%
% ID Number: 109802403
%
% Date Created: 11/9/21
%
% Date Modified: 11/19/21

%% Housekeeping
close all
clear
clc

%% Initial conditions
const = constants();

x0 = const.x0;
y0 = const.y0;
vx0 = const.v0(1);
vy0 = const.v0(2);
m_wat0 = const.V_water_i * const.rho_water;
m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
V_air0 = const.V_air_i; 

%% Values for ODE45
tspan = [0, 5];  % [s]
X0 = [x0; y0; vx0; vy0; m_wat0; m_air0; V_air0]; % State vector

%% ODE45 call
odend = odeset('Events', @myEvent);

[t, X] = ode45(@(t,X) rocketEOM(t, X), tspan, X0, odend);

%% Thrust calculation
% Found this on mathworks forum to access variables in ODE45 
[~, F_thrust] = cellfun(@(t, X) rocketEOM(t, X.'), num2cell(t), num2cell(X,2), 'uni', 0);
 
% convert F_thrust from cellfun from cell to matrix
F_thrust = cell2mat(F_thrust');

% Get into right dimensions
F_thrust = F_thrust';

%Create normalized thrust vector
F = zeros(length(F_thrust), 1);

    for i = 1:length(F_thrust)
        F(i) = norm(F_thrust(i, :));
    end

%% Plotting

figure
plot(t, F)
xlabel('Time [$s$]')
ylabel('Thrust [$N$]')
xlim([0 .45])
ylim([0 200])
title('Thrust vs Time')
legend('Thrust')
    
figure
subplot(2,2,1)
plot(X(:,1), X(:, 2));
xlabel('Distance [$m$]')
ylabel('Height [$m$]')
title('Position vs Height')
legend('XY-Position')

subplot(2,2,2)
plot(t, X(:,4))
xlabel('Time [$s$]')
ylabel('Y-Velocity [$m/s$]')
title('Y-Velocity vs Time')
legend('Y-Velocity')

subplot(2,2,3)
plot(t, X(:,3))
xlabel('Time [$s$]')
ylabel('X-Velocity [$m/s$]')
title('X-Velocity vs Time')
legend('X-Velocity')

subplot(2,2,4)
plot(t, X(:,7))
xlim([0 .25])
ylim([1 * 10^-3 2 * 10^-3])
xlabel('Time [$s$]')
ylabel('Volume of Air [$m^3$]')
title('Volume of Air vs Time')
legend('Volume of air')

set(0,'defaultTextInterpreter','latex')

height = max(X(:,2));
fprintf('Maximum height reached is: %f [m] \n', height);

distance = max(X(:,1));
fprintf('Maximum distance reached is: %f [m] \n', distance);

