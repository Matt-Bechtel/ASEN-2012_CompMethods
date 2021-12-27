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

%% Variation in theta

% IC
const = constants();

distances = ones(length(const.theta), 2);
figure();
for i = 1:length(const.theta)
    x0 = const.x0;
    y0 = const.y0;
    vx0 = const.v0(1);
    vy0 = const.v0(2);
    m_wat0 = const.V_water_i * const.rho_water;
    m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
    V_air0 = const.V_air_i;
    theta = const.theta(i);
    
    %% Values for ODE45
    tspan = [0, 5];  % [s]
    X0 = [x0; y0; vx0; vy0; m_wat0; m_air0; V_air0; theta]; % State vector

    %% ODE45 call
    odend = odeset('Events', @myEvent);

    [t, X] = ode45(@(t,X) rocketEOM(t, X), tspan, X0, odend);
    
    distances(i, 1) = max(X(:,1));
    distances(i, 2) = theta;
    
    plot(X(:,1), X(:, 2))
    hold on
end

fprintf('Max distance reached is: %f [m]  at angle %f [rads] \n', max(distances(:,1)), distances(find(distances(:,1) == max(distances(:,1))),2)) 

title('X vs Y Position ($\theta$ Varied from 0 to $\frac{\pi}{2}$)')
xlabel('X-Position [$m$]')
ylabel('Y-Position [$m$]')
hold off
set(0,'defaultTextInterpreter','latex')
    
%% variation in water volume




