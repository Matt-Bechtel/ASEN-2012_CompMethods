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
% Date Modified: 12/8/21

%% Housekeeping
close all
clear
clc

   %% Set up variation vectors
theta = pi/12:pi/48:(7*pi)/18;
drag = 0.2:.1:0.8;
P_gage = 0:1:80;
P_gage = P_gage .* 6895;
V_water_i = 0:.0005:(.002-.0008);

count = 0;

% Preallocate -- ran a few times so we knew had big these vctors would be
params = zeros(80, 5);
index = zeros(80, 5);

%% Variation
for i = 1:length(theta)
    for j = 1:length(drag)
        for k = 1:length(P_gage)
            for l = 1:length(V_water_i)
                % get constants for specific variation
                const = constants(i, j, k, l);
    
                % set values initial state vector
                x0 = const.x0;
                y0 = const.y0;
                vx0 = const.v0(1);
                vy0 = const.v0(2);
                m_wat0 = const.V_water_i * const.rho_water;
                m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
                V_air0 = const.V_air_i;
                
                %tracker = 1;
                
                %% Values for ODE45
                tspan = [0, 8];  % [s]
                X0 = [x0; y0; vx0; vy0; m_wat0; m_air0; V_air0]; % State vector

                %% ODE45 call
                odend = odeset('Events', @myEvent);

                [t, X, ~] = ode45(@(t,X) rocketEOM(t, X, i, j, k, l), tspan, X0, odend);
                
                % get all valid combinations
                if max(X(:,1)) >= 84.5 && max(X(:,1)) <= 85.5
                    count = count + 1;
                    fprintf('The max distance is %f\n', max(X(:,1)));
                    
                    % keeps track of value for paramter
                    params(count, 1) = theta(i);
                    params(count, 2) = drag(j);
                    params(count, 3) = P_gage(k);
                    params(count, 4) = V_water_i(l);
                    params(count, 5) = X(end, 1);
                    
                    % keeps track of index at which parameter occurs --
                    % bottom line keeps track of absolute value of
                    % difference from 85m
                    index(count, 1) = i;
                    index(count, 2) = j;
                    index(count, 3) = k;
                    index(count, 4) = l;
                    index(count, 5) = abs(85 - X(end, 1));
                    
                    % plot all trajectories
                    plot(X(:,1), X(:,2));
                    hold on
                    
                end
            end
        end
    end
end

set(0,'defaultTextInterpreter','latex')

best_index = find(index(:, 5) == min(index(:, 5)));
 fprintf('Theta = %f [degrees]\n Drag Coefficient = %f\n Initial Gage Pressure %f [Pa]\n Initial Water Volume %f [m^3]\n', ...
     params(best_index, 1) * (180 / pi), params(best_index, 2), params(best_index, 3), params(best_index, 4));
 
title('Distance vs Height for Successful Combinations')
xlabel('Distance [$m$]')
ylabel('Height [$m$]')

hold off
          

% %% Thrust calculation
% Found this on mathworks forum to access variables in ODE45 

figure();

for i = 1:length(params)
    [~, F_thrust] = cellfun(@(t, X) rocketEOM(t, X.', index(i, 1), index(i, 2), index(i, 3), index(i, 4)), num2cell(t), num2cell(X,2), 'uni', 0);

    % convert F_thrust from cellfun from cell to matrix
    F_thrust = cell2mat(F_thrust');

    % Get into right dimensions
    F_thrust = F_thrust';

    %Create normalized thrust vector
    F = zeros(length(F_thrust), 1);

    for j = 1:length(F_thrust)
        F(j) = norm(F_thrust(j, :));
    end
    
    plot(t, F)
    hold on
end

xlim([0 .45])
ylim([0 500])

set(0,'defaultTextInterpreter','latex')

title('Thrust vs Time for Successful Combinations')
xlabel('Time [$s$]')
ylabel('Thrust [$N$]')

hold off

%% ODE call for best case 

const = constants(index(best_index, 1), index(best_index, 2), index(best_index, 3), index(best_index, 4));

x0 = const.x0;
y0 = const.y0;
vx0 = const.v0(1);
vy0 = const.v0(2);
m_wat0 = const.V_water_i * const.rho_water;
m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
V_air0 = const.V_air_i;

% tracker = 1;

% Values for ODE45
tspan = [0, 8];  % [s]
X0 = [x0; y0; vx0; vy0; m_wat0; m_air0; V_air0]; % State vector

% ODE45 call
odend = odeset('Events', @myEvent);

[t, X] = ode45(@(t,X) rocketEOM(t, X, index(best_index, 1), index(best_index, 2), index(best_index, 3), index(best_index, 4)), tspan, X0, odend);

%% Position plots
figure ()

plot(X(:,1), X(:,2), 'LineWidth', 3);
hold on
plot(X(end, 1), X(end, 2), 'r.', 'MarkerSize', 20)

set(0,'defaultTextInterpreter','latex')

title('Distance vs Height')
xlabel('Distance [$m$]')
ylabel('Height [$m$]')
legend('Position Profile', 'Final Position')

%% Thrust 
[~, F_thrust, tracker] = cellfun(@(t, X) rocketEOM(t, X.', index(best_index, 1), index(best_index, 2), index(best_index, 3), index(best_index, 4)), num2cell(t), num2cell(X,2), 'uni', 0);

% convert F_thrust from cellfun from cell to matrix
F_thrust = cell2mat(F_thrust');
tracker = cell2mat(tracker');

% Get into right dimensions
F_thrust = F_thrust';
tracker = tracker';

%Create normalized thrust vector
F = zeros(length(F_thrust), 1);

for j = 1:length(F_thrust)
    F(j) = norm(F_thrust(j, :));
end

%% Thrust plot
figure ()

x = t';
y = F';
z = (zeros(size(x)));
c = tracker';
surface([x;x],[y;y],[z;z],[c;c],'facecolor','none','edgecolor','interp','linewidth',2);
title('Thrust vs Time')
xlabel('Time [$s$]')
ylabel('Thrust [$N$]')
legend('Thrust Profile')

xlim([0 .3])
ylim([0 250])

colorbar('Ticks', [1 2 3], 'TickLabels', {'Phase 1', 'Phase 2', 'Phase 3'})

beep on; beep;



