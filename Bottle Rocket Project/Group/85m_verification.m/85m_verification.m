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
% Date Modified: 11/19/21]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               NOTE                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Professor Anderson told us to publisht his code since our actual
% variation was too complex to publish. This is the published version using the parameters
% that land us between 84.5 and 85.5m. Our varied code contains plots with
% the phase changes but this one does not. Professor Anderson said 
% that this was okay. The code in the zip file is where the parameters
% are varied and plots with phase changes. 


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
tspan = [0, 8];  % [s]
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             FUNCTIONS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose: This function sets all constants used for the Bottle Rocket
%          project
% 
% Inputs: None
%
% Outputs: const -> struct of constants
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


function const = constants()

    const.g = [0; 9.81];                  % [m/s^2] gravity
    const.C_discharge = 0.8;             % coefficient of discharge
    const.rho_air = .961;                % [kg/m^3]
    const.V_bottle = .002;               % [m^3]
    const.P_amb = 12.1 * 6894.76;        % [Pa] ** remove 6894.76 to convert back to psi ambient
    const.gamma = 1.4;
    const.rho_water = 1000;              % [kg/m^2]
    const.D_throat = .021;                % [m]
    const.D_bottle = .105;               % [m]
    const.R = 287;                       % [J/kgK]
    const.M_bottle = .15;                 % [kg]
    const.C_drag = 0.2;
    const.P_gage = 427490;         % [Pa]
    const.V_water_i = .001;              % [m^3]
    const.T_air_i = 300;                 % [k]
    const.v0 = [0; 0];                   % [m/s]
    const.theta = 63.75 * (pi / 180);                  % [rad]
    const.x0 = 0;                        % [m]
    const.y0 = .25;                      % [m]
    const.length = 0.5;                  % [m]
    const.Ac = (pi / 4) * (const.D_bottle)^2;  % [m^2]
    const.At = (pi / 4) * (const.D_throat)^2; % [m^2]
    
    const.P_init = const.P_gage + const.P_amb; % [Pa]
    const.V_air_i = const.V_bottle - const.V_water_i;
    
    const.m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
    
end

% Purpose: The purpose of this script is to create the equaations of motion
%          and rates of change necessary for ODE45. These are primarily thrust,
%          gravity, and drag forces which are used to calculate acceleration.
%          Additionally, this scripts calculate sthe rate of change of the mass of
%          water, mass of air, and volume of air that ODE45 then numerically
%          integrates
% 
% Inputs: t -> time
%         X -> state vector containing position, velocity, mass of water,
%         mass of air, and volume of air
% Outputs: v -> velocity in x and y directions
%          a -> acceleration in x and y directions
%          m_water_dot -> rate of change of mass of water
%          m_air_dot -> rate of change of mass of air
%          dVdt -> rate of change of volume of air
% Assumptions: Air expansion in bottle is adiabatic and isentropic
% Author: Matt Bechtel
% ID Number: 109802403
% Date Created: 11/9/21
% Date Modified: 11/19/21

function [dX, F_thrust] = rocketEOM(t, X)
%% Extract from state vector
pos = [X(1); X(2)];     % position [x and y]
v = [X(3); X(4)];       % velocity [x and y]
m_water = X(5);
m_air = X(6);
V_air = X(7);

%% Get Constants
const = constants();

%% Heading
if pos(1) < const.length * cos(const.theta)
    heading = [cos(const.theta); sin(const.theta)];
else  
    heading = v / norm(v);
end

% Pressure used before water is exhausted
P_end = const.P_init * (const.V_air_i / V_air) ^ const.gamma;

P_end_2 = const.P_init * (const.V_air_i / const.V_bottle) ^ const.gamma;

% mass of rocket
mr = const.M_bottle + (const.rho_water * (const.V_bottle - V_air)) + m_air;

% Pressure used after water is exhausted
P = P_end_2 * (m_air / const.m_air0) ^ const.gamma;


%% Before water is exhausted
if V_air < const.V_bottle && m_water >= .002
    
    Ve = sqrt(2 * (P_end - const.P_amb) / const.rho_water);
    % mass flow rate for water
    m_water_dot = -const.C_discharge * const.rho_water * const.At * Ve;
    
   % Force calculations
    F_grav = -mr * const.g;
    F_drag = -heading .* .5 * const.rho_air * (norm(v)^2) * const.C_drag * const.Ac;
    F_thrust = -heading .* m_water_dot * Ve; 
    
    % rate of change of volume of air
    dVdt = const.C_discharge * const.At * Ve;
    
    %rate of change of mass of air
    m_air_dot = 0;
    
%% After water is exhausted
elseif P > const.P_amb
    
    % pressure, density, and temperature at any time
    rho = m_air / const.V_bottle;
    T = P / (rho * const.R);
    
    % Critical pressure
    P_crit = P * (2 / (const.gamma + 1)) ^ (const.gamma / (const.gamma - 1));
    
    % Choked
    if P_crit > const.P_amb
        
        Te = T * (2 / (const.gamma + 1));
        Pe = P_crit;
        rhoe = Pe / (const.R * Te);
        Ve = sqrt(const.gamma * const.R * ((2 / (const.gamma + 1)) * T));
        
    % Not choked
    else
  
       Me = sqrt(((P / const.P_amb) ^ ((const.gamma - 1) / const.gamma) - 1) / ((const.gamma - 1) / 2));
       Te = T / (1 + ((const.gamma - 1) / 2) * Me ^ 2);
       rhoe = const.P_amb / (const.R * Te);
       Pe = const.P_amb;
       Ve = Me * sqrt(const.gamma * const.R * Te); 
       
    end
    
    % rates of change
    m_air_dot = -const.C_discharge * rhoe * const.At * Ve;
    m_water_dot = 0;
    dVdt = 0;
    
    % total mass of rocket
    mr = const.M_bottle + m_air;
    
    % force calculations
    F_thrust = -heading .* (m_air_dot * Ve + (const.P_amb - Pe) * const.At); %change back to neg
    F_grav = -mr * const.g;
    F_drag = -heading .* (.5 * const.rho_air * (norm(v)^2) * const.C_drag * const.Ac);
    
%% Ballistic
else 
    % rates of change
    mr_dot = 0;
    m_air_dot = 0;
    dVdt = 0;
    mr = const.M_bottle + m_air;
    m_water_dot = 0;
    
    % Force calculations
    F_thrust = heading .* 0;
    F_drag = -heading .* (.5 * const.rho_air * (norm(v)^2) * const.C_drag * const.Ac); % can change norm(v) back to v.^2 if need be
    F_grav = -const.M_bottle  * const.g;
end

%% Solving for net force and acceleteration
fnet = F_thrust + F_drag + F_grav;

a = fnet / mr;

% output state vector
dX = [v(1); v(2); a(1); a(2); m_water_dot; m_air_dot; dVdt]; 

end

% Purpose: This function sets the limit on the ODE45 integration. Learned
%          this in Office Hours. Not quite sure how it works.
% 
% Inputs: NA
%
% Outputs: NA
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

function [val, ist, dur] = myEvent(~, X)
    val = (X(2) <= 0);
    ist = 1;
    dur = 0;
end
