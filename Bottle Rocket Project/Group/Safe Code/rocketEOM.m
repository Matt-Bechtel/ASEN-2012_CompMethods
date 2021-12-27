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
theta = X(8);

%% Get Constants
const = constants();

%% Heading
if pos(1) < const.length * cos(theta)
    heading = [cos(theta); sin(theta)];
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

dTheta = 0;

% output state vector
dX = [v(1); v(2); a(1); a(2); m_water_dot; m_air_dot; dVdt; dTheta]; 

end