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


function const = constants(n)

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
    const.C_drag = 0.5;
    const.P_gage = 50 * 6894.76;         % [Pa]
    const.V_water_i = .0001:.0001:.002;              % [m^3]    initial .001
    const.T_air_i = 300;                 % [k]
    const.v0 = [0; 0];                   % [m/s]
    const.theta = pi/4;                  % [rad]
    const.x0 = 0;                        % [m]
    const.y0 = .25;                      % [m]
    const.length = 0.5;                  % [m]
    const.Ac = (pi / 4) * (const.D_bottle)^2;  % [m^2]
    const.At = (pi / 4) * (const.D_throat)^2; % [m^2]
    
    const.P_init = const.P_gage + const.P_amb; % [Pa]
    const.V_air_i = const.V_bottle - const.V_water_i(n);
    
    const.m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
    
end
