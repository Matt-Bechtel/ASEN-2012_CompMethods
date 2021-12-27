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
% Date Modified: 12/8/21


function const = constants(i, j, k, l)
    
    theta = pi/12:pi/48:(7*pi)/18;
    drag = 0.2:.1:0.8;
    P_gage = 0:1:80;
    P_gage = P_gage .* 6895;
    V_water_i = 0:.0005:(.002-.0008);
    
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
    const.C_drag = drag(j);              %                          varies
    const.P_gage = P_gage(k);         % [Pa]
    const.V_water_i = V_water_i(l);      % [m^3]                    varies
    const.T_air_i = 300;                 % [k]
    const.v0 = [0; 0];                   % [m/s]
    const.theta = theta(i);         % [rad]                    varies
    const.x0 = 0;                        % [m]
    const.y0 = .25;                      % [m]
    const.length = 0.5;                  % [m]
    const.Ac = (pi / 4) * (const.D_bottle)^2;  % [m^2]
    const.At = (pi / 4) * (const.D_throat)^2; % [m^2]
    
    const.P_init = const.P_gage + const.P_amb; % [Pa]
    const.V_air_i = const.V_bottle - const.V_water_i;
    
    const.m_air0 = (const.P_init * (const.V_air_i)) / (const.R * const.T_air_i);
    
end
