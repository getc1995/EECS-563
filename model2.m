clear all
close all

% Assume straight lane profile
% Model 2
% Single integrator 3 states
% 	d | x_e - x_o |   | u_x - d_x |
% --- |    y_e    | = |    u_y    |
%  dt |    y_o    |   |    d_y    |

% Bounded u_x, u_y, d_x, d_y
% u_x, d_x in [0, u_max], u_y, d_y in [-vmax, vmax]
% Bounded y_o, y_o in [-ymax, ymax]
% (Other car keep on one lane [-ymax, ymax])






