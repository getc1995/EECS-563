clear all
close all
clear;clc;

% Assume straight lane profile
% Model 1
% Single integrator 2 states
% 	d | x_e-x_o |   | u_x - d_x |
% --- |   y_e   | = |    u_y    |
%  dt 
% Bounded u_x, u_y, d_x
% u_x, d_x in [u_min, u_max], u_y in [-vmax, vmax]
% Assume y_o = 0 all the time. 
% (Other car follows a straight line y = 0)

%%
opts = []
% opts = sdpsettings('solver','mosek','sos.model',1,...
%     'sos.scale',1,'verbose',1,'sos.newton',1,'sos.congruence',1);

% Define the system
A = zeros(2);
B = eye(2);
E = [-1; 0];

% Constants
Lx = 6; % meters
Ly = 3; % meters
umax = 20; % m/s, longitudinal velocity
umin = 10; % m/s
dumax = 16; % m/s
dumin = 14; % m/s
vmax = 1; % m/s, lateral velocity

% Variables
% let h = x_e - x_o, h is headway
sdpvar h ye ux uy dx k 
x = [h; ye];
% dxdt = [ux-dx; uy]; % dx/dt
% xd = [h; ye; dx]; % xd is x augmented with disturbance
[B_hat, coefB, mon] = polynomial(x,2); % B_hat is sos

% Input constraint

%% Initialize step
% initial_X is the polynomial used to initialize barrier
% such that initial_X in barrier 
initial_X = safe - 99; 
% Control gain
K = [-10, 0;
       0, -11];
gamma = 1;
eps = 1e-7;

% Non-empty set constraint
[s0, c0, m0] = polynomial(x,2);
initial_const = Barrier - s0*initial_X;

% Contained in safety set
[s1, c1, m1] = polynomial(x,2);
safe_const = -Barrier + s1*safe - eps; % >= 0, such that B in safe, or unsafe in neg(B)

% Control barrier constraint
[s2, c2, m2] = polynomial([x;dx],2);
[s3, c3, m3] = polynomial([x;dx],2);
% for ux, uy constraint
% [s4, c4, m4] = polynomial([x;ux;uy;dx],2);
% [s5, c5, m5] = polynomial([x;ux;uy;dx],2);

u_abs = (umax - umin)/2;
den1 = 1;% + ((-K(1,:)*x)/(2*u_abs))^2;
den2 = 1 + ((-K(2,:)*x)/(2*vmax))^2;
% 
control_const = dB * (A*x + B*[15;0] + E*dx) * den1*den2 ...
                + dB * B * [0; K(2,:)*(x-[0;-4])*den1] ...
                + gamma*Barrier * den1*den2 ...
                - s2*safe - s3*(dumax-dx)*(dx-dumin) - eps;
%                 - s4*(umax-ux)*(ux-umin) - s5*(vmax-uy)^2 - eps; % >= 0,
%             
% 
constraint = [ k>=0; sos(initial_const); sos(safe_const); sos(control_const);...
              sos(B_hat); sos(s0); sos(s1); sos(s2);sos(s3);];% sos(s3); sos(s4)];

obj = [];
variables = [k;coefB;c0;c1;c2;c3];%c4;c5];
[sol,v,Q,res] = solvesos(constraint, obj, opts, variables);


%% plotting
figure
range = [-10 10 -10 10];
fcontour(str2sym(sdisplay(safe)), range, '-r', 'LevelList', [0])
grid on;
axis equal;
hold on
% Substitute variables with values
B_plot = replace(Barrier, [k; coefB], [value(k); value(coefB)]);
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

