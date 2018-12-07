clear all
close all

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
opts = sdpsettings('solver','mosek','sos.model',1,...
    'sos.scale',1,'verbose',1,'sos.newton',1,'sos.congruence',1);

% Define the system
A = zeros(2);
B = eye(2);
E = [-1; 0];

% Constants
Lx = 6; % meters
Ly = 3; % meters
umax = 20; % m/s, longitudinal velocity
umin = 10; % m/s
vmax = 1; % m/s, lateral velocity

% Variables
% let h = x_e - x_o, h is headway
sdpvar h ye ux uy dx k 
x = [h; ye];
dxdt = [ux-dx; uy]; % dx/dt
% xd = [h; ye; dx]; % xd is x augmented with disturbance
[B_hat, coefB, mon] = polynomial(x,2); % B_hat is sos

% minimize k to maximize the volume of B >= 0.
Barrier = B_hat - k; % >= 0 means safe
dB = jacobian(Barrier,x);
% Safeset
safe = h^2/(Lx/2)^2 + ye^2/(Ly/2)^2 - 1; % >= 0
% Input constraint

%% Initialize step
% initial_X is the polynomial used to initialize barrier
% such that initial_X in barrier 
initial_X = safe - 99; 
% Control gain
K = [-10, 0;
       0, -10];
gamma = 1;

% Non-empty set constraint
[s0, c0, m0] = polynomial(x,2);
initial_const = Barrier - s0*initial_X;

% Contained in safety set
[s1, c1, m1] = polynomial(x,2);
safe_const = -Barrier + s1*safe ; % >= 0, such that B in safe, or unsafe in neg(B)

% Control barrier constraint
control_const = dB * dx + gamma * Barrier; % >= 0,

constraint = [sos(B_hat); sos(initial_const); sos(safe_const);...
              sos(s0); sos(s1)];

obj = [];
variables = [k;coefB;c0;c1];
[sol,v,Q,res] = solvesos(constraint, obj, opts, variables);


%% plotting
figure
range = [-50 50 -50 50];
fcontour(str2sym(sdisplay(safe)), range, '-r', 'LevelList', [0])
grid on;
axis equal;
hold on
% Substitute variables with values
B_plot = replace(Barrier, [k; coefB], [value(k); value(coefB)]);
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

