clear all
close all

% Assume straight lane profile
% Model 2
% double integrator 4 states
% 	d | x_e - x_o |   | u_x - d_x |
% --- | y_e - x_o | = |    u_y    |
%  dt |    y_o    |   |    d_y    |

% Bounded u_x, u_y, d_x, d_y
% u_x, d_x in [0, u_max], u_y, d_y in [-vmax, vmax]
% Bounded y_o, y_o in [-ymax, ymax]
% (Other car keep on one lane [-ymax, ymax])
opts = [];

% system
% x = Ax+B(u-Ed);
A = [zeros(2),eye(2);zeros(2),zeros(2)];
B = [zeros(2);eye(2)];
E = [0;-1];

% Constants
Lx = 6; % meters
Ly = 3; % meters
uxmax = 5; % m/s, longitudinal acceleration
uxmin = -5; % m/s
uymax = 1; % m/s, lateral acceleration
uymin = -1; % m/s
dumax = 0.1; % m/s
dumin = -0.1; % m/s
vxmax = 5; % m/s, longitudinal velocity
vymax = 1; % m/s, lateral velocity

% Variables
% let x = x_e - x_o, x is x diff
% let y = y_e - y_o, h is y diff
% u = -Kx
sdpvar px py vx vy ux uy d 
x = [px; py; vx; vy];
u = [ux;uy];
fx = A*x + B*(u+E*d);

% minimize k to maximize the volume of B >= 0.
sdpvar k
[B_hat, coefB, mon] = polynomial(x,2); % B_hat is sos
Barrier = k - B_hat; % >= 0 means safe
dB = jacobian(Barrier,x);
% Safeset
safe_set1 = px^2/(Lx/2)^2 + py^2/(Ly/2)^2 - 1; % >= 0 real safe
safe_set2 = 1 - (px+7)^2/(20)^2 - (py-1)^2/(8)^2; % >= 0 bounded region

% initial
init_set = 3^2 - (px+17)^2 - (py-1.5)^2;

initial_X = safe_set1 - 0.1;
K = lqr(A,B,eye(4),eye(2));
gamma = 10;
eps = 1e-7;

% Non-empty set constraint
[s0, c0, m0] = polynomial(x,2);
nept_const = -Barrier + s0*initial_X;

% initial set constraint
[s1, c1, m1] = polynomial(x,2);
init_const = Barrier - s1*init_set;

% Contained in safety set
[s21, c21, m21] = polynomial(x,2);
[s22, c22, m22] = polynomial(x,2);
safe_const1 = -(Barrier + eps) + s21*safe_set1; % >= 0, such that B in safe, or unsafe in neg(B)
safe_const2 = -(Barrier + eps) + s22*safe_set2; % >= 0, such that B in safe, or unsafe in neg(B)

% Control barrier constraint
[s3, c3, m3] = polynomial([x;d],2);
[s4, c4, m4] = polynomial([x;d],2);
xd = [5;5;0;0];
den1 = 1 + (-K(1,:)*(x-xd)/(2*uxmax))^2;
den2 = 1 + (-K(2,:)*(x-xd)/(2*uymax))^2;
control_const = -(dB * (A*x + B*E*d) * den1*den2 ...
                + dB * B * [K(1,:)*(x-xd)*den2; K(2,:)*(x-xd)*den1] ...
                + gamma*Barrier * den1*den2 + eps) ...
                + (s3*safe_set1 + s4*(dumax-d)*(d-dumin))*den1*den2;
            
% constraint = [ k>=0; sos(nept_const); sos(init_const); sos(safe_const); sos(control_const);...
%               sos(B_hat); sos(s0); sos(s1); sos(s2);sos(s3);sos(s4);];% sos(s3); sos(s4)];
          
constraint = [ k>=0; sos(init_const); sos(nept_const); sos(safe_const1); sos(safe_const2); sos(control_const);...
               sos(B_hat); sos(s1); sos(s0); sos(s21); sos(s22); sos(s3); sos(s4)];

obj = [];
variables = [k;coefB;c0;c1;c21;c22;c3;c4];%;c5];
%variables = [k;coefB];
[sol,v,Q,res] = solvesos(constraint, obj, opts, variables);

%% plotting
figure
range = [-35 20 -20 35];
hold on
fcontour(str2sym(sdisplay(safe_set1)), range, '-r', 'LevelList', [0])
fcontour(str2sym(sdisplay(safe_set2)), range, '-r', 'LevelList', [0])
fcontour(str2sym(sdisplay(init_set)), range, '-k', 'LevelList', [0])
fcontour(str2sym(sdisplay(initial_X)), range, '-g', 'LevelList', [0])
grid on;
axis equal;
% Substitute variables with values
B_plot = replace(Barrier, [k; coefB; vx; vy], [value(k); value(coefB); 0; 0]);
fb = matlabFunction(str2sym(sdisplay(B_plot)));
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])








