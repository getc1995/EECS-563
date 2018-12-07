clear all
close all
clc
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
% opts = sdpsettings('solver','sedumi','sos.model',2,'sos.scale',1,'verbose',1,...
%    'sos.newton',1,'sos.congruence',1);
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
gamma = 5;
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
% + dB * B * [K(1,:)*(x-[0;4])*den2; K(2,:)*(x-[0;4])*den1]
control_const = dB * (A*x + B*[(umin+umax)/2;0] + E*dx) * den1*den2 ...
                + dB * B * [0; K(2,:)*(x-[0;4])*den1] ...
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
B_plot = replace(Barrier, [k; coefB], [clean(value(k),1e-8); clean(value(coefB),1e-8)]);
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

%%
while true
% for idx = 1
    old_k = value(k);
    %% hold barrier
    old_Bhat = replace(B_hat, coefB, value(coefB));
    expand_const = old_k - k >= 0;
    
    [s0, c0, m0] = polynomial(x,2);
    safe_const = -(old_Bhat - k) + s0*safe - eps;
    
    dB = jacobian(old_Bhat - k, x);
    [ux_input, cx_input, mx_input] = polynomial(x,2);
    [uy_input, cy_input, my_input] = polynomial(x,2);
    
    [s11, c11, m11] = polynomial([x;dx],2);
    [s12, c12, m12] = polynomial([x;dx],2);
    control_const = dB * (A*x + B*[ux_input;uy_input] + E*dx)...
                    + gamma*Barrier...
                    - s11*safe - s12*(dumax-dx)*(dx-dumin) - eps;
    
    [s21, c21, m21] = polynomial([x;dx],2);
    [s22, c22, m22] = polynomial([x;dx],2);
    ux_const1 = ux_input - umin - s21*safe - s22*(dumax-dx)*(dx-dumin);
    
    [s31, c31, m31] = polynomial([x;dx],2);
    [s32, c32, m32] = polynomial([x;dx],2);
    ux_const2 = -ux_input + umax - s31*safe - s32*(dumax-dx)*(dx-dumin);
    
    [s41, c41, m41] = polynomial([x;dx],2);
    [s42, c42, m42] = polynomial([x;dx],2);
    uy_const1 = uy_input + vmax - s41*safe - s42*(dumax-dx)*(dx-dumin);
    
    [s51, c51, m51] = polynomial([x;dx],2);
    [s52, c52, m52] = polynomial([x;dx],2);
    uy_const2 = -uy_input + vmax - s51*safe - s52*(dumax-dx)*(dx-dumin);
    
    constraint = [k >= 0; expand_const; sos(safe_const); sos(control_const);...
                  sos(ux_const1);sos(ux_const2);...
                  sos(uy_const1);sos(uy_const2);...
                  sos(s0); sos(s11); sos(s12);...
                  sos(s21); sos(s22);...
                  sos(s31); sos(s32);...
                  sos(s41); sos(s42);...
                  sos(s51); sos(s52)];
    obj = [k];
    variables = [k;c0;c11;c12;c21;c22;c31;c32;c41;c42;c51;c52];
    [sol,v,Q,res] = solvesos(constraint, obj, opts, variables);
    %% Hold control u
    old_k = value(k);
    old_Bhat = replace(B_hat, coefB, value(coefB));
    [B_hat, coefB, mon] = polynomial(x,2);
    
    [s0, c0, m0] = polynomial(x,2);
    safe_const = -(old_Bhat - k) + s0*safe - eps;
    
    [s1, c1, m1] = polynomial(x,2);
    Barrier = B_hat - k;
    expand_const = Barrier - s1*(old_Bhat - old_k);
% %   hold control input
    ux_input = replace(ux_input, [cx_input], value(cx_input));
    uy_input = replace(uy_input, [cy_input], value(cy_input));
    [s11, c11, m11] = polynomial([x;dx],2);
    [s12, c12, m12] = polynomial([x;dx],2);
    
    dB = jacobian(B_hat - k,x);
    control_const = dB * (A*x + B*[ux_input;uy_input] + E*dx)...
                    + gamma*Barrier...
                    - s11*safe - s12*(dumax-dx)*(dx-dumin) - eps;
                
    variables = [k;coefB;c0;c11;c12;];
    
    constraint = [k >= 0; sos(expand_const); sos(safe_const); sos(control_const);...
                  sos(B_hat); sos(s0); sos(s1); sos(s11); sos(s12);...
                  sos(s51); sos(s52)];
    %%
    if abs(value(k) - old_k) < 1e-4
        break
    end
end
%%
figure
range = [-10 10 -10 10];
fcontour(str2sym(sdisplay(safe)), range, '-r', 'LevelList', [0])
grid on;
axis equal;
hold on
% Substitute variables with values
B_plot = replace(B_hat-k, [k; coefB], [clean(value(k),1e-8); clean(value(coefB),1e-8)]);
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])