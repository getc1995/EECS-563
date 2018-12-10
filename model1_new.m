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
opts = sdpsettings('solver','sedumi','sos.model',2,'sos.scale',1,'verbose',1,...
    'sos.newton',1,'sos.congruence',1);

% system
% x = Ax+B(u-Ed);
A = zeros(2);
B = eye(2);
E = [-1; 0];
% Constants
Lx = 6; % meters
Ly = 3; % meters
Lx = 6; % meters
Ly = 3; % meters
umax = 20; % m/s, longitudinal velocity
umin = 10; % m/s
um = 15;
dumax = 16; % m/s
dumin = 14; % m/s
vmax = 1; % m/s, lateral velocity

% Variables
% let x = x_e - x_o, x is x diff
% let y = y_e - y_o, h is y diff
% u = -Kx
sdpvar px py d 
x = [px; py];
%fx = A*x + B*(u+E*d);

% minimize k to maximize the volume of B >= 0.
sdpvar k
[B_hat, coefB, mon] = polynomial(x,2); % B_hat is sos
Barrier = k - B_hat; % >= 0 means safe
dB = jacobian(Barrier,x);
% Safeset
safe_set = px^2/(Lx/2)^2 + py^2/(Ly/2)^2 - 1;
safe_set1 = 1 - (px)^2/(10)^2 - (py-6)^2/(4)^2; % >= 0 real safe
%safe_set2 = 1 - (px+7)^2/(20)^2 - (py-1)^2/(8)^2; % >= 0 bounded region

% initial
init_set = 3^2 - (px+17)^2 - (py-1.5)^2;

initial_X = safe_set1 - 0.99;
K = -lqr(A,B,0.1*eye(2),100*eye(2));
gamma = 40;
eps = 1e-7;

% Non-empty set constraint
[s0, c0, m0] = polynomial(x,2);
nept_const = Barrier - s0*initial_X;

% initial set constraint
%[s1, c1, m1] = polynomial(x,2);
%init_const = Barrier - s1*init_set;

% Contained in safety set
[s21, c21, m21] = polynomial(x,2);
%[s22, c22, m22] = polynomial(x,2);
safe_const1 = -(Barrier + eps) + s21*safe_set1; % >= 0, such that B in safe, or unsafe in neg(B)
%safe_const2 = -(Barrier + eps) + s22*safe_set2; % >= 0, such that B in safe, or unsafe in neg(B)

% Control barrier constraint
[s3, c3, m3] = polynomial([x;d],2);
[s4, c4, m4] = polynomial([x;d],2);
xd = [0;6];
den1 = 1 + ((K(1,:)*(x-xd)-um)/(2*umax))^2;
den2 = 1 + ((K(2,:)*(x-xd))/(2*vmax))^2;
control_const = (dB * (A*x + B*(um+E*d)) * den1*den2 ...
                + dB * B * [(K(1,:)*(x-xd))*den2; K(2,:)*(x-xd)*den1] ...
                + gamma*Barrier * den1*den2 ) ...
                - (s3*safe_set1 + s4*(dumax-d)*(d-dumin));
            
% constraint = [ k>=0; sos(nept_const); sos(init_const); sos(safe_const); sos(control_const);...
%               sos(B_hat); sos(s0); sos(s1); sos(s2);sos(s3);sos(s4);];% sos(s3); sos(s4)];
          
constraint = [ k>=0; sos(nept_const); sos(safe_const1);sos(control_const)...
               sos(B_hat); sos(s0); sos(s21); sos(s3); sos(s4)];

obj = [];
variables = [k;coefB;c0;c21;c3;c4];%;c5];
%variables = [k;coefB];
[sol,v,Q,res] = solvesos(constraint, obj, opts, variables);

%% plotting
figure
range = [-35 20 -20 35];
hold on
fcontour(str2sym(sdisplay(safe_set)), range, '-r', 'LevelList', [0])
fcontour(str2sym(sdisplay(safe_set1)), range, '-k', 'LevelList', [0])
%fcontour(str2sym(sdisplay(init_set)), range, '-k', 'LevelList', [0])
fcontour(str2sym(sdisplay(initial_X)), range, '-g', 'LevelList', [0])
grid on;
axis equal;
% Substitute variables with values
B_plot = replace(Barrier, [k; coefB], [value(k); value(coefB)]);
fb = matlabFunction(str2sym(sdisplay(B_plot)));
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

while true
% for idx = 1
    old_k = value(k);
    %% hold barrier
    old_Bhat = replace(B_hat, coefB, value(coefB));
    expand_const = k - old_k >= 0;
    
    [s01, c01, m01] = polynomial(x,2);
    safe_const1 = -(k - old_Bhat + eps) + s01*safe_set1;

    
    dB = jacobian(k - old_Bhat, x);
    [ux_input, cx_input, mx_input] = polynomial([x;d],2);
    [uy_input, cy_input, my_input] = polynomial([x;d],2);
    
    [s11, c11, m11] = polynomial([x;d],2);
    [s12, c12, m12] = polynomial([x;d],2);
    control_const = (dB * (A*x + B*([ux_input;uy_input] + E*d))...
                    + gamma*(k-old_Bhat) + eps)...
                    - (s11*safe_set1 + s12*(dumax-d)*(d-dumin));
    
    [s21, c21, m21] = polynomial([x;d],2);
    [s22, c22, m22] = polynomial([x;d],2);
    ux_const1 = (ux_input - umin) - (s21*safe_set1 + s22*(dumax-d)*(d-dumin));
    
    [s31, c31, m31] = polynomial([x;d],2);
    [s32, c32, m32] = polynomial([x;d],2);
    ux_const2 = (-ux_input + umax) - (s31*safe_set1 + s32*(dumax-d)*(d-dumin));
    
    [s41, c41, m41] = polynomial([x;d],2);
    [s42, c42, m42] = polynomial([x;d],2);
    uy_const1 = (uy_input + vmax) - (s41*safe_set1 + s42*(dumax-d)*(d-dumin));
    
    [s51, c51, m51] = polynomial([x;d],2);
    [s52, c52, m52] = polynomial([x;d],2);
    uy_const2 = (-uy_input + vmax) - (s51*safe_set1 + s52*(dumax-d)*(d-dumin));
    
    constraint = [k >= 0; sos(safe_const1); sos(control_const);...
%                   sos(ux_const1);sos(ux_const2);...
%                   sos(uy_const1);sos(uy_const2);...
                  sos(s01); sos(s11); sos(s12)];%;...
%                   sos(s21); sos(s22);...
%                   sos(s31); sos(s32);...
%                   sos(s41); sos(s42);...
%                   sos(s51); sos(s52)];
    obj = [-k];
    variables = [k;c01;c11;c12;c21;c22;c31;c32;c41;c42;c51;c52];
    variables = [k;c01;c11;c12];
    [sol,v,Q,res] = solvesos(constraint, obj, opts, variables);
    
    figure
    range = [-35 20 -20 35];
    hold on
    fcontour(str2sym(sdisplay(safe_set)), range, '-r', 'LevelList', [0])
    fcontour(str2sym(sdisplay(safe_set1)), range, '-k', 'LevelList', [0])
    %fcontour(str2sym(sdisplay(safe_set2)), range, '-r', 'LevelList', [0])
    %fcontour(str2sym(sdisplay(init_set)), range, '-k', 'LevelList', [0])
    fcontour(str2sym(sdisplay(initial_X)), range, '-g', 'LevelList', [0])
    grid on;
    axis equal;
    % Substitute variables with values
    B_plot = replace(k-B_hat, [k; coefB], [value(k); value(coefB)]);
    fb = matlabFunction(str2sym(sdisplay(B_plot)));
    fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

    %% Hold control u
    old_k = value(k);
    old_Bhat = replace(B_hat, coefB, value(coefB));
    [B_hat, coefB, mon] = polynomial(x,2);
    
    [s0, c0, m0] = polynomial(x,2);
    safe_const = -(old_Bhat - k) + s0*safe_set - eps;
    
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






