clear all
close all

% any set function h(x) is assumed to be {x|h(x)>0}

% Assume straight lane profile
% Model 2
% double integrator 4 states
% 	d | px_e - px_o |   | vx_e - vx_o |
% --- | py_e - py_o |   | vy_e - vy_o |
%     | vx_e - vx_o | = |    ux - d   |
%  dt | vy_e - vy_o |   |      uy     |


% Bounded u1, u1, d
% ux in [-uxmax, uxmax] uxmin = -uxmax
% uy in [-uymax, uymax] uymin = -uymax
% d  in [-dumax, dumax] dumin = -dumax
% vx_e-vx_o in [-vxmax, vxmax] vxmin = -vxmax
% vy_e-vy_o in [-vymax, vymax] vymin = -vymax
% obstacle vehicle state in lane

% system
% x = Ax+Bu+Ed;
A = [zeros(2),eye(2);zeros(2),zeros(2)];
B = [zeros(2);eye(2)];
E = [0;0;-1;0];

opts = [];
opts = sdpsettings('solver','sedumi','sos.model',2,'sos.scale',1,'verbose',1,...
'sos.newton',1,'sos.congruence',1);

% Constants
Lx = 6; % meters
Ly = 3; % meters
uxmax = 50; % m/s^2, longitudinal acceleration
uxmin = -50; % m/s^2
uymax = 10; % m/s^2, lateral acceleration
uymin = -10; % m/s^2
dumax = 1; % m/s^2
dumin = -1; % m/s^2
vxmax = 5; % m/s, longitudinal velocity
vymax = 1; % m/s, lateral velocity

% Variables
% let px = px_e - px_o
% let py = py_e - py_o
% let vx = vx_e - vx_o
% let vy = vy_e - vy_o
% u = Kx
sdpvar px py vx vy ux uy d 
x = [px; py; vx; vy];


% maxmize k to maximize the volume of B >= 0.
sdpvar k
[B_hat, coefB, mon] = polynomial(x,4); % B_hat is sos
Barrier = k - B_hat; % >= 0 means invarance inside barrier
dB = jacobian(Barrier,x);

% Safeset
safe_set = px^2/(Lx/2)^2 + py^2/(Ly/2)^2 - 1; % >= 0 obstacle vehcile
safe_set1 = 1 - (px)^2/(10)^2 - (py-6)^2/(4)^2; % >= 0 real safe set
safe_set2 = (vxmax-vx)*(vx+vxmax); % >= 0 x velocity bound 
safe_set3 = (vymax-vy)*(vy+vymax); % >= 0 y velocity bound

% nonempty barrier
nept_set1 = safe_set1 - 0.99; % >= 0 position tiny set to remove trivial empty solution
nept_set2 = 0.1^2-vx^2; % x velocity
nept_set3 = 0.1^2-vy^2; % y velocity

% parameter
K = -lqr(A,B,eye(4),1*eye(2)); % lqr gain
gamma = 500; % CBF rate
eps = 1e-7; % samlm constant

% Non-empty set constraint
[s01, c01, m01] = polynomial(x,2);
[s02, c02, m02] = polynomial(x,2);
[s03, c03, m03] = polynomial(x,2);
nept_const = Barrier - (s01*nept_set1 + s02*nept_set2 + s03*nept_set1);

% Contained in safety set, (s1 is skipped)
[s21, c21, m21] = polynomial(x,2);
[s22, c22, m22] = polynomial(x,2);
[s23, c23, m23] = polynomial(x,2);
safe_const1 = -(Barrier + eps) + s21*safe_set1; % >= 0, such that B in safe1
safe_const2 = -(Barrier + eps) + s22*safe_set2; % >= 0, such that B in safe2
safe_const3 = -(Barrier + eps) + s23*safe_set3; % >= 0, such that B in safe3

% Control barrier constraint
[s31, c31, m31] = polynomial([x;d],2);
[s32, c32, m32] = polynomial([x;d],2);
[s33, c33, m33] = polynomial([x;d],2);
[s4, c4, m4] = polynomial([x;d],2);
xd = [0;6;0;0];
% % with control ocnstraints
% den1 = 1 + (K(1,:)*(x-xd)/(2*uxmax))^2;
% den2 = 1 + (K(2,:)*(x-xd)/(2*uymax))^2;
% u = K*(x-xd);
% control_const = (dB * (A*x + B*E*d)*den1*den2 ...
%                 +dB*B*diag([den2,den1])*u...
%             + gamma*Barrier * den1*den2 ) ...
%             - (s31*safe_set1 + s32*safe_set2 + s33*safe_set3 + s4*(dumax-d)*(d-dumin));

% without control constraints
u = K*(x-xd);
control_const = (dB*(A*x + B*u + E*d) + gamma*Barrier) ...
                - (s31*safe_set1 + s32*safe_set2 + s33*safe_set3 + s4*(dumax-d)*(d-dumin));

constraint = [ k>=0; sos(nept_const); sos(safe_const1); sos(safe_const2); sos(safe_const3);...
               sos(control_const);...
               sos(B_hat); sos(s01); sos(s02); sos(s03);...
               sos(s21); sos(s22); sos(s23);...
               sos(s31); sos(s32); sos(s33); sos(s4)];

obj = [];
variables = [k;coefB;c01;c02;c03;c21;c22;c23;c31;c32;c33;c4];
[sol,v,Q,res] = solvesos(constraint, obj, opts, variables);

%% plotting
figure
range = [-35 20 -20 35];
hold on
fcontour(str2sym(sdisplay(safe_set)), range, '-r', 'LevelList', [0])
fcontour(str2sym(sdisplay(safe_set1)), range, '-k', 'LevelList', [0])
fcontour(str2sym(sdisplay(nept_set1)), range, '-g', 'LevelList', [0])
grid on;
axis equal;
% Substitute variables with values
B_plot = replace(Barrier, [k; coefB; vx; vy], [value(k); value(coefB); 0; 0]);
fb = matlabFunction(str2sym(sdisplay(B_plot)));
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])


%% start iterating
old_k = 0;
while (value(k)-old_k)>1e-3
    % for idx = 1
    old_k = value(k);
    %% hold barrier
    old_Bhat = replace(B_hat, coefB, value(coefB));
    expand_const = k - old_k >= 0;
    
    Barrier2 = k - old_Bhat;

    [s01, c01, m01] = polynomial(x,2);
    [s02, c02, m02] = polynomial(x,2);
    [s03, c03, m03] = polynomial(x,2);
    safe_const1 = - (Barrier2 + eps) + s01*safe_set1;
    safe_const2 = - (Barrier2 + eps) + s02*safe_set2;
    safe_const3 = - (Barrier2 + eps) + s03*safe_set3;

    dB = jacobian(Barrier2, x);
    [ux_input, cx_input, mx_input] = polynomial(x,2);
    [uy_input, cy_input, my_input] = polynomial(x,2);

    [s11, c11, m11] = polynomial([x;d],2);
    [s12, c12, m12] = polynomial([x;d],2);
    [s13, c13, m13] = polynomial([x;d],2);
    u = [ux_input; uy_input];
    control_const = -(dB * (A*x + B*u + E*d) + gamma*Barrier)...
                    + (s11*safe_set1 + s12*safe_set2 + s13*safe_set3 + s12*(dumax-d)*(d-dumin));

    [s21, c21, m21] = polynomial([x;d],2);
    [s22, c22, m22] = polynomial([x;d],2);
    [s23, c23, m23] = polynomial([x;d],2);
    [s24, c24, m24] = polynomial([x;d],2);
    ux_const1 = (ux_input - uxmin) - (s21*safe_set1 + s22*safe_set2+ s23*safe_set3 + s24*(dumax-d)*(d-dumin));

    [s31, c31, m31] = polynomial([x;d],2);
    [s32, c32, m32] = polynomial([x;d],2);
    [s33, c33, m33] = polynomial([x;d],2);
    [s34, c34, m34] = polynomial([x;d],2);
    ux_const2 = (-ux_input + uxmax) - (s31*safe_set1 + s32*safe_set2 + s33*safe_set3 + s34*(dumax-d)*(d-dumin));

    [s41, c41, m41] = polynomial([x;d],2);
    [s42, c42, m42] = polynomial([x;d],2);
    [s43, c43, m43] = polynomial([x;d],2);
    [s44, c44, m44] = polynomial([x;d],2);
    uy_const1 = (uy_input - uymin) - (s41*safe_set1 + s42*safe_set2 + s43*safe_set3 + s44*(dumax-d)*(d-dumin));

    [s51, c51, m51] = polynomial([x;d],2);
    [s52, c52, m52] = polynomial([x;d],2);
    [s53, c53, m53] = polynomial([x;d],2);
    [s54, c54, m54] = polynomial([x;d],2);
    uy_const2 = (-uy_input + uymax) - (s51*safe_set1 + s52*safe_set2 + s53*safe_set3 + s54*(dumax-d)*(d-dumin));

    constraint = [k >= 0; expand_const; sos(safe_const1); sos(safe_const2); sos(safe_const3); sos(control_const);...
                  sos(ux_const1);sos(ux_const2);...
                  sos(uy_const1);sos(uy_const2);...
                  sos(s01); sos(s02); sos(s03)...
                  sos(s11); sos(s12); sos(s13)...
                  sos(s21); sos(s22); sos(s23); sos(s24)...
                  sos(s31); sos(s32); sos(s33); sos(s34)...
                  sos(s41); sos(s42); sos(s43); sos(s44)...
                  sos(s51); sos(s52); sos(s53); sos(s54)];
    % maximize k to enlarge area of Barrier
    obj = [-k];
    variables = [k;c01;c02;c03;c11;c12;c13;...
                 c21;c22;c23;c24;
                 c31;c32;c33;c34;
                 c41;c42;c43;c44;
                 c51;c52;c53;c54;];
    [sol,v,Q,res] = solvesos(constraint, obj, opts, variables);

    figure
    range = [-35 20 -20 35];
    hold on
    fcontour(str2sym(sdisplay(safe_set)), range, '-r', 'LevelList', [0])
    fcontour(str2sym(sdisplay(safe_set1)), range, '-k', 'LevelList', [0])
    fcontour(str2sym(sdisplay(nept_set1)), range, '-g', 'LevelList', [0])
    grid on;
    axis equal;
    % Substitute variables with values
    B_plot = replace(Barrier, [k], [value(k)]);
    fb = matlabFunction(str2sym(sdisplay(B_plot)));
    fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

    %% Hold control u
    old_k = value(k);
    old_Bhat = replace(B_hat, coefB, value(coefB));
    [B_hat, coefB, mon] = polynomial(x,2);
    
    Barrier = B_hat - k;
    dB = jacobian(Barrier,x);
    
    [s01, c01, m01] = polynomial(x,2);
    [s02, c02, m02] = polynomial(x,2);
    [s03, c03, m03] = polynomial(x,2);
    safe_const1 = - (Barrier2 + eps) + s01*safe_set1;
    safe_const2 = - (Barrier2 + eps) + s02*safe_set2;
    safe_const3 = - (Barrier2 + eps) + s03*safe_set3;

    [s1, c1, m1] = polynomial(x,2);
    expand_const = Barrier - s1*(old_Bhat - old_k);
    
    % hold control input
    ux_input = replace(ux_input, [cx_input], value(cx_input));
    uy_input = replace(uy_input, [cy_input], value(cy_input));
    [s11, c11, m11] = polynomial([x;d],4);
    [s12, c12, m12] = polynomial([x;d],4);
    [s13, c13, m13] = polynomial([x;d],4);
    [s2, c2, m2] = polynomial([x;d],2);
    u = [ux_input; uy_input];
    control_const = (dB*(A*x + Bu + E*d) + gamma*Barrier) ...
                    - (s11*safe_set1 + s12*safe_set2 + s13*safe_set3 + s2*(dumax-d)*(d-dumin));

    variables = [k;coefB;c01;c01;c02;c11;c12;c13;c2];

    constraint = [k >= 0; sos(expand_const); sos(safe_const1); sos(safe_const2); sos(safe_const3); sos(control_const);...
                  sos(B_hat); sos(s01); sos(s02); sos(s03);...
                  sos(s11); sos(s12); sos(s13); sos(s2)];
end
%% plotting
figure
range = [-35 20 -20 35];
hold on
fcontour(str2sym(sdisplay(safe_set)), range, '-r', 'LevelList', [0])
fcontour(str2sym(sdisplay(safe_set1)), range, '-k', 'LevelList', [0])
fcontour(str2sym(sdisplay(nept_set1)), range, '-g', 'LevelList', [0])
grid on;
axis equal;
% Substitute variables with values
B_plot = replace(Barrier, [k; coefB; vx; vy], [value(k); value(coefB); 0; 0]);
fb = matlabFunction(str2sym(sdisplay(B_plot)));
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])





