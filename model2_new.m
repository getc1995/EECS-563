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
A = [zeros(2),eye(2);zeros(2),zeros(2)];
B = [zeros(2);eye(2)];
E = [-1;0];

% Constants
Lx = 6; % meters
Ly = 3; % meters
uxmax = 50; % m/s, longitudinal acceleration
uxmin = -50; % m/s
uymax = 10; % m/s, lateral acceleration
uymin = -10; % m/s
dumax = 1; % m/s
dumin = -1; % m/s
vxmax = 5; % m/s, longitudinal velocity
vymax = 1; % m/s, lateral velocity

% Variables
% let x = x_e - x_o, x is x diff
% let y = y_e - y_o, h is y diff
% u = sKx
sdpvar px py vx vy ux uy d 
x = [px; py; vx; vy];
% u = [ux;uy];
% fx = A*x + B*(u+E*d);

% minimize k to maximize the volume of B >= 0.
sdpvar k
[B_hat, coefB, mon] = polynomial(x,2); % B_hat is sos
Barrier = k - B_hat; % >= 0 means safe
dB = jacobian(Barrier,x);
% Safeset
safe_set = px^2/(Lx/2)^2 + py^2/(Ly/2)^2 - 1;
%safe_set1 = 1 - (px)^2/(10)^2 - (py-6)^2/(4)^2; % >= 0 real safe
safe_set1 = 1 - (py-6)^2/(4)^2; % >= 0 real safe
safe_set2 = (vxmax-vx)*(vx+vxmax); % >= 0 velocity bound 
safe_set3= (vymax-vy)*(vy+vymax); % >= 0 velocity bound
%safe_set2 = 1 - (px+7)^2/(20)^2 - (py-1)^2/(8)^2; % >= 0 bounded region

% initial
init_set = 3^2 - (px+17)^2 - (py-1.5)^2;

initial_X = safe_set1 - 0.99;
K = -lqr(A,B,eye(4),1*eye(2));
gamma = 50;
eps = 1e-7;

% Non-empty set constraint
[s01, c01, m01] = polynomial(x,2);
[s02, c02, m02] = polynomial(x,2);
[s03, c03, m03] = polynomial(x,2);
nept_const = Barrier - (s01*initial_X + s02*(0.1^2-vx^2) + s03*(0.1^2-vy^2));

% initial set constraint
%[s1, c1, m1] = polynomial(x,2);
%init_const = Barrier - s1*init_set;

% Contained in safety set
[s21, c21, m21] = polynomial(x,2);
[s22, c22, m22] = polynomial(x,2);
[s23, c23, m23] = polynomial(x,2);
safe_const1 = -(Barrier + eps) + s21*safe_set1; % >= 0, such that B in safe, or unsafe in neg(B)
safe_const2 = -(Barrier + eps) + s22*safe_set2; % >= 0, such that B in safe, or unsafe in neg(B)
safe_const3 = -(Barrier + eps) + s23*safe_set3; % >= 0, such that B in safe, or unsafe in neg(B)

% Control barrier constraint
[s31, c31, m31] = polynomial([x;d],2);
[s32, c32, m32] = polynomial([x;d],2);
[s33, c33, m33] = polynomial([x;d],2);
[s4, c4, m4] = polynomial([x;d],2);
xd = [0;6;0;0];
den1 = 1 + (K(1,:)*(x-xd)/(2*uxmax))^2;
den2 = 1 + (K(2,:)*(x-xd)/(2*uymax))^2;
u = K*(x-xd);
% control_const = (dB * (A*x + B*E*d) * den1*den2 ...
%                 + dB * B * [K(1,:)*(x-xd)*den2; K(2,:)*(x-xd)*den1] ...
%                 + gamma*Barrier * den1*den2 ) ...
%                 - (s31*safe_set1 + s4*(dumax-d)*(d-dumin));
            
control_const = (dB * (A*x + B*(u./[den1;den2]+E*d)) * den1*den2 ...
                + gamma*Barrier * den1*den2 ) ...
                - (s31*safe_set1 + s32*safe_set2 + s33*safe_set3 + s4*(dumax-d)*(d-dumin));
            
% control_const = (dB * (A*x + B*(u./[den1;den2])) * den1*den2 ...
%                 + gamma*Barrier * den1*den2 ) ...
%                 - (s3*safe_set1);
            
% constraint = [ k>=0; sos(nept_const); sos(init_const); sos(safe_const); sos(control_const);...
%               sos(B_hat); sos(s0); sos(s1); sos(s2);sos(s3);sos(s4);];% sos(s3); sos(s4)];
          
constraint = [ k>=0; sos(nept_const); sos(safe_const1); sos(safe_const2); sos(safe_const3);...
               sos(control_const);...
               sos(B_hat); sos(s01); sos(s02); sos(s03); sos(s21); sos(s22); sos(s23);...
               sos(s31); sos(s32); sos(s33); sos(s4)];

% constraint = [ k>=0; sos(nept_const); sos(safe_const1); sos(control_const);...
%                sos(B_hat); sos(s0); sos(s21); sos(s3)];

% constraint = [ k>=0; sos(nept_const); sos(safe_const1);...
%                sos(B_hat); sos(s0); sos(s21)];


obj = [];
variables = [k;coefB;c01;c02;c03;c21;c22;c23;c31;c32;c33;c4];%;c5];
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
B_plot = replace(Barrier, [k; coefB; vx; vy], [value(k); value(coefB); 0; 0]);
fb = matlabFunction(str2sym(sdisplay(B_plot)));
fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

while true
% for idx = 1
    old_k = value(k);
    %% hold barrier
    old_Bhat = replace(B_hat, coefB, value(coefB));
    expand_const = k - old_k >= 0;
    
    [s01, c01, m01] = polynomial(x,2);
    [s02, c02, m02] = polynomial(x,2);
    safe_const1 = -(k - old_Bhat + eps) + s01*safe_set1;
    safe_const2 = -(k - old_Bhat + eps) + s02*safe_set2;
    
    dB = jacobian(k - old_Bhat, x);
    [ux_input, cx_input, mx_input] = polynomial(x,2);
    [uy_input, cy_input, my_input] = polynomial(x,2);
    
    [s11, c11, m11] = polynomial([x;d],2);
    [s12, c12, m12] = polynomial([x;d],2);
    control_const = -(dB * (A*x + B*([ux_input;uy_input] + E*d))...
                    + gamma*Barrier + eps)...
                    + (s11*safe_set1 + s12*(dumax-d)*(d-dumin));
    
    [s21, c21, m21] = polynomial([x;d],2);
    [s22, c22, m22] = polynomial([x;d],2);
    ux_const1 = (ux_input - uxmin) - (s21*safe_set2 + s22*(dumax-d)*(d-dumin));
    
    [s31, c31, m31] = polynomial([x;d],2);
    [s32, c32, m32] = polynomial([x;d],2);
    ux_const2 = (-ux_input + uxmax) - (s31*safe_set2 + s32*(dumax-d)*(d-dumin));
    
    [s41, c41, m41] = polynomial([x;d],2);
    [s42, c42, m42] = polynomial([x;d],2);
    uy_const1 = (uy_input - uymin) - (s41*safe_set2 + s42*(dumax-d)*(d-dumin));
    
    [s51, c51, m51] = polynomial([x;d],2);
    [s52, c52, m52] = polynomial([x;d],2);
    uy_const2 = (-uy_input + uymax) - (s51*safe_set2 - s52*(dumax-d)*(d-dumin));
    
    constraint = [k >= 0; expand_const; sos(safe_const1); sos(safe_const2); sos(control_const);...
                  sos(ux_const1);sos(ux_const2);...
                  sos(uy_const1);sos(uy_const2);...
                  sos(s01); sos(s02); sos(s11); sos(s12);...
                  sos(s21); sos(s22);...
                  sos(s31); sos(s32);...
                  sos(s41); sos(s42);...
                  sos(s51); sos(s52)];
    obj = [-k];
    variables = [k;c01;c02;c11;c12;c21;c22;c31;c32;c41;c42;c51;c52];
    [sol,v,Q,res] = solvesos(constraint, obj, opts, variables);
    
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
    B_plot = replace(k-B_hat, [k; coefB; vx; vy], [value(k); value(coefB); 0; 0]);
    fb = matlabFunction(str2sym(sdisplay(B_plot)));
    fcontour(str2sym(sdisplay(B_plot)), range, '-b', 'LevelList', [0])%, [-50 50 -50 50])

    %% Hold control u
    old_k = value(k);
    old_Bhat = replace(B_hat, coefB, value(coefB));
    [B_hat, coefB, mon] = polynomial(x,2);
    
    [s01, c01, m01] = polynomial(x,2);
    safe_const = -(old_Bhat - k) + s01*safe - eps;
    
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
                
    variables = [k;coefB;c01;c11;c12;];
    
    constraint = [k >= 0; sos(expand_const); sos(safe_const); sos(control_const);...
                  sos(B_hat); sos(s01); sos(s1); sos(s11); sos(s12);...
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





