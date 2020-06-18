function ds = quad3d_dynamics(x,u)

% Define inputs
T = u(1);
omega = u(2:4);

% Parameters
m = 0.38905;
beta = 0.5;
g = 9.81;

% Define states
% X = x(:,1:3); % Position vector
V = x(4:6); % Velocity vector
q0 = x(7);
q1 = x(8);
q2 = x(9);
q3 = x(10);

% Body to Earth matrix using quaternions:
c11 = q0.^2 + q1.^2 - q2.^2 - q3.^2;
c12 = 2*q1.*q2 - 2*q0.*q3;
c13 = 2*q0.*q2 + 2*q1.*q3;
c21 = 2*q0.*q3 + 2*q1.*q2;
c22 = q0.^2 - q1.^2 + q2.^2 - q3.^2;
c23 = 2*q2.*q3 - 2*q0.*q1;
c31 = 2*q1.*q3 - 2*q0.*q2;
c32 = 2*q2.*q3 + 2*q0.*q1;
c33 = q0.^2 - q1.^2 - q2.^2 + q3.^2;
C = [
    c11, c12, c13;
    c21, c22, c23;
    c31, c32, c33;
    ];

C_vec = [
        c13; 
        c23;
        c33;
        ];
    
r11 = 0; % = r22, r33 and r44
p = omega(1);
q = omega(2);
r = omega(3);
R = [
    r11, -p, -q, -r;
    p, r11, r, -q;
    q, -r, r11, p;
    r, q, -p, r11;
    ];

% Diff. equations: 
% -- Position
dX = V;
% -- Velocity
T_B = (T + m*g)/m;
% Drag
K_mat = [beta, 0, 0; 0, beta, 0; 0, 0, 0];
% C. aux and dV
T_E = C_vec*T_B;
D_E = C*(-K_mat)*abs(C*dX);
dV = [0; 0; -g] + T_E + D_E; 
% -- Angles
dQ = (1/2)*R*[q0, q1, q2, q3]';

%Define ODE right-hand side
ds = [dX; dV; dQ];

end