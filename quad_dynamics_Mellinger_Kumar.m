function ds = quad_dynamics_Mellinger_Kumar(states,controls,data)

% Parameters
% -- vehicle: 
m = data.m; % [kg] quadrotor mass
I = data.I; % [kg.m^2] quadrotor moment of inertia
L = data.L; % [m] (assuming + configuration) distance from axis of rotation of rotors to quad. CoM
k_F = data.k_F; % [?] rotor force coeficient 
k_M = data.k_M; % [?] rotor moment coeficient
% -- environment: 
g = data.g; % [m/s^2] gravity

% Define states
% X = states(:,1:3); % [m] Position vector
R = states(4:6); % [rad] Orientation vector (Euler)
V = states(7:9); % [m/s] Velocity vector
omega = states(10:12); % [rad/s] Angular velocity vector

% Calculate thrust and angular velocities:
G = [
    k_F, k_F, k_F, k_F;
    0, k_F*L, 0, -k_F*L;
    -k_K*L, 0, k_F*L, 0;
    k_M, -k_M, k_M, -k_M;
     ];
u = G*controls;

% Unit vectors:
z_W = [0, 0, 1]'; % z in the world frame. 

x_B = [ 
    cos(theta)*cos(psi); 
    sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi); 
    cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)
    ];
y_B = [
    cos(theta)*sin(psi);
    sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi);
    cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)
    ];
z_B = [
    -sin(theta);
    sin(phi)*cos(theta);
    cos(phi)*cos(theta);
    ];

% Direction cosine matrix (transformation matrix from Body to Euler)
R = [x_B, y_B, z_B];

% Angular velocity in skew-symmetric matrix form (usefull, for instance, for doing cross product)
omega_mat = [
    0, -omega(3), omega(2);
    omega(3), 0, -omega(1);
    -omega(2), omega(1), 0
    ];

% Diff. equations: 
% -- Positions
dX = V;
% -- Angles
dR = R*omega_mat; % R is the orientation
% -- Velocities
dV = (1/m)*(-m*g*z_W + u(1)*z_B);
% -- Angular velocities
domega = I^(-1)*(-omega_mat*I*omega + [u(2), u(3), u(4)]');

%Define ODE derivative vector:
ds = [dX; dR; dV; domega];

end