function dx = quad_3D_Dynamics_Kumar_Jelle(state,input,data)
% state x = [x,y,z,phi,theta,psi,xdot,ydot,zdot,omega_x,omega_y,omega_z]^T 
if not(isfield(data,'rotorcontrol'))
    data.rotorcontrol=0;
end

pos=state(1:3);
euler=state(4:6);
phi=euler(1);
theta=euler(2);
psi=euler(3);
vel=state(7:9);
omega=state(10:12); %rotational rate in body frame
Tbi=get_rotationmatrix([-1 -1 1].*euler','I2B');
R=get_rotationmatrix(euler','B2E'); %transformation body to world

% here i take a different take on transformations than the paper
% the paper doesn't use inertial reference frame and uses the earth to body
% transform instead
% Omega=Tbi*omega; %rotational rate in body rf (wx,wy,wz)

if data.rotorcontrol
    omega_rotor=input; 
    
    u=[data.k_F, data.k_F, data.k_F, data.k_F;
        0, data.k_F*data.L, 0, -data.k_F*data.L;
        -data.k_F*data.L, 0 data.k_F*data.L, 0 ;
        data.k_M, -data.k_M, data.k_M, -data.k_M]*omega_rotor.^2; 
else
    u=input;
end

I=data.I; 
domega=I^(-1)*(cross(-omega,I*omega)+u(2:4)); %rot acc in body frame
% dOmega=domega; %rot acc in inertial frame 

omega_skew=[0, -omega(3), omega(2);
            omega(3), 0, -omega(1);
            -omega(2), omega(1), 0];

% three different methods to go from body rotational rate to inertial rotational rate (i think dTheta3 is correct)
% the symmetric skew matrix method is adapted from http://arxiv.org/abs/1712.02402 (at one point they take the matrix multiplication which leaves 3x3 matrix and at another point they seem to take the dot product which stays 0 if starting from 0 euler
dTheta=Tbi*omega; %dphi dtheta dpsi (p,q,r)
dTheta2=(dot(R,omega_skew))';
dTheta3=R*omega_skew;
dTheta3=[dTheta3(3,2);dTheta3(1,3);dTheta3(2,1)];%(dot(R,omega_skew))';    
dTheta4=R*omega;
dpos=vel; 
dvel=(1/data.m)*([0;0;-data.m*data.g]+(R*[0;0;u(1)])); %acceleration in world frame (NWU)

dx=[dpos; dTheta3; dvel; domega];


end