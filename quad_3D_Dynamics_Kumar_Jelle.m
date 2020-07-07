function dx = quad_3D_Dynamics_Kumar_Jelle(state,input,data)
% state x = [x,y,z,phi,theta,psi,xdot,ydot,zdot,p,q,r]^T 
if not(isfield(data,'rotorcontrol'))
    data.rotorcontrol=0;
end

pos=state(1:3);
euler=state(4:6);
vel=state(7:9);
omega=state(10:12); %rotational rate in inertial rf (p,q,r)
Tbi=get_rotationmatrix(euler','I2B');
Teb=get_rotationmatrix(euler','B2E'); %transformation body to world

% here i take a different take on transformations than the paper
% the paper doesn't use inertial reference frame and uses the earth to body
% transform instead
Omega=Tbi*omega; %rotational rate in body rf (wx,wy,wz)

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

dOmega=I^(-1)*(cross(-Omega,I*Omega)+u(2:4)); %rot acc in body frame
domega=(Tbi')*dOmega; %rot acc in inertial frame 
dTheta=omega;

dpos=vel; 
dvel=(1/data.m)*([0;0;data.m*data.g]-(Teb*[0;0;u(1)])); %acceleration in world frame

dx=[dpos; dTheta; dvel; domega];


end