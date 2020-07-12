function [u, extra] = controller(state,reference,data)

%% get desired acceleration
a_ref=reference(1,1:3,1,3);
a_fb=data.K_pos*(reference(1,1:3,1,1)'-state(1:3))+data.K_vel*(reference(1,1:3,1,2)'-state(7:9));
a_des=(a_fb'+a_ref)'+[0;0;data.g];

%% get orientation
phi=state(4);
theta=state(5);
psi=state(6);
x_b = [ 
    cos(theta)*cos(psi); 
    cos(theta)*sin(psi);
    -sin(theta);
    ];

y_b = [
    sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi); 
    sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi);
    sin(phi)*cos(theta)
    ];

z_b = [
    cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
    cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
    cos(phi)*cos(theta);
    ];

Rb=[x_b,y_b,z_b];

%% Get desired orientation
if norm(a_des,2)==0    
    z_b_des = z_b;
else
    z_b_des=a_des./norm(a_des,2);
end
y_c=[-sin(psi),cos(psi),0]';
x_b_des=cross(y_c,z_b_des);
x_b_des=x_b_des./norm(x_b_des,2);
y_b_des=cross(z_b_des,x_b_des);
R_des=[x_b_des,y_b_des,z_b_des]; %desired orientation


% test2=[atan2(Rb(3,2),Rb(3,3)); atan2(-Rb(3,1),norm([Rb(3,2),Rb(3,3)],2)); atan2(Rb(2,1),Rb(1,1))];
%% Get normalized thrust command
    c_cmd=a_des'*z_b;
    extra.a_des=a_des;
%% Get desired body angular rate    
    omega_ref=data.omega;
    euler_des= [atan2(R_des(3,2),R_des(3,3)); atan2(-R_des(3,1),norm([R_des(3,2),R_des(3,3)],2)); atan2(R_des(2,1),R_des(1,1))];
    euler_error=(euler_des-state(4:6));
    if sum(euler_error(euler_error>pi))||(sum(euler_error(euler_error<-pi)))
        dummy=2;
    end
    euler_error(euler_error>pi)=euler_error(euler_error>pi)-2*pi;
    euler_error(euler_error<-pi)=euler_error(euler_error<-pi)+2*pi;
    
    %% using the inverse method of http://arxiv.org/abs/1712.02402 to convert from dTheta to domega
    d_Theta_fb=data.K_theta*euler_error; 
    d_Theta_skew=[0, -d_Theta_fb(3), d_Theta_fb(2);
                  d_Theta_fb(3), 0, -d_Theta_fb(1);
                  -d_Theta_fb(2), d_Theta_fb(1),0];
    omega_fb_skew=Rb'*d_Theta_skew;
    omega_fb=[omega_fb_skew(3,2);omega_fb_skew(1,3);omega_fb_skew(2,1)];
    omega_des=omega_ref+omega_fb;
    extra.omega_des=omega_des;
    extra.euler_des=euler_des;
%% Get desired body angular acceleration
    omega_dot_ref=data.omega_dot;
    omega_dot_fb=data.K_omega*(omega_des-state(10:12));
    omega_dot_des=omega_dot_fb+omega_dot_ref;
    extra.omega_dot_des=omega_dot_des;

%$ get desired torque:
tau=data.I*omega_dot_des+(cross(state(10:12),data.I*state(10:12)));
u=[c_cmd*data.m,tau'];
end