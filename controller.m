function [u, extra] = controller(state,R,reference,data)
phi=state(4);
theta=state(5);
psi=state(6);

%% get orientation
Rb=R;%[x_b,y_b,z_b];

% theta=-asin(Rb(3,1));
% psi=reference(1,4,1);%atan2(Rb(2,1)/cos(theta),Rb(1,1)/cos(theta));
z_b=Rb(:,3);
%% get desired acceleration
a_ref=reference(1,1:3,3);
a_fb=data.K_pos*(reference(1,1:3,1)'-state(1:3))+data.K_vel*(reference(1,1:3,2)'-state(7:9));
% a_fb=Rb'*a_fb;
a_des=(a_fb'+a_ref)'+[0;0;data.g];


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
    
    theta_des=-asin(R_des(3,1));
    phi_des=atan2(R_des(3,2)/cos(theta_des),R_des(3,3)/cos(theta_des));
    psi_des=atan2(R_des(2,1)/cos(theta_des),R_des(1,1)/cos(theta_des));
   
    euler_des=[phi_des;theta_des;psi_des];
%     euler_des= [atan2(R_des(3,2),R_des(3,3)); -asin(R_des(3,1)); atan2(R_des(2,1),R_des(1,1))];
%     euler_des2= [atan2(R_des(3,2),R_des(3,3)); atan2(-R_des(3,1),norm([R_des(3,2),R_des(3,3)],2)); atan2(R_des(2,1),R_des(1,1))];
    euler_error=(euler_des-state(4:6));
    
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
    extra.R_des=R_des;
%$ get desired torque:
tau=data.I*omega_dot_des+(cross(state(10:12),data.I*state(10:12)));
u=[c_cmd*data.m,tau'];
end