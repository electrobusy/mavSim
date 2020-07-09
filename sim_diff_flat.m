%% get trajectory
clear all

plotflag=0;
close all

if not(exist('poly_coeffs','var'))
    minimum_snap_3D_v0;
end

clearvars -except poly_coeffs t_m keyframes

 %% Drone Parameters
    data.m = 0.38905;
    data.beta = 0.5;
    data.g = 9.81;    
    data.L=0.2;
    data.I=diag([0.0049,0.0049,0.0069]); % kg m^2
    data.k_F=1.91e-6;
    data.k_M=2.6e-7;
    data.rotorcontrol=0;
    
    data.K_pos=eye(3).*[1;0;0]; %position to angle
    data.K_vel=eye(3).*[1;0;0]; %velocity to angle
%     data.K_acc=eye(3).*[0;0;0];
    data.K_theta=eye(3).*[10;10;10]; %angle to omega
    data.K_omega=eye(3).*[20;20;20]; %omega to omega_dot
    


%% Get input vaLues from differential flatness
T=[0,t_m]; %Timestamps
data.dt=1e-3;
[t,u,poly_coeffs,extra,disc]=get_inputs_Diff_Flatness(poly_coeffs,T,data);


%% Simulate 
state_0=[keyframes(1:3,1); 0;0;0; 0;0;0; extra.omega(1,:)']; % [pos, R,V, omega];  
% state x = [x,y,z,phi,theta,psi,xdot,ydot,zdot,omega_x,omega_y,omega_z]^T 
dt=mean(diff(t));
state=state_0;
statel=zeros(size(state,1),length(t)+1);
statel(:,1)=state_0;
dsl=[];
euler_des_l=[];
omega_des_l=[];
omega_dot_des_l=[];
for i=1:length(t)
    
%     ds=quad_dynamics_Mellinger_Kumar(state,u(i,:)',data);  
    data.omega=extra.omega(i,:)';
    data.omega_dot=extra.omega_dot(i,:)';
    uin=u(i,:);
    [u2,ctrlextra]=controller(state,disc(i,:,:,:),data); %disc= discretized polynomial
    euler_des_l=[euler_des_l,ctrlextra.euler_des];
    omega_des_l=[omega_des_l,ctrlextra.omega_des];
    omega_dot_des_l=[omega_dot_des_l,ctrlextra.omega_dot_des];
    
    
    uin=u2;
    ds= quad_3D_Dynamics_Kumar_Jelle(state,uin',data);
    dsl=[dsl,ds];
    state=state+(ds*dt);
    statel(:,i+1)=state;
    if i==400
        dummy=2
    end
    dymm=2;
end

r2d=180/pi;
figure()
plot3(statel(1,:),statel(2,:),statel(3,:))
xlabel('x')
ylabel('y')
zlabel('z')
title("NWU frame")
grid on  
axis equal

figure()
subplot(431)
plot(t,statel(4,1:end-1)*r2d)
hold on 
plot(t,extra.phi*r2d,'--');
hold on 
plot(t,euler_des_l(1,:)*r2d,'--');
title("Roll")
xlabel('t');
ylabel('phi')
legend("Simulated","From Flatness","Desired from controller",'Location','northwest');

subplot(432)
plot(t,statel(5,1:end-1)*r2d)
hold on 
plot(t,extra.theta*r2d,'--')
hold on 
plot(t,euler_des_l(2,:)*r2d,'--');
title("Pitch")
xlabel('t');
ylabel('theta')

subplot(433)
plot(t,statel(6,1:end-1)*r2d)
hold on 
plot(t,extra.psi*r2d,'--');
hold on 
plot(t,euler_des_l(3,:)*r2d,'--');
title("Yaw")
xlabel('t');
ylabel('Yaw')

subplot(434)
plot(t,statel(10,1:end-1))
hold on 
plot(t,extra.omega(:,1),'--')
hold on 
plot(t,omega_des_l(1,:),'--')
title('omega\_x')
subplot(435)
plot(t,statel(11,1:end-1))
hold on
plot(t,extra.omega(:,2),'--')
hold on 
plot(t,omega_des_l(2,:),'--')
title('omega\_y ')
subplot(436)
plot(t,statel(12,1:end-1))
hold on 
plot(t,extra.omega(:,3),'--');
hold on 
plot(t,omega_des_l(3,:),'--')
title('omega\_z');

subplot(437)
plot(t,dsl(10,:))
hold on 
plot(t,extra.omega_dot(:,1),'--')
hold on 
plot(t,omega_dot_des_l(1,:),'--')
title('omega\_x\_dot')
subplot(438)
plot(t,dsl(11,:));
hold on 
plot(t,extra.omega_dot(:,2),'--')
hold on 
plot(t,omega_dot_des_l(2,:),'--')
title('omega\_y\_dot');
subplot(439)
plot(t,dsl(12,:));
hold on 
plot(t,extra.omega_dot(:,3),'--')
hold on 
plot(t,omega_dot_des_l(3,:),'--')
title('omega\_z\_dot');

subplot(4,3,10)
plot(t,u(:,2))
title("Roll moment")
xlabel('t');
ylabel('Nm')
subplot(4,3,11)
plot(t,u(:,3))
title("Pitch moment")
xlabel('t');
ylabel('Nm')
subplot(4,3,12)
plot(t,u(:,4))
title("Yaw moment")
xlabel('t');
ylabel('Nm')

[tt,traj]=discretize_poly(poly_coeffs,T,mean(diff(t)));

figure()
plot(t,u(:,1))
title("Thrust");


omega_z2=[extra.omega(1,3)];
for i = 1:size(extra.omega_dot,1)
    az=extra.omega_dot(i,3);
    omega_z2=[omega_z2;omega_z2(end)+az*dt;];
    dummy=2;
end

