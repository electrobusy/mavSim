%% get trajectory
% clear all

plotflag=1;
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



%% Get input vaLues from differential flatness
T=[0,t_m]; %Timestamps
[t,u,poly_coeffs,extra]=get_inputs_Diff_Flatness(poly_coeffs,T,data);
omega=[extra.p,extra.q,extra.r];
u(:,2:4)=[-1,-1,1].*u(:,2:4);
%% Simulate 
state_0=[keyframes(1:3,1); 0;0;0; 0;0;0; -omega(1,:)']; % [pos, R,V, Omega];  
% state x = [x,y,z,phi,theta,psi,xdot,ydot,zdot,p,q,r]^T 
dt=mean(diff(t));
state=state_0;
statel=zeros(size(state,1),length(t)+1);
statel(:,1)=state_0;
dsl=[];
for i=1:length(t)
    
%     ds=quad_dynamics_Mellinger_Kumar(state,u(i,:)',data);  
    data.omega=omega(i,:)';
    data.alpha=extra.alpha(i,:)';
    ds= quad_3D_Dynamics_Kumar_Jelle(state,u(i,:)',data);
    dsl=[dsl,ds];
    state=state+(ds*dt);
    statel(:,i+1)=state;
end

r2d=180/pi;
figure()
plot3(statel(1,:),statel(2,:),-statel(3,:))
xlabel('x')
ylabel('y')
zlabel('-z')
grid on  
axis equal

figure()
subplot(431)
plot(t,statel(4,1:end-1)*r2d)
hold on 
plot(t,-extra.phi*r2d,'--');
title("Roll")
xlabel('t');
ylabel('phi')

subplot(432)
plot(t,statel(5,1:end-1)*r2d)
hold on 
plot(t,extra.theta*r2d,'--')
title("Pitch")
xlabel('t');
ylabel('theta')

subplot(433)
plot(t,statel(6,1:end-1)*r2d)
hold on 
plot(t,extra.psi*r2d,'--');

title("Yaw")
xlabel('t');
ylabel('Yaw')

subplot(434)
plot(t,statel(10,1:end-1))
hold on 
plot(t,-extra.p,'--')
title('Roll rate')
subplot(435)
plot(t,statel(11,1:end-1))
hold on
plot(t,-extra.q,'--')
title('pitch rate')
subplot(436)
plot(t,statel(12,1:end-1))
hold on 
plot(t,extra.r,'--');
title('yaw rate');

subplot(437)
plot(t,dsl(10,:))
hold on 
plot(t,-extra.alpha(:,1),'--')
title('roll acceleration')
subplot(438)
plot(t,dsl(11,:));
hold on 
plot(t,-extra.alpha(:,2),'--')
title('pitch acceleration');
subplot(439)
plot(t,dsl(12,:));
hold on 
plot(t,extra.alpha(:,3),'--')
title('yaw acceleration');

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


