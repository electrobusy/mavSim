%% get trajectory
clear all

plotflag=0;
close all

if not(exist('poly_coeffs','var'))
    minimum_snap_3D_v1;
end

clearvars -except poly_coeffs keyframes t 

 %% Drone Parameters
    data.m = 0.38905;
    data.beta = 0.5;
    data.g = 9.81;    
    data.L=0.2;
    data.I=diag([0.0049,0.0049,0.0069]); % kg m^2
    data.k_F=1.91e-6;
    data.k_M=2.6e-7;
    data.tau_g= [0;0;0]; %propeller gyroscopic torque coefficnet
    data.D=diag([0,0,0]); %mass normalized drag R*D*R'*V
    data.A=zeros(3);% rotational linear drag A*R'*v
    data.B=zeros(3);%   rotational rotation drag B*omega 
    
    data.rotorcontrol=0;    
    data.K_pos=eye(3).*[0.5,0.5,0.5]; %position error to desired acceleration  %set to zero for debugging the controller 
    data.K_vel=eye(3).*[0.3,0.3,0.3];  %velocity error to acceleration
%     data.K_acc=eye(3).*[0;0;0];
    data.K_theta=eye(3).*[10;10;10]; %angle to omega
    data.K_omega=eye(3).*[10;10;10]; %omega to omega_dot
    


%% Get input vaLues from differential flatness
T=t;%[0,t_m]; %Timestamps
data.dt=1e-2;
for i = 1:length(T)
    T(i)=T(i)-mod(T(i),data.dt);
end


[t,u,poly_coeffs,extra,disc]=get_inputs_Diff_Flatness(poly_coeffs,T,data);


%% Simulate 
state_0=[keyframes(1:3,1); 0;0;0; 0;0;0;extra.omega(1,:)']; % [pos, R,V, omega];  %for perfect feedforward control there should be a nonzero initial angular rate
% state x = [x,y,z,phi,theta,psi,xdot,ydot,zdot,omega_x,omega_y,omega_z]^T 
dt=mean(diff(t));
state=state_0;
statel=zeros(size(state,1),length(t)+1);
statel(:,1)=state_0;
dsl=[];
euler_des_l=[];
omega_des_l=[];
omega_dot_des_l=[];
u2_l=[];
Rb=get_rotationmatrix((state_0(4:6))','E2B');
R_list=[];
for i=1:length(t)
    

    data.omega=extra.omega(i,:)';
    data.omega_dot=extra.omega_dot(i,:)';
    uin=u(i,:);
    if(t(i)>10)
        dymm=2;
    end
    [u2,ctrlextra]=controller(state,Rb,disc(i,:,:),data); %disc= discretized polynomial

    euler_des_l=[euler_des_l,ctrlextra.euler_des];
    omega_des_l=[omega_des_l,ctrlextra.omega_des];
    omega_dot_des_l=[omega_dot_des_l,ctrlextra.omega_dot_des];
    
    u2_l=[u2_l;u2];
    uin=u2; %comment to disable controller and fly on feedforward (don't forget initial conditions should be perfect too)
%     ds=quad_dynamics_Mellinger_Kumar(state,u(i,:)',data);  
    [ds, Rbdot]= quad_3D_Dynamics_Kumar_Jelle(state,Rb,uin',data);
    dsl=[dsl,ds];
    state=state+(ds*dt);
    Rb=Rb+Rbdot*dt; %propagate orientation matrix
    Rb=Rb./vecnorm(Rb,2,1);%normalize columns
    if sum(state(4:6)>pi)
        state(find(state(4:6)>pi)+3)=state(find(state(4:6)>pi)+3)-2*pi;
    elseif sum(state(4:6)<-pi)
        state(find(state(4:6)<-pi)+3)=state(find(state(4:6)<-pi)+3)+2*pi;
    end
        
    statel(:,i+1)=state;
    R_list=cat(3,R_list,Rb);
  
end

%% plotting
% close all
r2d=180/pi;
figure()
plot3(statel(1,:),statel(2,:),statel(3,:))
% hold on

% patch(statel(1,2:end),statel(2,2:end),statel(3,2:end),vecnorm(dsl(1:3,:),2,1),'FaceColor','none','EdgeColor','interp')
% c = colorbar;  % Add a colorbar
% view(3)
% c.Label.String = 'V [m/s]';
xlabel('x [m]');
ylabel('y [m]');
grid on;
view([1,1,1]);

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
grid on
xlabel('t');
ylabel('phi')
legend("Simulated","From Flatness","Desired from controller",'Location','northwest');

subplot(432)
plot(t,statel(5,1:end-1)*r2d)
hold on 
plot(t,extra.theta*r2d,'--')
hold on 
plot(t,euler_des_l(2,:)*r2d,'--');
grid on
title("Pitch")
xlabel('t');
ylabel('theta')

subplot(433)
plot(t,statel(6,1:end-1)*r2d)
hold on 
plot(t,extra.psi*r2d,'--');
hold on 
plot(t,euler_des_l(3,:)*r2d,'--');
hold on
plot(t,disc(:,4,1)*r2d);
hold on 
 plot(T,keyframes(4,:)*r2d,'o');
grid on
title("Yaw")
xlabel('t');
ylabel('Yaw')

subplot(434)
plot(t,statel(10,1:end-1)*r2d)
hold on 
plot(t,extra.omega(:,1)*r2d,'--')
hold on 
plot(t,omega_des_l(1,:)*r2d,'--')

title('omega\_x')
grid on
subplot(435)
plot(t,statel(11,1:end-1)*r2d)
hold on
plot(t,extra.omega(:,2)*r2d,'--')
hold on 
plot(t,omega_des_l(2,:)*r2d,'--')
grid on
title('omega\_y ')
subplot(436)
plot(t,statel(12,1:end-1)*r2d)
hold on 
plot(t,extra.omega(:,3)*r2d,'--');
hold on 
plot(t,omega_des_l(3,:)*r2d,'--')
title('omega\_z');
grid on

subplot(437)
plot(t,dsl(10,:)*r2d)
hold on 
plot(t,extra.omega_dot(:,1)*r2d,'--')
hold on 
plot(t,omega_dot_des_l(1,:)*r2d,'--')
title('omega\_x\_dot')
grid on
subplot(438)
plot(t,dsl(11,:)*r2d);
hold on 
plot(t,extra.omega_dot(:,2)*r2d,'--')
hold on 
plot(t,omega_dot_des_l(2,:)*r2d,'--')
title('omega\_y\_dot');
grid on
subplot(439)
plot(t,dsl(12,:)*r2d);
hold on 
plot(t,extra.omega_dot(:,3)*r2d,'--')
hold on 
plot(t,omega_dot_des_l(3,:)*r2d,'--')
title('omega\_z\_dot');
grid on

subplot(4,3,10)
plot(t,u2_l(:,2),'--');
hold on 
plot(t,u(:,2),'--')
title("Roll moment")
xlabel('t');
ylabel('Nm')
grid on
subplot(4,3,11)
plot(t,u2_l(:,3),'--');
hold on 
plot(t,u(:,3),'--')
title("Pitch moment")
grid on
xlabel('t');
ylabel('Nm')
subplot(4,3,12)
plot(t,u2_l(:,4),'--');
hold on 
plot(t,u(:,4),'--')
title("Yaw moment")
xlabel('t');
ylabel('Nm')
grid on



figure()
subplot(2,3,1)
plot(t,statel(1,1:end-1));
hold on 
plot(t,disc(:,1,1),'--')
title('x')
grid on 

subplot(2,3,2)
plot(t,statel(2,1:end-1));
hold on 
plot(t,disc(:,2,1),'--')
title('y')
grid on 

subplot(2,3,3)
plot(t,statel(3,1:end-1));
hold on 
plot(t,disc(:,3,1),'--')
title('z')
grid on 

subplot(2,3,4)
plot(t,statel(7,1:end-1));
hold on 
plot(t,disc(:,1,2),'--')
title('vx')
grid on 

subplot(2,3,5)
plot(t,statel(8,1:end-1));
hold on 
plot(t,disc(:,2,2),'--')
title('vy')
grid on 

subplot(2,3,6)
plot(t,statel(9,1:end-1));
hold on 
plot(t,disc(:,3,2),'--')
title('vz')
grid on 



figure()
subplot(141)
plot(t,u(:,1));
grid on 
title("Thrust [N]") 

subplot(142)
plot(t,u(:,2));
grid on 
title("\tau\_x [Nm]")  

subplot(143)
plot(t,u(:,3));
grid on 
title("\tau\_y [Nm]") 

subplot(144)
plot(t,u(:,4));
grid on
title("\tau\_z [Nm]") 

omega_z2=[extra.omega(1,3)];
for i = 1:size(extra.omega_dot,1)
    az=extra.omega_dot(i,3);
    omega_z2=[omega_z2;omega_z2(end)+az*dt;];
    dummy=2;
end
