% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Comparison: Optimal Control vs LQR
% Model: quadrotor 3d

%% 
clear all; close all; clc;

% -- initial and final attitude: 
phi_i = 0*pi/180;
theta_i = 0*pi/180;
psi_i = 0*pi/180; 

[q0_i,q1_i,q2_i,q3_i] = Euler2Quat(phi_i,theta_i,psi_i);

phi_f = 0*pi/180;
theta_f = 0*pi/180;
psi_f = 0*pi/180;

[q0_f,q1_f,q2_f,q3_f] = Euler2Quat(phi_f,theta_f,psi_f);

% -- initial and final state: 
x_i = [0, 0, 0, 0, 0, 0, q0_i, q1_i, q2_i, q3_i]';
x_ref = [5, 0, 0, 0, 0, 0, q0_f, q1_f, q2_f, q3_f]';

% -- Operations: 
Q_cmd_mat = @(q0, q1, q2, q3) [q0, -q1, -q2, -q3; q1, q0, -q3, q2; q2, q3, q0, -q1; q3, -q2, q1, q0];

% -- parameters:
data.m = 0.38905; 
data.beta = 0.5; 
data.g = 9.81;

numStates = 10;
numInputs = 4;

% -- simulation time
T = 10; % [sec]
dt = 0.01; % [sec]
numPts = floor(T/dt); % [-]
t = linspace(0,T,numPts); % [sec]

% -- initialization:
x = zeros(numStates,numPts);
u = zeros(numInputs,numPts);

x(:,1) = x_i;
 
% LQR parameters: 
Q = diag([1/10,1/10,1/10,1/5^2,1/5^2,1/5^2,1/10^2,1/10^2,1/10^2,1/10^2]);
R = diag([1/(data.m*data.g)^2,1/10^2,1/10^2,1/10^2]);

% Reference trajectory (hovering case)
u_ref = [data.m*data.g, 0.01, 0.01, 0.01]';

%%
for i = 1:numPts-1
    
%     % Simulate instantaneous disturbances: 
%     if t(i) > 2 && t(i) < 2.2
%         x(1,i) = 2;
%     end
    
    % --> Hovering case (comes from the spatial tracking method -- which yields the desired values of state and input)
    u_0 = u_ref;
    x_0 = x_ref;
    
    % --> LQR command - to get the LQR gain
    [A,B] = quad_linearization(x_0,u_0,data);
    
    % print(rank(ctrb(A,B)));
    
    [K,~,~] = lqr(A,B,Q,R);
    
    % --> Control update
    u(:,i) = - K*(x(:,i) - x_0);
    
    % --> Control Saturation:
    if(u(1,i) > data.g*data.m)
        u(1,i) = data.g*data.m;
    elseif (u(1,i) < -data.g*data.m)
        u(1,i) = -data.g*data.m;
    end
    
    if(u(2,i) > 20)
        u(2,i) = 20;
    elseif (u(2,i) < -20)
        u(2,i) = -20;
    end   
    
    if(u(3,i) > 20)
        u(3,i) = 20;
    elseif (u(3,i) < -20)
        u(3,i) = -20;
    end 
    
    if(u(4,i) > 2)
        u(4,i) = 2;
    elseif (u(4,i) < -2)
        u(4,i) = -2;
    end 
    
    % --> Get state derivative (instantaneous variation in state):
    f = quad3d_dynamics(x(:,i),u(:,i),data);
    
    % --> Propagate state -- Euler integration:
    x(:,i+1) = x(:,i) + dt*f;
    
    % --> Quaternion normalization:
    q0 = x(7,i+1);
    q1 = x(8,i+1);
    q2 = x(9,i+1);
    q3 = x(10,i+1);
    
    x(7,i+1) = q0/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    x(8,i+1) = q1/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    x(9,i+1) = q2/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    x(10,i+1) = q3/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
     
    % 7) State Saturation:
    [phi, theta, psi] = Quat2Euler(x(7,i+1),x(8,i+1),x(9,i+1),x(10,i+1));
    if(phi > 20*pi/180)
        phi = 20*pi/180;
    elseif (phi < -20*pi/180)
        phi = -20*pi/180;
    end
    
    if(theta > 20*pi/180)
        theta = 20*pi/180;
    elseif (theta < -20*pi/180)
        theta = -20*pi/180;
    end
    
    [q0,q1,q2,q3] = Euler2Quat(phi,theta,psi);
    
    x(7,i+1) = q0;
    x(8,i+1) = q1;
    x(9,i+1) = q2;
    x(10,i+1) = q3;
end

%% Plots

% -- States
% figure(1);
% % Position
% plot3(x(1,:),x(2,:),x(3,:))
% if optimal
%     hold on;
%     plot3(solution.X(:,1),solution.X(:,2),solution.X(:,3),'*')
%     hold on;
%     plot3(xv(:,1),xv(:,2),xv(:,3))
%     title(['Final time = ' num2str(solution.tf)]);
% end
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% grid on

figure(2);
subplot(3,3,1);
plot(t,x(1,:));
line([t(1), t(end)], [x_ref(1) x_ref(1)], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('x [m]')
grid on

subplot(3,3,4);
plot(t,x(2,:));
line([t(1), t(end)], [x_ref(2) x_ref(2)], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('y [m]')
grid on

subplot(3,3,7);
plot(t,x(3,:))
line([t(1), t(end)], [x_ref(3) x_ref(3)], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('z [m]')
grid on

% -- Velocity
subplot(3,3,2);
plot(t,x(4,:))
line([t(1), t(end)], [x_ref(4) x_ref(4)], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('vx [m/s]')
grid on

subplot(3,3,5);
plot(t,x(5,:))
line([t(1), t(end)], [x_ref(5) x_ref(5)], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('vy [m/s]')
grid on

subplot(3,3,8);
plot(t,x(6,:))
line([t(1), t(end)], [x_ref(6) x_ref(6)], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('vz [m/s]')
grid on
% -- Angles -- Euler
% figure(4); % clf;
[phid, thetad, psid] = Quat2Euler(x(7,:),x(8,:),x(9,:),x(10,:));
[phiref, thetaref, psiref] = Quat2Euler(x_ref(7),x_ref(8),x_ref(9),x_ref(10));
subplot(3,3,3);
plot(t,phid*180/pi);
line([t(1), t(end)], [phiref phiref], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('\phi [deg]')
grid on

subplot(3,3,6);
plot(t,thetad*180/pi);
line([t(1), t(end)], [thetaref thetaref], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('\theta [deg]')
grid on

subplot(3,3,9);
plot(t,psid*180/pi);
line([t(1), t(end)], [psiref psiref], 'LineStyle', '-', 'Color', 'r');
xlabel('t [sec]')
ylabel('\psi [deg]')
grid on

% % -- Angles:
% figure(3); % clf;
% subplot(3,2,1);
% plot(t,x(7,:))
% xlabel('t [sec]')
% ylabel('q0 [-]')
% grid on
% 
% subplot(3,2,2);
% plot(t,x(8,:))
% xlabel('t [sec]')
% ylabel('q1 [-]')
% grid on
% 
% subplot(3,2,3);
% plot(t,x(9,:))
% xlabel('t [sec]')
% ylabel('q2 [-]')
% grid on
% 
% subplot(3,2,4);
% plot(t,x(10,:))
% xlabel('t [sec]')
% ylabel('q3 [-]')
% grid on
% 
% subplot(3,2,[5 6]);
% plot(t,sqrt(x(7,:).^2 + x(8,:).^2 + x(9,:).^2 + x(10,:).^2));
% xlabel('t [sec]');
% ylabel('|q|');
% grid on;

% -- Control inputs:
figure(4); % clf;
subplot(2,2,1);
plot(t,u(1,:));
xlabel('t [sec]')
ylabel('T [N]')
grid on

subplot(2,2,2);
plot(t,u(2,:));
xlabel('t [sec]')
ylabel('p [rad/s]')
grid on

subplot(2,2,3);
plot(t,u(3,:));
xlabel('t [sec]')
ylabel('q [rad/s]')
grid on

subplot(2,2,4);
plot(t,u(4,:));
xlabel('t [sec]')
ylabel('r [rad/s]')
grid on