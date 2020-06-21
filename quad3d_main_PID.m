% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Comparison: Optimal Control vs PID
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
x_ref = [-8, 0, 0, 0, 0, 0, q0_f, q1_f, q2_f, q3_f]'; % r_xy = [-4,4] m

% -- Operations: 
Q_cmd_mat = @(q0, q1, q2, q3) [q0, -q1, -q2, -q3; q1, q0, -q3, q2; q2, q3, q0, -q1; q3, -q2, q1, q0];

% -- parameters:
data.m = 0.38905;
data.beta = 0.5;
data.g = 9.81;

numStates = 10;
numInputs = 4;
d2r = pi/180;

% -- simulation time
T = 5; % [sec]
dt = 0.01; % [sec]
numPts = floor(T/dt); % [-]
t = linspace(0,T,numPts); % [sec]

% -- initialization:
x = zeros(numStates,numPts);
u = zeros(numInputs,numPts);

q_ref = zeros(4,numPts);
q_ref(:,1) = [q0_i, q1_i, q2_i, q3_i]';

u_T = zeros(1,numPts);
u_omega = zeros(3,numPts);

phi_theta_ref = zeros(2,numPts);
a_ref = zeros(2,numPts);

e_pos = zeros(3,numPts);
de_pos = zeros(3,numPts);
ie_pos = zeros(3,numPts);

v_e = zeros(2,numPts);
dv_e = zeros(2,numPts);
iv_e = zeros(2,numPts);

q_e = zeros(4,numPts);
dq_e = zeros(3,numPts);
iq_e = zeros(3,numPts);

% -- Gains:

% - altitude loop:
Kp_z = 4;
Kd_z = 3;
Ki_z = 0;

% - attitude loop:
Kp_quat = [20 20 5]';
Kd_quat = [0.2 0.2 0.01]';
Ki_quat = [0 0 0]';

% - velocity loop:
vel_loop = true;  
Kp_v = [50 20]'; % [3 3]'; % comment is without feedforward term
Kd_v = [-0.1 -0.1]'; % [-0.4 -0.4]';
Ki_v = [0 0]';

% position loop:
pos_loop = true; % NOTE: Set this to true after "vel_loop" is true. Otherwise, it will prompt an error!
Kp_xy = [0.5 0.5]';
Kd_xy = [0 0]';
Ki_xy = [0 0]';

% Initial states: 
x(:,1) = x_i;
x(:,2) = x_i;

if vel_loop
    v_ref = [5*ones(1,numPts); 0*ones(1,numPts)];
else
    v_ref = zeros(2,numPts); 
end

for i = 2:numPts-1
    
    pos = x(1:3,i);
    vel = x(4:6,i);
    q = x(7:end,i);
    
    [phi, theta, psi] = Quat2Euler(q(1),q(2),q(3),q(4));
    
    % --> Position x-y - PID loop -- Get commanded velocity
    if pos_loop
        % -- error computation:
        e_pos(1:2,i) = x_ref(1:2) - pos(1:2);
    
        % -- computation of error integration and derivative:
        de_pos(1:2,i) = vel(1:2); % from the system of equations of quadrotor: theta_dot = omega
        ie_pos(1:2,i) = ie_pos(1:2,i-1) + (e_pos(1:2,i))*dt;
    
        % -- compute controller:
        v_ref(:,i) = Kp_xy.*e_pos(1:2,i) + Kd_xy.*de_pos(1:2,i) + Ki_xy.*ie_pos(1:2,i);
    end
    
    % --> Velocity x-y - PID loop -- Pass to commanded angles
    if vel_loop
        % -- error computation:
        v_e(:,i) = v_ref(:,i) - vel(1:2);
        % -- computation of error integration and derivative:
        dv_e(:,i) = (x(4:5,i) - x(4:5,i-1))/dt; 
        iv_e(:,i) = iv_e(:,i-1) + v_e(:,i-1)*dt;
      
        % -- compute controller:
        a_ref(:,i) = Kp_v.*v_e(:,i) + Kd_v.*dv_e(:,i) + Ki_v.*iv_e(:,i) + v_e(:,i);
    
        % Now, linearization around hover yields
        % phi_theta_ref(:,i) = [cos(psi) -sin(psi); sin(psi) cos(psi)]'*a_ref(:,i);
        phi_theta_ref(:,i) = [sin(psi_f), -cos(psi_f); cos(psi_f), sin(psi_f)]*a_ref(:,i);
        phi_theta_ref(phi_theta_ref > 90) = 90;
        phi_theta_ref(phi_theta_ref < -90) = -90;
        phi_theta_ref = phi_theta_ref*d2r;
        % --- References: 
        % [paper - eq. 10] https://www.researchgate.net/publication/321448210_Quadrotor_trajectory_tracking_using_PID_cascade_control
        % [document - Page 10] https://repository.upenn.edu/cgi/viewcontent.cgi?article=1705&context=edissertations
        % ---------------
        [q_ref(1,i),q_ref(2,i),q_ref(3,i),q_ref(4,i)] = Euler2Quat(phi_theta_ref(1,i), phi_theta_ref(2,i), psi_f);
    else
        [q_ref(1,i),q_ref(2,i),q_ref(3,i),q_ref(4,i)] = Euler2Quat(phi_f, theta_f, psi_f);
    end
    
    [phi_ref_i, theta_ref_i, psi_ref_i] = Quat2Euler(q_ref(1,i),q_ref(2,i),q_ref(3,i),q_ref(4,i));
    
    plot(t(i), phi_theta_ref(2,i), 'gx');
    hold on;
    plot(t(i), theta_ref_i, 'r+')
    hold on;
    plot(t(i), theta, 'bo')
    hold on;
    
    % % --> 0) Angle - PID inner loop
    % -- error computation:
    q_e(:,i) = Q_cmd_mat(q_ref(1,i),q_ref(2,i),q_ref(3,i),q_ref(4,i))*[q(1), -q(2), -q(3), -q(4)]';
    
    % -- computation of error integration and derivative:
    dq_e(:,i) = u_omega(:,i-1); 
    iq_e(:,i) = iq_e(:,i-1) + q_e(2:end,i-1)*dt;
    
    % -- compute controller:
    u_omega(:,i) = Kp_quat.*q_e(2:end,i) + Kd_quat.*dq_e(:,i) + Ki_quat.*iq_e(:,i);
    
    % % --> a) Vertical position - PID loop
    % -- error computation:
    e_pos(3,i) = x_ref(3) - pos(3);
    
    % -- computation of error integration and derivative:
    de_pos(3,i) = -vel(3);
    ie_pos(3,i) = ie_pos(3,i-1) + (e_pos(3,i))*dt;
    
    % -- compute controller:
    u_T(i) = Kp_z*e_pos(3,i) + Kd_z*de_pos(3,i) + Ki_z*ie_pos(3,i);
    
    % --> Append controller:
    u(1,i) = u_T(i);
    u(2:end,i) = u_omega(:,i);
    
    % --> Control Saturation:
    if(u(1,i) > data.m*data.g)
        u(1,i) = data.m*data.g;
    elseif (u(1,i) < -data.m*data.g)
        u(1,i) = -data.m*data.g;
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
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% grid on

figure(2);
subplot(3,3,1);
plot(t,x(1,:));
xlabel('t [sec]')
ylabel('x [m]')
grid on

subplot(3,3,4);
plot(t,x(2,:))
xlabel('t [sec]')
ylabel('y [m]')
grid on

subplot(3,3,7);
plot(t,x(3,:));
line([t(1), t(end)], [x_ref(3) x_ref(3)], 'LineStyle', '-', 'Color', 'r');
axis([t(1), t(end), -0.5, x_ref(3) + 0.5]);
xlabel('t [sec]')
ylabel('z [m]')
grid on

% -- Velocity
subplot(3,3,2);
plot(t,x(4,:));
hold on;
if vel_loop || pos_loop
    plot(t,v_ref(1,:), 'LineStyle', '-', 'Color', 'r');
end
xlabel('t [sec]')
ylabel('vx [m/s]')
grid on

subplot(3,3,5);
plot(t,x(5,:));
hold on;
if vel_loop || pos_loop
    plot(t,v_ref(2,:), 'LineStyle', '-', 'Color', 'r');
end
xlabel('t [sec]')
ylabel('vy [m/s]')
grid on

subplot(3,3,8);
plot(t,x(6,:))
xlabel('t [sec]')
ylabel('vz [m/s]')
grid on
% -- Angles -- Euler
% figure(4); % clf;
[phid, thetad, psid] = Quat2Euler(x(7,:),x(8,:),x(9,:),x(10,:));
[phiref,thetaref,psiref] = Quat2Euler(q_ref(1,:),q_ref(2,:),q_ref(3,:),q_ref(4,:));
subplot(3,3,3);
plot(t,phid*180/pi);
hold on;
if vel_loop
    plot(t,phiref*180/pi, 'LineStyle', '-', 'Color', 'r');
else
    plot(t,phi_f*ones(numPts)*180/pi, 'LineStyle', '-', 'Color', 'r');
end
xlabel('t [sec]')
ylabel('\phi [deg]')
grid on

subplot(3,3,6);
plot(t,thetad*180/pi);
hold on;
if vel_loop
    plot(t,thetaref*180/pi, 'LineStyle', '-', 'Color', 'r');
else
    plot(t,theta_f*ones(numPts)*180/pi, 'LineStyle', '-', 'Color', 'r');
end
xlabel('t [sec]')
ylabel('\theta [deg]')
grid on

subplot(3,3,9);
plot(t,psid*180/pi);
hold on;
if vel_loop
    plot(t,psiref*180/pi, 'LineStyle', '-', 'Color', 'r');
else
    plot(t,psi_f*ones(numPts)*180/pi, 'LineStyle', '-', 'Color', 'r');
end
xlabel('t [sec]')
ylabel('\psi [deg]')
grid on
% 
% % -- Angles:
% figure(3); % clf;
% subplot(3,2,1);
% plot(t,x(7,:));
% hold on;
% plot(t,q_ref(1,:));
% xlabel('t [sec]')
% ylabel('q0 [-]')
% grid on
% 
% subplot(3,2,2);
% plot(t,x(8,:));
% hold on;
% plot(t,q_ref(2,:));
% xlabel('t [sec]')
% ylabel('q1 [-]')
% grid on
% 
% subplot(3,2,3);
% plot(t,x(9,:));
% hold on;
% plot(t,q_ref(3,:));
% xlabel('t [sec]')
% ylabel('q2 [-]')
% grid on
% 
% subplot(3,2,4);
% plot(t,x(10,:));
% hold on;
% plot(t,q_ref(4,:));
% xlabel('t [sec]')
% ylabel('q3 [-]')
% grid on
% 
% subplot(3,2,[5 6]);
% plot(t,sqrt(x(7,:).^2 + x(8,:).^2 + x(9,:).^2 + x(10,:).^2));
% hold on;
% plot(t,sqrt(q_ref(1,:).^2 + q_ref(2,:).^2 + q_ref(3,:).^2 + q_ref(4,:).^2));
% xlabel('t [sec]');
% ylabel('|q|');
% grid on;
% 
% % -- Control inputs:
% figure(4); % clf;
% subplot(2,2,1);
% plot(t,u(1,:));
% xlabel('t [sec]')
% ylabel('T [N]')
% grid on
% 
% subplot(2,2,2);
% plot(t,u(2,:));
% xlabel('t [sec]')
% ylabel('p [rad/s]')
% grid on
% 
% subplot(2,2,3);
% plot(t,u(3,:));
% xlabel('t [sec]')
% ylabel('q [rad/s]')
% grid on
% 
% subplot(2,2,4);
% plot(t,u(4,:));
% xlabel('t [sec]')
% ylabel('r [rad/s]')
% grid on
% 
% figure(5);
% subplot(2,1,1);
% plot(t,a_ref(1,:));
% xlabel('t [sec]')
% ylabel('a_x [m/s^2]')
% grid on
% subplot(2,1,2);
% plot(t,a_ref(2,:));
% xlabel('t [sec]')
% ylabel('a_y [m/s^2]')
% grid on