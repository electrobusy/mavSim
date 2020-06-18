% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Comparison: Optimal Control vs PID
% Model: quadrotor 3d

clear all; close all; 

% -- add path: 
p = genpath('C:\Users\Rohan\Desktop\thesiscode\direct\quadrotor3d\quad_quat_newtemplate');
addpath(p);

addpath('functions');

% -- controllers:
controller = 'PID';
optimal = false;

% -- initial and final States: 
phi_i = 0*pi/180;
theta_i = 0*pi/180;
psi_i = 0*pi/180;

[q0_i,q1_i,q2_i,q3_i] = Euler2Quat(phi_i,theta_i,psi_i);

phi_f = 0*pi/180;
theta_f = 0*pi/180;
psi_f = 0*pi/180;

[q0_f,q1_f,q2_f,q3_f] = Euler2Quat(phi_f,theta_f,psi_f);

x_i = [0, 0, 0, 0, 0, 0, q0_i, q1_i, q2_i, q3_i]'; % height should be between [-3,3]
x_ref = [0, 0, 0, 0, 0, 0, q0_f, q1_f, q2_f, q3_f]'; 

Q_cmd_mat = @(q0, q1, q2, q3) [q0, -q1, -q2, -q3; q1, q0, -q3, q2; q2, q3, q0, -q1; q3, -q2, q1, q0];

% -- parameters
m = 0.38905;
beta = 0.5;
g = 9.81;

numStates = 10;
numInputs = 4;

% -- simulation time
T = 5; % [sec]
dt = 0.01; % [sec]
numPts = floor(T/dt); % [-]
t = linspace(0,T,numPts); % [sec]

x = zeros(numStates,numPts);
u = zeros(numInputs,numPts);

x(:,1) = x_i;
x(:,2) = x_i;

% -- Gains:
Kp_xy = [0.4 0.4]';
Kd_xy = [0.5 0.5]';
Ki_xy = [0 0]';

Kp_v = [-0.5 -0.5]';
Kd_v = [-0.1 -0.1]';
Ki_v = [0 0]';

Kp_quat = [20 20 5]';
Kd_quat = [0.2 0.2 0.01]';
Ki_quat = [0 0 0]';

Kp_z = 4;
Kd_z = 3;
Ki_z = 0;

% Initialization variables:
u_T = zeros(1,numPts);
u_omega = zeros(3,numPts);

phi_theta_ref = zeros(2,numPts);

e_pos = zeros(3,numPts);
de_pos = zeros(3,numPts);
ie_pos = zeros(3,numPts);

v_e = zeros(2,numPts);
dv_e = zeros(2,numPts);
iv_e = zeros(2,numPts);

q_e = zeros(4,numPts);
dq_e = zeros(3,numPts);
iq_e = zeros(3,numPts);

for i = 2:numPts-1
    
    pos = x(1:3,i);
    vel = x(4:6,i);
    q = x(7:end,i);
    
    [phi, theta, psi] = Quat2Euler(q(1),q(2),q(3),q(4));
    
%     % --> 2) Position x-y - PID loop -- Pass to commanded velocity
%     % -- error computation:
%     e_pos(1:2,i) = x_ref(1:2) - pos(1:2);
%     
%     % -- computation of error integration and derivative:
%     de_pos(1:2,i) = vel(1:2); % from the system of equations of quadrotor: theta_dot = omega
%     ie_pos(1:2,i) = ie_pos(1:2,i-1) + (e_pos(1:2,i) - e_pos(1:2,i-1))*dt;
%     
%     % -- compute controller:
%     v_ref = Kp_xy.*e_pos(1:2,i) + Kd_xy.*de_pos(1:2,i) + Ki_xy.*ie_pos(1:2,i);
%     
%     % --> 1) Velocity x-y - PID loop -- Pass to commanded angles
%     % -- error computation:
%     v_ref = [1 -1]';
%     v_e(:,i) = v_ref - vel(1:2);
%     
%     % -- computation of error integration and derivative:
%     dv_e(:,i) = (x(4:5,i) - x(4:5,i-1))/dt; 
%     iv_e(:,i) = iv_e(:,i-1) + (v_e(:,i) - v_e(:,i-1))*dt;
%     
%     % -- compute controller:
%     a_ref = Kp_v.*v_e(:,i) + Kd_v.*dv_e(:,i) + Ki_v.*iv_e(:,i) + v_e(:,i);
%     
%     phi_theta_ref(:,i) = [cos(psi), -sin(psi); sin(psi), cos(psi)]'*a_ref;
%     
    [q0_ref,q1_ref,q2_ref,q3_ref] = Euler2Quat(phi_f, theta_f, psi_f);
    
    % % --> 0) Angle - PID inner loop
    % -- error computation:
    q_e(:,i) = Q_cmd_mat(q0_ref,q1_ref,q2_ref,q3_ref)*[q(1), -q(2), -q(3), -q(4)]';
    
    % -- computation of error integration and derivative:
    dq_e(:,i) = u_omega(:,i-1); % from the system of equations of quadrotor: theta_dot = omega
    iq_e(:,i) = iq_e(:,i-1) + (q_e(2:end,i) - q_e(2:end,i-1))*dt;
    
    % -- compute controller:
    u_omega(:,i) = Kp_quat.*q_e(2:end,i) + Kd_quat.*dq_e(:,i) + Ki_quat.*iq_e(:,i);
    
    % % --> a) Vertical position - PID loop
    % -- error computation:
    e_pos(3,i) = x_ref(3) - pos(3);
    
    % -- computation of error integration and derivative:
    de_pos(3,i) = -vel(3);
    ie_pos(3,i) = ie_pos(3,i-1) + (e_pos(3,i) - e_pos(3,i-1))*dt;
    
    % -- compute controller:
    u_T(i) = Kp_z*e_pos(3,i) + Kd_z*de_pos(3,i) + Ki_z*ie_pos(3,i);
    
    % --> Append controller:
    u(1,i) = u_T(i);
    u(2:end,i) = u_omega(:,i);
    
    % --> Control Saturation:
    if(u(1,i) > m*g)
        u(1,i) = m*g;
    elseif (u(1,i) < -m*g)
        u(1,i) = -m*g;
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
    
    % 3) Compute results obtained from dynamics:
    f = quad3d_dynamics(x(:,i),u(:,i));
    
    % 4) Propagate state -- Euler integration:
    x(:,i+1) = x(:,i) + dt*f;
    
    % 5) Quaternion normalization:
    q0 = x(7,i+1);
    q1 = x(8,i+1);
    q2 = x(9,i+1);
    q3 = x(10,i+1);
    
    x(7,i+1) = q0/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    x(8,i+1) = q1/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    x(9,i+1) = q2/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    x(10,i+1) = q3/sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    
%     % 7) State Saturation:
%     [phi, theta, psi] = Quat2Euler(x(7,i+1),x(8,i+1),x(9,i+1),x(10,i+1));
%     if(phi > 20*pi/180)
%         phi = 20*pi/180;
%     elseif (phi < -20*pi/180)
%         phi = -20*pi/180;
%     end
%     
%     if(theta > 20*pi/180)
%         theta = 20*pi/180;
%     elseif (theta < -20*pi/180)
%         theta = -20*pi/180;
%     end
%     
%     [q0,q1,q2,q3] = Euler2Quat(phi,theta,psi);
%     
%     x(7,i+1) = q0;
%     x(8,i+1) = q1;
%     x(9,i+1) = q2;
%     x(10,i+1) = q3;
end

%% Run Optimal Controller

if optimal
    % -- Number of Nodes
    numNodes = 100;
    
    % Homotopy parameter
    alpha = 0.8;
    
    % Parameters:
    in.beta = beta;          % [-]
    in.m =  m;          % [kg]
    in.gAcc = g;         % [m/s^2]
    % ----
    
    % -- Final time guess/limit:
    in.tGuess = 20;
    
    % -- State bounds
    in.angleLim = 20*pi/180;
    in.x_I = [0 0 1]';
    in.y_B = [0 0 1]';
    
    % -- Control Bounds:
    % - min
    in.TMin = -in.m*in.gAcc;    % [N]
    in.pMin = -20;              %  [rad/s]
    in.qMin = -20;
    in.rMin = -2;
    % - max
    in.TMax = in.m*in.gAcc;                 % [N]
    in.pMax = 20;                  % [rad/s]
    in.qMax = 20;
    in.rMax = 2;
    
    % Initial attitude
    phiBegin = phi_i;       % [rad]
    thetaBegin = theta_i;     % [rad]
    psiBegin = psi_i;       % [rad]
    
    % Final attitude:
    phiFinal = 0;
    thetaFinal = 0;
    psiFinal = 0;
    
    % -- Initial States
    x1Begin = x_i(1);          % [m]
    x2Begin = x_i(2);            % [m]
    x3Begin = x_i(3);            % [m]
    vx1Begin = x_i(4);           % [m/s]
    vx2Begin = x_i(5);           % [m/s]
    vx3Begin = x_i(6);           % [m/s]
    [q0Begin, q1Begin, q2Begin, q3Begin] = Euler2Quat(phiBegin,thetaBegin,psiBegin);
    
    stateInit = [x1Begin, x2Begin, x3Begin, vx1Begin, vx2Begin, vx3Begin, q0Begin, q1Begin, q2Begin, q3Begin]';
    
    % -- Target States
    in.x1Final = 0;
    in.x2Final = 0;
    in.x3Final = 0;
    in.vx1Final = 0;
    in.vx2Final = 0;
    in.vx3Final = 0;
    [in.q0Final, in.q1Final, in.q2Final, in.q3Final] = Euler2Quat(phiFinal,thetaFinal,psiFinal);
    
    % -- Fetch the problem definition
    [problem,guess] = myProblem(in,stateInit,alpha);
    % -- Get options and solver settings
    options = problem.settings(numNodes);
    [solution,MRHistory] = solveMyProblem(problem, guess, options);
    
    [ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01);
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
if optimal
    hold on;
    plot(solution.T(:,1),solution.X(:,1),'*')
    hold on;
    plot(tv,xv(:,1));
end
xlabel('t [sec]')
ylabel('x [m]')
grid on

subplot(3,3,4);
plot(t,x(2,:))
if optimal
    hold on;
    plot(solution.T(:,1),solution.X(:,2),'*');
    hold on;
    plot(tv,xv(:,2));
end
xlabel('t [sec]')
ylabel('y [m]')
grid on

subplot(3,3,7);
plot(t,x(3,:));
line([t(1), t(end)], [x_ref(3) x_ref(3)], 'LineStyle', '-', 'Color', 'r');
axis([t(1), t(end), -0.5, x_ref(3) + 0.5]);
if optimal
    hold on;
    plot(solution.T(:,1),solution.X(:,3),'*');
    hold on;
    plot(tv,xv(:,3));
end
xlabel('t [sec]')
ylabel('z [m]')
grid on

% -- Velocity
subplot(3,3,2);
plot(t,x(4,:));
% line([t(1), t(end)], [v_ref(1) v_ref(1)], 'LineStyle', '-', 'Color', 'r');
if optimal
    hold on;
    plot(solution.T(:,1),solution.X(:,4),'*');
    hold on;
    plot(tv,xv(:,4));
end
xlabel('t [sec]')
ylabel('vx [m/s]')
grid on

subplot(3,3,5);
plot(t,x(5,:));
% line([t(1), t(end)], [v_ref(2) v_ref(2)], 'LineStyle', '-', 'Color', 'r');
if optimal
    hold on;
    plot(solution.T(:,1),solution.X(:,5),'*');
    hold on;
    plot(tv,xv(:,5));
end
xlabel('t [sec]')
ylabel('vy [m/s]')
grid on

subplot(3,3,8);
plot(t,x(6,:))
if optimal
    hold on;
    plot(solution.T(:,1),solution.X(:,6),'*');
    hold on;
    plot(tv,xv(:,6));
end
xlabel('t [sec]')
ylabel('vz [m/s]')
grid on
% -- Angles -- Euler
% figure(4); % clf;
[phid, thetad, psid] = Quat2Euler(x(7,:),x(8,:),x(9,:),x(10,:));
subplot(3,3,3);
plot(t,phid*180/pi);
hold on;
plot(t,phi_f*180/pi, 'LineStyle', '-', 'Color', 'r');
if optimal
    [phi, theta, psi] = Quat2Euler(solution.X(:,7),solution.X(:,8),solution.X(:,9),solution.X(:,10));
    [phiv, thetav, psiv] = Quat2Euler(xv(:,7),xv(:,8),xv(:,9),xv(:,10));
    hold on;
    plot(solution.T(:,1),phi*180/pi,'*');
    hold on;
    plot(tv,phiv*180/pi);
end
xlabel('t [sec]')
ylabel('\phi [deg]')
grid on

subplot(3,3,6);
plot(t,thetad*180/pi);
hold on;
plot(t,phi_f*180/pi, 'LineStyle', '-', 'Color', 'r');
if optimal
    hold on;
    plot(solution.T(:,1),theta*180/pi,'*');
    hold on;
    plot(tv,thetav*180/pi);
end
xlabel('t [sec]')
ylabel('\theta [deg]')
grid on

subplot(3,3,9);
plot(t,psid*180/pi);
hold on; 
line([t(1), t(end)], [psi_f psi_f], 'LineStyle', '-', 'Color', 'r');
if optimal
    hold on;
    plot(solution.T(:,1),psi*180/pi,'*');
    hold on;
    plot(tv,psiv*180/pi);
end
xlabel('t [sec]')
ylabel('\psi [deg]')
grid on

% figure(3);
% plot(t,sqrt(x(4,:).^2 + x(5,:).^2 + x(6,:).^2))
if optimal
    % hold on;
    % plot(solution.T(:,1),sqrt(solution.X(:,4).^2 + solution.X(:,5).^2 + solution.X(:,6).^2),'*')
    % hold on;
    % plot(tv,sqrt(xv(:,4).^2 + xv(:,5).^2 + xv(:,6).^2));
end
% xlabel('t [sec]')
% ylabel('V_{tot} [m/s]')
% grid on

% % -- Angles:
% figure(3); % clf;
% subplot(3,2,1);
% plot(t,x(7,:));
% if optimal
%     hold on;
%     plot(solution.T(:,1),solution.X(:,7),'*');
%     hold on;
%     plot(tv,xv(:,7));
% end
% xlabel('t [sec]')
% ylabel('q0 [-]')
% grid on
% 
% subplot(3,2,2);
% plot(t,x(8,:));
% if optimal
%     hold on;
%     plot(solution.T(:,1),solution.X(:,8),'*');
%     hold on;
%     plot(tv,xv(:,8));
% end
% xlabel('t [sec]')
% ylabel('q1 [-]')
% grid on
% 
% subplot(3,2,3);
% plot(t,x(9,:));
% if optimal
%     hold on;
%     plot(solution.T(:,1),solution.X(:,9),'*');
%     hold on;
%     plot(tv,xv(:,9));
% end
% xlabel('t [sec]')
% ylabel('q2 [-]')
% grid on
% 
% subplot(3,2,4);
% plot(t,x(10,:));
% if optimal
%     hold on;
%     plot(solution.T(:,1),solution.X(:,10),'*');
%     hold on;
%     plot(tv,xv(:,10));
% end
% xlabel('t [sec]')
% ylabel('q3 [-]')
% grid on
% 
% subplot(3,2,[5 6]);
% plot(t,sqrt(x(7,:).^2 + x(8,:).^2 + x(9,:).^2 + x(10,:).^2));
% if optimal
%     hold on;
%     plot(solution.T(:,1),sqrt(solution.X(:,7).^2 + solution.X(:,8).^2 + solution.X(:,9).^2 + solution.X(:,10).^2),'*')
%     hold on;
%     plot(tv,sqrt(xv(:,7).^2 + xv(:,8).^2 + xv(:,9).^2 + xv(:,10).^2));
% end
% xlabel('t [sec]');
% ylabel('|q|');
% grid on;

% -- Control inputs:
figure(4); % clf;
subplot(2,2,1);
plot(t,u(1,:));
if optimal
    hold on;
    plot(solution.T(:,1),solution.U(:,1),'*');
    hold on;
    plot(tv,uv(:,1));
end
xlabel('t [sec]')
ylabel('T [N]')
grid on

subplot(2,2,2);
plot(t,u(2,:));
if optimal
    hold on;
    plot(solution.T(:,1),solution.U(:,2),'*');
    hold on;
    plot(tv,uv(:,2));
end
xlabel('t [sec]')
ylabel('p [rad/s]')
grid on

subplot(2,2,3);
plot(t,u(3,:));
if optimal
    hold on;
    plot(solution.T(:,1),solution.U(:,3),'*');
    hold on;
    plot(tv,uv(:,3));
end
xlabel('t [sec]')
ylabel('q [rad/s]')
grid on

subplot(2,2,4);
plot(t,u(4,:));
if optimal
    hold on;
    plot(solution.T(:,1),solution.U(:,4),'*');
    hold on;
    plot(tv,uv(:,4));
end
xlabel('t [sec]')
ylabel('r [rad/s]')
grid on
