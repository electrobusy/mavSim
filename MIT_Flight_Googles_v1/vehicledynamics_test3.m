%% Lockheed Martin Challenge: Alpha Pilot - TEST 3
% % Author: Rohan Chotalal
% % Supervisors: Guido de Croon & Christophe De Wagter 

% % Goal: "Reverse-engineer" the code used in the MIT - FlightGoggles 
% Simulator for the dynamics.

clc, clear all
close all

%% Parameters: 

% -- Conversion from degrees to radians:
deg2rad = pi/180;

% % Time step: 
dt = 1.0/960; % [sec] -- dt_secs

% % Vehicle properties:
vehicleMass = 1.0; % [kg]
vehicleInertia = [
                 0.0049 0 0;
                 0 0.0049 0;
                 0 0 0.0069;
                 ]; % [kg m^2]
motorTimeconstant = 0.02; % [s]
maxPropSpeed = 2200; % [rad/s] -- this information is per motor
momentArm = 0.08; % [m]

grav = 9.81; % [m/s^2]
thrustCoeff = 1.91e-6; % [N/(rad/s)^2] -- C_Thrust
torqueCoeff = 2.6e-7; % [N/(rad/s)^2]  -- C_torque
dragCoeff = 0.1; % [N/(m/s)^2]

% % UAV Dynamics flags: 
useAutoThrottle = false;
includeDrag = true;

% % Process Noise - Autocorrelation value:
angAccelProcessNoiseAutoCorrelation = 0.00025; % [rad^2/s^2]
linAccelProcessNoiseAutoCorrelation = 0.0005; % [m^2/s^3]   

% % Sensor Noise - Autocorrelation value:
accMeasNoiseVariance = 0.005; % [m^2/s^4]
gyroMeasNoiseVariance = 0.003; % [rad^2/s^2]

% % Filter the states - Low pass filter:
gainP = 35530.5758439217;
gainQ = 266.572976289502;

% % Inner loop - PID parameters: 
propGain = [9, 9, 9]';
intGain = [3, 3, 3]';
derGain = [0.3, 0.3, 0.3]'; 
intState = [0, 0, 0]';
intBound = [1000, 1000, 1000]';

thrust_z = vehicleMass*grav;

% Bounding parameters:
CTRL_MAX_SPEED = 2.5;
CTRL_MAX_PITCH = 20*deg2rad;
CTRL_MAX_ROLL = 17*deg2rad;
CTRL_MAX_R = 20*deg2rad;


%% Generate data:

% NOTE: Right now, only the inner loop rate controller is implemented. Thus,
% only angular velocity reference commands are used. 
% -- Total time:
T = 10; % [sec]
% -- Time:
t = 0:dt:T; % [sec]

% - Reference commands:
t1 = 0:dt:2-dt; % [sec]
omega_cmd_t1 = 0*deg2rad*ones(3,length(t1));
t2 = 2:dt:4-dt; % [sec]
omega_cmd_t2 = 0*deg2rad*ones(3,length(t2));
t3 = 4:dt:8-dt; % [sec]
omega_cmd_t3 = 0*deg2rad*ones(3,length(t3));
omega_cmd_t3(2,:) = 0*deg2rad*ones(1,length(t3));
t4 = 8:dt:T; % [sec]
omega_cmd_t4 = 0*deg2rad*ones(3,length(t4));
omega_cmd_t4(2,:) = 0*deg2rad*ones(1,length(t4));
% t5 = 10:dt:16-dt; % [sec]
% omega_cmd_t5 = 0*deg2rad*ones(3,length(t5));
% omega_cmd_t5(2,:) = 0*deg2rad*ones(1,length(t5));
% t6 = 16:dt:T; % [sec]
% omega_cmd_t6 = 0*deg2rad*ones(3,length(t6));
% omega_cmd_t6(2,:) = 0*deg2rad*ones(1,length(t6));

% - Concatenate all the reference commands: 
omega_cmd = [omega_cmd_t1 omega_cmd_t2 omega_cmd_t3 omega_cmd_t4];%  omega_cmd_t5 omega_cmd_t6];

figure();
plot(t,omega_cmd(1,:),t,omega_cmd(2,:),t,omega_cmd(3,:));
title('Angular velocity reference commands');
xlabel('t [sec]');
ylabel('p,q & r [rad/s]');
legend('p','q','r');

%% Simulation and control design:

% % State variables:
% -- Position:
position = zeros(3,length(t)); % [m]
% -- Velocity:
velocity = zeros(3,length(t)); % [m/s]
% -- Quaternions:
q = zeros(4,length(t)); % -- quaternions [q_x, q_y, q_z, q_w]
% -- Attitude -- Euler angles
att = zeros(3,length(t)); % -- Euler angles [phi, theta, psi]
% -- Angular velocity:
omega = zeros(3,length(t)); % [rad/s]
omega_filt = zeros(3,length(t)); % [rad/s]
omega_dot_filt = zeros(3,length(t)); % [rad/s^2]
% -- Angular acceleration:
angAccel = zeros(3,length(t)); % [rad/s^2]
% -- Specific Thrust:
specificThrust = zeros(1,length(t)); % [N]

% % Measured states/variables:
% -- Measured angular velocity:
m_angVel = zeros(3,length(t)); % [rad/s^2]
% -- Measured acceleration:
m_accel = zeros(3,length(t)); % [m/s^2]
% -- Specific force: 
specificForce = zeros(3,length(t)); % [N/kg] or [m/s^2]
% -- Specific force in the Body frame:
specificForceBodyFrame = zeros(3,length(t)); % [N/s] or [m/s^2]
% -- Propeller speed: 
propSpeed = zeros(4,length(t)); % [rad/s]
% -- Propeller speed command: 
propSpeedComand = zeros(4,length(t)); % [rad/s]
% -- Angular acceleration command (obtained from the PID):
angAcc_cmd = zeros(3,length(t)); % [rad/s^2]
% -- PID - error integrated:
int_error = zeros(3,length(t));

% % Initial Conditions:
% -- Initial position:
position(:,1) = [0 0 1.5]';

% -- Specific force -- acceleration
specificForce_fix = [0 0 -grav]'; 
specificForce(:,1) = specificForce_fix; 

% % - Compute propeller speed:
propSpeed(:,1) = sqrt(vehicleMass/4*grav/thrustCoeff); 

% % - Initial quaternion state:
theta = att(1,1);
phi = att(2,1);
psi = att(3,1);
q(:,1) = [
         sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2); % qx
         cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2); % qy
         cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2); % qz
         cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2); % qw
         ]; % qx, qy, qz, qw

for i = 1:length(t)-1
    % % ----------- Things provided by Lockheed Martin --------------
    % 1 - Filter the measurements, in this case, the angular velocity:
    [omega_filt(:,i),omega_dot_filt(:,i)] = lpfproceedState(omega(:,i),omega_filt(:,i),omega_dot_filt(:,i),gainP,gainQ,dt);
    
    % 2 - Use the PID control law, in order to obtain the motor inputs:
    [angAcc_cmd(:,i),int_error(:,i+1)] = controlUpdate(omega_cmd(:,i),omega_filt(:,i),omega_dot_filt(:,i),int_error(:,i),intBound,dt,propGain,intGain,derGain);
    
    % 3 - Compute the motor speed commands, in order to obtain the desired
    % response: 
    propSpeedComand(:,i) = computeMotorSpeedCommand(thrust_z,angAcc_cmd(:,i),thrustCoeff,torqueCoeff,maxPropSpeed,vehicleInertia,momentArm);
    % angAcc_cmd(:,i) = 0;
    
    % 4 - Compute the state update for position, velocity, etc.
    [position(:,i+1),velocity(:,i+1),q(:,i+1),omega(:,i+1),angAccel(:,i), ...
    propSpeed(:,i+1),specificForce(:,i+1),specificForceBodyFrame(:,i+1),dragForce,specificThrust(i)] = proceedState(position(:,i), ...
    velocity(:,i),q(:,i),omega(:,i),propSpeed(:,i),specificForce(:,i),propSpeedComand(:,i),...
    thrustCoeff,torqueCoeff,dragCoeff,momentArm,vehicleInertia,vehicleMass,...
    angAccelProcessNoiseAutoCorrelation,linAccelProcessNoiseAutoCorrelation,...
    motorTimeconstant,grav,dt,includeDrag);
    att(:,i+1) = quat2Euler(q(:,i+1)); % Conversion from Quaternions to Euler angles.
    
    % 5 - Obtain measurements from the IMU:
    [m_angVel(:,i+1),m_accel(:,i+1),~,~] = IMU_meas(omega(:,i+1),specificForceBodyFrame(:,i+1),gyroMeasNoiseVariance,accMeasNoiseVariance);
    % % ----------- Things provided by Lockheed Martin --------------
    
end

%% Time to plot everything: 

% % Low pass filter - Signals:
% -- Omega:
figure();
subplot(3,1,1);
plot(t,omega_cmd(1,:),t,omega(1,:));
xlabel('t [sec]');
ylabel('p [rad/s]'); 
legend('Command','Obtained');
subplot(3,1,2);
plot(t,omega_cmd(2,:),t,omega(2,:));
xlabel('t [sec]');
ylabel('q [rad/s]'); 
subplot(3,1,3);
plot(t,omega_cmd(3,:),t,omega(3,:));
xlabel('t [sec]');
ylabel('r [rad/s]'); 
suptitle('Angular velocities of the quadrotor');   

% -- Filtered omega:
figure();
subplot(3,1,1);
plot(t,omega_cmd(1,:),t,omega_filt(1,:));
xlabel('t [sec]');
ylabel('p_{filt} [rad/s]'); 
legend('Command','Obtained');
subplot(3,1,2);
plot(t,omega_cmd(2,:),t,omega_filt(2,:));
xlabel('t [sec]');
ylabel('q_{filt} [rad/s]'); 
subplot(3,1,3);
plot(t,omega_cmd(3,:),t,omega_filt(3,:));
xlabel('t [sec]');
ylabel('r_{filt} [rad/s]'); 
suptitle('Filtered angular velocity of the quadrotor');       

% -- Filtered omega dot:
figure();
subplot(3,1,1);
plot(t,omega_dot_filt(1,:));
xlabel('t [sec]');
ylabel('pdot_{filt} [rad/s]'); 
subplot(3,1,2);
plot(t,omega_dot_filt(2,:));
xlabel('t [sec]');
ylabel('qdot_{filt} [rad/s]'); 
subplot(3,1,3);
plot(t,omega_dot_filt(3,:));
xlabel('t [sec]');
ylabel('rdot_{filt} [rad/s]'); 
suptitle('Filtered angular acceleration of the quadrotor');   

% % Inner loop control signals:
% -- Commanded angular acceleration:
figure();
subplot(3,1,1);
plot(t,angAcc_cmd(1,:));
xlabel('t [sec]');
ylabel('pdot_{filt} [rad/s^2]'); 
subplot(3,1,2);
plot(t,angAcc_cmd(2,:));
xlabel('t [sec]');
ylabel('qdot_{filt} [rad/s^2]'); 
subplot(3,1,3);
plot(t,angAcc_cmd(3,:));
xlabel('t [sec]');
ylabel('rdot_{filt} [rad/s^2]'); 
suptitle('Commanded angular acceleration (after PID) for the quadrotor'); 

% % Commanded signals for the motor inputs: 
% -- Command propeller speeds:
figure();
subplot(2,2,1);
plot(t,propSpeedComand(1,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_1} (Cmd) [rad/s]'); 
subplot(2,2,2);
plot(t,propSpeedComand(2,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_2} (Cmd) [rad/s]'); 
subplot(2,2,3);
plot(t,propSpeedComand(3,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_3} (Cmd) [rad/s]'); 
subplot(2,2,4);
plot(t,propSpeedComand(4,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_4} (Cmd) [rad/s]'); 
suptitle('Propeller speed commands for the quadrotor');

% % State update - Drone dynamics and kinematics
% -- Position:
figure();
subplot(3,1,1);
plot(t,position(1,:));
xlabel('t [sec]');
ylabel('x [m]'); 
subplot(3,1,2);
plot(t,position(2,:));
xlabel('t [sec]');
ylabel('y [m]'); 
subplot(3,1,3);
plot(t,position(3,:));
xlabel('t [sec]');
ylabel('z [m]'); 
suptitle('Quadrotor position'); 

% -- Velocity:
figure();
subplot(3,1,1);
plot(t,velocity(1,:));
xlabel('t [sec]');
ylabel('v_x [m/s]'); 
subplot(3,1,2);
plot(t,velocity(2,:));
xlabel('t [sec]');
ylabel('v_y [m/s]'); 
subplot(3,1,3);
plot(t,velocity(3,:));
xlabel('t [sec]');
ylabel('v_z [m/s]'); 
suptitle('Quadrotor velocity'); 

% -- Quadrotor Euler Angles:
figure();
subplot(3,1,1);
plot(t,att(1,:));
xlabel('t [sec]');
ylabel('phi [rad]'); 
subplot(3,1,2);
plot(t,att(2,:));
xlabel('t [sec]');
ylabel('theta [rad]'); 
subplot(3,1,3);
plot(t,att(3,:));
xlabel('t [sec]');
ylabel('psi [rad]');  
suptitle('Quadrotor Euler angles');

% -- Quaternions:
figure();
subplot(2,2,1);
plot(t,q(1,:));
xlabel('t [sec]');
ylabel('q_x [-]'); 
subplot(2,2,2);
plot(t,q(2,:));
xlabel('t [sec]');
ylabel('q_y [-]'); 
subplot(2,2,3);
plot(t,q(3,:));
xlabel('t [sec]');
ylabel('q_z [-]'); 
subplot(2,2,4);
plot(t,q(4,:));
xlabel('t [sec]');
ylabel('q_w [-]'); 
suptitle('Quadrotor quaternions');

% -- Angular acceleration
figure();
subplot(3,1,1);
plot(t,angAccel(1,:),t,m_angVel(1,:));
xlabel('t [sec]');
ylabel('\alpha_x [rad/s^2]'); 
legend('True','Measured');
subplot(3,1,2);
plot(t,angAccel(2,:),t,m_angVel(2,:));
xlabel('t [sec]');
ylabel('\alpha_y [rad/s^2]'); 
subplot(3,1,3);
plot(t,angAccel(3,:),t,m_angVel(3,:));
xlabel('t [sec]');
ylabel('\alpha_z [rad/s^2]'); 
suptitle('Quadrotor angular acceleration'); 

% -- Specific Thrust:
figure();
plot(t,specificThrust);
xlabel('t [sec]');
ylabel('Specific Thrust [N]');
title('Specific Thrust');

% -- Propeller speeds:
figure();
subplot(2,2,1);
plot(t,propSpeed(1,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_1} [rad/s]'); 
% xL = get(gca,'XLim'); % maximum propeller speed
% p2= line(xL,[maxPropSpeed maxPropSpeed],'Linestyle','-.','Color','b',...
%    'Linewidth',1);

subplot(2,2,2);
plot(t,propSpeed(2,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_2} [rad/s]'); 
% p2= line(xL,[maxPropSpeed maxPropSpeed],'Linestyle','-.','Color','b',...
%    'Linewidth',1);

subplot(2,2,3);
plot(t,propSpeed(3,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_3} [rad/s]'); 
% p2= line(xL,[maxPropSpeed maxPropSpeed],'Linestyle','-.','Color','b',...
%    'Linewidth',1);

subplot(2,2,4);
plot(t,propSpeed(4,:));
xlabel('t [sec]');
ylabel('\omega_{{prop}_4} [rad/s]'); 
suptitle('Quadrotor propeller speed');
% p2= line(xL,[maxPropSpeed maxPropSpeed],'Linestyle','-.','Color','b',...
%    'Linewidth',1);

% Specific Force -- Acceleration:
figure();
subplot(3,1,1);
plot(t,specificForce(1,:));
xlabel('t [sec]');
ylabel('a_x [m/s^2]'); 
subplot(3,1,2);
plot(t,specificForce(2,:));
xlabel('t [sec]');
ylabel('a_y [m/s^2]'); 
subplot(3,1,3);
plot(t,specificForce(3,:));
xlabel('t [sec]');
ylabel('a_z [m/s^2]'); 
suptitle('Quadrotor specific force - acceleration'); 

% Specific Force -- Acceleration:
figure();
subplot(3,1,1);
plot(t,specificForceBodyFrame(1,:),t,m_accel(1,:));
xlabel('t [sec]');
ylabel('a_x [m/s^2]'); 
legend('True','Measured');
subplot(3,1,2);
plot(t,specificForceBodyFrame(2,:),t,m_accel(2,:));
xlabel('t [sec]');
ylabel('a_y [m/s^2]'); 
subplot(3,1,3);
plot(t,specificForceBodyFrame(3,:),t,m_accel(3,:));
xlabel('t [sec]');
ylabel('a_z [m/s^2]'); 
suptitle('Quadrotor specific force - acceleration (BODY FRAME)'); 
