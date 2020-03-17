function [position_k1,velocity_k1,attitude_k1,angVelocity_k1,angAccel_k0,...
    propSpeed_k1,specificForce_k1,specificForceBodyFrame_k1,dragForce,specificThrust] = proceedState(position_k0, ...
    velocity_k0,attitude_k0,angVelocity_k0,propSpeed_k0,specificForce_k0,propSpeedComand,...
    thrustCoeff,torqueCoeff,dragCoeff,momentArm,vehicleInertia,vehicleMass,...
    angAccelProcessNoiseAutoCorrelation,linAccelProcessNoiseAutoCorrelation,...
    motorTimeconstant,grav,dt_secs,includeDrag)

% % Explicit Euler integration:
% -- Position:
position_k1 = position_k0 + dt_secs*velocity_k0;
% -- Velocity:
velocity_k1 = velocity_k0 + dt_secs*specificForce_k0;
velocity_k1(3) = velocity_k1(3) - dt_secs*grav; % Gravity influence

% % Attitude derivative: 
% -- Rotational tranformation matrix
S_Q = [
        attitude_k0(4) -attitude_k0(3) attitude_k0(2);
        attitude_k0(3) attitude_k0(4) -attitude_k0(1);
        -attitude_k0(2) attitude_k0(1) attitude_k0(4);
        -attitude_k0(1) -attitude_k0(2) -attitude_k0(3);
      ];
% -- Compute attitude derivative in the quaternions coordinate frame:
attitudeDer_k0 = 0.5*S_Q*angVelocity_k0;

% % Explicit Euler integration:
% -- attitude:
attitude_k1 = attitude_k0 + dt_secs*attitudeDer_k0;

% % Quaternions normalization factor:
attNorm = sqrt(sum(attitude_k1.^2));

% % Normalize the quaternions:
if attNorm > 1  
    attitude_k1 = attitude_k1/attNorm;
end

% % Process noise:
% -- Angular Acceleration:
angAccelProcessNoise = sqrt(angAccelProcessNoiseAutoCorrelation/dt_secs)*normrnd(0,1,[3,1]);
% -- Linear Acceleration:
linAccelProcessNoise = sqrt(linAccelProcessNoiseAutoCorrelation/dt_secs)*normrnd(0,1,[3,1]);

% % Compute angular acceleration: -- use the Euler equation
angAccel_k0(1) = momentArm*thrustCoeff/vehicleInertia(1,1)*(propSpeed_k0(1)^2 + propSpeed_k0(2)^2 - propSpeed_k0(3)^2 - propSpeed_k0(4)^2) + angAccelProcessNoise(1);
angAccel_k0(2) = momentArm*thrustCoeff/vehicleInertia(2,2)*(-propSpeed_k0(1)^2 + propSpeed_k0(2)^2 + propSpeed_k0(3)^2 - propSpeed_k0(4)^2) + angAccelProcessNoise(2);
angAccel_k0(3) = torqueCoeff/vehicleInertia(3,3)*(-propSpeed_k0(1)^2 + propSpeed_k0(2)^2 - propSpeed_k0(3)^2 + propSpeed_k0(4)^2) + angAccelProcessNoise(3);

angAccel_k0 = angAccel_k0';

% % Explicit Euler integration:
% -- Angular velocity:
angVelocity_k1 = angVelocity_k0 + dt_secs*angAccel_k0;

% % Propeller speed:
propSpeed_k1 = propSpeed_k0 + min([dt_secs motorTimeconstant])*((propSpeedComand - propSpeed_k0)/motorTimeconstant);

% % Vehicle specific thrust:
specificThrust = thrustCoeff*(sum(propSpeed_k1.^2))/vehicleMass;

% % Vehicle total speed:
speed = sqrt(sum(velocity_k1.^2));

% % [Added this myself - won't influence anything in the code]
dragForce = zeros(3,1);

% % Include the drag model: 
if (includeDrag && speed > 0)
    % % Total drag:
    drag = dragCoeff*speed^2;
    % % Drag force:
    dragForce = drag*(-velocity_k1)/speed + linAccelProcessNoise;
    
    % % Express the drag in the body frame: 
    a = attitude_k1(4);
    b = attitude_k1(1);
    c = attitude_k1(2);
    d = attitude_k1(3);
    
    % Transformation matrix: 
    R = [
        (a^2 + b^2 - c^2 - d^2) (2*b*c + 2*a*d) (2*b*d - 2*a*c);
        (2*b*c - 2*a*d) (a^2 - b^2 + c^2 - d^2) (2*c*d + 2*a*b);
        (2*b*d + 2*a*c) (2*c*d - 2*a*b) (a^2 - b^2 - c^2 + d^2);
        ];
    specificForceBodyFrame_k1 = R*dragForce;
    specificForce_k1 = dragForce;
else
    % % Express the drag in the body frame: 
    a = attitude_k1(4);
    b = attitude_k1(1);
    c = attitude_k1(2);
    d = attitude_k1(3);
    
    % Transformation matrix: 
    R = [
        (a^2 + b^2 - c^2 - d^2) (2*b*c + 2*a*d) (2*b*d - 2*a*c);
        (2*b*c - 2*a*d) (a^2 - b^2 + c^2 - d^2) (2*c*d + 2*a*b);
        (2*b*d + 2*a*c) (2*c*d - 2*a*b) (a^2 - b^2 - c^2 + d^2);
        ];
    specificForceBodyFrame_k1 = R*linAccelProcessNoise;
    specificForce_k1 = linAccelProcessNoise;
end

specificForceBodyFrame_k1(3) = specificForceBodyFrame_k1(3) + specificThrust;

specificForce_k1 = specificForce_k1 + [
                2*attitude_k1(1)*attitude_k1(3) + 2*attitude_k1(2)*attitude_k1(4);
                2*attitude_k1(2)*attitude_k1(3) - 2*attitude_k1(1)*attitude_k1(4);
                -attitude_k1(1)^2 - attitude_k1(2)^2 + attitude_k1(3)^2 + attitude_k1(4)^2;
                ]*specificThrust;            
end