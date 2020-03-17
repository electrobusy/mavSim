function propSpeedCommand = computeMotorSpeedCommand(thrust_z,angAccComand,thrustCoeff,torqueCoeff,maxPropSpeed,vehicleInertia,momentArm)

% % - Compute  momentThrust:
momentThrust = [
                diag(vehicleInertia).*angAccComand;
                thrust_z;
               ];
           
% % - Compute motorSpeedsSquared:
aux = 1/(4*momentArm);
aux_C_T = aux/thrustCoeff;
aux_C_tau = aux/torqueCoeff;
G = [
    aux_C_T -aux_C_T -aux_C_tau*momentArm aux_C_T*momentArm;
    aux_C_T aux_C_T aux_C_tau*momentArm aux_C_T*momentArm;
    -aux_C_T aux_C_T -aux_C_tau*momentArm aux_C_T*momentArm;
    -aux_C_T -aux_C_T aux_C_tau*momentArm aux_C_T*momentArm;
    ];
 motorSpeedsSquared = G*momentThrust;
 
 % % - Prop speed commands (assuming that rotor speeds >= 0) 
 propSpeedCommand = min([sqrt(max([zeros(4,1) motorSpeedsSquared],[],2)) maxPropSpeed*ones(4,1)],[],2);

end