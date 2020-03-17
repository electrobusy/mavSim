% ---------
% This function yields the outputs of the IMU sensor already with noise (in 
% this case - AWGN - Artificial White Gaussian Noise)
% ---------
function [measured_angVel,measured_accel,measured_angVel_Cov,measured_acc_Cov] = IMU_meas(angVel,accel,gyroMeasNoiseVariance,accMeasNoiseVariance)

measured_angVel = angVel; % + sqrt(gyroMeasNoiseVariance)*normrnd(0,1,[3 1]);
measured_accel = accel; % + sqrt(accMeasNoiseVariance)*normrnd(0,1,[3 1]);

measured_angVel_Cov = zeros(9,1);
measured_acc_Cov = zeros(9,1);

for i = 1:9
    if i == 1 || i == 5 || i == 9
        measured_angVel_Cov(i) = gyroMeasNoiseVariance;
        measured_acc_Cov(i) = accMeasNoiseVariance;
    else
        measured_angVel_Cov(i) = 0;
        measured_acc_Cov(i) = 0;
    end
end

end