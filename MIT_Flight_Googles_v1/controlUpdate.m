function [out,intState_k1] = controlUpdate(command,curval,curder,intState_k0,intBound,dt_sec,propGain,intGain,derGain)
% - State deviation - error:
stateDev = command - curval;

% - Integration of the state deviation - integral of the error
intState_k1 = intState_k0 + dt_sec*stateDev;
intState_k1 = min([max([-intBound intState_k1],[],2) intBound],[],2);

% - PID control law:
out = propGain.*stateDev + intGain.*intState_k1 + derGain.*(-curder);
end