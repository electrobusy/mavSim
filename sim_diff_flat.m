%% get trajectory
plotflag=1;

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
    



%% Get input vaLues from differential flatness
T=[0,t_m]; %Timestamps
[t,u,poly_coeffs]=get_inputs_Diff_Flatness(poly_coeffs,T,data);

%% Simulate 
state_0=[keyframes(1:3,1); 0;0;0; 0;0;0; 0;0;0]; % [pos, R,V, Omega];  

dt=mean(diff(t));
state=state_0;
statel=zeros(size(state,1),length(t)+1);
statel(:,1)=state_0;
for i=1:length(t)
    
    ds=quad_dynamics_Mellinger_Kumar(state,u(i,:)',data);   
    state=state+(ds*dt);
    statel(:,i+1)=state;
end

figure()
plot3(statel(1,:),statel(2,:),-statel(3,:))
grid on 