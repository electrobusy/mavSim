close all
clear all
data=readmatrix("discretize_traject.txt");
inputs=readmatrix("u_diff.txt");
t=data(:,1);

traj.x=data(:,2);
traj.dx=data(:,3);
traj.ddx=data(:,4);
traj.dddx=data(:,5);
traj.ddddx=data(:,6);

traj.y=data(:,7);
traj.dy=data(:,8);
traj.ddy=data(:,9);
traj.dddy=data(:,10);
traj.ddddy=data(:,11);

traj.z=data(:,12);
traj.dz=data(:,13);
traj.ddz=data(:,14);
traj.dddz=data(:,15);
traj.ddddz=data(:,16);

traj.psi=data(:,17);
traj.dpsi=data(:,18);
traj.ddpsi=data(:,19);
traj.dddpsi=data(:,20);
traj.ddddpsi=data(:,21);

testy= [0,0,0,0.148148148223096,0.0148148148247847,0.000395061728794061,-2.94313291698298e-15];
y2=zeros(size(traj.z));
for i=1:length(testy)
    y2=y2+testy(i)*t.^(i-1);
end

figure()
plot3(traj.x,traj.y,traj.z)
xlabel('x')
ylabel('y')
zlabel('z')
grid on 

figure()
subplot(141)
plot(t,inputs(:,1));
grid on 
title("Thrust [N]") 

subplot(142)
plot(t,inputs(:,2));
grid on 
title("\tau\_x [Nm]")  

subplot(143)
plot(t,inputs(:,3));
grid on 
title("\tau\_y [Nm]") 

subplot(144)
plot(t,inputs(:,4));
grid on
title("\tau\_z [Nm]") 