% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Minimum Snap Algorithm (Mellinger and Kumar) - 3-D going from 1 Waypoint
% to another using a 6-th order polynomial

clc,
% clear all,
close all
if not(exist('plotflag','var'))
    plotflag = 1; 
end


% -- keyframes = [x y z psi]' (where x, y, z and psi are column vectors)
keyframes = [
        0,0,0,0     ;
        4,0,5,-0.5*pi;
        5,0,5,0.5*pi;
        7,2,5,pi
    ]';

% -- number of keyframes: 
[~,m] = size(keyframes);

% -- polynomial order: 
n = 6; % choose order 

% -- total time: 
% t_m = 5; % [sec]

% -- vector of times:
% t = [0, 3, t_m]; % [t_0, t_1, ..., t_m]
t=linspace(0,10,size(keyframes,2));
% t=[0,5,20,30];
t=[0,7,12,15];
t_m=t(end);
% t_vec = t(2:end-1); % [t_1, ..., t_{m-1}]

% -- derivatives to minimize: 
k_x = 4; % snap
k_psi = 2; % angular yaw acc. 

%% -----------------------------------------------------------------------

% Let's consider a polynomial: x(t) = t^n + t^{n-1} + ... + t^2 + t + 1 
x = ones(1,n+1);

% Same for y(t): 
y = ones(1,n+1);
% Same for z(t):
z = ones(1,n+1);
% Same for psi(t):
psi = ones(1,n+1);

% Think that in each term of the polynomial, the coeficient c_i is being
% multiplied, where i = {n, ... ,0}
% - for polynomial x
dx = polyder(x);
ddx = polyder(dx);
dddx = polyder(ddx);
ddddx = polyder(dddx);
% - for polynomial y: 
dy = polyder(y);
ddy = polyder(dy);
dddy = polyder(ddy);
ddddy = polyder(dddy);
% - for polynomial z: 
dz = polyder(z);
ddz = polyder(dz);
dddz = polyder(ddz);
ddddz = polyder(dddz);
% - for polynomial psi: 
dpsi = polyder(psi);
ddpsi = polyder(dpsi);


[A_x,b_x]=get_A_eq_b_eqxyz(keyframes,t,x,dx,ddx,dddx,ddddx,'x');
[A_y,b_y]=get_A_eq_b_eqxyz(keyframes,t,y,dy,ddy,dddy,ddddy,'y');
[A_z,b_z]=get_A_eq_b_eqxyz(keyframes,t,x,dz,ddz,dddz,ddddz,'z');
[A_psi,b_psi]=get_A_eq_b_eqxyz(keyframes,t,psi,dpsi,ddpsi,[],[],'psi');

A_eq = blkdiag(A_x, A_y, A_z, A_psi);
b_eq = [b_x; b_y; b_z; b_psi];

% -- Define A matrix and b vector 
% [NO INEQUALITY CONTRAINTS IN THIS CASE]

% -- Define Hessian matrix
H_x = zeros(n+1,n+1);
H_y = zeros(n+1,n+1);
H_z = zeros(n+1,n+1);
H_psi = zeros(n+1,n+1);

% - for x: 
for i = 1:length(ddddx)
    aux = zeros(1,length(ddddx));
    aux(i) = ddddx(i);
    conv_pol = conv(ddddx,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H_x(i,1:length(ddddx)) = res(res > 0);
end

% - for y:
for i = 1:length(ddddy)
    aux = zeros(1,length(ddddy));
    aux(i) = ddddy(i);
    conv_pol = conv(ddddy,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H_y(i,1:length(ddddy)) = res(res > 0);
end

% - for z:
for i = 1:length(ddddz)
    aux = zeros(1,length(ddddz));
    aux(i) = ddddz(i);
    conv_pol = conv(ddddz,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H_z(i,1:length(ddddz)) = res(res > 0);
end

% - for psi:
for i = 1:length(ddpsi)
    aux = zeros(1,length(ddpsi));
    aux(i) = ddpsi(i);
    conv_pol = conv(ddpsi,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H_z(i,1:length(ddpsi)) = res(res > 0);
end

H=H_x;
for i=2:(m-1)
    H=blkdiag(H,H_x);
end
for i=1:(m-1)
    H=blkdiag(H,H_y);
end

for i=1:(m-1)
    H=blkdiag(H,H_z);
end

for i=1:(m-1)
    H=blkdiag(H,H_psi);
end



% H = blkdiag(H_x,H_x,H_y,H_y,H_z,H_z,H_psi,H_psi);

% -- Solves quadratic programming problem: 
options=optimoptions('quadprog','Algorithm','interior-point-convex','ConstraintTolerance',1e-8,'MaxIterations',5000,'StepTolerance',1e-12);
f = zeros(size(A_eq,2),1);
[sol,cost,fl,output] = quadprog(H,f,[],[],A_eq,b_eq,[],[],[],options);

% -- Plots:
dt = 0.01;
t_1 = 0:dt:t(2);
t_2 = t(2):dt:t(end);

pol_x=[];
d_pol_x=[];
dd_pol_x=[];
ddd_pol_x=[];
dddd_pol_x=[];

pol_y=[];
d_pol_y=[];
dd_pol_y=[];
ddd_pol_y=[];
dddd_pol_y=[];

pol_z=[];
d_pol_z=[];
dd_pol_z=[];
ddd_pol_z=[];
dddd_pol_z=[];

pol_psi=[];
d_pol_psi=[];
dd_pol_psi=[];


poly_coeffs=zeros(max([length(pol_x),length(pol_y),length(pol_z),length(pol_psi)]),4,m-1,5); 

% Isolate coefficients from solution
t_s.t_0=0;
for i = 1:m-1
    pol_x=[pol_x ;x.*sol(1+(n+1)*(i-1):(n+1)*i)';];
    d_pol_x=[ d_pol_x ;dx.*sol(1+(n+1)*(i-1):(n+1)*i-1)';];
    dd_pol_x=[ dd_pol_x ;ddx.*sol(1+(n+1)*(i-1):(n+1)*i-2)';];
    ddd_pol_x=[ddd_pol_x ;dddx.*sol(1+(n+1)*(i-1):(n+1)*i-3)';];
    dddd_pol_x=[dddd_pol_x ;ddddx.*sol(1+(n+1)*(i-1):(n+1)*i-4)';];     
    
    pol_y=[pol_y ;y.*sol((m-1)*(n+1)+1+(n+1)*(i-1):((m-1)*(n+1)+(n+1)*i))';];
    d_pol_y=[ d_pol_y ;dy.*sol((m-1)*(n+1)+1+(n+1)*(i-1):((m-1)*(n+1)+(n+1)*i)-1)';];
    dd_pol_y=[ dd_pol_y ;ddy.*sol((m-1)*(n+1)+1+(n+1)*(i-1):((m-1)*(n+1)+(n+1)*i)-2)';];
    ddd_pol_y=[ddd_pol_y ;dddy.*sol((m-1)*(n+1)+1+(n+1)*(i-1):((m-1)*(n+1)+(n+1)*i)-3)';];
    dddd_pol_y=[dddd_pol_y ;ddddy.*sol((m-1)*(n+1)+1+(n+1)*(i-1):((m-1)*(n+1)+(n+1)*i)-4)';];  
%     
    pol_z=[pol_z ;z.*sol(2*(m-1)*(n+1)+1+(n+1)*(i-1):(2*(m-1)*(n+1)+(n+1)*i))';];
    d_pol_z=[ d_pol_z ;dz.*sol(2*(m-1)*(n+1)+1+(n+1)*(i-1):(2*(m-1)*(n+1)+(n+1)*i-1))';];
    dd_pol_z=[ dd_pol_z ;ddz.*sol(2*(m-1)*(n+1)+1+(n+1)*(i-1):(2*(m-1)*(n+1)+(n+1)*i-2))';];
    ddd_pol_z=[ddd_pol_z ;dddz.*sol(2*(m-1)*(n+1)+1+(n+1)*(i-1):(2*(m-1)*(n+1)+(n+1)*i-3))';];
    dddd_pol_z=[dddd_pol_z ;ddddz.*sol(2*(m-1)*(n+1)+1+(n+1)*(i-1):(2*(m-1)*(n+1)+(n+1)*i-4))';];  
    
    pol_psi=[pol_psi ;psi.*sol(3*(m-1)*(n+1)+1+(n+1)*(i-1):(3*(m-1)*(n+1)+(n+1)*i))';];
    d_pol_psi=[ d_pol_psi ;dpsi.*sol(3*(m-1)*(n+1)+1+(n+1)*(i-1):(3*(m-1)*(n+1)+(n+1)*i-1))';];
    dd_pol_psi=[ dd_pol_psi ;ddpsi.*sol(3*(m-1)*(n+1)+1+(n+1)*(i-1):(3*(m-1)*(n+1)+(n+1)*i-2))';];
    t0=t_s.(strcat("t_",num2str(i-1)));
    t0=t0(end);
    t_s.(strcat("t_",num2str(i)))= t0:dt:t(i+1);
    
    
    % put the coefficients in a format as read by the differential flatness code: 
    poly_coeffs=augment_arrays(poly_coeffs,flip(pol_x(i,:)),1,i,1);
    poly_coeffs=augment_arrays(poly_coeffs,flip(d_pol_x(i,:)),1,i,2);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dd_pol_x(i,:)),1,i,3);
    poly_coeffs=augment_arrays(poly_coeffs,flip(ddd_pol_x(i,:)),1,i,4);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dddd_pol_x(i,:)),1,i,5);
    
    poly_coeffs=augment_arrays(poly_coeffs,flip(pol_y(i,:)),2,i,1);
    poly_coeffs=augment_arrays(poly_coeffs,flip(d_pol_y(i,:)),2,i,2);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dd_pol_y(i,:)),2,i,3);
    poly_coeffs=augment_arrays(poly_coeffs,flip(ddd_pol_y(i,:)),2,i,4);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dddd_pol_y(i,:)),2,i,5);
    
    poly_coeffs=augment_arrays(poly_coeffs,flip(pol_z(i,:)),3,i,1);
    poly_coeffs=augment_arrays(poly_coeffs,flip(d_pol_z(i,:)),3,i,2);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dd_pol_z(i,:)),3,i,3);
    poly_coeffs=augment_arrays(poly_coeffs,flip(ddd_pol_z(i,:)),3,i,4);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dddd_pol_z(i,:)),3,i,5);
    
    poly_coeffs=augment_arrays(poly_coeffs,flip(pol_psi(i,:)),4,i,1);
    poly_coeffs=augment_arrays(poly_coeffs,flip(d_pol_psi(i,:)),4,i,2);
    poly_coeffs=augment_arrays(poly_coeffs,flip(dd_pol_psi(i,:)),4,i,3);
    
end







if plotflag
% -- plot x
figure();
subplot(3,2,[1 2])
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(pol_x(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('x [m]');
grid on;

subplot(3,2,3)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(d_pol_x(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('v_x [m/s]');
grid on;

subplot(3,2,4)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dd_pol_x(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('a_x [m/s^2]');
grid on;

subplot(3,2,5)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(ddd_pol_x(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('j_x [m/s^3]');
grid on;

subplot(3,2,6)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dddd_pol_x(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('s_x [m/s^4]');
grid on;

% -- plot y
figure();
subplot(3,2,[1 2])
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(pol_y(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('y [m]');
grid on;

subplot(3,2,3)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(d_pol_y(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('v_y [m/s]');
grid on;

subplot(3,2,4)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dd_pol_y(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('a_y [m/s^2]');
grid on;

subplot(3,2,5)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(ddd_pol_y(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('j_y [m/s^3]');
grid on;

subplot(3,2,6)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dddd_pol_y(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('s_y [m/s^4]');
grid on;

% -- plot z
figure();
subplot(3,2,[1 2])
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(pol_z(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('z [m]');
grid on;

subplot(3,2,3)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(d_pol_z(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('v_z [m/s]');
grid on;

subplot(3,2,4)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dd_pol_z(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('a_z [m/s^2]');
grid on;

subplot(3,2,5)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(ddd_pol_z(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('j_z [m/s^3]');
grid on;

subplot(3,2,6)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dddd_pol_z(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('s_z [m/s^4]');
grid on;

% -- plot psi
figure();
subplot(2,2,[1 2])
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(pol_psi(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('\psi [rad]');
grid on;

subplot(2,2,3)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(d_pol_psi(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('d\psi [rad/s]');
grid on;

subplot(2,2,4)
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
plot(ti,polyval(dd_pol_psi(i,:),ti));
hold on;
end
xlabel('t [sec]');
ylabel('dd\psi [rad/s^2]');
grid on;

%% 
figure();
V=[];
x_data=[];
y_data=[];
z_data=[];
for i=1:m-1
    ti=t_s.(strcat("t_",num2str(i)));
V=[V,sqrt(polyval(d_pol_x(i,:),ti).^2 + polyval(d_pol_y(i,:),ti).^2 + polyval(d_pol_z(i,:),ti).^2)];
x_data=[x_data,  polyval(pol_x(i,:),ti)];
y_data=[y_data,  polyval(pol_y(i,:),ti)];
z_data=[z_data,  polyval(pol_z(i,:),ti)];
hold on;
end


surf([x_data(:) x_data(:)], [y_data(:) y_data(:)], [z_data(:) z_data(:)], ...
    [V(:) V(:)], ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 2);            % Make a thicker line
% view(2);   % Default 2-D view
c = colorbar;  % Add a colorbar
c.Label.String = 'V [m/s]';
xlabel('x [m]');
ylabel('y [m]');
grid on;
for i = 1:m
    hold on;
    scatter3(keyframes(1,i),keyframes(2,i),keyframes(3,i),'ro');
end
view([1,1,1]);
axis equal
end

function coeffs = augment_arrays(coeffs,poly,statedim,sectiondim,derivativedim)
    coeffs(1:length(poly),statedim,sectiondim,derivativedim)=poly;
end