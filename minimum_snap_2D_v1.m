% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Minimum Snap Algorithm (Mellinger and Kumar) - 2-D using m Waypoints with
% a 6-th order polynomial

clc,
clear all,
close all

% -- keyframes = [x y]' (where x and y are column vectors)
keyframes = [
        0,0;
        4,5;
        10,6;
    ]';

% -- number of keyframes: 
[~,m] = size(keyframes);

% -- polynomial order: 
n = 6; % choose order 

% -- total time: 
t_m = 5; % [sec]

% -- vector of times:
t = [0, 3, t_m]; % [t_0, t_1, ..., t_m]
% t_vec = t(2:end-1); % [t_1, ..., t_{m-1}]

% -- derivatives to minimize: 
k_x = 4; % snap

%% -----------------------------------------------------------------------

% Let's consider a polynomial: x(t) = t^n + t^{n-1} + ... + t^2 + t + 1 
x = ones(1,n+1);
% Same for y(t): 
y = ones(1,n+1);

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

% -- Construct A_eq matrix and b_eq vector
% -> for x: 
% - polynomial 1 
A_x_p1 = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(x,t(1)); % t_0
    polyval_terms(x,t(2)); % t_1
    % derivatives (in the first waypoints) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dx,t(1)) 0;
    polyval_terms(ddx,t(1)) 0 0;
    % polyval_terms(dddx,t(1)) 0 0 0;
    % polyval_terms(ddddx,t(1)) 0 0 0 0; 
    ];

b_x_p1 = [
    keyframes(1,1);
    keyframes(1,2);
    0; 
    0; 
    % inf; 
    % inf;
    ];

% - polynomial 2 
A_x_p2 = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(x,t(2)); % t_1
    polyval_terms(x,t(end)); % t_end
    % derivatives (in the last waypoint) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dx,t(end)) 0;
    polyval_terms(ddx,t(end)) 0 0;
    % polyval_terms(dddx,t(end)) 0 0 0;
    % polyval_terms(ddddx,t(end)) 0 0 0 0; 
    ];

b_x_p2 = [
    keyframes(1,2);
    keyframes(1,3);
    0; 
    0; 
    % inf; 
    % inf;
    ];

% - continuity constraints: 
A_x_cont = [
    polyval_terms(x,t(2)), -polyval_terms(x,t(2));
    polyval_terms(dx,t(2)) 0, -polyval_terms(dx,t(2)) 0;
    polyval_terms(ddx,t(2)) 0 0, -polyval_terms(ddx,t(2)) 0 0;
    polyval_terms(dddx,t(2)) 0 0 0, -polyval_terms(dddx,t(2)) 0 0 0;
    polyval_terms(ddddx,t(2)) 0 0 0 0, -polyval_terms(ddddx,t(2)) 0 0 0 0;
    ];

b_x_cont = [
    0;
    0;
    0;
    0;
    0;
    ];

% - A_x:
A_x_aux = blkdiag(A_x_p1,A_x_p2);
A_x = [
    A_x_aux; 
    A_x_cont
    ];

b_x = [
    b_x_p1;
    b_x_p2;
    b_x_cont
    ];

% -> for y: 
% - polynomial 1 
A_y_p1 = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(y,t(1)); % t_0
    polyval_terms(y,t(2)); % t_1
    % derivatives (in the first waypoints) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dy,t(1)) 0;
    polyval_terms(ddy,t(1)) 0 0;
    % polyval_terms(dddy,t(1)) 0 0 0;
    % polyval_terms(ddddy,t(1)) 0 0 0 0; 
    ];

b_y_p1 = [
    keyframes(2,1);
    keyframes(2,2);
    0; 
    0; 
    % inf; 
    % inf;
    ];

% - polynomial 2 
A_y_p2 = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(y,t(2)); % t_1
    polyval_terms(y,t(end)); % t_end
    % derivatives (in the last waypoint) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dy,t(end)) 0;
    polyval_terms(ddy,t(end)) 0 0;
    % polyval_terms(dddy,t(end)) 0 0 0;
    % polyval_terms(ddddy,t(end)) 0 0 0 0; 
    ];

b_y_p2 = [
    keyframes(2,2);
    keyframes(2,end);
    0; 
    0; 
    % inf; 
    % inf;
    ];

% - continuity constraints: 
A_y_cont = [
    polyval_terms(y,t(2)), -polyval_terms(y,t(2));
    polyval_terms(dy,t(2)) 0, -polyval_terms(dy,t(2)) 0;
    polyval_terms(ddy,t(2)) 0 0, -polyval_terms(ddy,t(2)) 0 0;
    polyval_terms(dddy,t(2)) 0 0 0, -polyval_terms(dddy,t(2)) 0 0 0;
    polyval_terms(ddddy,t(2)) 0 0 0 0, -polyval_terms(ddddy,t(2)) 0 0 0 0;
    ];

b_y_cont = [
    0;
    0;
    0;
    0;
    0;
    ];

% - A_y:
A_y_aux = blkdiag(A_y_p1,A_y_p2);
A_y = [
    A_y_aux; 
    A_y_cont
    ];

b_y = [
    b_y_p1;
    b_y_p2;
    b_y_cont
    ];

% - Obtain A_eq using the block matrix form and concatenate b_x and b_y to
% get b_eq:
A_eq = blkdiag(A_x,A_y);
b_eq = [b_x; b_y];

% -- Define A matrix and b vector 
% [NO INEQUALITY CONTRAINTS IN THIS CASE]

% -- Define Hessian matrix
H_x = zeros(n+1,n+1);
H_y = zeros(n+1,n+1);

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

H = blkdiag(H_x,H_x,H_y,H_y);

% -- Solves quadratic programming problem: 
f = zeros((m-1)*2*(n+1),1);
sol = quadprog(H,f,[],[],A_eq,b_eq);

% -- Plots:
dt = 0.01;
t_1 = 0:dt:t(2);
t_2 = t(2):dt:t(end);

pol_x_1 = x.*sol(1:n+1)';
d_pol_x_1 = dx.*sol(1:n+1-1)';
dd_pol_x_1 = ddx.*sol(1:n+1-2)';
ddd_pol_x_1 = dddx.*sol(1:n+1-3)';
dddd_pol_x_1 = ddddx.*sol(1:n+1-4)';

pol_x_2 = y.*sol(n+2:2*(n+1))';
d_pol_x_2 = dy.*sol(n+2:2*(n+1)-1)';
dd_pol_x_2 = ddy.*sol(n+2:2*(n+1)-2)';
ddd_pol_x_2 = dddy.*sol(n+2:2*(n+1)-3)';
dddd_pol_x_2 = ddddy.*sol(n+2:2*(n+1)-4)';

pol_y_1 = y.*sol(2*(n+1)+1:3*(n+1))';
d_pol_y_1 = dy.*sol(2*(n+1)+1:3*(n+1)-1)';
dd_pol_y_1 = ddy.*sol(2*(n+1)+1:3*(n+1)-2)';
ddd_pol_y_1 = dddy.*sol(2*(n+1)+1:3*(n+1)-3)';
dddd_pol_y_1 = ddddy.*sol(2*(n+1)+1:3*(n+1)-4)';

pol_y_2 = y.*sol(3*(n+1)+1:4*(n+1))';
d_pol_y_2 = dy.*sol(3*(n+1)+1:4*(n+1)-1)';
dd_pol_y_2 = ddy.*sol(3*(n+1)+1:4*(n+1)-2)';
ddd_pol_y_2 = dddy.*sol(3*(n+1)+1:4*(n+1)-3)';
dddd_pol_y_2 = ddddy.*sol(3*(n+1)+1:4*(n+1)-4)';

%% 
figure(1);
subplot(3,2,[1 2])
plot(t_1,polyval(pol_x_1,t_1));
hold on; 
plot(t_2,polyval(pol_x_2,t_2));
xlabel('t [sec]');
ylabel('x [m]');
grid on;

subplot(3,2,3)
plot(t_1,polyval(d_pol_x_1,t_1));
hold on; 
plot(t_2,polyval(d_pol_x_2,t_2));
xlabel('t [sec]');
ylabel('v_x [m/s]');
grid on;

subplot(3,2,4)
plot(t_1,polyval(dd_pol_x_1,t_1));
hold on; 
plot(t_2,polyval(dd_pol_x_2,t_2));
xlabel('t [sec]');
ylabel('a_x [m/s^2]');
grid on;

subplot(3,2,5)
plot(t_1,polyval(ddd_pol_x_1,t_1));
hold on; 
plot(t_2,polyval(ddd_pol_x_2,t_2));
xlabel('t [sec]');
ylabel('j_x [m/s^3]');
grid on;

subplot(3,2,6)
plot(t_1,polyval(dddd_pol_x_1,t_1));
hold on; 
plot(t_2,polyval(dddd_pol_x_2,t_2));
xlabel('t [sec]');
ylabel('s_x [m/s^4]');
grid on;

figure(2);
subplot(3,2,[1 2])
plot(t_1,polyval(pol_y_1,t_1));
hold on; 
plot(t_2,polyval(pol_y_2,t_2));
xlabel('t [sec]');
ylabel('y [m]');
grid on;

subplot(3,2,3)
plot(t_1,polyval(d_pol_y_1,t_1));
hold on; 
plot(t_2,polyval(d_pol_y_2,t_2));
xlabel('t [sec]');
ylabel('v_y [m/s]');
grid on;

subplot(3,2,4)
plot(t_1,polyval(dd_pol_y_1,t_1));
hold on; 
plot(t_2,polyval(dd_pol_y_2,t_2));
xlabel('t [sec]');
ylabel('a_y [m/s^2]');
grid on;

subplot(3,2,5)
plot(t_1,polyval(ddd_pol_y_1,t_1));
hold on; 
plot(t_2,polyval(ddd_pol_y_2,t_2));
xlabel('t [sec]');
ylabel('j_y [m/s^3]');
grid on;

subplot(3,2,6)
plot(t_1,polyval(dddd_pol_y_1,t_1));
hold on; 
plot(t_2,polyval(dddd_pol_y_2,t_2));
xlabel('t [sec]');
ylabel('s_y [m/s^4]');
grid on;

%% 
figure(3);
V_1 = sqrt(polyval(d_pol_x_1,t_1).^2 + polyval(d_pol_y_1,t_1).^2);
V_2 = sqrt(polyval(d_pol_x_2,t_2).^2 + polyval(d_pol_y_2,t_2).^2);
V = [V_1 V_2];
x_1_data = polyval(pol_x_1,t_1);
y_1_data = polyval(pol_y_1,t_1);
x_2_data = polyval(pol_x_2,t_2);
y_2_data = polyval(pol_y_2,t_2);
x_data = [x_1_data x_2_data];
y_data = [y_1_data y_2_data];
surf([x_data(:) x_data(:)], [y_data(:) y_data(:)], [V(:) V(:)], ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 2);            % Make a thicker line
view(2);   % Default 2-D view
c = colorbar;  % Add a colorbar
c.Label.String = 'V [m/s]';
xlabel('x [m]');
ylabel('y [m]');
grid on;