% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Minimum Snap Algorithm (Mellinger and Kumar) - 1-D for m Waypoints with 
% a 6-th order polynomial

% clc,
% clear all,
close all

% -- keyframes = [x]'
keyframes = [0 5 10];

% -- number of keyframes: 
[~,m] = size(keyframes);

% -- polynomial order: 
n = 6; % choose order 

% -- total time: 
t_m = 4; % [sec]

% -- vector of times:
t = [0, 2, t_m]; % [t_0, t_1, ..., t_m]
% t_vec = t(2:end-1); % [t_1, ..., t_{m-1}]

% -- derivatives to minimize: 
k_x = 4; % snap

%% -----------------------------------------------------------------------

% Let's consider a polynomial: x(t) = t^n + t^{n-1} + ... + t^2 + t + 1 
x = ones(1,n+1);

% Think that in each term of the polynomial, the coeficient c_i is being
% multiplied, where i = {n, ... ,0}
dx = polyder(x);
ddx = polyder(dx);
dddx = polyder(ddx);
ddddx = polyder(dddx);

% -- Construct A_eq matrix and b_eq vector
% - polynomial 1 
A_p1 = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(x,t(1)); % t_0
    polyval_terms(x,t(2)); % t_1
    % derivatives (in the first waypoints) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dx,t(1)) 0;
    polyval_terms(ddx,t(1)) 0 0;
    % polyval_terms(dddx,t(1)) 0 0 0;
    % polyval_terms(ddddx,t(1)) 0 0 0 0; 
    ];

b_p1 = [
    keyframes(1);
    keyframes(2);
    0; 
    0; 
    % inf; 
    % inf;
    ];

% - polynomial 2 
A_p2 = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(x,t(2)); % t_1
    polyval_terms(x,t(end)); % t_end
    % derivatives (in the last waypoint) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dx,t(end)) 0;
    polyval_terms(ddx,t(end)) 0 0;
    % polyval_terms(dddx,t(end)) 0 0 0;
    % polyval_terms(ddddx,t(end)) 0 0 0 0; 
    ];

b_p2 = [
    keyframes(2);
    keyframes(end);
    0; 
    0; 
    % inf; 
    % inf;
    ];

% - continuity constraints: 
A_cont = [
    polyval_terms(x,t(2)), -polyval_terms(x,t(2));
    polyval_terms(dx,t(2)) 0, -polyval_terms(dx,t(2)) 0;
    polyval_terms(ddx,t(2)) 0 0, -polyval_terms(ddx,t(2)) 0 0;
    polyval_terms(dddx,t(2)) 0 0 0, -polyval_terms(dddx,t(2)) 0 0 0;
    polyval_terms(ddddx,t(2)) 0 0 0 0, -polyval_terms(ddddx,t(2)) 0 0 0 0;
    ];

% - get A_eq and b_eq
A_eq = blkdiag(A_p1,A_p2);
b_eq = [b_p1; b_p2];

A_eq = [
    A_eq; 
    A_cont
    ];

b_eq = [
    b_eq
    0;
    0;
    0;
    0;
    0
    ];

% -- Define A matrix and b vector 
% [NO INEQUALITY CONTRAINTS IN THIS CASE]

% -- Define Hessian matrix
% - polynomial 1
H_p1 = zeros(n+1,n+1);

for i = 1:length(ddddx)
    aux = zeros(1,length(ddddx));
    aux(i) = ddddx(i);
    conv_pol = conv(ddddx,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H_p1(i,1:length(ddddx)) = res(res > 0);
end

% - polynomial 2
H_p2 = zeros(n+1,n+1);

for i = 1:length(ddddx)
    aux = zeros(1,length(ddddx));
    aux(i) = ddddx(i);
    conv_pol = conv(ddddx,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H_p2(i,1:length(ddddx)) = res(res > 0);
end

% - final hessian
H = blkdiag(H_p1,H_p2);

% -- Solves quadratic programming problem: 
f = zeros((m-1)*(n+1),1);
sol = quadprog(H,f,[],[],A_eq,b_eq);

% -- Plots:
dt = 0.01;
t_1 = 0:dt:t(2);
t_2 = t(2):dt:t_m;

pol_x_1 = x.*sol(1:n+1)';
d_pol_x_1 = dx.*sol(1:n)';
dd_pol_x_1 = ddx.*sol(1:n-1)';
ddd_pol_x_1 = dddx.*sol(1:n-2)';
dddd_pol_x_1 = ddddx.*sol(1:n-3)';
pol_x_2 = x.*sol(n+2:end)';
d_pol_x_2 = dx.*sol(n+2:end-1)';
dd_pol_x_2 = ddx.*sol(n+2:end-2)';
ddd_pol_x_2 = dddx.*sol(n+2:end-3)';
dddd_pol_x_2 = ddddx.*sol(n+2:end-4)';

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