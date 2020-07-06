% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Minimum Snap Algorithm (Mellinger and Kumar) - 1-D going from 1 Waypoint
% to another using a 6-th order polynomial

% clc,
clear all,
close all

% -- keyframes = [x]'
keyframes = [0 10];

% -- number of keyframes: 
[~,m] = size(keyframes);

% -- polynomial order: 
n = 6; % choose order 

% -- total time: 
t_m = 2; % [sec]

% -- vector of times:
t = [0, t_m]; % [t_0, t_1, ..., t_m]
% t_vec = t(2:end-1); % [t_1, ..., t_{m-1}]

% -- derivatives to minimize: 
k_x = 4; % snap

%% -----------------------------------------------------------------------

% Let's consider a polynomial: x(t) = t^n + t^{n-1} + ... + t^2 + t + 1 
x = ones(m-1,n+1);

% Think that in each term of the polynomial, the coeficient c_i is being
% multiplied, where i = {n, ... ,0}
dx = polyder(x);
ddx = polyder(dx);
dddx = polyder(ddx);
ddddx = polyder(dddx);

% -- Construct A_eq matrix and b_eq vector
A_eq = [
    % waypoint constraints: sigma(t_i) = sigma_i
    polyval_terms(x,t(1)); % t_0
    polyval_terms(x,t(2)); % t_m
    % derivatives (in the first and last waypoints) d^{p}sigma/dt^{p} = 0 (or free)
    polyval_terms(dx,t(1)) 0;
    polyval_terms(ddx,t(1)) 0 0;
    % polyval_terms(dddx,t(1)) 0 0 0;
    % polyval_terms(ddddx,t(1)) 0 0 0 0; 
    polyval_terms(dx,t(2)) 0;
    polyval_terms(ddx,t(2)) 0 0;
    % polyval_terms(dddx,t(2)) 0 0 0;
    % polyval_terms(ddddx,t(2)) 0 0 0 0; 
    ];

b_eq = [
    keyframes(1);
    keyframes(2);
    0; 
    0; 
    % inf; 
    % inf;
    2; 
    0; 
    % inf; 
    % inf;
    ];

% -- Define A matrix and b vector 
% [NO INEQUALITY CONTRAINTS IN THIS CASE]

% -- Define Hessian matrix
H = zeros(n+1,n+1);

for i = 1:length(ddddx)
    aux = zeros(1,length(ddddx));
    aux(i) = ddddx(i);
    conv_pol = conv(ddddx,aux);
    aux_int_pol = polyint(conv_pol);
    res = polyval_terms(aux_int_pol,t_m);
    H(i,1:length(ddddx)) = res(res > 0);
end

% -- Solves quadratic programming problem: 
f = zeros(n+1,1);
sol = quadprog(H,f,[],[],A_eq,b_eq);

% -- Plots:
dt = 0.01;
t = 0:dt:t_m;

pol_x = x.*sol';
d_pol_x = dx.*sol(1:end-1)';
dd_pol_x = ddx.*sol(1:end-2)';
ddd_pol_x = dddx.*sol(1:end-3)';
dddd_pol_x = ddddx.*sol(1:end-4)';

subplot(3,2,[1 2])
plot(t,polyval(pol_x,t));
xlabel('t [sec]');
ylabel('x [m]');
grid on;

subplot(3,2,3)
plot(t,polyval(d_pol_x,t));
xlabel('t [sec]');
ylabel('v_x [m/s]');
grid on;

subplot(3,2,4)
plot(t,polyval(dd_pol_x,t));
xlabel('t [sec]');
ylabel('a_x [m/s^2]');
grid on;

subplot(3,2,5)
plot(t,polyval(ddd_pol_x,t));
xlabel('t [sec]');
ylabel('j_x [m/s^3]');
grid on;

subplot(3,2,6)
plot(t,polyval(dddd_pol_x,t));
xlabel('t [sec]');
ylabel('s_x [m/s^4]');
grid on;



