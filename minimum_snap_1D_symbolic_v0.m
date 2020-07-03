% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Version 0: Determine polynomial using minimum snap approach for 1-D case
% (trajectory between two waypoints)

clc
clear all
close all

%%

% -- Waypoints
wps = [0 4];

% -- Total time: 
t_m = 4;

% -- Time vector: 
t_vec = [0 t_m];

% -- Variable declarations: 
% - polynomial coefficients: 
syms c0 c1 c2 c3 c4 c5 c6 
% - time: 
syms t T 

% -- vector of coeficients
coef = [c6 c5 c4 c3 c2 c1 c0] ;
% -- vector of polynomial terms
tim = [t^6 t^5 t^4 t^3 t^2 t 1] ;
% -- polynomial: x = c_6.t^6 + c_5.t^5 + ... + c_1.t + c0
pol = coef*tim';

% -- 1rst derivative of the polynomial
d_pol = jacobian(pol,t);
% -- 2nd derivative of the polynomial
d2_pol = jacobian(d_pol,t);
% -- 3rd derivative of the polynomial
d3_pol = jacobian(d2_pol,t);
% -- 4th derivative of the polynomial 
d4_pol = jacobian(d3_pol,t);

% -- Minimizing snap squared:  
J = int(d4_pol^2,t,0,T);

% -- Determine hessian (d^2J/dc_ic_j):
H = hessian(J,coef)/2;

% -- Get matrix array: 
H_arr = double(subs(H,T,t_m));

% -- Define matrix A_eq and vector b_eq: 
d_tim = jacobian(tim,t)';
d2_tim = jacobian(d_tim,t)';
d3_tim = jacobian(d2_tim,t)';
d4_tim = jacobian(d3_tim,t)';

A_eq = [
    % waypoint constraints: sigma(t_i) = sigma_i
    double(subs(tim,t,t_vec(1)));
    double(subs(tim,t,t_vec(2)));
    % derivatives (in the first and last waypoints) d^{p}sigma/dt^{p} = 0 (or free)
    double(subs(d_tim,t,t_vec(1)));
    double(subs(d2_tim,t,t_vec(1)));
    % double(subs(d3_tim,t,t_vec(1)));
    % double(subs(d4_tim,t,t_vec(1)));
    double(subs(d_tim,t,t_vec(2)));
    double(subs(d2_tim,t,t_vec(2)));
    % double(subs(d3_tim,t,t_vec(2)));
    % double(subs(d4_tim,t,t_vec(2)));
    ];

b_eq = [
    wps(1);
    wps(2);
    0;
    0; 
    % 0;
    % 0;
    0;
    0;
    % 0;
    % 0;
    ];

% -- Solves quadratic programming problem: 
f = zeros(length(coef),1);
sol = quadprog(H_arr,f,[],[],A_eq,b_eq);

% -- Plots:
dt = 0.01;
t = 0:dt:t_m;

% Convert polynomials in symbolic form to matlab functions: 
tim_fun = matlabFunction(tim*sol);
d_tim_fun = matlabFunction(d_tim*sol);
d2_tim_fun = matlabFunction(d2_tim*sol);
d3_tim_fun = matlabFunction(d3_tim*sol);
d4_tim_fun = matlabFunction(d4_tim*sol);

subplot(3,2,[1 2])
plot(t,tim_fun(t));
xlabel('t [sec]');
ylabel('x [m]');
grid on;

subplot(3,2,3)
plot(t,d_tim_fun(t));
xlabel('t [sec]');
ylabel('v_x [m/s]');
grid on;

subplot(3,2,4)
plot(t,d2_tim_fun(t));
xlabel('t [sec]');
ylabel('a_x [m/s^2]');
grid on;

subplot(3,2,5)
plot(t,d3_tim_fun(t));
xlabel('t [sec]');
ylabel('j_x [m/s^3]');
grid on;

subplot(3,2,6)
plot(t,d4_tim_fun(t));
xlabel('t [sec]');
ylabel('s_x [m/s^4]');
grid on;