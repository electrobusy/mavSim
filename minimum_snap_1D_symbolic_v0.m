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
    double(subs(tim,t,t_vec(1)));
    double(subs(tim,t,t_vec(2)));
    double(subs(d_tim,t,t_vec(1)));
    double(subs(d2_tim,t,t_vec(1)));
    double(subs(d3_tim,t,t_vec(1)));
    double(subs(d4_tim,t,t_vec(1)));
    double(subs(d_tim,t,t_vec(2)));
    double(subs(d2_tim,t,t_vec(2)));
    double(subs(d3_tim,t,t_vec(2)));
    double(subs(d4_tim,t,t_vec(2)));
    ];

b_eq = [
    wps(1);
    wps(2);
    0;
    0; 
    0;
    0;
    0;
    0;
    0;
    0;
    ];

% -- Solves quadratic programming problem: 
f = zeros(length(coef),1);
sol = quadprog(H_arr,f,[],[],A_eq,b_eq);

% -- Plots:
dt = 0.01;
t = 0:dt:t_m;

% -- [THIS IS NOT WORKING, NEEDS TO BE CORRECTED]
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
ylabel('v_x [m]');
grid on;

subplot(3,2,4)
plot(t,polyval(dd_pol_x,t));
xlabel('t [sec]');
ylabel('a_x [m]');
grid on;

subplot(3,2,5)
plot(t,polyval(ddd_pol_x,t));
xlabel('t [sec]');
ylabel('j_x [m]');
grid on;

subplot(3,2,6)
plot(t,polyval(dddd_pol_x,t));
xlabel('t [sec]');
ylabel('s_x [m]');
grid on;