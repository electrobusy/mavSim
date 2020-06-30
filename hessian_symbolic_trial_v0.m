% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Version 0: Determine the Hessian matrix for the minimum snap approach 
% using one polynomial through symbolic variables (1-D)

clc
clear all
close all

%%

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
pol = coef*tim' ;

% -- 4th derivative of the polynomial 
d4_pol = jacobian(jacobian(jacobian(jacobian(pol,t),t),t),t) ;

% -- Minimizing snap squared:  
J = int(d4_pol^2,t,0,T);

% -- Determine hessian (d^2J/dc_ic_j):
H = hessian(J,coef)/2;

% -- Get matrix array: 
H_arr = double(subs(H,T,4));
