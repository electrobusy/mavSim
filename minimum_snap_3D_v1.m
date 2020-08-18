% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Minimum Snap Algorithm (Mellinger and Kumar) - 3-D going from 1 Waypoint
% to another using a 6-th order polynomial

clc,
clear all,
close all
if not(exist('plotflag','var'))
    plotflag = 1; 
end

global keyframes;
% -- keyframes = [x y z psi]' (where x, y, z and psi are column vectors)
keyframes = [
        0,0,0,0     ;
        4,0,5,0;
        5,0,5,0.5*pi;
        7,2,5,0
    ]';

% -- total time: 


% -- vector of times:
% t = [0, 3, t_m]; % [t_0, t_1, ..., t_m]
t=linspace(0,25,size(keyframes,2));
% t=[0,11,20,30];
t=[0,10,20,30];
T0=diff(t);
m=length(t);

optifunct=@(T)opti_function_section_time(T);

A1 = -eye(m-1,m-1) ;
B1 = zeros(m-1,1) ;
Aeq = [] ;
beq = [] ;
lb = [] ;
ub = [] ;
nonlcon = [] ;
maxiter = 10 ;
Tn=T0
% Matlab optimiser
options = optimoptions('fmincon','Display','iter',...
    'PlotFcn','optimplotfval','MaxIterations',maxiter) ;
% [Tn,fval,exitflag,output] = fmincon(optifunct,T0,A1,B1,Aeq,beq,lb,ub,nonlcon,options) ;
% 
t=0;
    for i = 1:length(Tn)
        t=[t,sum(Tn(1:i))];
    end

[cost,poly_coeffs]=generate_minimum_snap_trajectory(t ,keyframes,1);

function cost = opti_function_section_time(T)
    global keyframes
    t=0;
    for i = 1:length(T)
        t=[t,sum(T(1:i))];
    end
      [cost,~]=generate_minimum_snap_trajectory(t,keyframes,0);
end

% t_vec = t(2:end-1); % [t_1, ..., t_{m-1}]

