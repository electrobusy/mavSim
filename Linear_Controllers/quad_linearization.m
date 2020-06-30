function [A,B] = quad_linearization(x,u,data)

m = data.m;

% ----
% position
p = x(1:3);
% velocity 
v = x(4:6);
% attitude (quaternions)
q = x(7:end);

% angular velocities
omega = u(2:4);
c = u(1);
% ----

% vector norm:
nq = norm(q);

% partial derivative with respect to unit quaternions:
ddqu_dq = (eye(4) - q*q'*nq^(-2))/nq;

% partial derivatives: 
% -- position: 
ddp_dv = eye(3);
% -- orientation:
x_mat = @(x1,x2,x3) [
    0,  -x1, -x2, -x3;
    x1,   0,  x3, -x2;
    x2, -x3,   0,  x1;
    x3,  x2,  -x1,  0;
    ];
ddq_dq = 1/2*x_mat(omega(1),omega(2),omega(3))*ddqu_dq; 
ddq_domega = 1/2*[
    -q(2), -q(3), -q(4);
     q(1), -q(4), -q(3);
     q(4),  q(1),  q(2);
    -q(3),  q(2),  q(1);
    ];
% -- velocity:
aux_mat_ddv_dq = [
     q(3),  q(4),  q(1), q(2);
    -q(2), -q(1),  q(4), q(3);
     q(1), -q(2), -q(3), q(4);
    ];
ddv_dq = 2*c/m*aux_mat_ddv_dq*ddqu_dq;

ddv_dc = 1/m*[
    q(1)*q(3) + q(2)*q(4);
    q(3)*q(4) - q(1)*q(2);
    (q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2) ;
    ];

% Linearized matrices: 
A = [
        zeros(3),     ddp_dv, zeros(3,4);
        zeros(3),   zeros(3),     ddv_dq;
      zeros(4,3), zeros(4,3),     ddq_dq;
    ];

B = [
      zeros(3,1),   zeros(3); 
          ddv_dc,   zeros(3);
      zeros(4,1), ddq_domega;    
    ];
end