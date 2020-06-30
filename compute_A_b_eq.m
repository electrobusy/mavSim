% Based on the equality constraints: 

% Function definition: Compute matrices A_eq and b_eq.
% Inputs: 
% n - polynomial order 
% m - number of waypoints 
% k_r - position derivative
% k_psi - yaw derivative
% wp - vector of waypoints
% t_wp - vector of time at each waypoint 

function [A_eq,b_eq] = compute_A_b_eq(n,m,k_r,k_psi,wp,t_wp)

pol = ones(1,n+1);

% -- obtain k_r-th derivative of x/y/z polynomial 
dk_r_pol = pol;
for i = 1:k_r
    dk_r_pol = polyder(dk_r_pol); 
end

% -- obtain k_psi-th derivative of psi polynomial 
dk_psi_pol = pol;
for i = 1:k_psi
    dk_psi_pol = polyder(dk_psi_pol);
end 

% Construct A_eq matrix: 
A_eq = [];

A_x = zeros(n+1,n+1);
A_y = zeros(n+1,n+1);
A_z = zeros(n+1,n+1);
A_psi = zeros(n+1,n+1);

for i = 1:n+1 % rows
    for j = 1:n+1 % columns
        A_x(i,j) = 
    end
end

% Construct b_eq vector: 
b_eq = zeros(n+1,1);

end