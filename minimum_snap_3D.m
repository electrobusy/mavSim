% MSc Student: Rohan Chotalal
% Supervisors: Guido De Croon and Christophe D' Waghter
% Faculty: Aerospace Engineering, TU Delft.

% Minimum snap algorithm (Mellinger and Kumar) - Code 

% -- keyframes = [x y z psi]'
keyframes = [
    0, 0, 0, 0;
    0, 5, 0, pi/2;
]'; % <-- we transpose this vector, to make easier to port in the for loops

% -- number of keyframes: 
[~,m] = size(keyframes);

% -- polynomial order: 
n = 7; % choose 

% -- total time: 
t_m = 4; % [sec]

% -- vector of times:
t = [0, t_m]; % [t_0, t_1, ..., t_m]
t_vec = t(2:end-1); % [t_1, ..., t_{m-1}]

% -- derivatives to minimize: 
k_r = 4; % snap
k_psi = 2; % angular acceleration (yaw)

% -----------------------------------------------------------------------
% --> Variables that could not be used (maybe!)(randomly declared just to 
% understand the problem in hand):

% - polynomial coefficients (to be determined): 
x_T_ij = [ones(m,1) zeros(m,n)];
y_T_ij = [ones(m,1) zeros(m,n)];
z_T_ij = [ones(m,1) zeros(m,n)];
psi_T_ij = [ones(m,1) zeros(m,n)];

% - polynomial and its derivatives
% pol = cell(1,k+1);
% for idx_der=1:k % k represents the number of derivatives
%     pol_der{1,idx_der} = zeros(m,n-idx_der);
%     for idx_rows = 1:m
%         der = polyder();
%         if length(der) == n-k
%             pol_der{idx_der}(:,) = der(2:end);
%         else 
%         end    
%     end
% end

% -- Define H matrix 

% -- Define A matrix and b vector 

% -- Define A_eq matrix and b_eq vector 

% -- Solves quadratic programming problem: 
% x = quadprog(2*H,f,A,b,Aeq,beq)



