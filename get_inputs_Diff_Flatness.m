% function u = get_inputs_Diff_Flatness(coeffs,data)
data.m = 0.38905;
data.beta = 0.5;
data.g = 9.81;    
data.L=0.2;
data.I=diag([0.0049,0.0049,0.0069]); % kg m^2
data.k_F=1.91e-6;
data.k_M=2.6e-7;

% Parameters
% -- vehicle: 
m = data.m; % [kg] quadrotor mass
I = data.I; % [kg.m^2] quadrotor moment of inertia
L = data.L; % [m] (assuming + configuration) distance from axis of rotation of rotors to quad. CoM
k_F = data.k_F; % [?] rotor force coeficient 
k_M = data.k_M; % [?] rotor moment coeficient
% -- environment: 
g = data.g; % [m/s^2] gravity
T=linspace(0,12,nr_wp); %placeholder t0 and t_end of each segment
% placeholder polynomial
nr_wp=7; %nr waypoints
n_derivative=4;
coeffs=zeros(5,4,nr_wp-1,n_derivative+1);
coeffs(:,:,:,1)=randn(5,4,nr_wp-1,1); %polynomial coefficients [coefficientss,state (x,y,z,psi), nr polynomials,(nth-1 derivative)]


multiplier=repmat((1:1:size(coeffs,1)-1)',[1,size(coeffs,[2,3])]); % used to multiply the coefficient with the corresponding power when differentiating
for n = 1:4
coeffs([1:size(coeffs,1)-1],:,:,n+1)=coeffs(2:size(coeffs,1),:,:,n).*multiplier; 
end

coeffs_x=coeffs(:,1,:,:);
coeffs_y=coeffs(:,2,:,:);
coeffs_z=coeffs(:,3,:,:);
coeffs_psi=coeffs(:,4,:,:);

dummy=2;
% end
state