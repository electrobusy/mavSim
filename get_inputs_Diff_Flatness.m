function u = get_inputs_Diff_Flatness(coeffs,data)
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
    T=[0,2,3.4,5,9,10,11]; %placeholder t0 and t_end of each segment
    % placeholder polynomial
    nr_wp=length(T); %nr waypoints
    nr_derivative=4;
    nr_states=4;
    poly_order=5;

    %% build array with polynomial coefficients
    coeffs=zeros(poly_order,nr_states,nr_wp-1,nr_derivative+1);
    coeffs(:,:,:,1)=randn(poly_order,nr_states,nr_wp-1,1); %polynomial coefficients [coefficientss,state (x,y,z,psi), nr sections,(nth-1 derivative)]
    multiplier=repmat((1:1:size(coeffs,1)-1)',[1,size(coeffs,[2,3])]); % used to multiply the coefficient with the corresponding power when differentiating
    for n = 1:4
    coeffs([1:size(coeffs,1)-1],:,:,n+1)=coeffs(2:size(coeffs,1),:,:,n).*multiplier; %calculate coefficients of derivatives and shift so that the index corresponds with power-1
    end



    [t,out]=discretize_poly(coeffs,T,1e-2); % dimension out: [polynomial outputs for ti,state,section,derivative-1]
                                     %dimension t: [ti,power,section] 

    %% get u1 (collective thrust) 
    % thrust is simply derived from m*r_dot_dot=-m*g*z_w+u1*zb
    %zb is derived from the assumption thrust is aligned with zb and with
    %r_dot_dot * R_B|W
    t=out(:,1:3,:,3); %2nd derivative of [x,y,z];
    t(:,3,:)=t(:,3,:)+(ones(size(t,1),1,size(t,3))*data.g);
    z_b=t./vecnorm(t,2,2); %z_b is the unitvector of t (vecnorm takes the 2-norm for each x_dot_dot,y_dot_dot,z_dot_dot vector 
    z_w=[0,0,1];
    u1=data.m*vecnorm(t,2,2); %collective thrust input, goes together with u2...u4 in thrust mixing matrix to get motor speeds

    x_c=[cos(out(:,4,:,1)), sin(out(:,4,:,1)),zeros(size(out,1),1,size(out,3))]; % x unit vector of intermediate frame (only rotated for yaw)
    y_b=cross(z_b,x_c,2);
    y_b=y_b./vecnorm(y_b,2,2);
    x_b=cross(y_b,z_b,2);
 %TODO may have made errors when mutliplying with z_b 
    h_omega=(m./u1).*out(:,1:3,:,4)-(dot(z_b,out(:,1:3,:,4),2).*z_b);  %out(:,1:3,:,4) = a_dot = r_dot_dot_dot (jerk)
    p=dot(-h_omega,y_b,2); %roll rate 
    q=dot(h_omega,x_b,2); % pitch rate 
    r=dot(out(:,4,:,2).*z_w,z_b,2);      % yaw rate, not sure if z_w is used correctly
    omega_bw=[p,q,r];  %body angular rate 

    a_dot_dot=out(:,1:3,:,5);
    a_dot=out(:,1:3,:,4);
    %%TODO not sure whether the derivation for h_alpha is correct:
    h_alpha=m*a_dot_dot-(dot(z_b,m*a_dot_dot,2).*z_b)-cross(omega_bw,dot(z_b,m*a_dot,2).*z_b,2)...
        - cross(omega_bw,(dot(z_b,m*a_dot,2).*z_b+cross(omega_bw,u1.*z_b,2)),2);

    p_dot=dot(-h_alpha,y_b,2); %roll acceleration 
    q_dot=dot(h_alpha,x_b,2); % pitch acceleration
    r_dot=dot(out(:,4,:,3).*z_w,z_b,2); %yaw acceleration
    alpha_bw=[p_dot,q_dot,r_dot];
    %% Get u2...u4 from alpha_bw
    I=data.I'; %transpose because all vectors are horizontal (doesn't matter if I is diagonal ofcourse) 
    term_1=zeros(size(alpha_bw)); % assigned for matrix multiplication of I and alpha_bw
    term_2=zeros(size(alpha_bw)); % assigned for matrix multiplication of I and omega_bw
    for i=1:size(term_1,1)
        for s =1:size(term_1,3)
            term_1(i,:,s)=(I*alpha_bw(i,:,s)')';
            term_2(i,:,s)=(I*omega_bw(i,:,s)')';
        end
    end
    u24=term_1+cross(omega_bw,term_2,2);

    u=[u1,u24];
   
    function [t,poly] =  discretize_poly(coeffs,T,dt)
    nr_sections=size(coeffs,3);
    poly_order=size(coeffs,1);
    nr_states=size(coeffs,2);
    nr_derivative=size(coeffs,4)-1;

    t=nan((max(diff(T)))/dt+1,poly_order,nr_sections); %each section may have a different length, the ends are padded with nans
    % build time array to multiply coefficients with
    poly=nan(size(t,1),nr_states,nr_sections,nr_derivative+1);
    for s = 1:nr_sections
        ts=T(s):dt:T(s+1);
        for n=1:size(t,2)
            t(1:length(ts),n,s)=ts.^(n-1);
        end
        for d=1:nr_derivative+1
            poly(:,:,s,d)=t(:,:,s)*coeffs(:,:,s,d);
        end
    end

    end

end