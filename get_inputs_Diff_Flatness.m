function [time,u, coeffs,extra,out] = get_inputs_Diff_Flatness(coeffs,T,data)
    %$ this function calculates the required inputs to follow a polynomial
    %trajectory
    
    %inputs:
    %coeffs: polynomial coefficients
    %dimension: 
    %1=polynomial order 0-n
    %2=state [x,y,z,psi]
    %3=sections 
    %4=derivative 0-k (not required, is calculated here from 0th derivative)
    
    %data: struct containing system parameters, m,g,L,I,k_f,k_m 
    %---output: discretized time vector, inputs u, d
    %dimensions: 
    %1=input values,
    %2=u1 ... u4 (collective thrust, moment_x,moment_y,moment_z
    %3=section (valid for sections between waypoints)
 
   

    
    % -- vehicle: 
    m = data.m; % [kg] quadrotor mass
    I = data.I; % [kg.m^2] quadrotor moment of inertia
    L = data.L; % [m] (assuming + configuration) distance from axis of rotation of rotors to quad. CoM
    k_F = data.k_F; % [?] rotor force coeficient 
    k_M = data.k_M; % [?] rotor moment coeficient
    % -- environment: 
    g = data.g; % [m/s^2] gravity
%     T=[0,2,3.4,5,9,10,11]; %placeholder t0 and t_end of each segment
    % placeholder polynomial
    nr_wp=length(T); %nr waypoints
    nr_derivative=4; %how many derivatives (excluding 0th)
    nr_states=4;
    poly_order=5;

    %% build array with polynomial coefficients
    %don't forget to remove coeffs
%     coeffs=zeros(poly_order,nr_states,nr_wp-1,nr_derivative+1);
%     coeffs(:,:,:,1)=randn(poly_order,nr_states,nr_wp-1,1); %polynomial coefficients [coefficientss,state (x,y,z,psi), nr sections,(nth-1 derivative)]
    multiplier=repmat((1:1:size(coeffs,1)-1)',[1,size(coeffs,[2,3])]); % used to multiply the coefficient with the corresponding power when differentiating
    for n = 1:nr_derivative
    coeffs([1:size(coeffs,1)-1],:,:,n+1)=coeffs(2:size(coeffs,1),:,:,n).*multiplier; %calculate coefficients of derivatives and shift so that the index corresponds with power-1
    end



    [time,out]=discretize_poly(coeffs,T,data.dt); % dimension out: [polynomial outputs for ti,state,section,derivative-1]
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
    y_c=[-sin(out(:,4,:,1)), cos(out(:,4,:,1)),zeros(size(out,1),1,size(out,3))];
    y_b=cross(z_b,x_c,2);
    y_b=y_b./vecnorm(y_b,2,2);
    x_b=cross(y_b,z_b,2);
 %TODO may have made errors when mutliplying with z_b 
    h_omega=(m./u1).*out(:,1:3,:,4)-(dot(z_b,out(:,1:3,:,4),2).*z_b);  %out(:,1:3,:,4) = a_dot = r_dot_dot_dot (jerk)
    p=dot(-h_omega,y_b,2); %roll rate 
    q=dot(h_omega,x_b,2); % pitch rate 
    r=dot(out(:,4,:,2).*z_w,z_b,2);      % yaw rate,gives different values than other method below
    omega_bw=[p,q,r];  %body angular rate 
    
    a_dot_dot=out(:,1:3,:,5); %snap
    a_dot=out(:,1:3,:,4); %jerk 
    %% this is likely wrong: 
    h_alpha=m*a_dot_dot-(dot(z_b,m*a_dot_dot,2).*z_b)-cross(omega_bw,dot(z_b,m*a_dot,2).*z_b,2)...
        - cross(omega_bw,(dot(z_b,m*a_dot,2).*z_b+cross(omega_bw,u1.*z_b,2)),2);

    p_dot=dot(-h_alpha,y_b,2); %roll acceleration 
    q_dot=dot(h_alpha,x_b,2); % pitch acceleration
    r_dot=dot(out(:,4,:,3).*z_w,z_b,2); %yaw acceleration
    
    psi_dot=out(:,4,:,2); %heading rate 
    psi_dot_dot=out(:,4,:,3); %heading acc
    
    %% Different derivation of thrust (includes drag terms)
    alpha=out(:,1:3,:,3);

    %% Different derivation of angular acceleration: from http://arxiv.org/abs/1712.02402 without drag terms
    c=u1./m;
    gzw=zeros(size(z_b));
    gzw(:,3)=data.g;
    c2=dot(z_b,out(:,1:3,1,3)+gzw,2); %mindcheck
    c_dot=dot(z_b,a_dot,2);
    B1=c;
    C1=zeros(size(c));
    D1=dot(x_b,a_dot,2);
    A2=c; 
    C2=C1;
    D2=dot(-y_b,a_dot,2);
    B3=dot(-y_c,z_b,2);
    C3=vecnorm(cross(y_c,z_b),2,2);
    D3=psi_dot.*dot(x_c,x_b,2);
    
    omega_x=(-B1.*C2.*D3+B1.*C3.*D2-B3.*C1.*D2+B3.*C2.*D1)./(A2.*(B1.*C3-B3.*C1));
    omega_y=(-C1.*D3+C3.*D1)./(B1.*C3-B3.*C1);
    omega_z=(B1.*D3-B3.*D1)./(B1.*C3-B3.*C1); %values of omega_x and omega_y are identical to kumar's method, but omega_z is different
    E1=dot(x_b,a_dot_dot,2)-2*c_dot.*omega_y-c.*omega_x.*omega_z;
    E2=dot(-y_b,a_dot_dot,2)-2*c_dot.*omega_x+c.*omega_y.*omega_z;
    E3=psi_dot_dot.*dot(x_c,x_b,2)+2*psi_dot.*omega_z.*dot(x_c,y_b,2)-2*psi_dot.*omega_y.*dot(x_c,z_b,2)-omega_x.*omega_y.*dot(y_c,y_b,2)-omega_x.*omega_z.*dot(y_c,z_b,2);
    
    omega_x_dot=(-B1.*C2.*E3+B1.*C3.*E2-B3.*C1.*E2+B3.*C2.*E1)./(A2.*(B1.*C3-B3.*C1));
    omega_y_dot=(-C1.*E3+C3.*E1)./(B1.*C3-B3.*C1);
    omega_z_dot=(B1.*E3-B3.*E1)./(B1.*C3-B3.*C1);
    omega_dot=[omega_x_dot,omega_y_dot,omega_z_dot];
    %% Get u2...u4 from alpha_bw
    I=data.I'; %transpose because all vectors are horizontal (doesn't matter if I is diagonal ofcourse) 
    term_1=zeros(size(omega_dot)); % assigned for matrix multiplication of I and alpha_bw
    term_2=zeros(size(omega_dot)); % assigned for matrix multiplication of I and omega_bw
    for i=1:size(term_1,1)
        for s =1:size(term_1,3)
            term_1(i,:,s)=(I*omega_dot(i,:,s)')';
            term_2(i,:,s)=(I*omega_bw(i,:,s)')';
        end
    end
    u24=term_1+cross(omega_bw,term_2,2);

    u=[u1,u24];
   
    % calculate trajectory euler angles (for simulation initial state and debugging purposes)
    
    extra.Rb=cat(2,reshape(x_b',3,1,size(x_b,1)),reshape(y_b',3,1,size(y_b,1)),reshape(z_b',3,1,size(z_b,1)));%orientation matrix
    extra.phi=atan2(extra.Rb(3,2,:),extra.Rb(3,3,:)); %from https://www.geometrictools.com/Documentation/EulerAngles.pdf
    extra.phi=reshape(extra.phi,size(extra.phi,3),1);
    extra.theta=atan2(-extra.Rb(3,1,:),sqrt(extra.Rb(3,2,:).^2+extra.Rb(3,3,:).^2));
    extra.theta=reshape(extra.theta,size(extra.theta,3),1);
    extra.psi=atan2(extra.Rb(2,1,:),extra.Rb(1,1,:));
    extra.psi=reshape(extra.psi,size(extra.psi,3),1);
    extra.psi2=atan2(x_c(:,2),x_c(:,1));
    extra.omega_dot=omega_dot; %body angular rates.
    extra.omega=[omega_x,omega_y,omega_z];
    
    
    %% write coefficients to file
    nr_section=size(coeffs,3);
    nr_coeffs=size(coeffs,1);
    nr_states=size(coeffs,2);
    nr_derivatives=size(coeffs,4);
    coeffile=fopen("coeffcients.txt",'wt');
    
    for sect = 1:nr_section
        for state= 1:nr_states
            for d=1:nr_derivatives
                for c=1:nr_coeffs
                    if c<nr_coeffs
                        fprintf(coeffile,"%.9f, ",coeffs(c,state,sect,d));
                    else
                        fprintf(coeffile,"%.9f\n",coeffs(c,state,sect,d));
                    end
                end
            end
            fprintf(coeffile,"\n");
        end
    end
    
    
    
end