function T = get_rotationmatrix(euler,type)
phi=euler(:,1);
theta=euler(:,2);
psi=euler(:,3);
if strcmp(type,'E2B') || strcmp(type,'B2E')
c11 = reshape(cos(theta).*cos(psi),[1,1,length(theta)]);
c12 = reshape(cos(theta).*sin(psi),[1,1,length(theta)]);
c13 = reshape(-sin(theta),[1,1,length(theta)]);
c21 = reshape(sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi),[1,1,length(theta)]);
c22 = reshape(sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi),[1,1,length(theta)]);
c23 = reshape(sin(phi).*cos(theta),[1,1,length(theta)]);
c31 = reshape(cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi),[1,1,length(theta)]);
c32 = reshape(cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi),[1,1,length(theta)]);
c33 = reshape(cos(phi).*cos(theta),[1,1,length(theta)]);
end

if strcmp(type,'E2B')
    T=[
    c11, c12, c13; 
    c21, c22, c23;
    c31, c32, c33;
    ];
elseif strcmp(type,'B2E')
    T=[c11, c21, c31;
        c12, c22, c32;
        c13, c23, c33;
        ];
    
elseif strcmp(type,'B2I') % body to inertial ([omega_x ,omega_y,omega_z]-> [phi_dot, theta_dot, psi_dot])
    r11 = ones(1,1,length(theta));
    r12 = reshape(tan(theta).*sin(phi),[1,1,length(theta)]);
    r13 = reshape(tan(theta).*cos(phi),[1,1,length(theta)]);
    r21 = zeros(1,1,length(theta));
    r22 = reshape(cos(phi),[1,1,length(theta)]);
    r23 = reshape(-sin(phi),[1,1,length(theta)]);
    r31 = zeros(1,1,length(theta));
    r32 = reshape(sin(phi)./cos(theta),[1,1,length(theta)]);
    r33 = reshape(cos(phi)./cos(theta),[1,1,length(theta)]);
    T = [
    r11, r12, r13;
    r21, r22, r23;
    r31, r32, r33;
    ];
elseif strcmp(type,'I2B')
    r11 = ones(1,1,length(theta));
    r12 = reshape(tan(theta).*sin(phi),[1,1,length(theta)]);
    r13 = reshape(tan(theta).*cos(phi),[1,1,length(theta)]);
    r21 = zeros(1,1,length(theta));
    r22 = reshape(cos(phi),[1,1,length(theta)]);
    r23 = reshape(-sin(phi),[1,1,length(theta)]);
    r31 = zeros(1,1,length(theta));
    r32 = reshape(sin(phi)./cos(theta),[1,1,length(theta)]);
    r33 = reshape(cos(phi)./cos(theta),[1,1,length(theta)]);
    T = [
    r11, r21, r31;
    r12, r22, r32;
    r13, r23, r33;
    ];
elseif strcmp(type,'E2V') %earth to velocity frame 
    c11=reshape(cos(psi),[1,1,length(psi)]);
    c12=reshape(sin(psi),[1,1,length(psi)]);
    c13=zeros(1,1,length(psi));
    c21=reshape(-sin(psi),[1,1,length(psi)]);
    c22=reshape(cos(psi),[1,1,length(psi)]);
    c23=c13;
    c31=c13;
    c32=c13;
    c33=ones(1,1,length(psi));
    T=[
    c11, c12, c13;
    c21, c22, c23;
    c31, c32, c33;
    ];
elseif strcmp(type,'V2E') %velocity to earth frame 
    c11=reshape(cos(psi),[1,1,length(psi)]);
    c12=reshape(sin(psi),[1,1,length(psi)]);
    c13=zeros(1,1,length(psi));
    c21=reshape(-sin(psi),[1,1,length(psi)]);
    c22=reshape(cos(psi),[1,1,length(psi)]);
    c23=c13;
    c31=c13;
    c32=c13;
    c33=ones(1,1,length(psi));
    T=[
    c11, c21, c31;
    c12, c22, c32;
    c13, c23, c33;
    ];
    
else
    error("Unable to resolve rotationmatrix type");
end

end