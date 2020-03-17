function euler = quat2Euler(q)

% -- Quaternions
a = q(4);
b = q(1);
c = q(2);
d = q(3);

% -- Conversion from Quaternions to Euler:
phi = atan2(2*(a*b + c*d),1-2*(b^2 + c^2));
theta = asin(2*(a*c - d*b));
psi = atan2(2*(a*d + b*c), 1-2*(c^2 + d^2));

euler = [phi,theta,psi]';

end