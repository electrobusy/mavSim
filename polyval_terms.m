% vector with terms of the polynomial evaluated at point x
function [res] = polyval_terms(p,x)

m = length(p);

% order of the polynomial with coefficients 
n = m - 1; 

% vector of powers
pow = flip(0:1:n);

% x^n + x^{n-1} + ... + x^1 + 1
aux = x.^pow;

% c_nx^n + c_{n-1}x^{n-1} + ... + c_1x^1 + c_0
res = aux.*p;
    
end