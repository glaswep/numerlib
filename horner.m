function v = horner(c,xx,x)

%   INPUT:
%           c: vector of Newton coefficients
%           xx: vector of interpolation points
%           x: vector of points to be evaluated
%
%   OUTPUT:
%           v: vector of value at interpolation points. (should be f)

n = length(c); 
v = c(1)*ones(size(x));

for i=2:n
   v = (x-xx(n-i+1)).*v + c(i);
end


end