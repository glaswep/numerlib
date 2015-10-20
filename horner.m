function v = horner(c,xx,x)

%	HORNER evaluates a polynomial with Newton basis by a
%	Horner-like scheme
%
%   INPUT:
%           c: vector of Newton coefficients
%           xx: vector of interpolation points
%           x: vector of points to be evaluated
%
%   OUTPUT:
%           v: vector of value at interpolation points. (should be f)
%
%	See also:
%		Newtoninterpolate.m, newtonup.m
%
%	Yulun Zeng, Oct, 2015

n = length(c); 
v = c(1)*ones(size(x));

for i=2:n
   v = (x-xx(n-i+1)).*v + c(i);
end


end