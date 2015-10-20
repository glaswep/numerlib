function c = NewtonInterpolate(x,y)

% 	NEWTONINTERPOLATE find newton's interpolations points.
%
%   INPUT:
%           x: vector of x-coordinate of interpolation points
%           y: vector of y-coordinate of interpolation points
%
%   OUTPUT:
%           c: vector of Newton's interpolation points (reverse order)
%
%	See also:
%		horner.m, newtonup.m	
%
%	Yulun Zeng, Oct, 2015

n = length(x);
c = y;
for i = 2:n
    c(i:n) = (c(i:n) - c(i-1))./(x(i:n) - x(i-1));
end
c = fliplr(c);

end