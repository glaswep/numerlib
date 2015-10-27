function c = NewtonInterpolate(x,y)

%
%   INPUT:
%           x: vector of x-coordinate of interpolation points
%           y: vector of y-coordinate of interpolation points
%
%   OUTPUT:
%           c: vector of Newton's interpolation points (reverse order)
%

n = length(x);
c = y;
for i = 2:n
    c(i:n) = (c(i:n) - c(i-1))./(x(i:n) - x(i-1));
end
c = fliplr(c);

end