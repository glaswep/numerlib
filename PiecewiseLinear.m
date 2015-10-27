function [c0, c1] = PiecewiseLinear(f)

%   INPUT: 
%           f: vector of the value of f(x) at interpolation points
%
%   OUTPUT:
%           c0: vector of intercepts of each p(x)
%           c1: vector of slopes of each p(x)

if(~iscolumn(f))
    c0 = 0;
    c1 = 0;
    disp('Error: PiecewiseLinear: input should be column vector.')
else
    h = 1/(length(f) - 1);
    f1 = f(1:end-1);
    f2 = f(2:end);
    c1 = (f2 - f1)./h;
    x1 = linspace(0,1-h,length(f) - 1)';
    c0 = (f1 - c1.*x1);
end

end