function y = EvalPiecewiseLinear(c0, c1, x)

%   INPUT:
%           c0: vector of intercepts of each p(x)
%           c1: vector of slopes of each p(x)
%           x: vector of values to be evaluated
%
%   OUTPUT:
%           y: vector of values of p(x) at x


if( iscolumn(c0) && iscolumn(c1))
    
    h = 1/length(c0);
    ind = floor(x./h) + 1;
    ind(ind==length(c0) + 1) = length(c0);
    X = sparse(1:length(x), ind, x, length(x), length(c1));
    c0s = c0(ind);
    y = X*c1 + c0s;

else
    disp('Error: EvalPiecewiseLinear: c0 and c1 should be column vector.')
    y = 0;
end

end