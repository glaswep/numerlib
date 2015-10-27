function c = newtonup(c,x,xn,yn)

%   INPUT:
%           c: vector of newton coefficients
%           x: vector of interpolation points
%           xn: x-coordinate of the new point
%           yn: y-coordinate of the new point
%
%   OUTPUT:
%           c: vector of new interpolation point

if(~(any(x)==xn))
    
    c = fliplr(c);

    c = [c yn];
    x = [x xn];
    n = length(x);

    for i = 2:n
        c(n) = (c(n) - c(i-1))./(x(n) - x(i-1));
    end

    c = fliplr(c);
    
else
    disp('Error: newtonup: new interpolation point alreay used.')
end

end