function [x,t] = LMM(x0, t0, t1, f, alpha, beta, h)
   
%   LMM: give numerical solutions to first order explicit differential 
%   equations using linear multistep method.
%
%   (This scheme uses RK4.p to find initial values for multistep methods.)   
%
%   INPUT:      x0: initial value (scalar)
%               t0: initial time
%                f: function handle for dx/dt
%            alpha: coef's on the LHS of linear multistep method.
%             beta: coef's on the RHS of linear multistep method.
%                h: time step
%
%   OUTPUT:      x: solution vector at each time step
%                t: vector of time steps
%
%   Yulun Zeng, Nov, 2015
%

    if nargin ~= 7
        fprintf('LMM error: invalid number of arguments.\n')
        x = -1; t = -1;
        return
    end
    
    n = length(alpha) - 1;
    
    if length(beta) ~= n
        fprintf('LMM error: invalid alpha or beta.\n')
        x = -1; t = -1;
        return
    end
    if iscolumn(alpha) == 0, alpha = alpha';end
    if iscolumn(beta) == 0, beta = beta';end
    if h<0, fprintf('LMM error: negative time step.\n')
        return 
    end
    if nargin(f) ~= 2, fprintf('LMM error: invalid f.\n')
        return
    end
    
    steps = floor((t1 - t0)/h);
    
    x = zeros(steps + 1,1);
    fval = x;
    t = (h*(0:steps) + t0)';
    if steps + 1 < n
        x(1:end) = RK4(x0, t0, f, h, steps);
        fprintf('LMM warning: time step too large, result obtained by RK4.\n')
    else
        x(1:n) = RK4(x0, t0, f, h, n-1);
%         x(1:n) = [1, 1 + eps];
        fval(1:n) = f(x(1:n), t(1:n));

        for i = (n+1):steps+1
            rhs = (-alpha(1:(end-1)))'*x((i-n):(i-1)) + ...
                  h*(beta(1:end)'*fval((i-n):(i-1)));
            x(i) = rhs/alpha(end);
            fval(i) = f(x(i), t(i));
        end
    end

end