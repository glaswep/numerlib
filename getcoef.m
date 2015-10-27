function [alpha, beta] = getcoef(N, a, b, w)

%   GETCOEF finds the coefficient alpha and beta used to generate 
%   orthogonal polynomials by 3 term recurrence relationship:
%
%   phi_{k+1}(x) = (x-alpha_{k+1})phi_k(x) - beta_{k+1}phi_{k-1}(x) 
%
%       [alpha, beta] = getcoef(N, a, b, w) returns alpha and beta to
%       generate N orthogonal polynomials on interval [a,b] weighted by
%       function w.
% 
%   Examples:
%       w = @(x) sin(x);
%       [alpha, beta] = getcoef(10, -1, 1, w);
%
%   See also:
%       evalphik.m, lscoef.m
%
%   Yulun Zeng, Oct, 2015

phikm = @(x) 1;                         
phik = @(x) x - integral(@(x)x.*1.*w(x),a,b)./integral(@(x)1.*1.*w(x),a,b);
alpha = zeros(N, 1);                    
beta = zeros(N, 1);

for i = 2:N
    
    alpha(i) = integral(@(x) x.*phik(x).*phik(x).*w(x), a, b)...
        /integral(@(x) phik(x).*phik(x).*w(x), a, b);
    
    beta(i) = integral(@(x) x.*phik(x).*phikm(x).*w(x), a, b)...
        /integral(@(x) phikm(x).*phikm(x).*w(x), a, b);
    
    phikp = @(x) (x-alpha(i)).*phik(x) - beta(i).*phikm(x);
    phikm = phik;
    phik = phikp;
    
end

end