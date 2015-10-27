function c = lscoef(k, alpha, beta, f, w, a, b)

%   LSCOEF find the coef of orthogonal polynomials used to solve continuous
%   LS problem.
%
%   c = lscoef(k, alpha, beta, f, w, a, b) returns c as coefficient of k+1
%       polynomials.
%
%   See also:
%       getcoef.m, evalphik.m
%
%   Yulun Zeng, Oct, 2015

c = zeros(k+1,1);

phikm = @(x) 1;
phik = @(x) x - integral(@(x)x.*1.*w(x),a,b)./integral(@(x)1.*1.*w(x),a,b);

c(1) = integral(@(x) f(x).*phikm(x).*w(x), a, b)...
    ./integral(@(x) phikm(x).*phikm(x).*w(x), a, b);

c(2) = integral(@(x) f(x).*phik(x).*w(x), a, b)...
    ./integral(@(x) phik(x).*phik(x).*w(x), a, b);

for i = 3:k+1
    
    phikp = @(x) (x-alpha(i-1)).*phik(x) - beta(i-1).*phikm(x);
    phikm = phik;
    phik = phikp;
    
    c(i) = integral(@(x) f(x).*phik(x).*w(x), a, b)...
        ./integral(@(x) phik(x).*phik(x).*w(x), a, b);
end


end