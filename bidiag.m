function [B, U, V] = bidiag(A)

%   BIDIAG: find bidiagonal matrix B and unitary U and V such that
%
%                       B = U*A*V'
%
%   INPUT: 
%           A: complex or real m by n matrix
%
%   OUTPUT:
%           B: complex or real m by n bidiagonal matrix
%           U: complex or real m by m matrix
%           V: complex or real n by n matrix
%
%   Yulun Zeng
%   Oct, 18th 2015


[m,n] = size(A);                        % Initialize
U = eye(m);
V = eye(n);

for k=1:min(m-1,n)
   x  = A(k:m,k);                       % Build Householder reflector
                                        % for A
   uk  = x; 
   uk(1) = x(1) + exp(1i*angle(x(1)))*norm(x);
   uk  = uk/norm(uk);
   cu = 1;
   A(k:m,k:n) = A(k:m,k:n) - (1+cu)*uk*(uk'*A(k:m,k:n));
                                        % Apply Householder to A
   U(k:end,:) = U(k:end,:) - (1+cu)*uk*(uk'*U(k:end,:));
                                        % Apply Householder to U
   if k < n-1 
       x = A(k,k+1:n)';                 % Build Householder reflector 
                                        % for A^T
       vk = x;
       vk(1) = x(1) + exp(1i*angle(x(1)))*norm(x);
       vk = vk/norm(vk);
       cv = 1;
       A(k:m,k+1:n) = A(k:m,k+1:n) - (1+cv)*(A(k:m,k+1:n)*vk)*vk';
                                        % Apply Householder to A
       V(k+1:end,:) = V(k+1:end,:) - (1+cu)*vk*(vk'*V(k+1:end,:));
                                        % Apply Householder to U
   end
end
B = A;

end