function [Q, R] = QRHouseholder(A)

% 	QRHOUSEHOLDER uses Householder reflectors to do QR decomposition on
%	m by n matrix A
%
%	[Q, R] = QRHouseholder(A) where A is m by n matrix
%
%	See also:
%		InvertUpperTriangular.m, Invert.m
%
%	Yulun Zeng, Oct, 2015

[m,n] = size(A);                            % Find size of A
Q = eye(m);                                 % Intialize Q

for k=1:min(m-1,n)
   ak  = A(k:m,k);
   vk  = ak; vk(1) = vk(1) + exp(1i*angle(ak(1)))*norm(ak);
                                            % Find v for householder
   vk  = vk/norm(vk);                       % Normalize
    %    c = (ak'*vk)/(vk'*ak);             % Alternative c value
   c = 1;
   A(k:m,k:n) = A(k:m,k:n) - (1+c)*vk*(vk'*A(k:m,k:n));
                                            % Apply Householder on subset
                                            % of matrix A
   Q(:,k:end) = Q(:,k:end) - (1+c)*Q(:,k:end)*(vk*vk');
                                            % Apply Householder on subset
                                            % of matrix Q                                           
end
R = A;                                      % Now A is just R

end