function Z = InvertUpperTriangular(R)

    if(sum(abs(diag(R)) < 1e-10)>0)               % Detect singular matrix
        disp('Matrix is singular to working precision.')
        return
    else
        
        n = size(R,1);
        Z = diag(1./diag(R));                     % Solve diagonal values
        for j = 2:n
            for i = j-1:-1:1
                Z(i,j) = -(R(i,i+1:end)*Z(i+1:end,j))/R(i,i);
                                                  % Solve entries on the 
                                                  % upper triangular part.
            end
        end
    end
end