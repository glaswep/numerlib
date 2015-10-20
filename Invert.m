function Y = Invert(A)
    [Q,R] = QRHouseholder(A);
    Y = InvertUpperTriangular(R)*Q';
end