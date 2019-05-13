function D = conDst(X1, X2, bAngle)

n1 = size(X1, 2);
n2 = size(X2, 2);
if size(X1, 1) == 1
    X1 = [X1; zeros(1, n1)]; 
    X2 = [X2; zeros(1, n2)]; 
end
XX1 = sum(X1 .* X1); XX2 = sum(X2 .* X2); 

X12 = X1' * X2;
D = repmat(XX1', [1, n2]) + repmat(XX2, [n1, 1]) - 2 * X12;
if bAngle
    idx = find (D>1);
    D(idx) = (sqrt(D(idx))-2).^2;
end