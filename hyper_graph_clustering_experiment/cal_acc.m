function acc = cal_acc(P,nOutlier,GT)
%P: Nx*Ny, GT: Nx*Ny
% nOutlier = 0;
[Nx Ny]=size(P);
if nargin<3
    GT = zeros(Nx,Ny);
    GT(:,1:Nx) = diag(ones(1,Nx));
end
if Ny==Nx
    P = P(1:end-nOutlier,1:end-nOutlier);
    Nx = size(P,1);
    GT = GT(1:end-nOutlier,1:end-nOutlier);
end

acc = sum(sum(abs(P-GT),2)==0)/Nx;