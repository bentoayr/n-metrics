function Xd=discretisationMatching_hungarian(X,E12,varargin);
X=-X;
X=X-min(X(:));
X(E12==0)=Inf;
[Xd,score]=hungarian(X);
