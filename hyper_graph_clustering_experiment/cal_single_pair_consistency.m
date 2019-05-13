function pairConst = cal_single_pair_consistency(X,Xij,i,j,nodeCnt,graphCnt,massOutlierMode,inlierMask)

if nargin<7
    massOutlierMode = 0;
end
is = (i-1)*nodeCnt+1:i*nodeCnt;
js = (j-1)*nodeCnt+1:j*nodeCnt;
errX = 0;
if massOutlierMode
    estInCnt = sum(inlierMask(:,1));

    b = mat2vec(repmat(inlierMask(:,i)',nodeCnt,1));
end
for k=1:graphCnt
    if k==i || k==j, continue;end
    ks = (k-1)*nodeCnt+1:k*nodeCnt;

    aggX = X(is,ks)*X(ks,js);
    if massOutlierMode
        errX = errX + sum(abs(Xij(:) - aggX(:)).*b);
    else
        errX = errX + sum(abs(Xij(:) - aggX(:)));
    end
end
if massOutlierMode
    pairConst = 1-errX/(graphCnt*2*estInCnt);
else
    pairConst = 1-errX/(graphCnt*2*nodeCnt);
end