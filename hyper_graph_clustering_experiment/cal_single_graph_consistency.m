function cstAdj = cal_single_graph_consistency(X,nodeCnt,graphCnt,massOutlierMode,inlierMask)
if nargin<4
    massOutlierMode = 0;
end

if massOutlierMode
    estInCnt = sum(inlierMask(:,1));
    b = zeros(nodeCnt^2,graphCnt);
    for i=1:graphCnt
        b(:,i) = mat2vec(repmat(inlierMask(:,i)',nodeCnt,1));
    end
end
cstAdj = zeros(graphCnt,1);
for ref = 1:graphCnt 
    viewk = 1; err = zeros((graphCnt+1)*(graphCnt-2)/2,1);
    rscope = (ref-1)*nodeCnt+1:ref*nodeCnt;
    for i = 1:graphCnt
        iscope = (i-1)*nodeCnt+1:i*nodeCnt;
        for j = i+1:graphCnt
            jscope = (j-1)*nodeCnt+1:j*nodeCnt;
            Xirj=X(iscope,rscope)*X(rscope,jscope);           
             if massOutlierMode
                 err(viewk) = sum(abs(Xirj(:) - mat2vec(X(iscope,jscope))).*b(:,i));
             else
                 err(viewk) = sum(sum(abs(Xirj-X(iscope,jscope)),2));
             end
            viewk = viewk + 1;
        end
    end 
    if massOutlierMode
        cstAdj(ref) = 1-mean(err)/(2*estInCnt);
    else
        cstAdj(ref) = 1-mean(err)/(2*nodeCnt);
    end
end