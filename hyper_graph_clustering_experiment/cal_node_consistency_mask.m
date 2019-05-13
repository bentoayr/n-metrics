function nodeCstPerGraph = cal_node_consistency_mask(X,nodeCnt,graphCnt,inCnt)
nodeCstPerGraph = zeros(nodeCnt,graphCnt);
for ref = 1:graphCnt
    rscope = (ref-1)*nodeCnt+1:ref*nodeCnt;
    for i = 1:graphCnt-1
        iscope = (i-1)*nodeCnt+1:i*nodeCnt;
        for j = i+1:graphCnt
            jscope = (j-1)*nodeCnt+1:j*nodeCnt;
            Xrij=X(rscope,iscope)*X(iscope,jscope);
            nodeCstPerGraph(:,ref) = nodeCstPerGraph(:,ref) + (1-sum(abs(Xrij-X(rscope,jscope)),2)/2);
        end
    end
end
nodeCstPerGraph = nodeCstPerGraph/(graphCnt*(graphCnt-1)/2);
[~,IX] = sort(nodeCstPerGraph,1,'descend');
nodeCstPerGraph2 = zeros(nodeCnt,graphCnt);
for ref = 1:graphCnt
    nodeCstPerGraph2(IX(1:inCnt,ref),ref) = 1;
end
nodeCstPerGraph = nodeCstPerGraph2;