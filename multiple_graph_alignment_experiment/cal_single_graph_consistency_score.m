function [singleGraphConstList, singleGraphScoreList] = cal_single_graph_consistency_score(X,nodeCnt,graphCnt)
global affinity
% ��һ��ͼk�ĽǶ����㣬 cstAdj��һ��graphCnt*1�ľ���ÿ��Ԫ�ش���һ��graph��consistency�̶�
singleGraphConstList = zeros(graphCnt,1);
singleGraphScoreList = zeros(graphCnt,1);
for ref = 1:graphCnt 
    viewk = 1; 
    err = zeros((graphCnt+1)*(graphCnt-2)/2,1);
    scr = zeros((graphCnt+1)*(graphCnt-2)/2,1);
    rscope = (ref-1)*nodeCnt+1:ref*nodeCnt;
    for i = 1:graphCnt
        iscope = (i-1)*nodeCnt+1:i*nodeCnt;
        for j = i+1:graphCnt
            jscope = (j-1)*nodeCnt+1:j*nodeCnt;
            Xirj=X(iscope,rscope)*X(rscope,jscope);
            err(viewk) = sum(sum(abs(Xirj-X(iscope,jscope))))/2/nodeCnt;
            p = mat2vec(Xirj);
            gt = mat2vec(affinity.GT(iscope,jscope));
            scr(viewk) = (p'*affinity.K{i,j}*p)/(gt'*affinity.K{i,j}*gt);
            viewk = viewk + 1;
        end
    end 
    singleGraphConstList(ref) = 1-mean(err);
    singleGraphScoreList(ref) = mean(scr);
end