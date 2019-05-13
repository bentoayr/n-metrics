
function affinity = generateRandomAffinity_for_direct_comparison(list_of_graphs, nodeCnt,graphCnt)

nInlier = nodeCnt;
testk = 1;

global target
affinity.BiDir = target.config.affinityBiDir;
affinity.edgeAffinityWeight = target.config.edgeAffinityWeight;
affinity.angleAffinityWeight = target.config.angleAffinityWeight;

bGraphMatch = target.config.bGraphMatch;
bUnaryEnable = target.config.bUnaryEnable;
bEdgeEnable = target.config.bEdgeEnable;
graphCnt = target.config.graphCnt;
Sacle_2D = target.config.Sacle_2D;
deform = target.config.deform;
density = target.config.density;
nOutlier = target.config.nOutlier;
complete = target.config.complete;
affinity.graphCnt = graphCnt;
affinity.nodeCnt = nodeCnt;
affinity.EG = cell(graphCnt,1);
target.pairwiseMask{1} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);

basePoint = rand(nInlier,1);
if bGraphMatch
    baseEdge = tril(rand(nInlier),-1); baseEdge = baseEdge+baseEdge';
else
    basePointSet = randn(2, nInlier);
end
if complete<1
    target.pairwiseMask{testk} = generateRandomCompleteMask(nodeCnt,graphCnt,target.config.complete);
else 
    target.pairwiseMask{testk} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);
end


for gc=1:length(list_of_graphs) 
    affinity.nP{gc} = nodeCnt;
    if bGraphMatch
        randomEdge = tril(rand(nodeCnt),-1);
        N = deform*tril(randn(nodeCnt),-1);N=N+N';
        affinity.edge{gc} = randomEdge+randomEdge';
        affinity.edge{gc}(1:nInlier,1:nInlier) = baseEdge;

        affinity.edge{gc} = affinity.edge{gc} + N;    
        affinity.edgeRaw{gc} = affinity.edge{gc};     
        
        affinity.edge{gc}  = affinity.edge{gc} *0;
        affinity.edgeRaw{gc}  = affinity.edgeRaw{gc} *0;
        
        P = tril(rand(nodeCnt),-1); P = P+P';
        
        A_pert = list_of_graphs{gc};
                
        affinity.adj{gc} = logical(A_pert);
        affinity.edge{gc}(~A_pert) = NaN;
       
        
        randomPoint = rand(1,nodeCnt);
        PN = deform*randn(1,nodeCnt);
        affinity.pointFeat{gc} = randomPoint;
        affinity.pointFeat{gc}(1:nInlier) = basePoint;
        affinity.pointFeat{gc} = affinity.pointFeat{gc} + PN;
        
    else
        pointSet = [basePointSet + deform*randn(2,nInlier) randn(2,nOutlier)];
        pointSet = pointSet';
       [X,Y] = meshgrid(1:nodeCnt,1:nodeCnt);
       LL = [X(:),Y(:)];   
        G = pointSet(LL(:,1),:)-pointSet(LL(:,2),:);
        G = sqrt(G(:,1).^2+G(:,2).^2);
        G = reshape(G, [nodeCnt nodeCnt]);
        affinity.edge{gc} = G;
        affinity.edgeRaw{gc} = G; 
        affinity.pointFeat{gc} = rand(1,nodeCnt);
        affinity.pointFeat{gc}(1:nInlier) = basePoint;
        affinity.pointFeat{gc} = affinity.pointFeat{gc} + deform*randn(1,nodeCnt);
    end
    
    
    
    affinity.nE{gc} = sum(sum(~isnan(affinity.edge{gc})));
    [r,c] = find(~isnan(affinity.edge{gc}));
    affinity.EG{gc}=[r,c]';
    affinity.edgeFeat{gc} = affinity.edge{gc}(~isnan(affinity.edge{gc}))';
    
    
    affinity.pointFeat{gc} = affinity.pointFeat{gc}*0;
    affinity.edgeFeat{gc} = affinity.edgeFeat{gc}*0;
    
    affinity.G{gc} = zeros(nodeCnt,affinity.nE{gc});
    for c = 1 : affinity.nE{gc}
        affinity.G{gc}(affinity.EG{gc}(:, c), c) = 1;
    end
    affinity.H{gc} = [affinity.G{gc}, eye(nodeCnt)];
end
for xview = 1:graphCnt
    if affinity.BiDir
        yviewSet = [1:xview-1,xview+1:graphCnt];
    else
        yviewSet = xview+1:graphCnt;
    end
    for yview = yviewSet
        if bUnaryEnable
            affinity.KP{xview,yview} = exp(-conDst(affinity.pointFeat{xview}, affinity.pointFeat{yview},0) / Sacle_2D);
        else
            affinity.KP{xview,yview} = zeros(nodeCnt,nodeCnt);
        end
        if bEdgeEnable
            affinity.KQ{xview,yview} = exp(-conDst(affinity.edgeFeat{xview}, affinity.edgeFeat{yview},0) / Sacle_2D);%���ߺͱ�֮����affinity Kq����ָ����ʽ�Ƚ�ƽ�������Ż�,�γ�ָ������ı߱߾���n1*n2
        else
            affinity.KQ{xview,yview} = zeros(nodeCnt^2,nodeCnt^2);
        end
        affinity.K{xview,yview} = conKnlGphKU(affinity.KP{xview,yview}, affinity.KQ{xview,yview}, affinity.EG{xview},affinity.EG{yview});%EG�ǱߵĶ˵��������2*n1,2*n
        affinity.K{xview,yview} = full(affinity.K{xview,yview});
    end
end