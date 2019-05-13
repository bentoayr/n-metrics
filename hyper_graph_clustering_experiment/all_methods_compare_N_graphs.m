function [X, g_align_unrounded_score] = all_methods_compare_N_graphs(list_of_graphs, num_nodes,num_graphs)

g_align_unrounded_score = nan;
nInlier = num_nodes;

global affinity target

varyMinGrhCnt=num_graphs;varyMaxGrhCnt=num_graphs;grhTestCnt = 1;
testCnt = 1;

setPlotColor;
setObsoleteVariables;
target.config.testType = 'formal';

algpar = setPairwiseSolver();

target.config.database = 'synthetic';
target.config.Sacle_2D = 0.05;
iterRange = 6;
graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
if strcmp(target.config.testType,'massOutlier')
    algSet.algNameSet = {'mpm','rrwm','cao_s','cao_','cao_c_s','cao_c_','cao_uc_s','cao_uc_','cao_pc_s','cao_pc_','mOpt','mSync'};
    algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1,1];
    target.config.bGraphMatch = 0;
    target.config.category = 'outlier';
    target.config.inCntType = 'exact';
    nInlier = 6;target.config.nOutlier = 12;target.config.deform = .05;
    target.config.density = 1;target.config.complete = 1;
    graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
else
    algSet.algNameSet = {'rrwm','cao_','cao','cao_c_','cao_c','cao_uc_','cao_uc','cao_pc_','cao_pc','mOpt','mSync'};
    algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1];
    
    target.config.bGraphMatch = 1;
    target.config.inCntType = 'all';
    target.config.category = 'deform';
    switch target.config.category
        case 'deform'
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.15;
            target.config.density = .9;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    end
end

graphRange = graphMinCnt:1:graphMaxCnt;
target.config.initConstWeight = .2;
target.config.constStep = 1.1;
target.config.constWeightMax = 1;
target.config.constIterImmune = 2;
target.config.edgeAffinityWeight = 1;
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:graphMaxCnt;
paraCnt=length(graphRange);

[~,rrwmIdx] = ismember('rrwm',algSet.algNameSet);
[~,cao_Idx] = ismember('cao_',algSet.algNameSet);
[~,caoIdx] = ismember('cao',algSet.algNameSet);

algCnt = length(algSet.algEnable);
X=cell(algCnt,1);
target.config.nodeCnt = length(target.config.selectNodeMask);
target.config.graphCnt = min(max(graphRange),length(target.config.selectGraphMask{1}));
nodeCnt = num_nodes;
graphCnt = num_graphs;

% now start to test

for testk = 1:testCnt
    
    affinity = generateRandomAffinity_for_direct_comparison(list_of_graphs, num_nodes, num_graphs);
    affinity.GT = repmat(eye(nodeCnt,nodeCnt),graphCnt,graphCnt);%just use identity matrix as grund truth matchings
    
    
    rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,testk);
    
    switch target.config.inCntType
        case 'exact'
            target.config.inCnt = nodeCnt - target.config.nOutlier;
        case 'all'
            target.config.inCnt = nodeCnt;
        case 'spec'
            target.config.inCnt = specNodeCnt;
    end
    scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    for parak=1:paraCnt
        viewCnt=graphRange(parak);
        offsetStep = graphRange(end);
        affinity.viewCnt = viewCnt;
        for grhOffset = 1:offsetStep:graphCnt-viewCnt+1
            rawscopeX = (grhOffset-1)*nodeCnt+1:(grhOffset+viewCnt-1)*nodeCnt;
            rawScopeY = (grhOffset-1)*nodeCnt+1:(grhOffset+viewCnt-1)*nodeCnt;
            affinity.Xgt = affinity.GT(rawscopeX,rawScopeY);
            baseMat = rawMat(rawscopeX,rawScopeY);
            
            X{rrwmIdx} = baseMat;
            scrDenom = max(max(scrDenomMatInCnt(1:viewCnt,1:viewCnt)));
            
            singleGraphConstList = cal_single_graph_consistency_score(baseMat,nodeCnt,viewCnt);
            [~,refConstGraph] = max(singleGraphConstList);
            
            cstGrhList = rankGrhByConsistencyToRefGrh(baseMat,refConstGraph,nodeCnt,viewCnt);
            updGrhList = [cstGrhList,refConstGraph];
            
            % nipsMethod
            [~,nipsIdx] = ismember('mSync',algSet.algNameSet);
            if nipsIdx>0&&algSet.algEnable(nipsIdx)
                X{nipsIdx} = SynchronizePermute(baseMat,nodeCnt,viewCnt,'sync');
            end
            [~,iccvIdx] = ismember('mOpt',algSet.algNameSet);
            if iccvIdx>0&&algSet.algEnable(iccvIdx)
                algpar.bPathSelect = 1;
                X{iccvIdx} = ConsistMultiMatch(updGrhList,nodeCnt,viewCnt,algpar,baseMat);
            end
            
            if cao_Idx>0&&algSet.algEnable(cao_Idx)
                X{cao_Idx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'afnty',1);
                [X{cao_Idx}, ~, g_align_unrounded_score] = g_align_distance(affinity, nodeCnt, viewCnt, 1, 1);
                X{caoIdx} = X{cao_Idx};
                if caoIdx&&algSet.algEnable(caoIdx)
                    X{caoIdx} = SynchronizePermute(X{cao_Idx},nodeCnt,viewCnt,'sync');
                end
            end
            
        end
        
    end
end

end

