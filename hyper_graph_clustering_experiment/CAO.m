
function P = CAO(rawMat,nodeCnt,graphCnt,iterMax,scrDenom,optType,useCstInlier)
    if iterMax==0,P = rawMat;return;end
    global  matchMat 
    global target 
    inCnt = target.config.inCnt;
    inlierMask =  zeros(nodeCnt,graphCnt); 
    if strcmp(target.config.testType,'massOutlier'),massOutlierMode = 1;else massOutlierMode = 0;end
    constIterImmune = target.config.constIterImmune;
    curIdx = 1;lastIdx = 2;
    constStep = target.config.constStep;
    initConstWeight = target.config.initConstWeight;
    constWeightMax = target.config.constWeightMax;
    if massOutlierMode 
        if useCstInlier
            inlierMask = cal_node_consistency_mask(rawMat,nodeCnt,graphCnt,inCnt);
        else
            inlierMask = cal_node_affinity_mask(rawMat,nodeCnt,graphCnt,inCnt);
        end
    end

    iter = 0;
    err = ones(iterMax,1);
    I = eye(nodeCnt*graphCnt,nodeCnt*graphCnt);
    matchMat = rawMat;
    constWeight = initConstWeight;
    unaryListConsistency{curIdx} = cal_single_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode,inlierMask);
    unaryListConsistency{lastIdx} = unaryListConsistency{curIdx};
    pairListConsistency{curIdx} = cal_pair_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode,inlierMask);
    pairListConsistency{lastIdx} = pairListConsistency{curIdx};

    [X,Y] = meshgrid(1:graphCnt,1:graphCnt);
    X = X(:);Y = Y(:);
    while iter<iterMax 
        tempMat = zeros(nodeCnt*graphCnt,nodeCnt*graphCnt);
        for vk = 1:graphCnt^2
            xview = X(vk); yview = Y(vk);
            if mod(vk-1,graphCnt) == 0,xscope = (xview-1)*nodeCnt+1:xview*nodeCnt;end
            if X(vk)>= Y(vk),continue;end
            yscope = (yview-1)*nodeCnt+1:yview*nodeCnt;	
            if iter<constIterImmune&&(strcmp(optType,'exact')||strcmp(optType,'pair')||strcmp(optType,'unary'))% for CAO-C, CAO-PC, CAO-UC, in early immune stage, only affinity is boosted
                tempMat(xscope,yscope) = find1stOrderPathByScoreUnaryPairConsistency(...
                    xview,yview,unaryListConsistency{curIdx},pairListConsistency{curIdx},nodeCnt,graphCnt,constWeight,scrDenom,'afnty',massOutlierMode,inlierMask);
                tempMat(xscope,yscope) = find1stOrderPathByScoreUnaryPairConsistency(...
                    xview,yview,unaryListConsistency{curIdx},pairListConsistency{curIdx},nodeCnt,graphCnt,constWeight,scrDenom,optType,massOutlierMode,inlierMask);
            end
        end
        matchMatOld = matchMat;
        matchMat = tempMat+tempMat'+I;
        if iter>=constIterImmune
            constWeight = min([constWeightMax,constStep*constWeight]);
        else
            if massOutlierMode
                if useCstInlier
                    inlierMask = cal_node_consistency_mask(matchMat,nodeCnt,graphCnt,inCnt);
                else
                    inlierMask = cal_node_affinity_mask(matchMat,nodeCnt,graphCnt,inCnt);
                end
            end
        end
        unaryListConsistency{lastIdx} = unaryListConsistency{curIdx};
        unaryListConsistency{curIdx} = cal_single_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode,inlierMask);%0.001秒的运行时间
        pairListConsistency{lastIdx} = pairListConsistency{curIdx};
        pairListConsistency{curIdx} = cal_pair_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode,inlierMask);

        iter = iter + 1;
        err(iter) = sum(abs(matchMatOld(:) - matchMat(:)));
        if err(iter)==0
            break;
        end
    end

    P = matchMat;

function P = find1stOrderPathByScoreUnaryPairConsistency(xview,yview,unaryListConsistency,pairListConsistency,...
    nodeCnt,graphCnt,constWeight,scrDenom,metricType,massOutlierMode,inlierMask)
    global affinity 
    global matchMat
    pairCon = zeros(graphCnt,1);
    xscope = (xview-1)*nodeCnt+1:xview*nodeCnt;
    yscope = (yview-1)*nodeCnt+1:yview*nodeCnt;
    Y = zeros(nodeCnt*nodeCnt,graphCnt);
    for anchor=1:graphCnt
        ascope = (anchor-1)*nodeCnt+1:anchor*nodeCnt;
        P1 = matchMat(xscope,ascope);P2 = matchMat(ascope,yscope);
        Y(:,anchor) = mat2vec(P1*P2);
        if strcmp(metricType,'pair')
            pairCon(anchor) = sqrt(pairListConsistency(xview,anchor)*pairListConsistency(anchor,yview));
        end
    end
    [~,m1,n1] = unique(Y','rows','first');
    uniLen = length(m1);
    uniCon = zeros(uniLen,1);uniAfty = zeros(uniLen,1);
    for i=1:uniLen
        a = m1(i);
        p = sparse(Y(:,a));
        if massOutlierMode,b=repmat(inlierMask(:,a)',nodeCnt,1);p = p.*b(:);end
        if ~strcmp(metricType,'cstcy')
            uniAfty(i) = p'*affinity.K{xview,yview}*p/scrDenom;
        end

        if strcmp(metricType,'exact')||strcmp(metricType,'cstcy')
			uniCon(i) = cal_single_pair_consistency(matchMat,vec2mat(Y(:,a),nodeCnt,nodeCnt),xview,yview,nodeCnt,graphCnt,massOutlierMode,inlierMask);
		end
    end
    afntyScr = uniAfty(n1);
    if strcmp(metricType,'exact')||strcmp(metricType,'cstcy')
        exaCon = uniCon(n1);
    end
    
    switch metricType
        case 'unary'% CAO-UC
            fitness = (1-constWeight)*afntyScr + constWeight*unaryListConsistency;
        case 'pair'% CAO-PC
            fitness = (1-constWeight)*afntyScr + constWeight*pairCon;
        case 'exact'% CAO-C
            fitness = (1-constWeight)*afntyScr + constWeight*exaCon;
        case 'afnty'% CAO
            fitness = afntyScr;
        case 'cstcy'
            fitness = exaCon;
    end
 
	[~, idx] = max(fitness);
	P = vec2mat(Y(:,idx(1)),nodeCnt,nodeCnt);
