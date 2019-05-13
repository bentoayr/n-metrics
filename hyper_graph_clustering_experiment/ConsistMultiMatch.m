
function X = ConsistMultiMatch(updtOrder,nodeCnt,graphCnt,algpar,rawMat)
global affinity
bPathSelect = algpar.bPathSelect;
iterMax = algpar.iccvIterMax;
refGrh = updtOrder(end);
P = cell(graphCnt,graphCnt,2);
F = cell(graphCnt,graphCnt);
J = zeros(graphCnt,2);
bestIdx = 1;
currIdx = 2;
rscope = (refGrh-1)*nodeCnt+1:refGrh*nodeCnt;

for x=[1:refGrh-1,refGrh+1:graphCnt]
    if nargin == 5 
        xscope = (x-1)*nodeCnt+1:x*nodeCnt;
        P{x,refGrh,bestIdx} = rawMat(xscope,rscope);
    else % otherwise compute the pairwise matching
        P{x,refGrh,bestIdx} = vec2mat(pairMatchSolver(x,refGrh,algpar),nodeCnt,nodeCnt);
    end
    P{x,refGrh,currIdx} = P{x,refGrh,bestIdx};
end

P{refGrh,refGrh} = eye(nodeCnt,nodeCnt);
I = sparse(P{refGrh,refGrh});

for iterk=1:iterMax
    for gk = 1:(graphCnt-1) 
        uptGrh = updtOrder(gk);
        totalK = affinity.K{uptGrh,refGrh};
        for fixGrh = 1:graphCnt
            if fixGrh==uptGrh||fixGrh==refGrh, continue;end
            F{fixGrh,refGrh} = kron(P{fixGrh,refGrh,bestIdx},I);
            totalK = totalK + F{fixGrh,refGrh}'*affinity.K{uptGrh,fixGrh}*F{fixGrh,refGrh};
        end
        totalK = totalK./(graphCnt-1);
        if strcmpi(algpar.algMethod,'GAGM'),totalK=sparse(totalK);end
        exeString = ['p=',algpar.algMethod,'(totalK,nodeCnt,nodeCnt,algpar);'];
        eval(exeString);
         if algpar.bDisc
            E12 = ones(affinity.nP{uptGrh},affinity.nP{refGrh});
            p = convert2Discrete(E12,p);
         end
         P{uptGrh,refGrh,currIdx} = vec2mat(p,nodeCnt,nodeCnt);
         if bPathSelect
             J(uptGrh,currIdx) = mat2vec(P{uptGrh,refGrh,currIdx})'*totalK*mat2vec(P{uptGrh,refGrh,currIdx});%本次迭代的优化函数分数
             J(uptGrh,bestIdx) = mat2vec(P{uptGrh,refGrh,bestIdx})'*totalK*mat2vec(P{uptGrh,refGrh,bestIdx});%本次迭代的优化函数分数
             if J(uptGrh,currIdx)>J(uptGrh,bestIdx)
                P{uptGrh,refGrh,bestIdx} = P{uptGrh,refGrh,currIdx};
            end
         else 
             P{uptGrh,refGrh,bestIdx} = P{uptGrh,refGrh,currIdx};
         end
    end
end
X = zeros(nodeCnt*graphCnt,nodeCnt*graphCnt);
for x = 1:graphCnt
    xscope = (x-1)*nodeCnt+1:x*nodeCnt;
    for y = x+1:graphCnt
        yscope = (y-1)*nodeCnt+1:y*nodeCnt;
        X(xscope,yscope) = P{x,refGrh,bestIdx}*P{y,refGrh,bestIdx}';
    end
end
X = X + X' + eye(nodeCnt*graphCnt,nodeCnt*graphCnt);
