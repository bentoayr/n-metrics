function [ X ] = RRWM( M, nP1, nP2, algpar)



param = struct( ...
    'c', 0.2, ...                   
    'amp_max', 30, ...              
    'thresConvergence', 1e-25, ...  
    'tolC', 1e-3 ...                
);
E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);

strField = fieldnames(param);
for i = 1:length(strField), eval([strField{i} '=param.' strField{i} ';']); end
[idx1 ID1] = make_groups(group1);
[idx2 ID2] = make_groups(group2);
% 
if ID1(end) < ID2(end)
    [idx1 ID1 idx2 ID2 dumVal dumSize] = make_groups_slack(idx1, ID1, idx2, ID2);
    dumDim = 1;
elseif ID1(end) > ID2(end)
    [idx2 ID2 idx1 ID1 dumVal dumSize] = make_groups_slack(idx2, ID2, idx1, ID1);
    dumDim = 2;
else
    dumDim = 0; dumVal = 0; dumSize = 0;
end
idx1 = idx1-1; idx2 = idx2-1;


d = sum(M, 1); 
maxD = max(d);
Mo = M ./ maxD; 

nMatch = length(M);
prev_score = ones(nMatch,1)/nMatch; 
prev_score2 = prev_score;         
prev_assign = ones(nMatch,1)/nMatch; 

bCont = 1;  iter_i = 0;

while bCont && iter_i < algpar.iterMax1
    
    iter_i = iter_i + 1;
    
    cur_score = Mo * ( c*prev_score + (1-c)*prev_assign );
    
    sumCurScore = sum(cur_score); 
    if sumCurScore>0, cur_score = cur_score./sumCurScore; end
    
    cur_assign = cur_score;
    amp_value = amp_max/ max(cur_assign);  
    cur_assign = exp( amp_value*cur_assign );  
    
    X_slack = [cur_assign; dumVal*ones(dumSize,1)];
    X_slack = mexBistocNormalize_match_slack(X_slack, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC, dumDim, dumVal, int32(1000));
    cur_assign = X_slack(1:nMatch);   
    
    sumCurAssign = sum(cur_assign); 
    if sumCurAssign>0, cur_assign = cur_assign./sumCurAssign; end
    
    if 1
        diff1 = sum((cur_score-prev_score).^2);
        diff2 = sum((cur_score-prev_score2).^2); 
        diff_min = min(diff1, diff2);
        if diff_min < thresConvergence
            bCont = 0;
        end
    else
        normed_cur_score = cur_score/norm(cur_score);
        if norm(M*normed_cur_score - la*normed_cur_score,1) < thresConvergence2
            bCont = 0;
        end
        la = normed_cur_score'*M*normed_cur_score;
    end

    prev_score2 = prev_score;
    prev_score = cur_score;
    prev_assign = cur_assign;
 
end 

X = cur_score;
end