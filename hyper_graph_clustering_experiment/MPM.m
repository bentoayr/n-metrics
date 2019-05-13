function [ X ] = MPM( M, nP1, nP2, algpar)

param = struct( ...
    'amp_max', 30, ...              
    'thresConvergence', 1e-15 ...  
);

E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);
strField = fieldnames(param);
for i = 1:length(strField), eval([strField{i} '=param.' strField{i} ';']); end


nMatch = length(M);
prev_score = ones(nMatch,1)/nMatch; 
prev_score2 = prev_score;        

bCont = 1;  iter_i = 0;



cur_score = zeros(nMatch,1);

while bCont && iter_i < algpar.iterMax1
    
    iter_i = iter_i + 1;
    cur_score = RMP_mult(full(M), prev_score, group1);
    sumCurScore = sqrt(sum(cur_score.^2)); 
    if sumCurScore>0, cur_score = cur_score./sumCurScore; end
    
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
 
end 

X = cur_score;