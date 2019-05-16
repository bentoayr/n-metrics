function [ave_best_val, std_val, min_val, max_val, best_thr_ix, all_vals] = find_best_thres_for_clust_50_50_split(clust_for_alg, num_reps)

    num_thres = length(clust_for_alg{1});

    ave_best_val = inf;
    best_thr_ix = nan;
    
    for i = 1:num_thres
       
        tmp = [];
        for j = 1:num_reps
            val = compute_clustering_score_for_two_clust_50_50_split(clust_for_alg{j}{i});
            tmp = [tmp, val];
        end
        av = mean(tmp);
        if (av < ave_best_val)
            ave_best_val = av;
            best_thr_ix = i;
            std_val = std(tmp)/sqrt(num_reps);
            min_val = min(tmp);
            max_val = max(tmp);
            all_vals = tmp;
        end
        
    end
    
end