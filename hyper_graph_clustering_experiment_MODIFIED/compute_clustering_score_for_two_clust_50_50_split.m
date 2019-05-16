function val = compute_clustering_score_for_two_clust_50_50_split(clust)

            val = min([norm(clust - kron([1,0;0,1],ones(10,1)),1),norm( clust - kron([0,1;1,0],ones(10,1)),1)]) / 20;



end