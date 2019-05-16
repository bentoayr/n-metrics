n = 10;
k = 4;
p = 0.5;
frac_pert = 0.1;

A_ref = rand(n) < p;
A_ref = triu(A_ref,1) +triu(A_ref,1)';

for i = 1:k
    affinity.adj{i} = A_ref;
    edges_to_flip = randi(n,round(frac_pert*p*n*(n-1)*0.5),2);
    for j = 1:length(edges_to_flip)
        affinity.adj{i}(edges_to_flip(j,1),edges_to_flip(j,2)) = 1 - affinity.adj{i}(edges_to_flip(j,1),edges_to_flip(j,2));
    end
end

[P, score] = g_align_distance(affinity, n, k, 0, 1);