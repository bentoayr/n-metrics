function G = generate_G(nodes,num_rand_graphs, alpha, beta)
G_sparse = [];
G_dense = [];

for i = 1:num_rand_graphs/2
    G_sparse{i} = randomAdj(nodes,alpha);
end

for j = 1:num_rand_graphs/2
    G_dense{j} = randomAdj(nodes,beta);
end

G = [G_sparse G_dense];

end
