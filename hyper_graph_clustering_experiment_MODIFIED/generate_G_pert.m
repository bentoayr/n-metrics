function G = generate_G_pert(nodes,num_rand_graphs, alpha, beta)
    G_A1 = cell(num_rand_graphs/2,1);
    G_A2 = cell(num_rand_graphs/2,1);

    A_1 = rand(nodes) < alpha;
    A_1 = triu(A_1,1) +triu(A_1,1)';

    A_2 = rand(nodes) < alpha;
    A_2 = triu(A_2,1) +triu(A_2,1)';

    for i = 1:num_rand_graphs/2
        A_pert = A_1;
        frac_pert = beta;
        edges_to_flip = randi(nodes,ceil(frac_pert*alpha*nodes*(nodes-1)*0.5),2);
        for j = 1:length(edges_to_flip)
            A_pert(edges_to_flip(j,1),edges_to_flip(j,2)) = 1 - A_pert(edges_to_flip(j,1),edges_to_flip(j,2));
        end
        G_A1{i} = A_pert;
    end

    for i = 1:num_rand_graphs/2
        A_pert = A_2;
        frac_pert = beta;
        edges_to_flip = randi(nodes,ceil(frac_pert*alpha*nodes*(nodes-1)*0.5),2);
        for j = 1:length(edges_to_flip)
            A_pert(edges_to_flip(j,1),edges_to_flip(j,2)) = 1 - A_pert(edges_to_flip(j,1),edges_to_flip(j,2));
        end
        G_A2{i} = A_pert;
    end

    G = [G_A1 G_A2];

end
