function [ new_X ] = fix_X_matrix(X, perm_mat_list)

    new_X = X;
    
    graph_size = size(perm_mat_list{1},1);

    num_graphs = size(X,1) / graph_size;
    
    for i = 1:num_graphs
        for j = 1:num_graphs
            new_X( (1:graph_size)   +   (i-1)*graph_size, (1:graph_size)   +   (j-1)*graph_size ) = inv(perm_mat_list{i})*X( (1:graph_size)   +   (i-1)*graph_size, (1:graph_size)   +   (j-1)*graph_size )*(perm_mat_list{j});
        end
    end

end