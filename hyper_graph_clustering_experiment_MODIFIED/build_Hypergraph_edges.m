function [hyper_edges] = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,distance_func)

%G = generate_G(nodes,num_rand_graphs);
%Example: Hypergraph_clustering(G,7,3,11,20,3,'gAlign')

hyper_edges = zeros(num_hyp_edges,graphs_in_hyper_edge+1); %graph indecies and the distance

for t=1:num_hyp_edges
    
    clear global;
    clearvars -except t hyper_edges G nodes graphs_in_hyper_edge num_hyp_edges num_rand_graphs distance_func;
    
    Result = 1;
    while Result
        g = randperm(num_rand_graphs,graphs_in_hyper_edge); % randi([1 20],1,3) with no repeated values
        [Result,~] = ismember(g,unique(hyper_edges(:,1:graphs_in_hyper_edge),'rows'),'rows'); %LocResult
    end
    
    hyper_edges(t,1:graphs_in_hyper_edge) = g;
    
    for ii = 1:graphs_in_hyper_edge
        all_graphs{ii} = G{g(ii)};
    end
    
    g_align_score = nan;
    
    [X,g_align_score]  = all_other_methods_compare_N_graphs(all_graphs, nodes,graphs_in_hyper_edge);
    
    if strcmp(distance_func,'gAlign')
        Aligns = X{2};
        Aligns_1 = X{3};
    else if strcmp(distance_func,'rrwm')
            Aligns = X{1};
        else if strcmp(distance_func,'mOpt')
                Aligns = X{10};
            else if strcmp(distance_func,'mSync')
                    Aligns = X{11};
                end
            end
        end
    end
    
    hyper_edges(t,4) = compute_distance_from_alignment(Aligns, all_graphs,graphs_in_hyper_edge, nodes);
    
    if strcmp(distance_func,'gAlign')
        hyper_edges(t,5) = compute_distance_from_alignment(Aligns_1, all_graphs,graphs_in_hyper_edge, nodes);
        hyper_edges(t,6) = g_align_score;
    end
    
end



end

