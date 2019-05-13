function [hyper_edges] = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,distance_func)

hyper_edges = zeros(num_hyp_edges,graphs_in_hyper_edge+1); %graph indecies and the distance

for t=1:num_hyp_edges
    
    clear global;
    clearvars -except t hyper_edges G nodes graphs_in_hyper_edge num_hyp_edges num_rand_graphs distance_func;
    
    Result = 1;
    while Result
        g = randperm(num_rand_graphs,graphs_in_hyper_edge);
        [Result,~] = ismember(g,unique(hyper_edges(:,1:graphs_in_hyper_edge),'rows'),'rows');
    end
    
    hyper_edges(t,1:graphs_in_hyper_edge) = g;
    
    for ii = 1:graphs_in_hyper_edge
        all_graphs{ii} = G{g(ii)};
    end
    
    g_align_score = nan;
    
    [X,g_align_score]  = all_methods_compare_N_graphs(all_graphs, nodes,graphs_in_hyper_edge);
    
    if strcmp(distance_func,'gAlign')
        Aligns = X{2};
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
        hyper_edges(t,5) = g_align_score; %g_align unrounded score
    end
    
end

fclose all;
end

