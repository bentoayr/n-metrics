function [clusters_1, clusters_2] = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges,threshold,distance_func)
clusters_1 = nan;
clusters_2 = nan;

if (1)
    good_hyper_edges = hyper_edges((hyper_edges(:,graphs_in_hyper_edge+1)< threshold),:);
    good_hyper_edges_inx = good_hyper_edges(:,1:graphs_in_hyper_edge);
    nn = num_rand_graphs; %size(unique(good_hyper_edges(:,1:3)),1); %number of vertices (graphs)
    mm = size(good_hyper_edges_inx,1); % number of edges (hyper-edges)
    msize = graphs_in_hyper_edge*ones(1,size(good_hyper_edges_inx,1));
    mlist = good_hyper_edges_inx;
    ng = 2;	%maximum number of clusters
    NR = 5;	%number of sampling initial conditions
    clusters_1 = hgc(nn,mm,msize,mlist,ng,NR);
end

if strcmp(distance_func,'gAlign')
        
    good_hyper_edges = hyper_edges((hyper_edges(:,graphs_in_hyper_edge+2)< threshold),:);
    good_hyper_edges_inx = good_hyper_edges(:,1:graphs_in_hyper_edge);
    nn = num_rand_graphs; %size(unique(good_hyper_edges(:,1:3)),1); %number of vertices (graphs)
    mm = size(good_hyper_edges_inx,1); % number of edges (hyper-edges)
    msize = graphs_in_hyper_edge*ones(1,size(good_hyper_edges_inx,1));
    mlist = good_hyper_edges_inx;
    ng = 2;	%maximum number of clusters
    NR = 5;	%number of sampling initial conditions
    clusters_2 = hgc(nn,mm,msize,mlist,ng,NR);
    
       
end

