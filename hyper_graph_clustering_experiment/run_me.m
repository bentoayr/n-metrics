%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code was used to produce the results included in Section 8.2 of
% ``Tractable $n$-metrics for multiple graphs'', published in ICML 2019

% The code takes about 1h to run on an intel i7 processor for one repetition, i.e. num_repts = 1.
% To fully replicate the experiment in the paper, num_repts = 50, the code should run for about two days

clc

clear

save_workspace_flag = 0; % set to 1 to save workspace after each repetition
workspacefilename = ['Workspace_', datestr(datetime)];

generate_plots_flag = 0; % set to 1 to generate plots

nodes = 7; %number of nodes in each graph
num_rand_graphs = 20; % total number of graphs in two clusters
graphs_in_hyper_edge = 3; %number of nodes included in a hyper edge
num_hyp_edges = 100; % number of hyper edges
thresholds = 0.05:0.025:0.5; %thresholds to decide wheather to create a hyper edge based on the distance among correspinding graphs
num_repts = 2;


%% Generate hyperedges values

for i = 1:num_repts
    
    G = generate_G_pert(nodes,num_rand_graphs, 0.7 ,  0.1);
    Graphs{i} = G;
    
    hyper_edges_gAlign{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'gAlign');
    hyper_edges_rrwm{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'rrwm');
    hyper_edges_mOpt{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'mOpt');
    hyper_edges_mSync{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'mSync');
    
    if save_workspace_flag == 1
        save(workspacefilename);
    end
    
end


%% Compute clusters

all_res ={};

for i = 1:num_repts
    
    for th_ix = 1:length(thresholds)
        
        thres_gAlign = min(hyper_edges_gAlign{i}(:,4)) + (max(hyper_edges_gAlign{i}(:,4)) - min(hyper_edges_gAlign{i}(:,4)))*thresholds(th_ix);
        thres_rrwm = min(hyper_edges_rrwm{i}(:,4)) + (max(hyper_edges_rrwm{i}(:,4)) - min(hyper_edges_rrwm{i}(:,4)))*thresholds(th_ix);
        thres_mOpt = min(hyper_edges_mOpt{i}(:,4)) + (max(hyper_edges_mOpt{i}(:,4)) - min(hyper_edges_mOpt{i}(:,4)))*thresholds(th_ix);
        thres_mSync = min(hyper_edges_mSync{i}(:,4)) + (max(hyper_edges_mSync{i}(:,4)) - min(hyper_edges_mSync{i}(:,4)))*thresholds(th_ix);
        
        [clusters_gAlign_1, clusters_gAlign_2] = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_gAlign{i}, thres_gAlign,'N');
        clusters_rrwm = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_rrwm{i}, thres_rrwm,'N');
        clusters_mOpt = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_mOpt{i}, thres_mOpt,'N');
        clusters_mSync = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_mSync{i}, thres_mSync,'N');
        
        all_res{1}{1}{i}{th_ix} = clusters_gAlign_1;
        all_res{1}{2}{i}{th_ix} = clusters_gAlign_2;
        
        all_res{2}{i}{th_ix} = clusters_rrwm;
        all_res{3}{i}{th_ix} = clusters_mOpt;
        all_res{4}{i}{th_ix} = clusters_mSync;
        
        if save_workspace_flag == 1
            save(workspacefilename);
        end
        
    end
end



%% Compute errors vs alpha

all_errors = {};

% First calulate errors for G-align
for tt=1:2 %loop over Ours and Ours*
    all_err = nan(num_repts,length(thresholds));
    alg_ix = 1;
    for j = 1:num_repts
        for i = 1:length(thresholds)
            val = min([norm(all_res{alg_ix}{tt}{j}{i} - kron([1,0;0,1],ones(10,1)),'fro')^2,norm(all_res{alg_ix}{tt}{j}{i} - kron([0,1;1,0],ones(10,1)),'fro')^2]);
            all_err(j,i) =  val;
        end
    end
    all_errors{1}{tt} = all_err;
end

% Now calulate errors for other algorithms
for t = 2:4 %loop of Pairwise, mOpt, matchSync
    all_err =  nan(num_repts,length(thresholds));
    
    alg_ix = t;
    for j = 1:num_repts
        for i = 1:length(thresholds)
            val = min([norm(all_res{alg_ix}{j}{i} - kron([1,0;0,1],ones(10,1)),'fro')^2,norm(all_res{alg_ix}{j}{i} - kron([0,1;1,0],ones(10,1)),'fro')^2]);
            all_err(j,i) =  val;
        end
    end
    all_errors{t} = all_err;
end

%% Generate plots

if generate_plots_flag == 1
    hold on
    
    % Error vs Threshold
    errorbar(thresholds(1:length(thresholds)),mean(all_errors{1}{1}(:,1:length(thresholds))',2),std(all_errors{1}{1}(:,1:length(thresholds))',1,2),'r')
    errorbar(thresholds(1:length(thresholds)),mean(all_errors{1}{2}(:,1:length(thresholds))',2),std(all_errors{1}{2}(:,1:length(thresholds))',1,2),'b')
    
    errorbar(thresholds(1:length(thresholds)),mean(all_errors{2}(:,1:length(thresholds))',2),std(all_errors{2}(:,1:length(thresholds))',1,2),'m')
    errorbar(thresholds(1:length(thresholds)),mean(all_errors{3}(:,1:length(thresholds))',2),std(all_errors{3}(:,1:length(thresholds))',1,2),'c')
    errorbar(thresholds(1:length(thresholds)),mean(all_errors{4}(:,1:length(thresholds))',2),std(all_errors{4}(:,1:length(thresholds))',1,2),'k')
    
end
