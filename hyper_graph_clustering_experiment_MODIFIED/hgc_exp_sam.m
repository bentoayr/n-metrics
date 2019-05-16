%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

nodes = 7;
num_rand_graphs = 20;
graphs_in_hyper_edge = 3;
num_hyp_edges = 100;
thresholds = [0.05:0.025:0.5]; 
num_repts = 100;

workspacefilename = ['Workspace_hgc_A_pert_clustering_glighn_rrwm_mOpt_mSync_', datestr(datetime)];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

for i = 1:num_repts

    G = generate_G_pert(nodes,num_rand_graphs, 0.7 ,  0.1);
    Graphs{i} = G;

    hyper_edges_gAlign{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'gAlign');
    hyper_edges_rrwm{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'rrwm');
    hyper_edges_mOpt{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'mOpt');
    hyper_edges_mSync{i} = build_Hypergraph_edges(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'mSync');
    
    save(workspacefilename);

end

toc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%all_res ={};

for i = 1:50%num_repts

    for th_ix = 1:length(thresholds)
        
        thres_gAlign = min(hyper_edges_gAlign{i}(:,4)) + (max(hyper_edges_gAlign{i}(:,4)) - min(hyper_edges_gAlign{i}(:,4)))*thresholds(th_ix);
        thres_rrwm = min(hyper_edges_rrwm{i}(:,4)) + (max(hyper_edges_rrwm{i}(:,4)) - min(hyper_edges_rrwm{i}(:,4)))*thresholds(th_ix);
        thres_mOpt = min(hyper_edges_mOpt{i}(:,4)) + (max(hyper_edges_mOpt{i}(:,4)) - min(hyper_edges_mOpt{i}(:,4)))*thresholds(th_ix);
        thres_mSync = min(hyper_edges_mSync{i}(:,4)) + (max(hyper_edges_mSync{i}(:,4)) - min(hyper_edges_mSync{i}(:,4)))*thresholds(th_ix);

        [clusters_gAlign_1, clusters_gAlign_2, clusters_gAlign_3] = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_gAlign{i}, thres_gAlign,'N');
        clusters_rrwm = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_rrwm{i}, thres_rrwm,'N');
        clusters_mOpt = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_mOpt{i}, thres_mOpt,'N');
        clusters_mSync = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges_mSync{i}, thres_mSync,'N');

        all_res{1}{1}{i}{th_ix} = clusters_gAlign_1;
        %all_res{1}{2}{i}{th_ix} = clusters_gAlign_2;
        %all_res{1}{3}{i}{th_ix} = clusters_gAlign_3;

        all_res{2}{i}{th_ix} = clusters_rrwm;
        all_res{3}{i}{th_ix} = clusters_mOpt;
        all_res{4}{i}{th_ix} = clusters_mSync;
        
        %save(workspacefilename);

    end
end



%% generate plots

% computing errors vs alpha
all_errors = {}
for t = 2:4
all_err = nan(10,10);

alg_ix = t;
for j = 4:10
    for i = 1:10
        val = min([norm(all_res{alg_ix}{j}{i} - kron([1,0;0,1],ones(10,1)),'fro')^2,norm(all_res{alg_ix}{j}{i} - kron([0,1;1,0],ones(10,1)),'fro')^2]);
        all_err(i,j) =  val;
    end
end
all_errors{t} = all_err;
end


for tt=1:3
all_err = nan(10,10);
alg_ix = 1;
for j = 4:10
    for i = 1:1
        val = min([norm(all_res{alg_ix}{tt}{j}{i} - kron([1,0;0,1],ones(10,1)),'fro')^2,norm(all_res{alg_ix}{tt}{j}{i} - kron([0,1;1,0],ones(10,1)),'fro')^2]);
        all_err(i,j) =  val;
    end
end
all_errors{1}{tt} = all_err;
end




hold on
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{1}{1}(:,4:10)',2),std(all_errors{1}{1}(:,4:10)',1,2),'r')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{1}{2}(:,4:10)',2),std(all_errors{1}{2}(:,4:10)',1,2),'b')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{1}{3}(:,4:10)',2),std(all_errors{1}{3}(:,4:10)',1,2),'g')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{2}(:,4:10)',2),std(all_errors{2}(:,4:10)',1,2),'m')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{3}(:,4:10)',2),std(all_errors{3}(:,4:10)',1,2),'c')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{4}(:,4:10)',2),std(all_errors{4}(:,4:10)',1,2),'k')