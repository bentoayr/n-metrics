clc
clear
tic
nodes = 7;
num_rand_graphs = 20;
graphs_in_hyper_edge = 3;
num_hyp_edges = 100;
threshold = [0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3];
num_repts = 10;



for alphaix = 5:5
    alpha = 0.02 + (alphaix - 1 )*0.05;
    for i = 1:num_repts
        
        G = generate_G(nodes,num_rand_graphs,alpha,1 - alpha);
        Graphs{alphaix}{i} = G;
        hyper_edges = Hypergraph_clustering(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,'gAlign');
    end
end
toc

res ={};

for alphaix = 5:5
    alpha = 0.02 + (alphaix - 1 )*0.05;
    for i = 1:num_repts
        
        for th = 1:8%2:8
            res{alphaix}{i}{th} = Hypergraph_clustering_th(graphs_in_hyper_edge,num_rand_graphs,hyper_edges, threshold(th)); %clusters_gAlign_3
            
            %[clusters_gAlign_1, clusters_gAlign_2, clusters_gAlign_3] = Hypergraph_clustering(G,nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,threshold,'gAlign');
            
            %all_res{1}{1}{alphaix}{i} = clusters_gAlign_1;
            %all_res{1}{2}{alphaix}{i} = clusters_gAlign_2;
            %all_res{1}{3}{alphaix}{i} = clusters_gAlign_3;
            
            
            %         clusters_rrwm = Hypergraph_clustering(Graphs{alphaix}{i},nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,threshold,'rrwm');
            %         all_res{2}{alphaix}{i} = clusters_rrwm;
            %         clusters_mOpt = Hypergraph_clustering(Graphs{alphaix}{i},nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,threshold,'mOpt');
            %         all_res{3}{alphaix}{i} = clusters_mOpt;
            %         clusters_mSync = Hypergraph_clustering(Graphs{alphaix}{i},nodes,graphs_in_hyper_edge,num_hyp_edges,num_rand_graphs,threshold,'mSync');
            %         all_res{4}{alphaix}{i} = clusters_mSync;
        end
    end
end

save(['Workspace_hgc_10_reps_1_alphas_glighn_rrwm_mOpt_mSync_', datestr(datetime)]);
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


for tt=3:3
all_err = nan(10,10);
alg_ix = 1;
for j = 4:4
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

plot([0.02+0.05*6],all_errors{1}{3}(1,4:4),'*y','MarkerSize',10,'LineWidth',2)

errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{2}(:,4:10)',2),std(all_errors{2}(:,4:10)',1,2),'m')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{3}(:,4:10)',2),std(all_errors{3}(:,4:10)',1,2),'c')
errorbar([0.02+0.05*3:0.05:0.48],mean(all_errors{4}(:,4:10)',2),std(all_errors{4}(:,4:10)',1,2),'k')
%%
b = 10;% number of repetitions
c = 8;% number of thresholds
all_errors = {};
err = nan(b,c);%nan(1,8);
for j = 5:5 %alpha
    for i = 1:10 %rep
        for th=1:8
        val = min([norm(res{j}{i}{th} - kron([1,0;0,1],ones(10,1)),'fro')^2,norm(res{j}{i}{th} - kron([0,1;1,0],ones(10,1)),'fro')^2]);
        err(i,th) =  val;
        all_errors{j} = err;
        end
    end  
end