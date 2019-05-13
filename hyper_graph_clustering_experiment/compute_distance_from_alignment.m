function value = compute_distance_from_alignment(Aligns, Graphs,numgraphs, numnodes)
 
    value = 0;
    for i = 1:numgraphs
        for j = 1:numgraphs
         
            Pij = Aligns((i-1)*numnodes+ [1:numnodes],(j-1)*numnodes+ [1:numnodes]);
            Ai = Graphs{i};
            Aj = Graphs{j};
             
            value = value + norm(Ai*Pij - Pij*Aj,'fro');
             
        end
    end
 
    value = value / (numgraphs*numgraphs);
     
end