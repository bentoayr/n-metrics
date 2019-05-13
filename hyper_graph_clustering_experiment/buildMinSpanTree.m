function treeAdj = buildMinSpanTree(scoreAdj)
[treeAdj, ~] =  UndirectedMaximumSpanningTree(scoreAdj);
