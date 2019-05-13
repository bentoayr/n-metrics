function [ Tree,Cost ] =  UndirectedMaximumSpanningTree (CostMatrix)


n = size (CostMatrix,1); 
EdgeWeights = 0;         
EdgeWeightsCounter = 0;
for i = 1:n
    for j = (i+1):n
        if ((CostMatrix(i,j))~=0)
            EdgeWeightsCounter = EdgeWeightsCounter + 1;
            EdgeWeights(EdgeWeightsCounter,1) = CostMatrix(i,j);
            EdgeWeights(EdgeWeightsCounter,2) = i;
            EdgeWeights(EdgeWeightsCounter,3) = j;
        end
    end
end

SortedEdgeWeights = 0;
SortedEdgeWeights = sortrows(EdgeWeights);

m = size(SortedEdgeWeights,1); 

global ParentPointer ;
ParentPointer = 0;
ParentPointer(1:n) = 1:n;


TreeRank = 0;
TreeRank(1:n) = 0;

MSTreeEdges = 0;
MSTreeEdgesCounter = 0; i = m;
while ((MSTreeEdgesCounter < (n-1)) && (i>=1))

    root1=0; root2=0; temproot=0;
    temproot = SortedEdgeWeights(i,2);
    root1 = FIND_PathCompression(temproot);
  
    temproot = SortedEdgeWeights(i,3);
    root2 = FIND_PathCompression(temproot);
    
    if (root1 ~= root2)
        MSTreeEdgesCounter = MSTreeEdgesCounter + 1;
        MSTreeEdges(MSTreeEdgesCounter,1:3) = SortedEdgeWeights(i,:);
        if (TreeRank(root1)>TreeRank(root2))
            ParentPointer(root2)=root1;
        else
            if (TreeRank(root1)==TreeRank(root2))
               TreeRank(root2)=TreeRank(root2) + 1;
            end
            ParentPointer(root1)=root2;
        end
    end
    i = i - 1;
end

MSTreeEdgesCounter = 0;
Tree = 0;
Tree(1:n,1:n)=0;
while (MSTreeEdgesCounter < (n-1))
    MSTreeEdgesCounter = MSTreeEdgesCounter + 1;
    Tree(MSTreeEdges(MSTreeEdgesCounter,2),MSTreeEdges(MSTreeEdgesCounter,3))=1;
    Tree(MSTreeEdges(MSTreeEdgesCounter,3),MSTreeEdges(MSTreeEdgesCounter,2))=1;
end


Cost = 0;
for p = 1:n
    for q = p+1:n
       if Tree( p,q ) == 1 
          Cost = Cost + CostMatrix( p,q );
       end
    end
end

end

function [parent] = FIND_PathCompression(temproot)

global ParentPointer;
ParentPointer(temproot);
if (ParentPointer(temproot)~=temproot)
    ParentPointer(temproot) = FIND_PathCompression(ParentPointer(temproot));
end
parent = ParentPointer(temproot);
end