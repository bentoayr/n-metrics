function [G] = randomAdj(n,p)
G = rand(n) < p;
G = triu(G,1) + triu(G,1)';
end

