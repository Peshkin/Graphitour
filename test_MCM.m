function test_MCM(filename);
% test_MCM(filename);    test Maximal Cardinality Matching
% 
% INPUT:  a name of file containing one line per edge(I,J) in the format 
%         "I J  1", e.g. "5 14   1" means there is and edge (5,14)=(14,5)
%         [we assume no isolated nodes and undirected edges]
% OUTPUT:  a graph represented by ADJ-acency matrix with edge IDs as elements.
%
% by Dr. Leonid Peshkin MIT AI Lab Dec 2003 
% http://www.ai.mit.edu/~pesha
%
if exist(filename, 'file')
    edges = load(filename);
    if (max(size(edges)) < 2)
        warning('Degenerate or erroneous graph in file! \n');    
    end
else
    warning('File does not exist! \n');
    edges = 0; 
end
Nnds = max([edges(:,1)' edges(:,2)']);    % Num of nodes in the orig. graph
if (Nnds == 0)
    return
end
Nedgs = length(edges);                    % Num of edges in the orig. graph
adj = sparse(edges(:,1), edges(:,2), (Nnds+1):(Nedgs + Nnds));  % load adjacency matrix 
adj(Nnds,Nnds) = 0;                       % ADJacency matrix with edge IDs as elements
adj = adj - diag(diag(adj)) + adj';           % make it symmetric  

close all, ttl1 = sprintf('Original graph ''%s'' ', filename); 
[x, y, labels] = draw_dot(adj);  title(ttl1); 

[mate] = card_match(adj), adj 
matched = find(mate);
score = length(matched)/2;
adj_match = sparse(matched, mate(matched), ones(1,score*2));

ttl2 = sprintf('Matching of cardinality %d', score); figure 
graph_draw(adj_match, 'node_labels', labels, 'node_shapes', zeros(size(x,2),1), 'X', x, 'Y', y);
title(ttl2); axis square 