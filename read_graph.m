function [edges] = read_graph(filename);
% [adj, edges, Nnds, Nedgs] = read_graph(filename) 
% INPUT:  a name of file containing one line per edge(I,J) in the format 
%         "I J  1", e.g. "5 14   1" means there is and edge (5,14)=(14,5)
%         [we assume no isolated nodes and undirected edges]
% OUTPUT:  a graph represented by ADJ-acency matrix with edge IDs as elements.
%
% by Dr. Leonid Peshkin MIT AI Lab Dec 2003 
% http://www.ai.mit.edu/~pesha
%
%      1 2 3                                                   4-edge(1,2)
%  1   0 4 5   sample  ADJ matrix for the following graph: (2)---(1)---(3)
%  2   4 0 0                                                         5-edge(1,3)
%  3   5 0 0
%   we assume no isolated vertices and un-directed graph (symmetric ADJ) 

if exist(filename, 'file')
    edges = load(filename);
    if (max(size(edges)) < 2)
        warning('Degenerate or erroneous graph in file! \n');    
    end
else
    warning('File ''%s'' does not exist! \n', filename);
    edges = 0; 
end