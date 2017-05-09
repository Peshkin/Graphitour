# Graphitour
Satelite MATLAB implementation of lossless graph compression algorithm Graphitour as described in https://arxiv.org/abs/cs/0703132
Please cite:  "Structure induction by lossless graph compression". L. Peshkin. DCC, page 53-62. IEEE Computer Society, (2007)

This work is motivated by the necessity to automate the discovery of structure in vast and evergrowing collection of relational data commonly represented as graphs, for example genomic networks. A novel algorithm, dubbed Graphitour, for structure induction by lossless graph compression is presented and illustrated by a clear and broadly known case of nested structure in a DNA molecule. This work extends to graphs some well established approaches to grammatical inference previously applied only to strings. The bottom-up graph compression problem is related to the maximum cardinality (non-bipartite) maximum cardinality matching problem. The algorithm accepts a variety of graph types including directed graphs and graphs with labeled nodes and arcs. The resulting structure could be used for representation and classification of graphs. 

Note - this code is old, created in 2003/2004 since MATLAB is not backward compatible it might need some debugging. line cell2num function got replaced by cell2mat but probably more fixes needed. 

It depends on "neato" graph layout library http://www.graphviz.org/category/graphviz-terms/neato

Can be run like this 
    [graph, lex, adj, edges] = parseG('graphfile', 'graphs/dna')
     when 'graphfile'  is 'X', edge list is assumed to be in 'X.dat' and node labels list in 'X.index'
 thre are sample files provided : dna.index and dna.dat as used in the paper. 
 
 There are two structures GRAPH and LEXICON, the first contains all elements of the graph and refers to TYPE in a LEXICON for descripion.
  GRAPH contains all nodes followed by all edges first.
 
The parse.m routine relies on card_match.m which is my implementation of maximal cardinality matching Gabow's N^3 algorithm.
Depends also on my_intersect, my_setdiff - optimized versions of MATLAB routines. 
