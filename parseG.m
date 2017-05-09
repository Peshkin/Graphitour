function [graph, lex, adj, edges] = parse(varargin)
%
% Graph compression / grammar induction 
% July 17th' 2003  last change   Feb 08, 2004
% Example: 
% [graph, lex, adj, edges] = parse('graphfile','graphs/XXX');  'integrated.dat'
% 
%   when 'graphfile'  is 'X', edge list is assumed to be in 'X.dat' and 
%                             node labels list in 'X.index'
% There are two structures GRAPH and LEXICON, the first contains all
% elements of the graph and refers to TYPE in a LEXICON for descripion.
% GRAPH contains all nodes followed by all edges first.
%
% This routine relies on card_match.m which is my implementation
% of maximal cardinality matching Gabow's N^3 algorithm.
% also on my_intersect, my_setdiff - optimized versions of MATLAB routines
% as a result of streamlining.
% also I introduced more flags to account for edge type between two nodes,
% and ignore cardinality of nodes as node type 
global flag lex 

flag.nodenames = 0;            % take nodenames into account when assigning node types
flag.ignorecardinality = 0;    % ignore cardinality defining initial node types  
flag.ignoredgetype = 0;        % ignore what edge connects two nodes ? 
flag.debug = 0;                             % debug = 2 the most detailed printouts 
flag.break_ties = 1;                        % flag indicating whether to break ties 
flag.size_bias  = -1;           % do we merge first large[1], small[-1] or any[0] elements ?
flag.assume_directed = 1;            % do we assume that input file contans DIRECTED graph ? 
filename = 'datag/stories_graph'; 
for i = 1:2:nargin                           %  INITIALIZATION 
    switch varargin{i}
        case 'graphfile', filename = varargin{i+1}; 
        case 'ratio', ratio = varargin{i+1}; 
        otherwise, error(['unrecognized argument' varargin{i}])
    end
end
[edges, labels] = read_graph(filename); 
[lex, graph, adj, MCM, adj0] = init_structs(edges, labels, flag); % INITIAL PROCESSING of NODES and EDGES  
[MCM] = update_MCM(adj, graph, MCM, 1); 
best_score = max(MCM.score); 
while (best_score > 1) % the most popular edgetype has more than one instance  
    rplc_type = find(MCM.score == best_score);
    rplc_type = pick_rplc_type(rplc_type);  % in case there are several
    [matched, mate, rplc_edges] = find(MCM.mate{rplc_type});
    rplc_edges = my_unique(rplc_edges); 
    if flag.debug >= 1, 
        fprintf('\n\t\t\tReplacing %d elements of type %d, size %d \t\t\t vvvvvvv\n', ...
        MCM.score(rplc_type), rplc_type, lex.weight(rplc_type)); end 
                % CYCLE THROUGH EDGEs_2B_REPLACED 
    while ~isempty(rplc_edges)         % -- have to do it an ugly way instead of FOR loop 
        edge = rplc_edges(1);          % -- since RPLC_EDGES changes within the loop 
               % a pair of vertises making replaced edge  
        [pair(1) pair(2)] = find(triu(adj) == edge);      % TriU because ADJ symmetric !!!!
        [adj, lex, MCM, rplc_edges] = remove_edge(pair, graph, adj, lex, MCM, rplc_edges);
        [graph] = create_elem(edge, pair, graph);                           %  
        if flag.debug == 2, fprintf('Merging elements %d & %d into %d. Family ties:',pair(1),pair(2),edge); end 
              % MAKE a LIST of NODES LINKED TO the (replaced/collapsed) edge = FAMILY 
        family = my_unique([find(adj(pair(2),:)) find(adj(pair(1),:))]);
        family = my_setdiff(family, pair);
        for node = family             % LOOP THROUGH the LIST of nodes adjacent to the edge 
            [adj, lex, MCM, rplc_edges, dead_type(1)] = remove_edge([node pair(1)], graph, adj, lex, MCM, rplc_edges);
            [adj, lex, MCM, rplc_edges, dead_type(2)] = remove_edge([node pair(2)], graph, adj, lex, MCM, rplc_edges);
            [graph, lex, MCM, adj] = create_edge(node, edge, pair, dead_type, graph, lex, adj, MCM, rplc_type, flag);     
        end 
        if flag.debug == 2, fprintf('... \n'); end 
    end
    if flag.debug >= 1, fprintf('\t\t\tDone replacing    elements \t\t\t\t^^^^^^^ \n'); 
        input('Hit <any> key to continue ....'); end
    [MCM] = update_MCM(adj, graph, MCM, 1); 
    best_score = max(MCM.score); 
end  
if flag.debug >= 1, fprintf('    GRAPH COMPRESSION is FINISHED \n'); end
display_results(adj0, adj, graph, lex, labels);

             % ------------------------------  PICK  REPLACE-TYPE  --------------   NOT USED 

function [rplc_type] = pick_rplc_type(rplc_type)    % lex, flag 
global lex flag 
              % PICK THE  (hyper)EDGE TYPE to COLLAPSE (HEURISTICALY) 
              % WE MIGHT PREFER SIZE to MCM Score ALLTOGETHER !!!!
count_ties = 0; count_tiesN = 0; 
LnNdxS = length(rplc_type);
if (LnNdxS > 1)  % there are several identically good matchings
    if (flag.size_bias == 1)          % prefer merging larger components first
        count_ties = count_ties + 1, 
        best_weight = max(lex.weight(rplc_type));  % only among top scored MCMs
        ndx_w = find(lex.weight(rplc_type) == best_weight);
        rplc_type = rplc_type(ndx_w);                   % find out which ones were these
        LnNdxW = length(rplc_type);  % how many of highest weights among highest scores ?
        if (LnNdxW > 1) & flag.break_ties 
            [trash, which] = max(rand(1, LnNdxW));  % uniformly random among best_scored ones
            rplc_type = rplc_type(which);
        else
            rplc_type = rplc_type(1); % always the 1st for reproductable results and debugging
        end
    elseif (flag.size_bias == -1)     % prefer merging smallest components first 
        count_tiesN = count_tiesN + 1, 
        best_weight = min(lex.weight(rplc_type));  % only among top scored MCMs
        ndx_w = find(lex.weight(rplc_type) == best_weight);
        rplc_type = rplc_type(ndx_w);                   % find out which ones were these
        LnNdxW = length(rplc_type);  % how many of highest weights among highest scores ?
        if (LnNdxW > 1) & flag.break_ties 
            [trash, which] = max(rand(1, LnNdxW));  % uniformly random among best_scored ones
            rplc_type = rplc_type(which);
        else
            rplc_type = rplc_type(1); % always the 1st for reproductable results and debugging
        end
    elseif flag.break_ties   
        [trash, which] = max(rand(1, LnNdxS));  % uniformly random among best_scored ones 
        rplc_type = rplc_type(which);
    else
        rplc_type = rplc_type(1);   % always the 1st for reproductable results and debugging 
    end
end

             % ------------------------------  UPDATE  Max Card Matching  ---------------  NOT USED! 

function [MCM] = update_MCM(adj, graph, MCM, best_score)
global lex flag 
adj_type = adj;                % replacing edge IDs with their EdgeType IDs
adj_type(find(adj)) = graph.type(adj(find(adj)));   % workaround Matlab indexing
if flag.debug >= 1, fprintf('\n Max. Cardinality Matching: \n');  end
types = my_intersect(find(MCM.score == 0), find(lex.count(lex.NndTyps + 1 : lex.size) > 1) + lex.NndTyps);
for type = types %(marked for updating) and (at lest two edges of that type)
     sub_adj = (adj_type == type);          % sub graph on edges of a given type
     sub_Nedgs = full(sum(sum(sub_adj))/2);  % # of edges in sub-graph 
     if (sub_Nedgs >= best_score)        % potentially better MCM score (Lazy evaluation)    
         if flag.debug >= 1, fprintf('on %d edges(type %d) ', sub_Nedgs, type);  end
         mate = card_match(adj.*sub_adj); if flag.debug >= 1, fprintf(' is'); end 
         matched = find(mate);
         score = length(matched)/2;  if flag.debug >= 1, fprintf(' %d elems, ', score); end
         MCM.score(type) = score;    if flag.debug == 2, [matched; mate(matched)], end  
         MCM.mate{type} = sparse(matched,mate(matched),ones(1,score*2));  % save the result of expensive MCM calculation          
         N = size(adj); MCM.mate{type}(N,N) = 0;    % change the "official" matrix size
         MCM.mate{type} = MCM.mate{type}.*adj; 
     end
end

              % ------------------------------ READ  GRAPH  ----------------------------
              
function [edges, labels] = read_graph(filename);
% [adj, edges, Nnds, Nedgs] = read_graph(filename) 
% INPUT:  a name of file containing one line per edge(I,J) in the format 
%         "I J  1", e.g. "5 14   1" means there is and edge (5,14)=(14,5)
%         [we assume no isolated nodes and undirected edges]
% OUTPUT:  a graph represented by ADJ-acency matrix with edge IDs as elements.
%      1 2 3                                                   4-edge(1,2)
%  1   0 4 5   sample  ADJ matrix for the following graph: (2)---(1)---(3)
%  2   4 0 0                                                         5-edge(1,3)
%  3   5 0 0
%   we assume no isolated vertices and un-directed graph (symmetric ADJ) 
edgefile = strcat(filename,'.dat');
if exist(edgefile, 'file')
    edges = load(edgefile);
    if (max(size(edges)) < 2)
        warning('Degenerate or erroneous graph in file! \n');    
    end
else
    warning('File does not exist! \n');
    edges = 0; 
end
lblsfile = strcat(filename,'.index');
if exist(lblsfile, 'file')
    [trash, labels] = textread(lblsfile, '%d %s'); 
else
    labels=[];
end

               % ------------------------------ INIT STRUCTS ------------------------

function [lex, graph, adj, MCM, adj0] = init_structs(edges, labels, flag)
graph.Nedgs = length(edges);                    % Num of edges in the orig. graph
graph.Nnds = max([edges(:,1)' edges(:,2)']);    % Num of nodes in the orig. graph
adj = sparse(edges(:,1), edges(:,2), (graph.Nnds+1):(graph.Nedgs + graph.Nnds));  % load adjacency matrix 
adj(graph.Nnds, graph.Nnds) = 0;                       % ADJacency matrix with edge IDs as elements
adj0 = adj';
adj = adj - diag(diag(adj)) + adj';        % make it symmetric  
graph.size = graph.Nedgs + graph.Nnds; 
if flag.assume_directed == 0, adj0 = adj; end   %  FORCE UN-DIRECTED for now 

       % PROCESS NODES into GRAPH and LEXICON structure     
if isempty(labels), flag.nodenames = 0; end 
if ~flag.ignorecardinality      % node degrees are used insted of labels here 
    degrees = sum(adj > 0, 1);              % array of degrees 
    [node_types, trash, type_ndx] = unique(degrees); 
    if flag.nodenames   % need to account for BOTH nodenames and cardinality 
        [lbl_types,trash, lbl_ndx] = unique(labels); Nlbls = length(lbl_types);
        complex = (degrees(:)-1)*(Nlbls+1)+lbl_ndx; 
        [node_types, trash, type_ndx] = unique(complex);    
    end
elseif flag.nodenames 
    [lbl_types, trash, type_ndx] = unique(labels); 
    lex.NndTyps = length(lbl_types);      % number of unique node types    
    node_types = 1:lex.NndTyps; 
else
    node_types = 1;       %  all nodes are the same type initially, makes all edges same  
    type_ndx   = ones(1, graph.Nnds); 
end           %% ^^^ alternatively we could have loaded node types ^^^^^           
lex.NndTyps = length(node_types);      % number of unique node types    
lex.junior = node_types(:)';           % description of nodes is their degree
lex.senior = zeros(1, lex.NndTyps);
lex.weight = ones(1, lex.NndTyps);     % these elements have one member
for type = 1:lex.NndTyps            % loop through node types 
    lex.count(type) = sum(type_ndx == type); % number of nodes of this type
end
graph.type(1:graph.Nnds) = type_ndx(1:graph.Nnds);
graph.members(1:graph.Nnds) = deal(num2cell(1:graph.Nnds));
graph.depth(1:graph.Nnds) = zeros(1,graph.Nnds);
          % NOW  PROCESSING  INITIAL EDGES ... 
for i = 1:graph.Nedgs       
    type1 = type_ndx(edges(i,1));     % here we order types of nodes
    type2 = type_ndx(edges(i,2));     % connected by edge to count
    jun_type(i) = min(type1, type2);  % edges of the same type
    sen_type(i) = max(type2, type1);
    graph.members(i + graph.Nnds) = {edges(i,1:2)};
end
graph.depth(graph.Nnds+1:graph.size) = ones(1,graph.Nedgs);
[edge_types, killme, type_ndx] = unique([jun_type' sen_type'], 'rows');
lex.NedgTyps = size(edge_types,1);    % number of edge types  
for type = 1:lex.NedgTyps 
    lex.count(lex.NndTyps + type) = sum(type_ndx == type); 
end
lex.size = lex.NndTyps + lex.NedgTyps;
lex.junior(lex.NndTyps + 1 : lex.size) = edge_types(:,1);
lex.senior(lex.NndTyps + 1 : lex.size) = edge_types(:,2);
lex.weight(lex.NndTyps + 1 : lex.size) = 2*ones(1,lex.NedgTyps);  % two nodes in each
graph.type((graph.Nnds+1):graph.size) = lex.NndTyps + type_ndx(1:graph.Nedgs);

MCM.score = zeros(1, lex.size);           % Maximum Cardinality matching structures 
MCM.mate = {}; 
lex.junE_type = 0;
lex.senE_type = 0;
lex.countC = lex.count;    % Cumulative count 
return

               % ------------------------------ CREATE  EDGE  ------------------------

function [graph, lex, MCM, adj] = create_edge(node, new_node, pair, dead_type, graph, lex, adj, MCM, rplc_type, flag) 
n_type = graph.type(node);    % FIND the type OF the family node
jun_type = min(n_type, rplc_type);  % edges of the same type
sen_type = max(n_type, rplc_type);
                           %% !!! Find here what kind of link - 
                  % CHECK (jun, sen) EDGE_TYPE in a LEXICON
new_type = my_intersect(find(lex.senior == sen_type), find(lex.junior == jun_type)); 
if ~flag.ignoredgetype
  extra_ = my_intersect(find(lex.junE_type == min(dead_type)), find(lex.senE_type == max(dead_type))); 
  new_type = my_intersect(new_type, extra_);
end
if isempty(new_type)     % IF does NOT exist, create a new TYPE entry in a LEXICON
    lex.size = lex.size + 1;
    new_type = lex.size; 
    lex.junior(new_type) = jun_type;    
    lex.senior(new_type) = sen_type;
    lex.count(new_type) = 0;   
    lex.countC(new_type) = 0;
    lex.weight(new_type) = lex.weight(jun_type) + lex.weight(sen_type); 
    lex.junE_type(new_type) = min(dead_type);
    lex.senE_type(new_type) = max(dead_type);
end
lex.count(new_type) = lex.count(new_type) + 1; 
lex.countC(new_type) = lex.countC(new_type) + 1; 
MCM.score(new_type) = 0;   % force to re-calculate Max Card matching on edges of this type 
graph.size = graph.size + 1;           % ADD a new (hyper)EDGE to a GRAPH 
graph.type(graph.size) = new_type; 
adj(new_node, node) = graph.size; 
adj(node, new_node) = graph.size; 
if flag.debug == 2, fprintf('link %d-%d(type %d)  ', node, new_node, new_type); end    

function [graph] = create_elem(edge, pair, graph)
graph.members(edge) = {[cell2mat(graph.members(pair(1))) cell2mat(graph.members(pair(2)))]}; 
graph.children(edge) = {pair};
graph.depth(edge) = max(graph.depth(pair(1)), graph.depth(pair(2))) + 1; 

               % ------------------------------ REMOVE  EDGE  ------------------------

function [adj, lex, MCM, rplc_edges, edge_type] = remove_edge(pair, graph, adj, lex, MCM, rplc_edges)  
edge = adj(pair(2), pair(1));
edge_type = 0;
if edge ~= 0
    edge_type = graph.type(edge);
    lex.count(edge_type) = lex.count(edge_type) - 1;
    MCM.score(edge_type) = 0;    % force to re-calculate Max Card matching on edges of this type 
    MCM.mate{edge_type} = {};
    adj(pair(2), pair(1)) = 0; adj(pair(1), pair(2)) = 0;     % Wipe involved edges out of ADJ 
    rplc_edges = my_setdiff(rplc_edges, edge); % REMOVE involved EDGEs of a given TYPE from the LIST of EDGEs_2B_REPLACED   
end
    
              % ------------------------------  DISPLAY  RESULTS  ---------------
              
function display_results(adj0, adj, graph, lex, labels);
   % let's draw some results
types = find(lex.count((graph.Nnds+1):lex.size)) + graph.Nnds;  % which types are present in final graph
elements = find(ismember(graph.type, types));           % elemnt IDs of the final graph
elements2 = find(ismember(graph.type, [lex.junior(types) lex.senior(types)]));
[ndx1 ndx2] = find(adj);  ndx1 = my_unique(ndx1)'; ndx2 = my_unique(ndx2)';
fprintf('\n   Ultimate graph has %d distinct elements. \n', length(ndx1));
close all, x = [];
if graph.Nnds < 300, 
    [x, y, labels] = draw_dot(adj0 > 0);  title('Original graph'); end 
MYlabels = cellstr(strcat(num2str(graph.type(ndx2)'), '.',int2str(ndx2')));
draw_dot(adj(ndx1,ndx2) > 0, MYlabels);
title('Top-level resulting graph. Labels are non-terminal types. Click to examine.');
text_elms = findall(gca,'Type','text');  
for ndx = 1:length(text_elms)
    callbk = 'my_call(str2num(get(gcbo,''String'')))'; 
    set(text_elms(ndx), 'ButtonDownFcn', callbk);  % assume the node label is a number
end
disj = zeros(graph.Nnds); 
for i = ndx1
    elems = graph.members{i};
    ttl = sprintf('Structure %4d of type %3d, depth %3d', i, graph.type(i), graph.depth(i));
    fprintf(ttl); fprintf(': { '); 
    fprintf('%d ',sort(elems)); fprintf('} \n');
    if (length(elems) > 5)
        if isempty(labels)
            draw_dot(adj0(elems,elems)> 0); %, labels(elems)); 
        else
            draw_dot(adj0(elems,elems)> 0, labels(elems)); 
        end
        title(ttl);
    elseif length(elems) > 3
        figure; axis square; 
        if (~isempty(x) & ~isempty(labels)) 
            graph_draw(adj0(elems,elems)> 0, 'node_labels', labels(elems), ...   % cellstr(int2str(elems'))
               'node_shapes', zeros(length(elems),1), 'X', x(elems), 'Y', y(elems));
        elseif ~isempty(labels) 
            graph_draw(adj0(elems,elems)> 0, 'node_labels', labels(elems));
        end
        title(ttl);
    end
    disj(elems,elems) = disj(elems,elems) + adj0(elems,elems); 
end
if graph.Nnds < 100, 
    figure; axis square; graph_draw(double(disj > 0), 'node_labels',labels, ...
    'node_shapes', zeros(1, graph.Nnds), 'X', x, 'Y', y); title('Disjoint components'); 
end 

 % funcion flatten_hypernode (hypernode, hyperedge, node_id, depth)
   % reconstructs the nested structure of a hypernode for a hypernode specified by NODE_ID
   % recursively iterates for the DEPTH levels or to the basic level if DEPTH = 0
