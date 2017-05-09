From pesha@bistromath0.ai.mit.edu Thu Feb 12 17:47:21 2004
Date: Thu, 12 Feb 2004 00:16:21 -0500 (EST)
From: Leon Peshkin <pesha@bistromath0.ai.mit.edu>
To: savova@jhu.edu
Subject: graph_to_dot.m

function graph_to_dot(adj, varargin)

% dag_to_dot(adj, ...)  Makes a GraphViz (AT&T) format file representing 
%                     a graph given by an adjacency matrix.
% Optional arguments should be passed as name/value pairs [default]
%
% 'filename' - if omitted, writes to 'tmp.dot'
% 'arc_label' - arc_label{i,j} is a string attached to the i-j arc [""]
% 'node_label' - node_label{i} is a string attached to the node i ["i"]
% 'width'     - width in inches [10]
% 'height'    - height in inches [10]
% 'leftright' - 1 means layout left-to-right, 0 means top-to-bottom [0]
% 'directed'  - 1 means use directed arcs, 0 means undirected [1]
%
% For details on dotty, See http://www.research.att.com/sw/tools/graphviz
%
% by Leon Peshkin, Jan 2004 modified from Kevin Murphy's BNT
                   
node_label = [];   arc_label = [];   % set default args
width = 10;        height = 20;
leftright = 0;     directed = 1;     filename = 'tmp.dot';
           
for i = 1:2:nargin-1                    % get optional args
    switch varargin{i}
        case 'filename', filename = varargin{i+1};
        case 'node_label', node_label = varargin{i+1};
        case 'arc_label', arc_label = varargin{i+1};
        case 'width', width = varargin{i+1};
        case 'height', height = varargin{i+1};
        case 'leftright', leftright = varargin{i+1};
        case 'directed', directed = varargin{i+1};
    end
end

fid = fopen(filename, 'w');
if directed
    fprintf(fid, 'digraph G {\n');
    arctxt = '->'; 
    if isempty(arc_label)
        labeltxt = '';
    else
        labeltxt = '[label="%s"]';
    end
else
    fprintf(fid, 'graph G {\n');
    arctxt = '--'; 
    if isempty(arc_label)
        labeltxt = '[dir=none]';
    else
        labeltext = '[label="%s",dir=none]';
    end
end
fprintf(fid, 'center = 1;\n');
fprintf(fid, 'size=\"%d,%d\";\n', width, height);
if leftright
    fprintf(fid, 'rankdir=LR;\n');
end
Nnds = length(adj);
for node = 1:Nnds               %  process nodes 
    if isempty(node_label)
        fprintf(fid, '%d;\n', node);
    else
        fprintf(fid, '%d [ label = "%s" ];\n', node, node_label{node});
    end
end
edgeformat = strcat(['%d ',arctxt,' %d ',labeltxt,';\n']);
for node1 = 1:Nnds   % process edges
    if directed
        arcs = find(adj(node1,:));         % children(adj, node);
    else
        arcs = find(adj(node1,node1+1:Nnds)) + node1; % remove duplicate arcs
    end
    for node2 = arcs
        fprintf(fid, edgeformat, node1, node2);
    end
end
fprintf(fid, '}');
fclose(fid);
