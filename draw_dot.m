function [x, y, labels] = draw_dot(adj, drawfil, labels);
%
% [x, y] = draw_dot(adj)   draw a graph defined by adjacency matrix 
%  
% Sample code illustrating use of dot_to_graph.m function
% Leon Peshkin  pesha @ ai.mit.edu  /~pesha     Jan 2004
[n,m] = size(adj);
if n ~= m, warning('not a square adjacency matrix!'); end
if ~isequal(diag(adj),zeros(n,1)), warning('Self-loops in adjacency matrix!');end
if isequal(triu(adj,1),tril(adj,-1)'), directed = 0; else, directed = 1; end 
adj = double(adj > 0);   % make sure it is a binary matrix cast to double type
      % to be platform independant no use of directories in temp. filenames
tmpDOTfile = '_GtDout.dot';           tmpLAYOUT  = '_LAYout.dot'; 
         %  labels have to be UNIQUE, strict syntaxis, no funny chars .... 
graph_2_dot(adj, 'directed', directed, 'filename', tmpDOTfile); %'node_label', labels);  % save in file
if ispc, shell = 'dos'; else, shell = 'unix'; end  %  Which OS ?
%cmnd = strcat(shell,'(''neato -V'')');    % request version to check NEATO is there
%status = eval(cmnd);
%if status == 1,  warning('DOT/NEATO not accessible'); end
neato = '(''neato -Tdot -Gcenter -Goverlap=false -Gmaxiter=500 -o drawfil'; % -Gstart="regular" -Gregular  
cmnd = strcat([shell neato tmpLAYOUT ' ' tmpDOTfile ''')']); % -x compact
status = eval(cmnd);              %  get NEATO to layout

[trash, names, x, y] = dot_to_graph(tmpLAYOUT);   % load NEATO layout
[ignore,lbl_ndx] = sort(str2num(char(names))');  % recover from dot_to_graph permutation 
x = x(lbl_ndx); y = y(lbl_ndx);  
if nargin == 2                                    % if no labels were provided 
    labels = names(lbl_ndx);
end

if n > 40, fontsz = 7; elseif n < 12, fontsz = 12; else fontsz = 9; end 
figure; clf; axis square      %  now plot 
[x, y, h] = graph_draw(double(adj>0), 'node_labels', labels, 'fontsize', fontsz, ...
    'node_shapes', zeros(size(x,2),1), 'X', x, 'Y', y);
delete(tmpLAYOUT); delete(tmpDOTfile);     % clean up temporary files
