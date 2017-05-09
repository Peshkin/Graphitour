function [x, y, h] = graph_plot(adj, str, max1)
% Alexi Savov 

N = size(adj,1);   % This is the number of nodes in the graph, as defined by the adjacency matrix.
str = str';        % These are the names of the nodes.

if ~isstr(str)
    str = int2str(str);
end;
if ~iscellstr(str) 
    str = cellstr(str);
end;

edges = [];           % This section creates a matrix of the edges of the graph, where each row represents
for i=1:N-1           % an edge as defined by the two nodes that make up that edge.
    for j=i+1:N
        if adj(i,j)==1 | adj(j,i)==1
            edge_ij = [i j];
            edges = [edges; edge_ij];
        end;
    end;
end;

[Y, s] = radial_layout(adj);          % This is the initial layout, constructed as described below.
K = Y;                                % K will always be our best (with fewest intersections) layout.
                                      % If an improved layout is ever found, it will be named K.
IXN = edge_intersections(Y,edges);    % This analyzes the current layout for edge intersections.
best = IXN(N+1,1);                    % IXN(N+1,1) gives the total number of intersections in the layout. It serves as a score.
k = 1;                                % 'best' will always designate the lowest score yet obtained (namely the score of K).

while k <= max1                      % This is the main loop of the program, it is executed max1 times, as specified.
    old = IXN(N+1, 1);               % The score before each optimization.
    if (old == 0)                    % If 'old' ever equals 0, no improvement can be made and the program is complete.
        break;
    end;
    if rand > .5                          % There are two optimization algorithms (see below). 
        Z = optimize_plot1(Y, adj, IXN);  % The program chooses at random.
    else
        Z = optimize_plot2(Y, IXN);
    end;
    IXN = edge_intersections(Z, edges);  % After the optimization, the new layout is analyzed for intersections.
    new = IXN(N+1,1);                    % The score of the optimized layout.
    if new < best                        % If the new score is better than the current best score, it becomes 'best' and
        best = new;                      % the layout it represents becomes our best layout K.
        K = Z;   
    end;
    if (new<old) | (rand-(new-old)/100>.5)  % Here, the program keeps the optimized layout if it is better than the old layout OR
        Y = Z;                              % if a random variable is sufficiently large. The chance of keeping the optimization
                                            % decreases as 'new' becomes much larger (worse) than 'old'. If it is larger by 50 or
        clf; axis off;                      % more intersections, there is no chance of keeping the optimization.
        for i = 1:size(edges,1)             % If one or both conditions are met, the optimization is plotted on the figure screen.
           line([Y(edges(i,1),1) Y(edges(i,2),1)], [Y(edges(i,1),2) Y(edges(i,2),2)]);
        end;
        for i=1:N
            g = extent(text(Y(i,1), Y(i,2), str(i), 'HorizontalAlignment', 'center'));
            rectangle('Position', g, 'Curvature', [1 1], 'FaceColor','w');
            text(Y(i,1), Y(i,2), str(i), 'HorizontalAlignment', 'center');
        end;
        pause(.1);
    end;
    k = k+1;
end;

X = zeros(N,2);                     % The best layout obtained through the optimization is checked for visually confusing
while any(X(:)-K(:))                % properties: (1) whether 3 or more nodes lie on the same line and are connected in
    X = K;                          % such a way that it is no longer clear which node is connected to which, and (2) whether
    K = colinear(X, edges, adj, s); % two or more nodes have been plotted too close to each other. These problems are remedied
    K = too_close(K, str, s);       % as described below in the respective algorithms.
end;

if ~isempty(s)                  % s is the set of nodes of degree 0. These nodes (if any) are plotted on the bottom of the screen.
    bottom = min(K(:,2));
    left = min(K(:,1));
    right = max(K(:,1));
    for i=1:size(s,2)
        K(s(i),1) = left+i*(right-left)/size(s,2);
    end;
    K(s,2) = bottom-.1;
end;

clf; axis off;                  % Finally, the best layout, K, is plotted in the figure window.
for i = 1:size(edges,1)
    line([K(edges(i,1),1) K(edges(i,2),1)], [K(edges(i,1),2) K(edges(i,2),2)]);
end;
for i=1:N
    g = extent(text(K(i,1), K(i,2), str(i), 'HorizontalAlignment', 'center'));
    rectangle('Position', g, 'Curvature', [1 1], 'FaceColor','w');
    text(K(i,1), K(i,2), str(i), 'HorizontalAlignment', 'center');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, s] = radial_layout(adj)    % This function produces the initial layout

N = size(adj,1);
X = zeros(N,2);
unplotted = 1:N;                        % 'unplotted' will be the list of nodes for whom coordinates have not yet been determined.
angl = zeros(N,1);                      % This matrix records the angle the position vecor of each node makes with the X-axis.

for i=1:N                               % This section determines the degree of each node.
    adj(i,i) = 0;
    out = find(adj(i,:));
    into = find(adj(:,i))';
    deg(i) = size(union(out,into), 2);
end;

s=find(deg==0);                         % This section determines the nodes of 0 degree. Their coordinates are at the origin for now.
if ~isempty(s)
    unplotted(s) = 0;
end;
n = max(deg);                       % All nodes that have the maximum degree are identified.
row = find(deg==n);                 % They become the 1st row.
next_level = [];                    % This will eventually be filled with the nodes that the nodes of 'row' are connected to.
n = size(row,2);
p = n;

for i = 1:n                           % The nodes of the 1st row create a "ring of power," a circle centered at the origin.
    angl(row(i)) = 2*pi*i/n;
    X(row(i),1) = (p-1)*cos(angl(row(i)));
    X(row(i),2) = (p-1)*sin(angl(row(i)));
end;

unplotted(row) = 0;                 % The plotted nodes are removed from 'unplotted'.
unplotted = nonzeros(unplotted)';
r = 1;                                % r is the level of each row, it is 1 for the "ring of power".

while ~isempty(unplotted)           % This loop runs for as long as there are unplotted nodes.
    for i=1:n
        out = find(adj(row(i),unplotted));
        into = find(adj(unplotted,row(i)));
        nxt = union(unplotted(out),unplotted(into));  % For each node row(i) of the current row, a set of offspring ('nxt') is made
        m = size(nxt,2);                            % up of the nodes that are connected to row(i) and are as of yet unplotted.
        for j=1:m
            if r==1                                 % In the first run-through of the loop, each ring node's
                angl(nxt(j)) = (-1)^r*pi/12+2*pi*i/n+2*pi*j/m; % offspring are plotted in a circle around their parent.
                X(nxt(j),1) = X(row(i),1)+p*cos(angl(nxt(j))); % (-1)^r*pi/12 offsets each node so as to prevent it from falling in
                X(nxt(j),2) = X(row(i),2)+p*sin(angl(nxt(j))); % line with a previous edge.
            else                                               % On subsequent levels, offspring are plotted away
                angl(nxt(j)) = (-1)^r*pi/12+angl(row(i))-pi/2+j*pi/(m+1);   % from the center of the graph, in a halfcircle around
                X(nxt(j),1) = X(row(i),1)+p*cos(angl(nxt(j))); % their parent, once again they are offset by a small
                X(nxt(j),2) = X(row(i),2)+p*sin(angl(nxt(j))); % angle.
            end;
        end;
        for j = nxt                   % The nodes just plotted are removed from the list of unplotted nodes.
            u = find(unplotted==j);

            unplotted(u) = 0;
        end;
        unplotted = nonzeros(unplotted)';
        next_level = [next_level nxt];     % The offspring of all nodes in a row create the next level of nodes.
    end;
    row = next_level;                      % This next level then becomes the row and the loop is re-executed.
    if isempty(row)                        % An empty row would mean that there are still unplotted nodes, whose degree is greater
                                           % than 0, but nopne of them are connected to any of the nodes we've already plotted.
        text(0, 0, 'This graph is not connected', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); break; 
    end;
    n = size(row,2);
    next_level = [];
    r = r+1;
    p = 2*p/3;                             % The distance between a parent and its offpsring decreases with every level.
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = colinear(X, edges, adj, s)    % This function checks if any 3 or more nodes are plotted ambiguously in a
                                           % straight line. Such nodes are then shifted slightly.
N = size(X,1);
cnt = 1:N;
cnt(s) = 0;                                % Nodes of 0 degree are ignored.
cnt = nonzeros(cnt)';

for i=cnt                                  % For each node, the edges in which it does not participate are found.
    a = find(edges(:,1) ~= i);
    b = find(edges(:,2) ~= i);
    chk = intersect(a, b)';                % These edges will be checked to see if the node lies on any of them.
    for j=chk
        k = X(edges(j,1),:)-X(edges(j,2),:); % This is the vector of the edge.
        l = X(i,:) - X(edges(j,1),:);        % Ths is a vector from the node to one of the nodes of the edge.
        m = subspace(k', l');                 % If the two are nearly parallel (linearly dependent), then the 3 nodes almost
        if m < .04                             % lie on the same line.
            if (adj(i,edges(j,1))==1 | adj(edges(j,1),i)==1) % If this is true, node i is connected to one of the two nodes
                X(i,:) = X(i,:)+.1*[-k(2) k(1)];             % of the edge, which makes visually unclear which ones are connected.
                                                             % It is therefore shifted perpendicularly to the edge.
            elseif (X(i,1)>min(X(edges(j,[1 2]),1)) & X(i,1)<max(X(edges(j,[1 2]),1))) & (X(i,2)>min(X(edges(j,[1 2]),2)) & X(i,2)>min(X(edges(j,[1 2]),2)))
                X(i,:) = X(i,:)+.1*[-k(2) k(1)];             %In this case node i lies directly on the edge. It is again shifted.
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = too_close(X, str, s)             % This function determines if any two nodes have been plotted on top of each other.
N = size(X,1);
g = 0;
cnt = 1:N;
cnt(s) = 0;
cnt = nonzeros(cnt)';

for i=cnt
    h = extent(text( 0, 0, str(i), 'visible', 'off'));
    g1 = max(h([3 4])); % This is the size of the bubble around each node.
    if g1 > g
        g = g1;         % The largest bubble is taken as the standard for separation.
    end;
end;
    
for i=1:size(cnt,2)-1
    for j=i+1:size(cnt,2)
        d = X(cnt(i),:)-X(cnt(j),:);
        e = norm(d);
        if e <= 1.5*g     % If any two nodes are plotted within 1.5 times the standard for separation, they are broken apart.
            X(cnt(j),:) = X(cnt(j),:)+3*g*d/e;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IXN = edge_intersections(X,edges)     % This algorithm studies the layout for intersections.
                                               % by taking the cross product of every two sides of the
N = size(X,1);                                 % quadrangle given by any two edges. If the cross products
E = size(edges,1);                             % at all four corners of a quadrangle have the same sign, 
IXN = zeros(N+1,N+1);                          % then the edges that make up its diagonal intersect.
sides = zeros(4,4);
vertex_angles = zeros(4,3);

for i=1:E-1
    for j = i+1:E
        skip = intersect(edges(i,[1 2]),edges(j,[1 2]));
        if isempty(skip)
            sides(1,[3 4]) = X(edges(j,2),:)-X(edges(i,1),:);
            sides(4,[1 2]) = -sides(1,[3 4]);
            
            sides(1,[1 2]) = X(edges(j,1),:)-X(edges(i,1),:);
            sides(3,[3 4]) = -sides(1,[1 2]);
            
            sides(2,[3 4]) = X(edges(j,1),:)-X(edges(i,2),:);
            sides(3,[1 2]) = -sides(2,[3 4]);
            
            sides(2,[1 2]) = X(edges(j,2),:)-X(edges(i,2),:);
            sides(4,[3 4]) = -sides(2,[1 2]);
 
            for p = 1:4
                vertex_angles(p,1:3) = cross([sides(p,[1 2]) 0],[sides(p,[3 4]) 0]);
            end;
            
            vertex_angles = vertex_angles(:,3);
            
            if all(vertex_angles < 0) | all(vertex_angles > 0)
                IXN(edges(i,[1 2]),1) = IXN(edges(i,[1 2]),1)+1;
                IXN(edges(i,[1 2]),edges(j,1)+1) = IXN(edges(i,[1 2]),edges(j,1)+1)+1;
                IXN(edges(i,[1 2]),edges(j,2)+1) = IXN(edges(i,[1 2]),edges(j,2)+1)+1;
                IXN(edges(j,[1 2]),1) = IXN(edges(j,[1 2]),1)+1;
                IXN(edges(j,[1 2]),edges(i,1)+1) = IXN(edges(j,[1 2]),edges(i,1)+1)+1;
                IXN(edges(j,[1 2]),edges(i,2)+1) = IXN(edges(j,[1 2]),edges(i,2)+1)+1;
            end;
        end;    
    end;
end;

IXN(N+1) = sum(IXN(1:N))/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = optimize_plot1( X, adj, IXN)   % This algorithm moves the node whose edges are involved
                                            % in the most intersections towards its center of gravity
N = size(adj, 1);                           % by a random amount.
m = max(IXN(1:N,1));
trbl = find(IXN(1:N,1)==m);
Z = [0 0];

if IXN(N+1,1)==factorial(N)
    return;
end;

for i=trbl'
    adj(i,i) = 0;
    out = find(adj(i,:));
    into = find(adj(:,i));
    csns = union(out,into);
    Z(1) = sum(X(csns,1))/size(csns,2);
    Z(2) = sum(X(csns,2))/size(csns,2);
    X(i,:) = 2*rand*([Z(1) Z(2)]-X(i,:));
    pause(.1);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = optimize_plot2( X, IXN)               % This part swaps two troublesome nodes
                                                   % to try and untangle them.
N = size(X, 1);

if (IXN(N+1,1)==factorial(N)) | (IXN(N+1,1)==0)
    return;
end;

m = max(IXN(1:N,1));
tbl = find(IXN(1:N,1)==m)';
d = 0;

for i = tbl
    s = max(IXN(i,2:N+1));
    s = find(IXN(i,2:N+1)==s);
    r = randperm(size(s,2));
    if (IXN(N+1,1)==factorial(N)) | (IXN(N+1,1)==0)
       return;
    end
end;

m = max(IXN(1:N,1));
tbl = find(IXN(1:N,1)==m)';
d = 0;

for i = tbl
    s = max(IXN(i,2:N+1));
    s = find(IXN(i,2:N+1)==s);
    r = randperm(size(s,2));
    s = s(r(1));
    if IXN(i,s+1)>d
        tsp = i;
        d = IXN(i,s+1);
        swp = s;
    end;
end;

q = X(tsp,[1 2]);
X(tsp,[1 2]) = X(swp,[1 2]);
X(swp,[1 2]) = q;
