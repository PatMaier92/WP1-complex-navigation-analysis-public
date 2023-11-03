function plotTestFigure(polyshape,G,graph_x,graph_y)
% plotTestFigure: Creates a test figure for ComplexMaze WP1.
% 
% Input: 
% polyshape is array of polyshapes
% graph is graph with nodes and edges
% graph_x, graph_y are vectors with x-/y-coordinates of nodes (float)
%
% Returns: a nice figure

figure('Position',[500 200 580 500]); 
plot(polyshape); hold on; 
set(gca,'xtick',[0 1],'ytick',[0 1]);
axis([0 1 0 1]);
p=plot(G,'XData',graph_x,'YData',graph_y);
labelnode(p,[1:numel(G.Nodes.Names)],G.Nodes.Names);
hold off; 

end