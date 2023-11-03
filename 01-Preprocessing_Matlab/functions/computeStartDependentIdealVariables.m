function [x_line, y_line, o_x_line, o_y_line, x_line_chosen, y_line_chosen,...
    x_line_ego, y_line_ego, goal_x_ego, goal_y_ego, ego_alley_i, ego_zone_i,...
    ideal_path, ideal_path_chosen, ideal_ego_path]=computeStartDependentIdealVariables(...
    Graph, graph_x, graph_y, start, goal,...
    chosen_x, chosen_y, goal_x_in_alleys, goal_y_in_alleys,...
    alley_out_x, alley_out_y, alley_in_x, alley_in_y,...
    rec_x, rec_y, cp_polyshape, full_poly)
% computeStartDependentIdealVariables: Function for determining starting point dependent
% variables in Starmaze WP1. Requires Matlab 2021a for 'shortestpath' function and 
% 'interparc' function by John D'Errico (Matlab File Exchanger).
%
% Input:
% Graph represents all start-goal connections (graph)
% graph_x, graph_y are x-/y-coordinates of graph nodes (float)
% start, goal are identifiers (integer)
% chosen_x, chosen_y are final x-/y-coordinates (float)
% alley_*, rec_* are x-/y-coordinate vectors in polyshape-ready form
% (i.e. repeating initial x/y for closed shape) (float)
% cp_polyshape is inner ring as one combined polyshape (polyshape)
% full_poly is array of all polyshapes (polyshape)
%
% Returns:
% *x_line*, *y_line* are vectors with ideal x-/y-coordinates (float)
% goal_x_ego, goal_y_ego are x-/y-coordinates of egocentric goal (float)
% ego_alley [1-14], ego_zone [1-28] are identifiers for hypothetical egocentric goal alley (integer)
% ideal_path* are ideal path length values (float)

%% shortest path from original start to goal
start_node=1; 
goal_node=size(Graph.Nodes,1)-size(goal_x_in_alleys,1)+goal;
[o_x_line, o_y_line, ~]=computeShortestPath(true,...
    cp_polyshape, full_poly, Graph, graph_x, graph_y, start_node, goal_node,...
    0, 0);

% % test plot
% figure; plot(full_poly); hold on;
% plot(o_x_line, o_y_line, 'k-');
% xlim([0 1]); ylim([0 1]); hold off;

% % test plot
% [path_nodes,~]=shortestpath(Graph, start_node, goal_node);
% figure; plot(full_poly); hold on; 
% pl = plot(Graph,'XData',graph_x,'YData',graph_y);
% highlight(pl,path_nodes,'EdgeColor','r');
% xlim([0 1]); ylim([0 1]); hold off; 
% clear path_nodes; 

%% shortest path from actual start to egocentric goal (using a rotation matrix) 
% define x- and y-data for original line
v = [o_x_line ; o_y_line];
% define center of rotation
x_center = 0.5; y_center = 0.5;
% create a matrix
center = repmat([x_center; y_center], 1, length(v));
% set number of arms (rotation steps) 
n_arms = size(alley_out_x,2); 

% define rotation matrix
theta=-360/n_arms*(start-1); % to rotate 360/7° clockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

% do rotation
vo = R*(v - center) + center;

% get rotated x- and y-data
r_x_line = vo(1,:);
r_y_line = vo(2,:);

% correct rotated x- and y-data to account for measurement/rotation errors
% i.e., find nearest vertex in actual starmaze boundaries
% otherwise, your egocentric paths might be slightly off/outside the maze. 
[vid,~,~] = nearestvertex(cp_polyshape,r_x_line(2:end-1),r_y_line(2:end-1));

% save ego path and goal location
x_line_ego=[r_x_line(1); cp_polyshape.Vertices(vid,1); r_x_line(end)];
y_line_ego=[r_y_line(1); cp_polyshape.Vertices(vid,2); r_y_line(end)];
goal_x_ego=r_x_line(end); goal_y_ego=r_y_line(end);
clear vid; 

% % test plot
% figure; plot(full_poly); hold on;
% plot(o_x_line, o_y_line, 'b-', x_line_ego, y_line_ego,...
%     'rx', x_line_ego, y_line_ego, 'k-', x_center, y_center, 'bo');
% xlim([0 1]); ylim([0 1]); hold off;

% get egocentric alley & zone integer
ego_alley_i=0; ego_zone_i=0;
for c=1:n_arms
    if inpolygon(goal_x_ego,goal_y_ego,alley_out_x(:,c),alley_out_y(:,c))
        ego_alley_i=c*2-1;
        ego_zone_i=c*4-3; 
        break
    elseif inpolygon(goal_x_ego,goal_y_ego,alley_in_x(:,c),alley_in_y(:,c))
        ego_alley_i=c*2-1;
        ego_zone_i=c*4-2; 
        break
    elseif inpolygon(goal_x_ego,goal_y_ego,rec_x(:,c),rec_y(:,c))
        ego_alley_i=c*2;
        ego_zone_i=c*4; 
        break
    end
end

% calculate ideal ego path length value (external function)
ideal_ego_path=computePathLength(x_line_ego, y_line_ego);

%% shortest path from actual start to goal
goal_node=size(Graph.Nodes,1)-size(goal_x_in_alleys,1)+goal;
[x_line, y_line, ideal_path]=computeShortestPath(true,...
    cp_polyshape, full_poly, Graph, graph_x, graph_y, start, goal_node,...
    0, 0);

% % test plot
% figure; plot(full_poly); hold on;
% plot(x_line, y_line, 'k-');
% xlim([0 1]); ylim([0 1]); hold off;

%% shortest path from actual start to chosen goal 
[x_line_chosen, y_line_chosen, ideal_path_chosen]=computeShortestPath(false,...
    cp_polyshape, full_poly, Graph, graph_x, graph_y, start, 0,...
    chosen_x, chosen_y);

% % test plot
% figure; plot(full_poly); hold on; 
% plot(x_line_chosen, y_line_chosen, 'p--');
% xlim([0 1]); ylim([0 1]); hold off; 

end