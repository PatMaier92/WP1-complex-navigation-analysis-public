function [goal_x,goal_y]=setGoalValues(goal,xmin,xmax,ymin,ymax)
% setGoalValues:  Normalizes x-/y-coordinates of goal positions. 
%
% Input: 
% goal is x-/y-coordinates of goal location (float).
% xmin, xmax, ymin, ymax are minimum & maximum x-/y-coordinates (float).
%
% Returns: goal_x, goal_y are normalized x-/y-coordinates (float).

[row,col]=size(goal);
for r=1:row
    goal_x(r,1)=spatialNormalization(goal(r,1),xmin,xmax);
    goal_y(r,1)=spatialNormalization(goal(r,2),ymin,ymax);
end

end