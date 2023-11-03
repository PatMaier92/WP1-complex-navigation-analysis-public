function ED=computeDistance(x1,x2,y1,y2)
% computeDistance Calculates Euclidean distance between two points.
%
% Input: 
% y1, x1 coordinates (float).
% y2, x2 coordinates (float).
%
% Returns: D is Euclidean distance (float). 

ED=sqrt((x2-x1)^2+(y2-y1)^2); 