function [start_x,start_y]=setStartValues(start,xmin,xmax,ymin,ymax)
% setStartValues: Normalizes x-/y-coordinates of start positions. 
%
% Input: 
% start are x-/y-coordinates of start positions (float).
% xmin, xmax, ymin, ymax are minimum & maximum x-/y-coordinates (float).
%
% Returns: start_x,start_y are normalized x-/y-coordinates (float).

[row,col]=size(start);
for r=1:row
    start_x(r,1)=spatialNormalization(start(r,1),xmin,xmax);
    start_y(r,1)=spatialNormalization(start(r,2),ymin,ymax);
end

end