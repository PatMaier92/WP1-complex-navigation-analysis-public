function path_length=computePathLength(x_line, y_line)
% computePathLength Calculates path length as sum of Euclidean distance.
%
% Input: 
% x_line, y_line are vectors with x-/y-coordinate points (float). 
%
% Returns:  
% path length as sum of Euclidean distance (float).

path_length=0; 
for i=1:length(x_line)-1
    path_length=path_length+computeDistance(x_line(i),x_line(i+1), y_line(i),y_line(i+1)); 
end