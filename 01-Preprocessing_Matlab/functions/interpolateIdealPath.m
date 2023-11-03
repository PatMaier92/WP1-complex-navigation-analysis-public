function [xi, yi]=interpolateIdealPath(x_line, y_line, ideal_distance)
% interpolateIdealPath: Interpolating data for ideal path vectors. 
% This method takes into account distance values to create equally spaced 
% interpolated values using the 'interparc' function by John D'Errico (Matlab File Exchanger). 
%
% Input:
% x_line, y_line are vectors with x-/y-coordinates (float)
% ideal_distance is ideal path distance values (float)
%
% Returns:
% xi, yi are vectors with interpolated x-/y-coordinates (float)

% interpolate
crit=round(ideal_distance*1000);
xy=interparc(crit,x_line,y_line,'linear');
xi=xy(:,1); yi=xy(:,2);

end