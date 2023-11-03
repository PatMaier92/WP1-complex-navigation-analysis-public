function [cp_x,cp_y,cp,pent_x,pent_y]=setPentagonValues(alley_x,alley_y,pentagon_x,pentagon_y,xmin,xmax,ymin,ymax)
% setPentagonValues Normalizes x-/y-coordinates for corners of inner pentagon  
% and creates pentagon polyshape. 
%
% Input: 
% alley_x, alley_y are x-/y-coordinates of alley corners (float vector) 
% pentagon_x, pentagon_y are x-/y-coordinates of inner pentagon (float vector) 
% xmin, xmax, ymin, ymax are minimum & maximum x-/y-coordinates (float).
%
% Returns: 
% pent_x, pent_y are normalized x-/y-coordinates of central polygon (float vector)
% cP_x, cP_y are x-/y-coordinates of central polygon in polyshape-ready form (float vector)
% cP is a polyshape of the central pentagon (polyshape)

[~,n] = size(pentagon_x);
pent_x = zeros(1,n);
for p=1:n
    pent_x(1,p)=spatialNormalization(pentagon_x(1,p),xmin,xmax);
end

pent_y = zeros(1,n);
for p=1:n
    pent_y(1,p)=spatialNormalization(pentagon_y(1,p),ymin,ymax);
end

cp_x=[alley_x(4,1);alley_x(3,1);alley_x(4,2);alley_x(3,2);alley_x(4,3);alley_x(3,3);...
    alley_x(4,4);alley_x(3,4);alley_x(4,5);alley_x(3,5);alley_x(4,6);alley_x(3,6);alley_x(4,7);alley_x(3,7);NaN;...
    pent_x(1,1);pent_x(1,7);pent_x(1,6);pent_x(1,5);...
    pent_x(1,4);pent_x(1,3);pent_x(1,2);pent_x(1,1)];

cp_y=[alley_y(4,1);alley_y(3,1);alley_y(4,2);alley_y(3,2);alley_y(4,3);alley_y(3,3);...
    alley_y(4,4);alley_y(3,4);alley_y(4,5);alley_y(3,5);alley_y(4,6);alley_y(3,6);alley_y(4,7);alley_y(3,7);NaN;...
    pent_y(1,1);pent_y(1,7);pent_y(1,6);pent_y(1,5);...
    pent_y(1,4);pent_y(1,3);pent_y(1,2);pent_y(1,1)];

cp=polyshape(cp_x,cp_y);

end