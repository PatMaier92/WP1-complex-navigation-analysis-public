function [tri_x,tri_y,tri]=setTriPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y)
% setTriPolyshape: Creates a polyshape for triangles. 
%
% Input: 
% n_alleys is number of available alleys (integer)
% alley_x, alley_y are x-/y-coordinates of alley corners (float vector) 
% pentagon_x, pentagon_y are x-/y-coordinates of inner pentagon (float vector) 
%
% Returns: 
% tri_x, tri_y are are x-/y-coordinates of triangles (float vector)
% tri is array of triangle polyshapes (polyshape) 

for a=1:n_alleys
    tri_x(:,a)=[alley_x(4,a);alley_x(3,a);pentagon_x(1,a);alley_x(4,a)];
    tri_y(:,a)=[alley_y(4,a);alley_y(3,a);pentagon_y(1,a);alley_y(4,a)];
    tri{a}=polyshape(tri_x(:,a),tri_y(:,a));
end

end
