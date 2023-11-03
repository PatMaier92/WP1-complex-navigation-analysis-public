function [abs_presence, rel_presence, abs_time_in_zone, index]=computePresence(zone_i,...
    x, y, joint_poly, time)
% computePresence: Compute absolute and relative presence in zone. 
% 
% Input: 
% zone_i is zone identifier (integer) 
% x, y are x-/y-coordinates (float) 
% joint_poly (polyshape)
% time is total time (float)
% 
% Returns:
% abs_presence (integer), rel_presence (float), 
% abs_time_in_zone(float), index (boolean vector)

% get index(es) in zone 
target_poly=joint_poly(zone_i);
index=inpolygon(x,y,target_poly.Vertices(:,1),target_poly.Vertices(:,2));

% % figure; plot(sm.coord.full_poly); hold on;
% plot(target_poly, 'FaceColor', 'k');
% plot(x(index), y(index), 'rx');

% compute values
abs_presence=numel(x(index));
rel_presence=abs_presence/numel(x);
abs_time_in_zone=time*rel_presence;
                
end
