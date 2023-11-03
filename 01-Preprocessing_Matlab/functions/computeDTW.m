function [dtw_distance]=computeDTW(xy_ref, xy_trial,...
    n_arms, original_start, allo_start)
% computeDTW Calculates dynamic time warping (DTW) distance between two paths.
% DTW distance is calculated according to Tao et al. (2021) 
% as root sum of squared distance normalized by max(m,n)
% (max(m,n) is smallest number of cells that need to be visited). 
%
% Input:
% xy_ref are x-/y-coordinates for reference path (float)
% xy_trial are x-/y-coordinates for trial path (float)
% n_arms is number of possible rotations (integer)
% *_start are identifier for trial start (integer) 
%
% Returns:
% dtw_distance

% rotate xy_ref to current start orientation 
v=xy_ref';
% define center of rotation
x_center=0.5; y_center=0.5;
% create a matrix
center=repmat([x_center; y_center], 1, length(v));
% define rotation matrix
theta=-360/n_arms*(allo_start-original_start); % to rotate 360/n_arms clockwise
R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];

% do rotation
vo=R*(v - center) + center;
% % helper plot
% plot(xy_ref(:,1), xy_ref(:,2), 'r+'); hold on; 
% plot(vo(1,:), vo(2,:), 'b*'); xlim([0 1]); ylim([0 1]); hold off; 

% compute DTW distance between reference path and trial path 
% root sum of squared distance normalized by max(m,n) (Tao et al., 2021)
[dtw_d,i1,i2]=dtw(vo, xy_trial', 'squared');
dtw_distance=sqrt(dtw_d/max([i1; i2])); 

% % helper plot
% figure; plot(poly); hold on;
% plot(xy_trial(:,1), xy_trial(:,2), 'k-', vo(1,:), vo(2,:), 'r-', xy_ref(:,1), xy_ref(:,2), 'b-');
% xlim([0 1]); ylim([0 1]); hold off;

end