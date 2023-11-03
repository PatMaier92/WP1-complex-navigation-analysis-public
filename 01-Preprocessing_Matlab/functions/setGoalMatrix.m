function [goal_x_in_alleys,goal_y_in_alleys]=setGoalMatrix(n_arms,alley_int,goal_x,goal_y)
% setGoalMatrix: Create matrix with hypothetical goal locations in each
% alley (by rotating goal locations)
% 
% Input: 
% goal_x, goal_y are x-/y-coordinates of all goal locations (float).
%
% Returns: 
% goal_x_in_alleys, goal_y_in_alleys are x-/y-coordinates of rotated locations (float).
% Impossible rotations are coded with zero. 

n_goals=length(goal_x); 
goal_x_in_alleys=zeros(n_goals,n_arms*2);
goal_y_in_alleys=zeros(n_goals,n_arms*2);
x_center = 0.5; y_center = 0.5;

% for each rotation theta
for i=0:(n_arms-1)
    % set rotation matrix
    theta = -360/n_arms*i;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    
    % for each goal
    for j=1:n_goals
        % set alley integer
        alley=alley_int(j); 
        
        % rotate
        v = [goal_x(j) ; goal_y(j)];
        center = repmat([x_center; y_center], 1, length(v));
        vo = R*(v - center) + center;
        
        % store
        integer = alley + i*2;
        if integer > n_arms*2
            integer = integer - n_arms*2;
        end
        goal_x_in_alleys(j,integer) = vo(1,1);
        goal_y_in_alleys(j,integer) = vo(2,1);
    end
end
    
% % test plot
% figure; plot(sm.coord.full_poly); hold on;
% plot(goal_x_in_alleys(1,1:6), goal_y_in_alleys(1,1:6), 'rx', ...
%     goal_x_in_alleys(2,1:6), goal_y_in_alleys(2,1:6), 'bx', ...
%     goal_x_in_alleys(3,1:6), goal_y_in_alleys(3,1:6), 'gx', ...
%     x_center, y_center, 'bo');
% xlim([0 1]); ylim([0 1]); hold off;

end