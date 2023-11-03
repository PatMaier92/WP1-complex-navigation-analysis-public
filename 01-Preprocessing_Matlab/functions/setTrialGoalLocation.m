function [goal_x,goal_y,goal_i,goal_s,alley_i,zone_i]=setTrialGoalLocation(goal_info,...
    all_goal_x,all_goal_y,goal_locs,alley_rec_locs,zone_locs)
% setTrialGoalLocation Return goal information for this trial.  
%
% Input: 
% goal_info (char)
% all_goal_x, all_goal_y are x-/y-coordinates of the goals (float)
% goal_locs are ordered goal names [1-16] (string)
% alley_rec_locs are ordered alley & rectangle names [1-14] (string)
% zone_locs are ordered zone names [1-28] (string)
%
% Returns: goal x-/y-coordinates (float), 
% goal_i (integer), goal_s (string), alley_i (integer), zone_i (integer)

% strings & integer 
len=length(goal_info); 
goal_s=string(goal_info(len-1:len));
alley_s=string(goal_info(len-1));
goal_i=find(contains(goal_locs, goal_s)); 
alley_i=find(contains(alley_rec_locs, alley_s)); 
zone_i=find(contains(zone_locs, goal_s)); 

% coordinates
goal_x=all_goal_x(goal_i); 
goal_y=all_goal_y(goal_i);

% correct for false input 
if isempty(goal_i) || isempty(alley_i) || isempty(zone_i)
    goal_i=999; 
    alley_i=999;
    zone_i=999;
    goal_x=0; 
    goal_y=0; 
end

end