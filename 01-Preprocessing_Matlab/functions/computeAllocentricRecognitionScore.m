function [chosen_obj_i, chosen_goal_i, chosen_alley_i,...
    accuracy_loc, accuracy_alley, memory_score]=computeAllocentricRecognitionScore(...
    rand_dict, goal_s, chosen_obj_s, goal_locs, alley_locs, goal_i,...
    goal_x, goal_y, distribution, poly)
% computeAllocentricRecognitionScore Analysis of allocentric recognition task. 
%
% Input: 
% rand_dict contains randomization info (dictionary structure)
% goal_s, chosen_obj_s are this trial's (chosen) goal info (char)
% goal_locs, alley_locs are all goal info (strings)
% goal_i is this trial's goal info (integer)
% goal_x, goal_y are x-/y-coordinates of all goals (float)
% distribution are Euclidean distance values between all goals (float)
% poly is maze polyshape (for plotting only)
%
% Returns: 
% chosen_obj_i, chosen_goal_i, chosen_alley_i (integer), 
% accuracy_loc, accuracy_alley (integer)
% memory_score (float) 

% get integer for chosen object
chosen_obj_i=chosen_obj_s(1:2); 

% get chosen location from rand dict
fields = fieldnames(rand_dict); 
for i=1:length(fields)
    if string(chosen_obj_s)==string(rand_dict.(fields{i}).object)
        chosen_loc=char(fields(i));
        break
    end
end
chosen_goal_i=find(contains(goal_locs, chosen_loc));
chosen_alley_i=find(contains(alley_locs, chosen_loc(1))); 

% compute accuracy 
accuracy_loc=double(isequal(goal_s, chosen_loc));
accuracy_alley=double(isequal(goal_s(1), chosen_loc(1)));

% compute memory score 
final_distance=computeDistance(goal_x(goal_i), goal_x(chosen_goal_i), goal_y(goal_i), goal_y(chosen_goal_i)); 
memory_score=computeMemoryScore(distribution, final_distance, goal_i, 1); 

% % helper plot 
% figure; plot(poly); hold on;
% plot(goal_x(goal_i), goal_y(goal_i), 'k+');
% plot(goal_x(chosen_goal_i), goal_y(chosen_goal_i), 'r+');
% xlim([0 1]); ylim([0 1]); hold off; 

end