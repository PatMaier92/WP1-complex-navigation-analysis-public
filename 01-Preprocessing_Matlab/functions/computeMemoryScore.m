function memory_score=computeMemoryScore(distribution, final_distance,...
    goal_i, alley_i)
% computeMemoryScore: Computes memory score, i.e. Euclidean distance to goal 
% is set in relation to randomly distributed points in the WP1 Starmaze. 
%
% Input: 
% distribution is matrix of Euclidean distance values to random points (float)
% final_distance is Euclidean distance value for this trial (float)
% goal_i indicates goal location for this comparison (integer) 
% alley_i indicates desired alley for this comparison (integer)
%
% Returns: standardized memory_score (0=low, 1=high) (float) 

% compute score 
percentage_below = sum(distribution{goal_i,alley_i} < final_distance) / size(distribution{goal_i,alley_i},1);
memory_score = 1 - percentage_below;
  
end

