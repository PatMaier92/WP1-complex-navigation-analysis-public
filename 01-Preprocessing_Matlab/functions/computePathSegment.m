function [seg_index]=computePathSegment(xy_ref, xy_trial)
% computePathSegment Get segment of trial path with certain path length. 
% If trial path is shorter than reference path, the segement consists of
% the whole trial path. 
% 
% Input:
% xy_ref are x-/y-coordinates for reference path (float)
% xy_trial are x-/y-coordinates for trial path (float)
%
% Returns:
% seg_index is index for path segmentation (integer) 

% compute path lengths
path_length_ref=computePathLength(xy_ref(:,1),xy_ref(:,2));
path_length_trial=computePathLength(xy_trial(:,1),xy_trial(:,2));

% set default of seg_index to whole trial 
seg_index=length(xy_trial); 

% if trial path is longer than reference path 
% find index for segmentation
if path_length_trial > path_length_ref
    temp_pl=0;
    for seg_index=1:length(xy_trial)-1
        temp_pl=temp_pl+computeDistance(xy_trial(seg_index,1),xy_trial(seg_index+1,1),...
            xy_trial(seg_index,2),xy_trial(seg_index+1,2));
        if temp_pl >= path_length_ref
            break;
        end
    end   
end

end