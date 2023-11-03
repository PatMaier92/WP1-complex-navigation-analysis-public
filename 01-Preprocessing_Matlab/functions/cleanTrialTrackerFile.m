function data=cleanTrialTrackerFile(data,s,k)
% cleanTrialTrackerFile Data cleaning of trial tracker files. 
%
% Input: 
% data is input trial tracker data (table).
% 
% Returns: data as cleaned trial tracker data (table). 

% clean rows
data=data(data.trialEvent==1,:); % remove rows before/after trial % alternative: keep trialEvent==2 (after goal was triggered) 
data(strcmp(data.gameIsPaused,'True'),:)=[]; % remove rows when gameIsPaused
[~,i,~]=unique(data.time); % remove rows with non-unique time info
if length(unique(i))~=length(data.time)
    fprintf('Non-unique time row(s) in tracker files was/were removed in session %d, trial %d.\n', s, k); 
end
data=data(i,:); 

% clean columns 
data.pos_y=[]; 
data.rot_x=[]; 
data.rot_z=[]; 
data.gameIsPaused=[]; 
data.trialEvent=[]; 

end
    