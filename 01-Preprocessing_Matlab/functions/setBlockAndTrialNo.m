function [corr_block, corr_trial_in_block]=setBlockAndTrialNo(...
    session, block, trial_in_block)
% setBlockAndTrialNo Corrects block and trial in block information. 
%
% Input: 
% session (integer)
% block (integer) 
% trial_in_block (integer)
%
% Returns: corr_block (integer), corr_trial_in_block (integer) 

corr_block = block;
corr_trial_in_block=trial_in_block;
    
if session==1 
    corr_block = block-1;
elseif session==2
    if mod(block,2)==0
        corr_block = block/2; 
        corr_trial_in_block = trial_in_block+8;
    elseif mod(block,2)==1
        corr_block = (block+1)/2; 
    end  
end

end 
