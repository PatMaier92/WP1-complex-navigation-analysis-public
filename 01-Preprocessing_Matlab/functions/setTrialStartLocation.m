function [start_i]=setTrialStartLocation(start_s,start_locs)
% setTrialStartLocation Return start information for this trial. 
%
% Input: 
% start_s is trial start position (string)
% start_locs are all start positions(string vector)
%
% Returns: start_i(integer)

if start_s=="X"  % X is egocentric player from alley A  
    start_i=1;
elseif start_s=="e" % last letter of 'none' for forced allocentric recognition
    start_i=999;  
else
    start_i=find(contains(start_locs, start_s)); 
end

end