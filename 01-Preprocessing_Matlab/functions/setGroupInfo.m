function [group]=setGroupInfo(id)
% setGroupInfo Assigns group information based on participant id. 
%
% Input: 
% participant id (integer) 
%
% Returns: group(string) 

id_string = num2str(id);

% assign group based on 3rd last element of string 
if id_string(numel(id_string)-2) == '1'
    group='Delay1H';
elseif id_string(numel(id_string)-2) == '2'
    group='Delay1D';
elseif id_string(numel(id_string)-2) == '3' 
    group='Delay2W';
else
    group='F';
    disp('The ID is out of limits. Group is set to "F"')
end
    
end
