function [rand_dict]=setRandomizationDict(log_data, subject)
% setRandomizationDict Preprocessing of log data with information on
% randomization of goals and starts.
%
% Input: 
% cleaned log_data from csv. file (cell structure)
% subject (integer)
%
% Returns: 
% rand_dict (structure) contains randomization info for each
% participant, session and goal location. 

if size(log_data,1) ~= 16
    disp('Error: Log file is too short. Please check the file.\n');
else
    for i=1:16
        line=split(log_data(i),';');
        line=line(~cellfun('isempty',line));
        
        % id and key
        t=split(line(1));
        id=str2double(t{3});
        if id ~= subject
            disp('Error: ID mismatch detected during log data randomization check. Please check the raw data.\n');
            % break 
        end
        key=t{end}(2:3);
               
        % rest of information
        inner_dict={};
        for n=2:size(line,1)
            kv=split(line(n));
            kv=replace(kv, {'[',']'}, ''); 
            kv=replace(kv, '.', '_');
            kv=kv(~cellfun('isempty',kv));
            k=kv{1}; v=kv{2};
            inner_dict.(k)=v;
        end
        
        rand_dict.(key)=inner_dict;
    end
end

end