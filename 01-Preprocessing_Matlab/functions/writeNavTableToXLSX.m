function writeNavTableToXLSX(data_folder)
% writeNavTableToXLSX Write data summary table in xlsx format. 
%
% Input: path to data folder 
%
% Returns: data summary table as xlsx file 

% load data 
file_path=[data_folder '\WP1_results']; 
file_name='\wp1_results_navigation.mat';
full_file=fullfile(file_path, file_name); 
if isfile(full_file)
    load(full_file, 'sm'); 
else 
    disp('Your input data does not exist. Please check your data folder.'); 
    return; 
end 
data=sm.participant; clear sm; 

% process data (from structure to table) 
p=size(data,2);
for i=1:p
    for s=1:length(data(i).session)
        % save temporary session data 
        session_data=data(i).session(s).trial;
        
        % get general information 
        r=size(session_data,2);
        g_data=table(repmat(data(i).id,r,1), string(repmat(data(i).group,r,1)),...
            repmat(data(i).session(s).session_num,r,1), repmat(data(i).session(s).session_duration,r,1),...
            'VariableNames',{'id' 'group' 'session' 'session_duration'}); 
 
        % remove empty rows in data (3 exploration trials in session 1: 1, 2, 67)
        if s==1
            session_data(67)=[]; session_data(2)=[]; session_data(1)=[]; 
            g_data(1:3,:)=[];        
        end
        
        % normalize missing values in data (otherwise struct2table creates messy 0x0 cells!)
        fields = fieldnames(session_data);
        for j=1:size(session_data, 2)
            for k=1:numel(fields)
                if isempty(session_data(j).(fields{k}))
                    session_data(j).(fields{k})=999;
                end
            end
        end
        
        % covert data to table 
        t_data=struct2table(session_data); 
        gt_data=[g_data t_data];
        
        % merge data
        if i==1 && s==1
            data_table=gt_data;
        else
            data_table=outerjoin(data_table,gt_data,'MergeKeys',true);
        end
        
        clear *_data r j k;
    end
end

% sort order 
data_table=sortrows(data_table,{'id','session','trial'});

% write data 
format='yymmdd'; date=datestr(now, format); 
file_name2=['wp1_navigation_data_' date '.xlsx']; 
writetable(data_table,fullfile(file_path,file_name2));

end