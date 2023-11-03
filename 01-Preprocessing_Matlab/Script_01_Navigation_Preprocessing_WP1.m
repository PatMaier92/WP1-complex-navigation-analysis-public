clear; close all; clc; format compact;
addpath(genpath(pwd)) % add subfolder functions to path

%% Starmaze Data Processing
% @ date 2019-10-01 @ author Deetje Iggena
% @ date 2020-09-25 update by @ Patrizia Maier & now tracked by git
% for complex maze version WP1 201022
% requires Matlab 2021a or later (for shortestpath())
% requires Signal Processing Toolbox for downsample() and dtw()T:\Analysis\WP1\WP1_data

% Input (.csv files)
% One tracker file per trial with timestamp, x- and y-coordinates for movement
% and z-coordinates for rotation.
% One trial_results file with general experimental information.

% Block 1: Set up input and output folders, Starmaze and Practise environment
% Block 2: Data preparation 
% Block 3: Variable computation 
% Block 4: Write data to xlsx file

%% Block 1: Set up input and output folders, maze environment 
%% data input folder and participant information
[data_folder] = setInputPath(); % provide folder with all raw data
[participant_start,participant_end] = setParticipants(); % provide participant range
n_sessions = 3; 

%% result folder
result_folder=[data_folder '\WP1_results'];
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
    disp('Your output folder did not exist, a new folder was created.')
end

%% load data table or create new one 
% load existing data
file_name         = '\wp1_results_navigation.mat';
file_path         = fullfile(result_folder, file_name);
if isfile(file_path)
    load(file_path)
end

% initialize if non-existing
if ~exist('sm','var')
    sm = []; 
    
    %% set up maze environment 
    % coordinates of min/max values
    values=table2array(readtable('wp1_values.csv'));
    [sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax]=setMinMaxValues(values);

    % coordinates of start positions (normalized!)
    start=table2array(readtable('wp1_start.csv'));
    [sm.coord.start_x,sm.coord.start_y]=setStartValues(start,sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax); 
    
    % coordinates of alley corners (normalized!)
    alley_x=table2array(readtable('wp1_alley_x.csv'));
    [n_corners,n_alleys] = size(alley_x);
    for i_alley=1:n_alleys
        for i_corner=1:n_corners
            alley_x(i_corner,i_alley)=spatialNormalization(alley_x(i_corner,i_alley),sm.coord.xmin,sm.coord.xmax);
        end
    end
    alley_y=table2array(readtable('wp1_alley_y.csv'));
    for i_alley=1:n_alleys
        for i_corner=1:n_corners
            alley_y(i_corner,i_alley)=spatialNormalization(alley_y(i_corner,i_alley),sm.coord.ymin,sm.coord.ymax);
        end
    end

    % coordinates of combined pentagon (normalized!) and central polyshape
    pentagon_x=table2array(readtable('wp1_pentagon_x.csv'));
    pentagon_y=table2array(readtable('wp1_pentagon_y.csv'));
    [sm.coord.central_x,sm.coord.central_y,sm.coord.central_poly,pentagon_x,pentagon_y]=setPentagonValues(alley_x,alley_y,pentagon_x,pentagon_y,...
        sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax);

    % create other polyshapes
    % alley polyshape
    [sm.coord.alley_full_x,sm.coord.alley_full_y,sm.coord.alley_poly,...
        sm.coord.alley_half_out_x,sm.coord.alley_half_out_y,sm.coord.alley_poly_out,...
        sm.coord.alley_half_in_x,sm.coord.alley_half_in_y,sm.coord.alley_poly_in]=setAlleyPolyshape(alley_x,alley_y);
    % rectangle polyshape
    [sm.coord.rec_x,sm.coord.rec_y,sm.coord.rec_poly]=setRectPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y);
    % triangle polyshape
    [sm.coord.tri_x,sm.coord.tri_y,sm.coord.tri_poly]=setTriPolyshape(n_alleys,alley_x,alley_y,pentagon_x,pentagon_y);
    
    % create joint polyshapes 
    % joint alley & triangle (unionized) polyshape
    alley_tri_polyshape=union(sm.coord.alley_poly{1}, sm.coord.tri_poly{1}); 
    for i=2:numel(sm.coord.alley_poly)
        alley_tri_polyshape=[alley_tri_polyshape union(sm.coord.alley_poly{i}, sm.coord.tri_poly{i})];
    end
    sm.coord.alley_tri_poly=alley_tri_polyshape; clear alley_tri_polyshape;
    % joint zone polyshape
    zone_poly=[sm.coord.alley_poly_out{1}, sm.coord.alley_poly_in{1}, sm.coord.tri_poly{1}, sm.coord.rec_poly{1}];
    for i=2:numel(sm.coord.alley_poly)
        zone_poly=[zone_poly sm.coord.alley_poly_out{i}, sm.coord.alley_poly_in{i}, sm.coord.tri_poly{i}, sm.coord.rec_poly{i}];
    end
    sm.coord.zone_poly=zone_poly; clear zone_poly; 
    % joint full polyshape
    sm.coord.full_poly=[sm.coord.alley_poly_out{1,1} sm.coord.alley_poly_in{1,1}...
        sm.coord.alley_poly_out{1,2} sm.coord.alley_poly_in{1,2}...
        sm.coord.alley_poly_out{1,3} sm.coord.alley_poly_in{1,3}...
        sm.coord.alley_poly_out{1,4} sm.coord.alley_poly_in{1,4}...
        sm.coord.alley_poly_out{1,5} sm.coord.alley_poly_in{1,5}...
        sm.coord.alley_poly_out{1,6} sm.coord.alley_poly_in{1,6}...
        sm.coord.alley_poly_out{1,7} sm.coord.alley_poly_in{1,7}...
        sm.coord.central_poly];
    
    % coordinates of goal positions (normalized!)
    goal=table2array(readtable('wp1_goal.csv'));
    [sm.coord.goal_x,sm.coord.goal_y]=setGoalValues(goal,sm.coord.xmin,sm.coord.xmax,sm.coord.ymin,sm.coord.ymax);
    sm.coord.goal_alley_i=setAlleyInteger(sm.coord.goal_x,sm.coord.goal_y,sm.coord.alley_poly,sm.coord.rec_poly);
    sm.coord.start_alley_i=setAlleyInteger(sm.coord.start_x,sm.coord.start_y,sm.coord.alley_poly,sm.coord.rec_poly);
    [sm.coord.goal_x_in_alleys,sm.coord.goal_y_in_alleys]=setGoalMatrix(numel(sm.coord.alley_poly),sm.coord.goal_alley_i,...
        sm.coord.goal_x,sm.coord.goal_y);

    % compute random x-/y-coordinate distribution in Starmaze
    [sm.coord.random_x_1000, sm.coord.random_y_1000]=computeRandomLocations(sm.coord.full_poly, 1000);
    
    % compute Euclidean distance between locations and random locations
    sm.coord.final_distance_distribution_1000=computeDistanceDistribution(sm.coord.goal_x_in_alleys,...
        sm.coord.goal_y_in_alleys, sm.coord.random_x_1000, sm.coord.random_y_1000);
    
    % compute Euclidean distance between all goal locations 
    sm.coord.final_distance_distribution_goals=computeDistanceDistribution(sm.coord.goal_x,...
        sm.coord.goal_y,sm.coord.goal_x,sm.coord.goal_y); 
    
    % create graph
    % for automated shortest path calculation (requires Matlab 2021a)
    [sm.coord.graph,sm.coord.graph_x,sm.coord.graph_y]=setGraph(sm.coord.start_x,sm.coord.start_y,...
        sm.coord.tri_x,sm.coord.tri_y,sm.coord.goal_x,sm.coord.goal_y);

     % information (ordered)
     sm.coord.goal_names=["C1" "C2" "D1" "E1" "E2" "F1" "G1" "G2" "I1" "I2" "J1" "K1" "K2" "L1" "M1" "M2"];
     sm.coord.start_names=["A" "C" "E" "G" "I" "K" "M"];
     sm.coord.alley_rec_names=["A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N"];
     sm.coord.zone_names=["A2" "A1" "Tab0" "B1" "C2" "C1" "Tcd0" "D1" "E2" "E1" "Tef0" "F1" ...
         "G2" "G1" "Tgh0" "H1" "I2" "I1" "Tij0" "J1" "K2" "K1" "Tkl0" "L1" "M2" "M1" "Tmn0" "N1"];
     
%     % create test figure plot
%     plotTestFigure(sm.coord.full_poly,sm.coord.graph,sm.coord.graph_x,sm.coord.graph_y);

    clear i j m distribution alley_x alley_y pentagon_x pentagon_y i_alley n_alleys i_corner n_corners...
        start goal values alley_rec_poly alley_tri_polyshape; 
    
    %% initialize participant index
    p=1; 
else
    p=0; % default value 
end

%% Block 2: Data preprocessing
for id=participant_start:participant_end
tic;     
    % set participant index
    if p~=1 
        % check if ID exists in data 
        p_ind = find([sm.participant.id]==id); 
        if isempty(p_ind) % if not: append participant data    
            [~,n]=size(sm.participant);
            p=n+1; clear n p_ind; 
        else % if yes: overwrite participant data
            p=p_ind; clear p_ind;
            fprintf('Data for participant %d already existed and is overwritten now.\n', id);
        end
    end 
    
    % loop over sessions 
    for s=1:n_sessions        
        %% set individual input and output folder
        % input folder
        input_folder=[data_folder '\' num2str(id) '\S00' num2str(s)]; 
        if ~exist(input_folder, 'dir')
            fprintf('Folder for participant %d, session %d does not exist. Continue with next iteration.\n', id, s);
            continue;
        end
        
        % output folder (only for trial plots)
        output_folder=[data_folder '\' num2str(id) '\plots'];
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
            fprintf('Your output folder for participant %d did not exist, a new folder was created.\n', id);
        end
        
        %% read-in trial results file
        opts=detectImportOptions([input_folder, '\trial_results.csv']);
        opts=setvaropts(opts,'timestamp','InputFormat','MM/dd/uuuu hh:mm:ss aa');
        trial_data=readtable([input_folder, '\trial_results.csv'], opts); clear opts; 
        
        % store participant information
        sm.participant(p).id=id;
        sm.participant(p).group=setGroupInfo(sm.participant(p).id);
        sm.participant(p).session(s).session_num=trial_data.session_num(1,1);
        
        %% read-in log file (only once in session 1)
        if s==1
            opts=detectImportOptions([input_folder, '\log.csv'], 'VariableNamesLine', 1, 'Delimiter', ',');
            opts.DataLines=[2,Inf]; 
            opts.SelectedVariableNames = {'message'};
            log_data=table2cell(readtable([input_folder, '\log.csv'], opts));
            log_data=log_data(contains(log_data,'ID is'));
            [sm.participant(p).rand_dict]=setRandomizationDict(log_data, id);
            clear log_data opts; 
        end 
        
        %% get individual trial tracker file names 
        d = dir(fullfile(input_folder, 'tracker_movement*.csv'));
        files = {d.name}; clear d;
    
        % loop over trials 
        for k=1:numel(files)
            name=files{k};
            
            % skip exploration trials 
            if s==1
                pattern=("_T001"|"_T002"|"_T067");
                if contains(name, pattern)
                    continue
                end
            end
            
            % read-in and clean trial tracker data
            data=readtable(fullfile(input_folder, name));
            data=cleanTrialTrackerFile(data,s,k);
                      
            % extract data 
            t=data.time; % time
            x=data.pos_x; % coordinates
            y=data.pos_z; % coordinates
            r=data.rot_y; R=deg2rad(r); Q=unwrap(R); r=Q; % yaw rotations
            clear data R Q;
            
            % spatial normalization
            x=spatialNormalization(x,sm.coord.xmin,sm.coord.xmax); y=spatialNormalization(y,sm.coord.ymin,sm.coord.ymax);
            % save end points 
            sm.participant(p).session(s).trial(k).x_n=x(end); sm.participant(p).session(s).trial(k).y_n=y(end);
            
            % temporal normalization 
            % get original sampling rate 
            original_sampling_rate=(t(end)-t(1))/length(t); 
            % set new sampling rate (default 0.05 corresponds to 20 frames per second)    
            new_sampling_rate=0.05;
            [ti, xi, yi, ri]=temporalNormalization(new_sampling_rate, t, x, y, r); 
            clear t x y r *sampling_rate; 
            
            %% get single trial info from trial_results
            sm.participant(p).session(s).session_duration=round(minutes(trial_data.timestamp(numel(files),1) - trial_data.timestamp(1,1))); 
            
            sm.participant(p).session(s).trial(k).trial=trial_data.trial_num(k,1);
            [sm.participant(p).session(s).trial(k).block, ...
                 sm.participant(p).session(s).trial(k).trial_in_block]=setBlockAndTrialNo(...
                 sm.participant(p).session(s).session_num, trial_data.block_num(k,1), trial_data.trial_num_in_block(k,1)); 
                     
            sm.participant(p).session(s).trial(k).goal_vis=string(trial_data.trial_goal_vis(k,1));
            sm.participant(p).session(s).trial(k).arrow_vis=string(trial_data.trial_path_vis(k,1)); 
            
            sm.participant(p).session(s).trial(k).trial_type=string(trial_data.trial_type(k,1));
            
            sm.participant(p).session(s).trial(k).goal_identity=string(trial_data.trial_goal_identity(k,1));
            
            [sm.participant(p).session(s).trial(k).goal_x,sm.participant(p).session(s).trial(k).goal_y,...
                sm.participant(p).session(s).trial(k).goal_i, goal_s,...
                goal_alley,goal_zone]=setTrialGoalLocation(char(string(trial_data.trial_goal(k,1))),...
                sm.coord.goal_x,sm.coord.goal_y,sm.coord.goal_names,sm.coord.alley_rec_names,sm.coord.zone_names);
            
            start_name=char(trial_data.trial_player(k,1)); start_s=string(start_name(end)); clear start_name; 
            [sm.participant(p).session(s).trial(k).start_i]=setTrialStartLocation(start_s,sm.coord.start_names);
 
            %% For allocentric / place ability trial (no navigation)
            if sm.participant(p).session(s).trial(k).trial_type=="testA"
                [chosen_object_i, ~, ~, sm.participant(p).session(s).trial(k).correct_object, ~,...
                    sm.participant(p).session(s).trial(k).memory_score]=computeAllocentricRecognitionScore(...
                    sm.participant(p).rand_dict, char(goal_s),...
                    char(trial_data.chosen_goal(k,1)), sm.coord.goal_names, sm.coord.alley_rec_names,...
                    sm.participant(p).session(s).trial(k).goal_i, sm.coord.goal_x, sm.coord.goal_y,...
                    sm.coord.final_distance_distribution_goals, sm.coord.full_poly);
                % fprintf('Allocentric / place ability analysis done for %d, session %d, file no %d.\n', id, s, k);
            else
                sm.participant(p).session(s).trial(k).correct_object=999; 
                
             %% For all navigation trials               
                %% compute support variables depending on this trial's settings
                % ideal path coordinates & ideal path length
                % Requires Matlab 2021a for 'shortestpath' function and
                % 'interparc' function by John D'Errico (Matlab File Exchanger)
                [x_line, y_line, ~, ~, x_line_chosen, y_line_chosen, x_line_ego, y_line_ego,...
                    sm.participant(p).session(s).trial(k).goal_x_ego, sm.participant(p).session(s).trial(k).goal_y_ego,...
                    ego_alley, ego_zone, ~, ideal_chosen_path_length, ~]=computeStartDependentIdealVariables(...
                    sm.coord.graph, sm.coord.graph_x, sm.coord.graph_y,...
                    sm.participant(p).session(s).trial(k).start_i, sm.participant(p).session(s).trial(k).goal_i,...
                    sm.participant(p).session(s).trial(k).x_n, sm.participant(p).session(s).trial(k).y_n,...
                    sm.coord.goal_x_in_alleys, sm.coord.goal_y_in_alleys,...
                    sm.coord.alley_half_out_x, sm.coord.alley_half_out_y, sm.coord.alley_half_in_x, sm.coord.alley_half_in_y,...
                    sm.coord.rec_x, sm.coord.rec_y, sm.coord.central_poly, sm.coord.full_poly);           
                
%                 % test plot
%                 figure; plot(sm.coord.full_poly); hold on; 
%                 plot(xi, yi, 'k-', x_line, y_line, 'b*',... 
%                     x_line_ego, y_line_ego, 'r+',...
%                     x_line_chosen, y_line_chosen, 'gx');
%                 xlim([0 1]); ylim([0 1]); hold off; 
                               
                %% Block 3: Variable computation 
                %% accuracy analysis
                % compute chosen goal location
                [~, ~, chosen_zone_i, ~, chosen_alley_i, ~]=computeChosenGoals(...
                    sm.participant(p).rand_dict, char(trial_data.chosen_goal(k,1)), ...
                    sm.coord.alley_poly_out, sm.coord.alley_poly_in, sm.coord.rec_poly, sm.coord.tri_poly,...
                    sm.coord.zone_names, sm.coord.alley_rec_names, sm.coord.goal_names,...
                    sm.participant(p).session(s).trial(k).x_n, sm.participant(p).session(s).trial(k).y_n);
                
                % evaluate chosen goal location
                if sm.participant(p).session(s).trial(k).goal_vis=="false"  
                    % EUCLIDEAN DISTANCE to ORIGINAL & ALLOCENTRIC / PLACE goal
                    final_distance=computeDistance(...
                        sm.participant(p).session(s).trial(k).goal_x, sm.participant(p).session(s).trial(k).x_n, ...
                        sm.participant(p).session(s).trial(k).goal_y, sm.participant(p).session(s).trial(k).y_n);
                    % EUCLIDEAN DISTANCE to EGOCENTRIC / RESPONSE goal
                    final_distance_ego=computeDistance(...
                        sm.participant(p).session(s).trial(k).goal_x_ego, sm.participant(p).session(s).trial(k).x_n,...
                        sm.participant(p).session(s).trial(k).goal_y_ego, sm.participant(p).session(s).trial(k).y_n);
                    
                    if sm.participant(p).session(s).session_num~=2
                        % CORRECT ALLEY (for baseline learning criterion)
                        sm.participant(p).session(s).trial(k).correct_alley=goal_alley==chosen_alley_i;
                        
                        % MEMORY SCORE to goal 
                        sm.participant(p).session(s).trial(k).memory_score=computeMemoryScore(sm.coord.final_distance_distribution_1000,...
                            final_distance, sm.participant(p).session(s).trial(k).goal_i, goal_alley);
                        
                        % default values
                        sm.participant(p).session(s).trial(k).strategy_score_allo=999;
                        sm.participant(p).session(s).trial(k).strategy_score_ego=999;
                    elseif sm.participant(p).session(s).session_num==2
                        % STRATEGY SCORE to ALLOCENTRIC / PLACE goal
                        sm.participant(p).session(s).trial(k).strategy_score_allo=computeMemoryScore(sm.coord.final_distance_distribution_1000,...
                            final_distance, sm.participant(p).session(s).trial(k).goal_i, goal_alley);
                        
                        % STRATEGY SCORE to EGOCENTRIC / PLACE goal
                        sm.participant(p).session(s).trial(k).strategy_score_ego=computeMemoryScore(sm.coord.final_distance_distribution_1000,...
                            final_distance_ego, sm.participant(p).session(s).trial(k).goal_i, ego_alley);
                        
                        % default values
                        sm.participant(p).session(s).trial(k).correct_alley=999;
                        sm.participant(p).session(s).trial(k).memory_score=999;
                    end
                else
                    % default values
                    sm.participant(p).session(s).trial(k).correct_alley=999; 
                    sm.participant(p).session(s).trial(k).memory_score=999; 
                    sm.participant(p).session(s).trial(k).strategy_score_allo=999; 
                    sm.participant(p).session(s).trial(k).strategy_score_ego=999; 
                end          
                % fprintf('Accuracy analysis done for %d, session %d, file no %d.\n', id, s, k);
                
                %% time analysis
                % TIME
                sm.participant(p).session(s).trial(k).time=ti(end)-ti(1);
                
                %% standard coordinate analysis using x-/y-coordinates
                % PATH LENGTH 
                path_length=computePathLength(xi,yi); 
                
                % EXCESS PATH LENGTH (to chosen goal) 
                sm.participant(p).session(s).trial(k).excess_path_length=path_length-ideal_chosen_path_length;

                %% rotation analysis
                % TOTAL ROTATION
                % calculate total rotation as cumulative absolute change in yaw rotation (r)
                % this value includes rotation due to x-/y-trajectory (i.e. left-forward movement)
                [total_rotation]=computeRotation(ri);
                    
                % INITIAL ROTATION in this trial's START AREA (alley plus triangle)
                % same method as above 
                % get rotation index
                [~, ~, ~, rot_index]=computePresence(sm.participant(p).session(s).trial(k).start_i,xi,yi,...
                    sm.coord.alley_tri_poly, sm.participant(p).session(s).trial(k).time);
                [rot_index]=computeFirstSegment(rot_index); 
                % compute initial rotation 
                [initial_rotation]=computeRotation(ri(rot_index));

                % INITIAL ROTATION VELOCITY in this trial's START AREA (alley plus triangle)
                % also called angular velocity (IdPhi) 
                time_index=ti(rot_index); initial_time=time_index(end)-time_index(1);
                sm.participant(p).session(s).trial(k).initial_rotation_velocity=...
                    initial_rotation/initial_time;
                clear rot_index time_index initial_time;
                % fprintf('Standard time, path and rotation analysis done for %d, session %d, file no %d.\n', id, s, k);
                
                %% additional behavioral analyses for session 2
                % save path trajectories from egocentric / response ability (testE) in session 1 
                % for path similarity analysis with dynamic time warping (DTW)
                if sm.participant(p).session(s).session_num==1 && sm.participant(p).session(s).trial(k).trial_type=="testE" 
                    dtw_dict_testE.(goal_s)=[xi yi];
                end 
                
                % analyse probe trials in session 2 
                if sm.participant(p).session(s).session_num==2 && sm.participant(p).session(s).trial(k).trial_type=="testN"     
                    % DYNAMIC TIME WARPING DISTANCE to EGO / RESPONSE ABILITY (testE) in session 1
                    [sm.participant(p).session(s).trial(k).dtw_to_testE]=computeDTW(...
                        dtw_dict_testE.(goal_s),[xi yi],...
                        size(sm.coord.alley_full_x,2), 1, sm.participant(p).session(s).trial(k).start_i);  
                else 
                    % default values
                    sm.participant(p).session(s).trial(k).dtw_to_testE=999;
                end
                % fprintf('Additional analysis for session 2 done for %d, session %d, file no %d.\n', id, s, k);
                                        
            end
            %% set marker for excluded trials
            % criteria: timeout, any start-goal combination where direct path is ambiguous (S2-G11, S3-G14, S6-G3, S7-G6)
            % or no movement/very short trial time (i.e. path_length=0, rotation=0, or time < 3)
            sm.participant(p).session(s).trial(k).exclude_trial_matlab=0;
            if sm.participant(p).session(s).trial(k).trial_type~="testA" % only for navigation trials (not testA)
                if char(trial_data.chosen_goal(k,1))=="timeout"
                    sm.participant(p).session(s).trial(k).exclude_trial_matlab=1;
                    fprintf('Session %d, trial %d marked for exclusion due to timeout.\n',s,k);
                elseif (sm.participant(p).session(s).trial(k).start_i==2 && sm.participant(p).session(s).trial(k).goal_i==11) ...
                        || (sm.participant(p).session(s).trial(k).start_i==3 && sm.participant(p).session(s).trial(k).goal_i==14) ...
                        || (sm.participant(p).session(s).trial(k).start_i==6 && sm.participant(p).session(s).trial(k).goal_i==3) ...
                        || (sm.participant(p).session(s).trial(k).start_i==7 && sm.participant(p).session(s).trial(k).goal_i==6)
                    sm.participant(p).session(s).trial(k).exclude_trial_matlab=1;
                    fprintf('Session %d, trial %d marked for exclusion due to ambiguous shortest path between start %d and goal %d.\n',...
                        s,k,sm.participant(p).session(s).trial(k).start_i,sm.participant(p).session(s).trial(k).goal_i);
                elseif (path_length<=0.1 || total_rotation==0 || sm.participant(p).session(s).trial(k).time<3)
                    sm.participant(p).session(s).trial(k).exclude_trial_matlab=1;
                    fprintf('Session %d, trial %d marked for exclusion due lack of movement or trial time < 3 sec.\n',s,k);
                end
                
                %% create trial plot               
%                 plotTrialTrack(sm.participant(p).id, sm.participant(p).session(s).session_num,...
%                     sm.participant(p).session(s).trial(k).trial, sm.participant(p).session(s).trial(k).trial_type, sm.participant(p).session(s).trial(k).memory_score,...
%                      sm.participant(p).session(s).trial(k).strategy_score_allo, sm.participant(p).session(s).trial(k).strategy_score_ego,...
%                      sm.participant(p).session(s).trial(k).time, sm.participant(p).session(s).trial(k).excess_path_length,...
%                      sm.participant(p).session(s).trial(k).initial_rotation_velocity, sm.participant(p).session(s).trial(k).dtw_to_testE,...
%                      sm.coord.full_poly, xi, yi, x_line, y_line, x_line_ego, y_line_ego, x_line_chosen, y_line_chosen,...
%                      sm.participant(p).session(s).trial(k).goal_x, sm.participant(p).session(s).trial(k).goal_y, output_folder); 
            end
            
            clear x* y* ti ri origin* goal_* ego_* start_* *chosen* name final_distance* path_length total_rotation initial_rotation;  
        end
        
        clear files trial_data input_folder;
        fprintf('Analysis done for %d, session %d.\n', id, s);
    end
    
    p=0; % reset value
    save(file_path, 'sm');
    clear dtw_dict*; 
    toc;    
end

%% Block 4: Write data to xlsx file
% [data_folder]  = setInputPath();
% writeNavTableToXLSX(data_folder); 

clear; 