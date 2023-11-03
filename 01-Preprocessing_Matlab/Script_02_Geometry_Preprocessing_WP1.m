clear; clc;
addpath(genpath(pwd)) % add subfolder functions to path

%% Cognitive Map Drawing Analysis %%
% Author: Patrizia Maier

% The script requires black and white .png images as input. 

%% General settings
% define data input folder 
[input_folder] = setInputPath(); 

% load template image 
in_file = fullfile(input_folder, 'WP1_template.svg'); 
t = loadsvg(in_file, 0.5, false); 
dataT = normalize(vertcat(t{:}),1);
clear t in_file; 

% define ID range 
[id_s, id_e] = setParticipants(); 

% define output folder 
result_folder = fullfile('..', 'WP1_data', 'WP1_results'); 
if ~exist(result_folder, 'dir')
    mkdir(result_folder);
    disp('Your output folder did not exist, a new folder was created.')
end

% set up results data structure
% out_file = fullfile(result_folder, '\wp1_results_map_scoring.mat');
out_file2 = fullfile(result_folder, '\wp1_results_map_scoring.txt');
if isfile(out_file)
    load(out_file)
    p = 0; 
else
    results = [];
    p = 1;
end

%% Load data image 
for id = id_s:id_e   
    % check if ID is valid
    file = fullfile(input_folder, [num2str(id), '_geometry.svg']); 
    if ~isfile(file)
        continue;
    end
    
    % get participant ID index 
    if p ~= 1
        % check if ID already exists
        p_i = find([results.id]==id);
        if isempty(p_i) % if not: append participant data
            n = numel([results.id]);
            p = n+1; clear n p_i;
        else % if yes: overwrite participant data
            p = p_i; clear p_i;
        end
    end 
     
    % load image
    i = loadsvg(file, 0.5, false);
    dataI = normalize(vertcat(i{:}),1); 
    clear i;
    %     figure; plot(dataT(:,1), dataT(:,2), 'k+-'); hold on;
    %     plot(dataI(:,1), dataI(:,2), 'rx-');
    
    %% Use dynamic time warping for interpolation
    % get indizes
    [~,indT,indI] = dtw(dataT', dataI');
    
    % interpolate
    dataTi = dataT(indT,:); clear indT;
    dataIi = dataI(indI,:); clear indI;
    %     figure; plot(dataTi(:,1), dataTi(:,2), 'k+-'); hold on;
    %     plot(dataIi(:,1), dataIi(:,2), 'rx-');
    
    %% Compute prokrustes distance
    % a measure of dissimilarity between two shapes, returned as a numeric
    % scalar in the range [0,1].
    
    % compute prokrustes distance
    [d,Z] = procrustes(dataTi, dataIi);
     
%     figure;
%     plot(dataTi(:,1),dataTi(:,2),'k+-'); hold on;
%     plot(dataIi(:,1),dataIi(:,2),'r+-');
%     plot(Z(:,1),Z(:,2),'s-')
%     legend('Target shape (X)','Comparison shape (Y)', 'Transformed shape (Z)');
%     hold off;
        
    %% Add to results
    results(p).id = id;
    results(p).prokrustes_distance = d;
    
    p = p+1;
    
    clear id d Z dataI*; 
end

%% Save results 
% % as .mat
% save(out_file, 'results');

% as .txt
r=struct2table(results);
writetable(r, out_file2);
   
clear all; 

%% Tests for Prokrustes method: Influence of rotation, translation, scaling and interpolation 
% % load template image
% t = loadsvg('WP1_template.svg', 0.5, false); 
% dataT = normalize([t{1,1}; t{1,2}],1); clear t; 
% % figure; plot(dataT(:,1), dataT(:,2), 'k+-'); 

%% Test 1: Prokrustes distance to identical image is zero 
% procrustes(dataT, dataT) < 0.001
 
%% Test 2: Prokrustes distance to rotated and translated original image is zero 
% % define rotation matrix
% x_center = 0; y_center = 0;
% center = repmat([x_center; y_center], 1, length(dataT));
% theta = 30; 
% R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% 
% % rotate image 
% vo = R*(dataT - center')' + center;
% % figure; plot(vo(1,:), vo(2,:), 'x-'); hold on; 
% % plot(dataT(:,1), dataT(:,2), 'k+-');
% 
% % Prokrustes distance to rotated and translated original image is zero
% procrustes(dataT, [vo(1,:)' vo(2,:)']) < 0.001
% 
% clear *center theta R vo;

%% Test 3: Prokrustes distance to scaled original image is zero
% % scale image 
% dataTs = dataT * 1.5;   
% % figure; plot(dataTs(:,1), dataTs(:,2),'k+-'); hold on;
% % plot(dataT(:,1), dataT(:,2), 'k+-'); 
% 
% % Prokrustes distance to scaled original image is zero
% [d,Z] = procrustes(dataT, dataTs); 
% d < 0.001
% 
% % figure;
% % plot(dataTs(:,1), dataTs(:,2), 'k+-'); hold on;
% % plot(dataT(:,1), dataT(:,2),'k+-');
% % plot(Z(:,1),Z(:,2),'s-')
% % legend('Target shape (X)','Comparison shape (Y)', 'Transformed shape (Z)');
% % hold off;
% 
% clear dataTs d Z; 
 
%% Test 4: Prokrustes distance to interpolated original image is very close to zero
% i = loadsvg('WP1_template_ap.svg', 0.5, false); 
% dataI = normalize([i{1,1}; i{1,2}],1); clear i; 
% % figure; plot(dataI(:,1), dataI(:,2), 'k+-'); 
% 
% % DTW interpolation
% [~,indT,indI] = dtw(dataT', dataI');
% dataTi = dataT(indT,:); clear indT;
% dataIi = dataI(indI,:); clear indI; 
% % figure; plot(dataTi(:,1), dataTi(:,2), 'k-'); hold on; 
% % plot(dataIi(:,1), dataIi(:,2), 'r-');
% 
% % Prokrustes distance to scaled and interpolated original image is very close to zero
% [d,Z] = procrustes(dataTi, dataIi); 
% d < 0.001 
% 
% % figure;
% % plot(dataTi(:,1),dataTi(:,2),'k+-'); hold on;
% % plot(dataIi(:,1),dataIi(:,2),'r+-');
% % plot(Z(:,1),Z(:,2),'s-')
% % legend('Target shape (X)','Comparison shape (Y)', 'Transformed shape (Z)');
% % hold off;
% 
% clear d Z dataI* dataT*;