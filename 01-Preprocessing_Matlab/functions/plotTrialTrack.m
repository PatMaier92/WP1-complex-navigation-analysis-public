function plotTrialTrack(id, session, trial, type, memory_score, ...
    strategy_score_allo, strategy_score_ego, time, excess_path, rotation, dtw, ...
    polyshape, x, y, x_line, y_line, x_line_ego, y_line_ego, ...
    x_line_chosen, y_line_chosen, goal_x, goal_y, folder)
% plotTrialTrack Creates track plots for each individual trial.
%
% Input: Information for creating and naming the plot
%
% Returns: A nice trial track plot.

ID=num2str(id);
Session=num2str(session);
Trial=int2str(trial);
MS=num2str(round(memory_score,2)); 
SSa=num2str(round(strategy_score_allo,2)); 
SSe=num2str(round(strategy_score_ego,2)); 
TI=num2str(round(time,1)); 
EP=num2str(round(excess_path,2)); 
RP=num2str(round(rotation,2)); 
DTW=num2str(round(dtw,2));

wfig=figure('visible','off');
plot(polyshape,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.1);
axis([0 1 0 1]); xticks(0:0.1:1); yticks(0:0.1:1); 
hold on;

viscircles([goal_x goal_y],0.01); 

if type=="learn" || type=="testE" || (type=="testA" && start==1)
    line1=plot(x,y,'k -', 'LineWidth', 1);
    line2=plot(x_line,y_line,'r -.', 'LineWidth', 1);
    line3=plot(x_line_chosen,y_line_chosen,'b .:', 'LineWidth', 1);
    legend([line1 line2 line3],{'path','ideal path','ideal path chosen'});
else % "testN" & new starts 
    line1=plot(x,y,'k -', 'LineWidth', 1);
    line2=plot(x_line,y_line,'r -.', 'LineWidth', 1);
    line3=plot(x_line_chosen,y_line_chosen,'b .:', 'LineWidth', 1);
    line4=plot(x_line_ego,y_line_ego,'g .:', 'LineWidth', 1);
    legend([line1 line2 line3 line4],{'path','ideal path','ideal path chosen','ideal path ego'});
end

if type=="learn"
    Type = 'Training'; 
elseif type=="testN"
    Type = 'Probe';
elseif type=="testE"
    Type = 'Forced egocentric';
else 
    Type = ' (XXXXX)'; 
end
title({[ID ', Session: ' Session ', Trial: ' Trial ', Condition: ' Type]; ...
    ['time: ' TI, ', ex. path: ' EP, ', memory: ' MS]; ...
    ['allo: ', SSa, ', ego: ', SSe, ', rot. vel: ' RP, ', dtw: ', DTW]});
hold off; 

% save plot
file_name = ['Plot_' ID '_' Session '_' Trial '.jpeg'];
file_path = fullfile(folder, file_name);
saveas(wfig,file_path);

end