function [goal_i, zone_s, zone_i, alley_s, alley_i, obj_at_location]=computeChosenGoals(rand_dict,...
    goal_s, alley_poly_out, alley_poly_in, rec_poly, tri_poly,...
    zone_names, alley_names, goal_names, xend, yend)
% computeChosenGoals: Returns chosen goal location as integer and string.
% The goal string in trial_data is error-prone for the triangle
% intersections and thus recalculated by this function. 
%
% Input: 
% rand_dict with randomization info (dictionary)
% goal_s [only relevant for "timeout", else recalculated] (string) 
% alley_poly, rec_poly, tri_poly (polyshapes)
% alley_names, goal_names (string)
% xend, yend are chosen final x-/y-coordinates (float)
%
% Returns: 
% goal_i (integer), alley_s (string), alley_i (integer), obj_at_location (string) 

if goal_s == "timeout"
    zone_s = "timeout"; 
    alley_s = "timeout"; 
    zone_i = 999;
    alley_i = 999; 
    goal_i = 999; 
    obj_at_location = "999";  
else 
    % recalculate goal location 
    for c=1:size(alley_poly_out,2)
        % chosen location in alley (outer half of outer arm)
        if inpolygon(xend, yend, alley_poly_out{c}.Vertices(:,1), alley_poly_out{c}.Vertices(:,2))
            zone_i = c*4-3; % [ 1 5 9 13 17 21 25 ]
            break; 
        % chosen location in alley (inner half of outer arm)
        elseif inpolygon(xend, yend, alley_poly_in{c}.Vertices(:,1), alley_poly_in{c}.Vertices(:,2))
            zone_i = c*4-2; % [ 2 6 10 14 18 22 26 ]
            break;   
        % chosen location in triangle (intersection)
        elseif inpolygon(xend, yend, tri_poly{c}.Vertices(:,1), tri_poly{c}.Vertices(:,2))
            zone_i = c*4-1; % [ 3 7 11 15 19 23 27 ]
            break;
        % chosen location in rectangle (inner arm)
        elseif inpolygon(xend, yend, rec_poly{c}.Vertices(:,1), rec_poly{c}.Vertices(:,2)) 
            zone_i = c*4; % [ 4 8 12 16 20 24 28 ]
            break;
        end
    end
    
    % chosen zone & alley string
    zone_s = zone_names(zone_i);
    alley_s = extractBefore(zone_s, digitsPattern); 
    alley_i = find(contains(alley_names, alley_s)); 
    if isempty(alley_i)
        alley_i = 999;
    end
    
    % chosen goal string
    goal_i = find(contains(goal_names, zone_s));
    if isempty(goal_i)
        goal_i = 999;
    end
    
    % correct object at chosen location
    obj_at_location="999";
    fields = fieldnames(rand_dict);
    for i=1:length(fields)
        if zone_s == string(fields{i})
            key = fields{i};
            obj_at_location = string(rand_dict.(key).object);
        end
    end
end 

end
