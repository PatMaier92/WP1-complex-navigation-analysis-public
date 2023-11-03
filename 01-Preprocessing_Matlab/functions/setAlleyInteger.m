function [alley_int]=setAlleyInteger(loc_x,loc_y,alley_poly,rec_poly)
% setAlleyInteger: Determine alley integers based on x-/y-coordinates
% of goal & start locations and polyshapes. 
%
% Input: 
% loc_x, loc_y are x-/y-coordinates of goal or start locations (float vector)
% alley_poly, rec_poly are alley and rectangle polyshapes (array of polyshapes)
%
% Returns: alley_int (integer vector) 

int=0;
alley_int=zeros(1,numel(loc_x));

for i=1:numel(loc_x)
    for j=1:numel(alley_poly)
        if inpolygon(loc_x(i),loc_y(i),...
                alley_poly{j}.Vertices(:,1), alley_poly{j}.Vertices(:,2))
            int=j*2-1;
            break; 
        end
    end
    if int==0
        for j=1:numel(rec_poly)
            if inpolygon(loc_x(i),loc_y(i),...
                    rec_poly{j}.Vertices(:,1), rec_poly{j}.Vertices(:,2))
                int=j*2;
                break; 
            end
        end
    end
    alley_int(i)=int;
    int=0;
end
    
end