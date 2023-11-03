function [G,x,y]=setGraph(start_x, start_y, tri_x, tri_y, goal_x, goal_y)
% setGraph: Define graph between all locations in Starmaze WP1. 
% 
% Input: start_x,start_y,tri_x,tri_y,goal_x,goal_y are vectors of x-/y-coordinates
%
% Returns: G (graph), x and y are x-/y-coordinate vectors (float).

% start nodes
st=[ 1 1 ... % start connections (to corners)
   2 2 ... % 2 2 ... % optional: to goals in same area
   3 3 ... % 3 3 ... 
   4 4 ... % 4 4 ...
   5 5 ... % 5 5 ...
   6 6 ... % 6 6 ...
   7 7 ... % 7 7 ...
   29 29 ... % goal connections (to corners)
   30 30 ...
   31 31 31 31 ... 
   32 32 ...
   33 33 ...
   34 34 34 34 ...
   35 35 ...
   36 36 ...
   37 37 ...
   38 38 ... 
   39 39 39 39 ...
   40 40 ...
   41 41 ...
   42 42 42 42 ...
   43 43 ...
   44 44 ...
   8 11 14 17 20 23 26 ... % outer line (between corners)
   9 12 15 18 21 24 27 ... % inner line (between corners)
   8 9 ... % cross-connections (between corners)
   11 12 ...
   14 15 ...
   17 18 ...
   20 21 ...
   23 24 ...
   26 27 ]; 

tar=[ 8 10 ... % start connections (to corners) 
    11 13 ... % 29 30 ... % optional: to goals in same area
    14 16 ... % 32 33 ...
    17 19 ... % 35 36 ...
    20 22 ... % 37 38 ...
    23 25 ... % 40 41 ...
    26 28 ... % 43 44 ...
    11 13 ... % goal connections (to corners)
    11 13 ...
    16 15 12 11 ...
    14 16 ...
    14 16 ...
    19 18 15 14 ...
    17 19 ...
    17 19 ...
    20 22 ...
    20 22 ...
    25 24 21 20 ...
    23 25 ...
    23 25 ...
    28 27 24 23 ...
    26 28 ... 
    26 28 ... 
    13 16 19 22 25 28 10 ... % outer line (between corners)
    12 15 18 21 24 27 9 ... % inner line (between corners)
    12 13 ... % cross-connections (between corners)
    15 16 ...
    18 19 ...
    21 22 ...
    24 25 ...
    27 28 ...
    9 10 ]; 

% node names 
nodenames={ 's1(A)','s2(C)','s3(E)','s4(G)','s5(I)','s6(K)', 's7(M)', ...
    'tri_1_r','tri_1_b','tri_1_l','tri_2_r','tri_2_b','tri_2_l','tri_3_r','tri_3_b','tri_3_l',...
    'tri_4_r','tri_4_b','tri_4_l','tri_5_r','tri_5_b','tri_5_l','tri_6_r','tri_6_b','tri_6_l','tri_7_r','tri_7_b','tri_7_l',...
    'g1(C)','g2(C)','g3(D)','g4(E)','g5(E)','g6(F)','g7(G)','g8(G)',...
    'g9(I)','g10(I)','g11(J)','g12(K)','g13(K)','g14(L)','g15(M)','g16(M)' };

% create graph
G=graph(st, tar);

% graph coordinates (must be in correct numerical order) 
x = [transpose(start_x) ...
    tri_x(2,1) tri_x(3,1) tri_x(4,1) tri_x(2,2) tri_x(3,2) tri_x(4,2) ...
    tri_x(2,3) tri_x(3,3) tri_x(4,3) tri_x(2,4) tri_x(3,4) tri_x(4,4) ...
    tri_x(2,5) tri_x(3,5) tri_x(4,5) tri_x(2,6) tri_x(3,6) tri_x(4,6) tri_x(2,7) tri_x(3,7) tri_x(4,7) ...
    transpose(goal_x)];
y = [transpose(start_y) ...
    tri_y(2,1) tri_y(3,1) tri_y(4,1) tri_y(2,2) tri_y(3,2) tri_y(4,2) ...
    tri_y(2,3) tri_y(3,3) tri_y(4,3) tri_y(2,4) tri_y(3,4) tri_y(4,4) ...
    tri_y(2,5) tri_y(3,5) tri_y(4,5) tri_y(2,6) tri_y(3,6) tri_y(4,6) tri_y(2,7) tri_y(3,7) tri_y(4,7) ...
    transpose(goal_y)];

% calculate distance for edge weights  
[sn,tn] = findedge(G);
dx = x(sn) - x(tn);
dy = y(sn) - y(tn);
D = hypot(dx,dy);

% add edge weights 
G.Edges.Weight = D';

% add names 
G.Nodes.Names=nodenames'; % access with e.g., G.Nodes.Names{1} = 's1'

% % plot 
% figure; 
% p = plot(G,'XData',x,'YData',y);
% % p = plot(G,'XData',x,'YData',y,'EdgeLabel',G.Edges.Weight);
% labelnode(p,[1:numel(nodenames)],nodenames);
% xlim([0 1]);
% ylim([0 1]);

end