% Nikhil Jayswal
% MATH 3890
% Machine Problem 6
% 02 Mar 2021

clc; clear; close all

%% Get Triangulation Data

% get triangulation data
% n = number of points/vertices
% nt = number of triangles
% (x, y) = coordinates
% TRI = triangle-vertex matrix
[n,x,y,nt,TRI] = readtri;

%% Plot the Triangulation

triplot(TRI,x,y);

%% Create lists

[v1,v2,v3,e1,e2,e3,ie1,ie2,area] = mylists(x,y,TRI);

%% Print tables
%
tbl = table;
tbl.i = [1:nt]';
tbl.v1 = v1;
tbl.v2 = v2;
tbl.v3 = v3;
tbl.e1 = e1;
tbl.e2 = e2;
tbl.e3 = e3;
tbl.area = area;
% 
tbl2 = table;
tbl2.i = [1:length(ie1)]';
tbl2.ie1 = ie1;
tbl2.ie2 = ie2;
% 
fprintf('\n\n Triangles \n\n')
disp(tbl)
fprintf('\n\n Edges \n\n')
disp(tbl2)

%% Label edges, vertices, and triangles

% label vertices
hold on
a = [1:n]'; b = num2str(a); c = cellstr(b);
text(x, y, c, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'Red');
% label edges
a = [1:length(ie1)]'; b = num2str(a); c = cellstr(b);
text(0.5*(x(ie1) + x(ie2)), 0.5*(y(ie1) + y(ie2)), c, ....
    'FontSize', 11, 'Color', 'Blue');
% label triangles
a = [1:nt]'; b = num2str(a); c = cellstr(b);
xt = (1/3)*(x(v1) + x(v2) + x(v3));
yt = (1/3)*(y(v1) + y(v2) + y(v3));
text(xt, yt, c, 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'Green');






