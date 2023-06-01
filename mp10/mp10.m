% Nikhil Jayswal
% MATH 3890
% Machine Problem 10
% 30 Mar 2021

clc; clear; close all

%% Get Triangulation Data

[n, x, y, ~, TRI] = readtri;

%% Set up lists

[nb, ne, nt, v1, v2, v3, e1, e2, e3, ie1, ie2, tril, trir, bdy, ...
    vadj, eadj, adjstart, tadj, tstart, area, TRI] = trilists(x, y, TRI);

%% Prompt for d

% d = input('Enter the value of d: ');
d = 3;

%% Franke's function

f = @(x, y) franke2(x, y);

%% Compute coefficients

z = f(x, y);
c = scat0d(d,x,y,z,v1,v2,v3,e1,e2,e3,ie1,ie2);

%% Evaluate the spline

ng = 71;
a = min(x); b = max(x); aw = min(y); bw = max(y);
[xg, yg, g] = valspgrid(d, x, y, v1, v2, v3, e1, e2, e3, ie1, c, ng, ...
    a, b, aw, bw);

%% Plot the spline

figure; surfl(xg,yg,g'); colormap(copper);
% titlestring = ['d = ', num2str(d), ' | file = type2.', num2str(289)];
% title(titlestring)

%% Compute errors

e = errg(xg,yg,g,@franke2);
fprintf('emax = %5.2e, RMS = %5.2e\n',norm(e,inf),erms(e));