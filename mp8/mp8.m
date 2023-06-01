% Nikhil Jayswal
% MATH 3890
% Machine Problem 8
% 16 April 2021

clc; clear; close all

%% Get Triangulation Data

[n, x, y, nt, TRI] = readtri;

%% Set up lists

[nb, ne, nt, v1, v2, v3, e1, e2, e3, ie1, ie2, tril, trir, bdy, ...
    vadj, eadj, adjstart, tadj, tstart, area, TRI] = trilists(x, y, TRI);

%% Prompt for d and lambda

d = input('Enter the value of d: ');
lambda = input('Enter the value of lambda: ');

%% Franke's function

f = @(x, y) franke2(x, y);

%% Setup matrix E

E = c1smooth(d,x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir);

%% Compute coefficients

tic
c = intDP1(d,x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,f,E,lambda);
tcomp = toc;

%% Evaluate the spline

ng = 51;
a = min(x); b = max(x); aw = min(y); bw = max(y);
[xg, yg, g] = valspgrid(d, x, y, v1, v2, v3, e1, e2, e3, ie1, c, ng, ...
    a, b, aw, bw);

%% Plot the spline

figure; surfl(xg,yg,g'); colormap(copper);
titlestring = ['Interpolating Spline, d = ', num2str(d), ...
                ', \lambda = ', num2str(lambda)];
title(titlestring)

%% Compute errors

e = errg(xg,yg,g,@franke2);
fprintf('emax =%5.2e, RMS = %5.2e\n',norm(e,inf),erms(e));

%% Report computation time and c1ck value

fprintf('Computation time = %f seconds\n', tcomp);
c1ck(d,x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,c);





