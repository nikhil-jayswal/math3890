% Nikhil Jayswal
% MATH 3890
% Machine Problem 12
% 12 April 2021

clc; clear; close all

%% Get Triangulation Data

[n, xo, yo, ~, TRI] = readtri;

%% Set up lists

[nb, ne, nt, v1o, v2o, v3o, e1o, e2o, e3o, ie1o, ie2o, trilo, triro, ...
    bdy, vadj, eadj, adjstart, tadj, tstart, area, TRI] = trilists(xo, yo, TRI);

%% Set zo

zo = sigmoid(xo, yo);

%% Compute coefficients

% refine the triangulation
[x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,A] = nmdsps(xo,yo,v1o,v2o,v3o,...
    e1o,e2o,e3o,ie1o,ie2o,trilo,triro);

% coefficients (minimal energy interpolation)
[c,M22,t1,t2] = menps(v1o,v2o,v3o,x,y,zo,v1,v2,v3,e1,e2,e3,ie1,A);

%% Evaluate spline on a grid

ng = 71; d = 2;
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
[xg,yg,g] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c,ng,xmin,xmax,ymin,ymax);
figure; surfl(xg,yg,g'); colormap(copper); title('Interpolating Spline')

%% Compute the max and RMS errors

e = errg(xg,yg,g,@hill);
fprintf('emax = %5.2e, RMS = %5.2e\n',norm(e,inf),erms(e));

%% Plot derivative in northeast direction

% Evaluate the directional derivative on the grid
u = [1,1];
[xg,yg,gu] = valspdergrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c,ng,u,...
    xmin,xmax,ymin,ymax);

% Plot the directional derivative
figure; surfl(xg,yg,gu'); colormap(copper); title('Directional derivative')

%% Find and print minimal value of directional derivative on grid

fprintf('The minimum value of the directional derivative of the spline = %5.2e\n', min(min(gu)));






