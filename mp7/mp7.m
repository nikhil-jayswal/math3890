% Nikhil Jayswal
% MATH 3890
% Machine Problem 7
% 16 April 2021

clc; clear; close all

%% Read Triangulation

[nv, x, y, nt, TRI] = readtri;

%% Set up lists

[nb,ne,nt,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,bdy,vadj,eadj, ...
    adjstart,tadj,tstart,area,TRI] = trilists(x,y,TRI);

%% Prompt for d

d = input('Enter the value of d: ');

%% Function to be interpolated

f = @(x, y) franke2(x, y);

%% Compute coefficients

c = intDP(d, x, y, v1, v2, v3, e1, e2, e3, ie1, ie2, f);

%% Evaluate spline

ng = 51;
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
[xg,yg,g] = valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c,ng,xmin,xmax,ymin,ymax);

%% Plot spline

figure; surfl(xg,yg,g'); colormap(copper);
titlestring = ['Interpolating Spline, d = ', num2str(d)];
title(titlestring)

%% Compute & print max and rms errors

e = errg(xg,yg,g,f);
fprintf('emax = %5.2e, RMS = %5.2e\n',norm(e,inf),erms(e));



