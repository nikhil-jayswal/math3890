%
% Nikhil Jayswal
% MATH 3890
% Machine Problem 9
% 24 Mar 2021
%

clc; clear; close all

%% Read Triangulation from a file

[n, x, y, nt, TRI] = readtri;

%% Call refinect() and print [x, y]

[xr, yr, TRIr] = refinect(x, y, TRI);
fprintf('The [x, y] pairs (vertex coordinates) are listed below:\n');
for i = 1:length(xr)
    fprintf('[%f, %f]\n', xr(i), yr(i));
end
fprintf('\n\n')

%% Plot refined triangulation

figure(1)
triplot(TRIr, xr, yr);
title('CT Refinement')

%% Call refinecti and print [x, y]

[xr, yr, TRIr] = refinecti(x, y, TRI);
fprintf('The [x, y] pairs (vertex coordinates) are listed below:\n');
for i = 1:length(xr)
    fprintf('[%f, %f]\n', xr(i), yr(i));
end
fprintf('\n\n')

%% Plot refined triangulation

figure(2)
triplot(TRIr, xr, yr);
title('Refinement using Incenters')
