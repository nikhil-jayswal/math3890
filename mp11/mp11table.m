% Nikhil Jayswal
% MATH 3890
% Machine Problem 11
% 06 April 2021

clc; clear; close all

%% Read Grid Points

[npts, x, y] = readxy;

%% Franke's function

f = @(x, y) franke2(x, y);
z = f(x, y);

%% Radial basic function
rbf = @(eps, r) exp(-(eps*r).^2);

%% eps values

epslist = 1:7;

%% Construct table

tbl = table;
tbl.eps = epslist';

maxerrors = zeros(size(epslist));
rmserrors = zeros(size(epslist));
cnum = zeros(size(epslist));

for m = 1:length(epslist)
    
    eps = epslist(m);
    % Compute coefficients

    [c, M] = scatrbf(x, y, z, eps, rbf);

    % Compute difference between RBF interpolant and f

    % Evaluate the RBF interpolant and Franke's function on a grid
    ng = 71;  
    xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);
    xg = linspace(xmin,xmax,ng); yg = linspace(ymin,ymax,ng);
    interp_value = zeros(ng, ng);
    exact_value = zeros(ng, ng);
    for i = 1:ng
        for j = 1:ng
            interp_value(i, j) = 0;
            for k = 1:length(c)
                r = sqrt((xg(i) - x(k))^2 + (yg(j) - y(k))^2);
                interp_value(i, j) = interp_value(i, j) + c(k)*rbf(eps, r);
                exact_value(i, j) = franke2(xg(i), yg(j));
            end
        end
    end

    % Compute difference
    err = exact_value - interp_value;
    err = reshape(err, ng*ng, 1);

    % compute errors and condition number
    maxerrors(m) = norm(err, inf);
    rmserrors(m) = erms(err);
    cnum(m) = cond(M);
    
end

tbl.Max_Error = maxerrors';
tbl.RMS_Error = rmserrors';
tbl.Condition_Number = cnum';

disp(tbl);