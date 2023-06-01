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

%% Prompt for eps

eps = input('Enter the value of eps: ');

%% Compute coefficients

[c, M] = scatrbf(x, y, z, eps, rbf);

%% Compute difference between RBF interpolant and f

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

%% Plot the interpolant

figure; surfl(xg',yg',interp_value'); colormap(copper);

%% Print max and RMS errors

fprintf('Maximum error: %5.2e,   RMS = %5.2e\n', norm(err,inf),erms(err))

%% Print condition number of M

fprintf('Condition number of M = %g\n', cond(M))





