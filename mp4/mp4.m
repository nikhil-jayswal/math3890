% Nikhil Jayswal
% MATH 3890
% Machine Problem 4
% 15/2/2021

clc; clear
close all

% anonymous testing function
f = @(x) exp(x).*sin(2*pi*x);

% interval of interest
a = -1;
b = 1;

% get inputs from user
% N = input('Enter the value of N: ');
% eps = input('Enter the value of epsilon: ');
% d = input('Enter the degree of spline: ');
% k = input('Enter the number of knots: ');
N = 201;
d = 3;
k = 15;
eps = 0;

% knot vector
knots = linspace(a, b, k+2);
% compute the extended knot vector
% dimension of spline space
% length of knot vector = k+2 -> spline consists of (k+1) pieces
% (k) interior points
dim = (k+1)*(d+1) - k*d;

% construct extended knot vector
y = zeros(1, dim+d+1);
% first (d+1) points
y(1:d+1) = knots(1); 
% last (d+1) points
y(dim+1:dim+d+1) = knots(end);
% points in between
y(d+2:dim) = knots(2:end-1);

% vector of sample points
t = linspace(a, b, N);
% w = noise uniformly distributed in [-eps, eps]
w = eps*(-1 + (1+1)*rand(N,1));
% vector of data values
z = f(t) + w';
 
% compute coefficient vector
c = lsqsplo(d, y, t, z);

% print extended knot vector and coefficient vector
fprintf('The extended knot vector is: \n\n')
disp(y')
fprintf('\n\n')
fprintf('The coefficients are: \n\n')
disp(c)
 
% evaluate and compute errors
t = linspace(a, b, 501);
% evaluate
val = sval2(d, y, c, t);

% print maxmimum norm of (f-s)
e = f(t)-val;
fprintf('\nThe maximum error is: %f\n', norm(e, inf));

% plot f and s
plot(t, f(t), 'g-', 'LineWidth', 3);
hold on
plot(t, val, 'b--', 'LineWidth', 2);
xlabel('t')
legend('Function f(t)', 'Spline s(t)', 'Location', 'best')
titlestring = ['\epsilon', ' = ', num2str(eps), ' | ', ...
    'max error = ', num2str(norm(e, inf))];
title(titlestring)
