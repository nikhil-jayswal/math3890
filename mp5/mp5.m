% Nikhil Jayswal
% MATH 3890
% Machine Problem 5
% 22 Feb 2021

%%
clc; clear; close all

% interval limits
a = 0;
b = 1;

% u'' + pu' + qu = f
p = @(x) exp(x);
q = @(x) sin(x.*x);

u = @(x) cos(x) + sin(3*x);
ux = @(x) -sin(x) + 3*cos(3*x);
uxx = @(x) -cos(x) - 9*sin(3*x);

f = @(x) uxx(x) + p(x).*ux(x) + q(x).*u(x);

% get inputs from user
% d = input('Enter the value of d: ');
% k = input('Enter the value of k: ');
% m = input('Enter the value of m: ');
d = 3;
k = 17;
m = 10;

% knot vector
knots = linspace(a, b, k+2);

%% compute the extended knot vector
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

%% compute coefficient vector for spline that approximates u
ua = u(a); ub = u(b);
c = lsc(d, y, m, p, q, f, ua, ub);

%% print extended knot vector and coefficient vector
fprintf('The extended knot vector is: \n\n')
disp(y')
fprintf('\n\n')
fprintf('The coefficients are: \n\n')
disp(c)

%% compute errors
t = linspace(a, b, 301);
% evaluate
val = sval2(d, y, c, t);
% print maxmimum norm of (u-s)
e = u(t)-val;
fprintf('\nThe maximum error is: %f\n', norm(e, inf));
% print RMS error
fprintf('\nThe RMS error is: %f\n', erms(e));

%% plot s
plot(t, val, 'LineWidth', 1);
hold on
plot(t, u(t), 'LineWidth', 2);
xlabel('t')
legend('Spline s', 'Function u', 'Location', 'best')