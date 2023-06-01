% Nikhil Jayswal
% MATH 3890
% Machine Problem 1, Question 1
% 2/1/2021

clc; clear
close all

% endpoints
a = -5;
b = 5;

% n = # of intervals
n = input('Enter an integer (n): ');
% h = length of interval
h = (b - a)/n;
% x = vector of interpolation points
x = a:h:b;

% function to be interpolated
f = @(x) 1./(1 + x.^2);

% function values for interpolation
p = f(x);

% find coeffs. of interpolating polynomial
A = zeros(n+1, n+1);
for i = 1:(n+1)
    for j = 1:(n+1)
        A(i, j) = x(i)^(j-1);
    end
end
c = A\p'; % c = coeffs. vector

% display interpolating polynomial
fprintf('\n\n')
fprintf('The coefficients of interpolating polynomial are listed below.')
tbl = table;
term = ["constant"];
for i = 2:n+1
    s = ['x^', num2str(i-1)];
    s = string(s);
    term = [term, s];
end
tbl.term = term';
tbl.coeffs = c;
fprintf('\n\n')
disp(tbl)

% use horner() to evaluate polynomial
% N = # of evaluation points
N = 201;
h = (b-a)/(N-1);
t = a:h:b;
v = horner(c, t);

% make plots
figure()
plot(t, v, 'r--', 'LineWidth', 2)
hold on
plot(t, f(t), 'b-', 'LineWidth', 1)
xlabel('x')
ylabel('y')
legendstr = ['Interpolating polynomial of degree ', num2str(n)];
legend(legendstr, 'Function f(x)', 'Location', 'best')

% compute error
error = 0;
for i = 1:N
    e = abs(v(i) - f(t(i)));
    if e > error
        error = e;
    end
end

% annotate plot with error
titlestr = ['n = ', num2str(n), ' | max. error = ', num2str(error)];
title(titlestr)
