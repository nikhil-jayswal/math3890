% Nikhil Jayswal
% MATH 3890
% Machine Problem 3
% 7/2/2021

clc; clear
close all

% define function to be interpolated (anonymous)
f = @(x) exp(x).*sin(2*pi*x);

% input number of sample points
% n = input('Number of sample points, n = ');

% for testing
n = 6;

% endpoints
a = 0;
b = 1;

% vector of sample-points
t = linspace(a, b, n);

% vector of data values
z = f(t);

% compute extended knot vector and coefficient vector for natural spline
[y, c] = cubnat(t, z);

% print extended knot vector and coefficient vector
fprintf('The extended knot vector is: \n\n')
disp(y')
fprintf('\n\n')
fprintf('The coefficients are: \n\n')
disp(c)

% evaluate and compute errors
t = linspace(a, b, 501);
d = 3;
% evaluate
val = sval2(d, y, c, t);
% error
e = abs(f(t) - val);

% print maxmimum and RMS error
fprintf('The maximum error is: %f\n', max(e));
fprintf('The rms error is %f\n', erms(e));

% plot spline and interpolated function
plot(t, f(t), 'g-', 'LineWidth', 3);
hold on
plot(t, val, 'b--', 'LineWidth', 2);
xlabel('t')
legend('Function f(t)', 'Spline s(t)', 'Location', 'best')

% create table similar to the one in Example 1.10
nlist = [5 9 17 33 65];
emax = zeros(size(nlist));
rms = zeros(size(nlist));

for i = 1:length(nlist)
    n = nlist(i);
    t = linspace(a, b, n);
    z = f(t); 
    [y, c] = cubnat(t, z);
    t = linspace(a, b, 501);
    d = 3;
    val = sval2(d, y, c, t);
    e = abs(f(t) - val);
    emax(i) = max(e);
    rms(i) = erms(e);
end

tbl1 = table;
tbl1.m = nlist';
tbl1.emax = emax';
tbl1.rms = rms';
fprintf('\n\n')
disp(tbl1)

ratios_emax = zeros(1, length(nlist)-1);
ratios_rms = zeros(1, length(nlist)-1);
for i = 2:length(nlist)
    ratios_emax(i-1) = emax(i-1)/emax(i);
    ratios_rms(i-1) = rms(i-1)/rms(i);
end
tbl2 = table;
tbl2.ratios_emax = ratios_emax';
tbl2.ratios_rms = ratios_rms';
fprintf('\n\n')
disp(tbl2)