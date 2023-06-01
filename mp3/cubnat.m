function [y, c] = cubnat(t,z)
    % Nikhil Jayswal
    % MATH 3890
    % Machine Problem 3, Question 1
    % 7/2/2021
    
    % z = [z1 z2 ... zn] = data values
    % t = [a=t1 t2 ... tn=b] = points
    % y = extened knot vector
    % c = coefficient vector
    
    % degree of polynomial = 3, i.e. cubic spline
    d = 3; 
    % number of sample points/knots
    n = length(t);
    % find dimension of spline space
    dim = n+2;
    
    % construct extended knot vector
    y = zeros(1, dim+d+1);
    % first (d+1) points
    y(1:d+1) = t(1); 
    % last (d+1) points
    y(dim+1:dim+d+1) = t(end);
    % points in between
    y(d+2:dim) = t(2:end-1);
    
    % construct basis functions (B-splines)
    % construct matrix M
    M = zeros(dim, dim);
    % value constraints
    % first 'n' rows
    for i = 1:n
%         % since s(t) = \Sigma c_i N_i(t)
%         % N_i(t) = s(t) when c_j = 1 iff j = i
%         coeff = zeros(1, dim);
%         coeff(i) = 1;
%         M(1:n, i) = sval2(d, y, coeff, t); 
        % find interval containing t(i)
        l = findinterval(dim, y, t(i));
        % get all b-splines with non-zero value at the point t(i)
        b = bspl(d, l, t(i), y);
        % assemble row #i of [M]
        for j = 1:dim
            if j < (l-d)
                M(i, j) = 0;
            end
            if j > l
                M(i, j) = 0;
            end
            M(i, l-d:l) = b;
        end
    end   
    % natural spline, double derivative at end points = 0
    % left end
% WRONG (RMS ERRORS ARE OFF!)    
%    for i = 1:dim
%         coeff = zeros(1, dim);
%         coeff(i) = 1;
%         [yd, cd] = derspl(d, y, coeff);
%         [yd2, cd2] = derspl(d-1, yd, cd);
%         l = findinterval(dim-2, yd, t(1));
%         M(n+1, i) = sval(d-2, yd2, cd2, t(1), l);
%    end
    l = findinterval(dim, y, t(1));
    % get d2/dt2 of all b-splines with non-zero value at 'a'
    b = bsplder(d, y, l, t(1), 2);
    % assemble #(n+1) row
    for j = 1:dim
        if j < (l-d)
            M(n+1, j) = 0;
        end
        if j > l
            M(n+1, j) = 0;
        end
    end
    M(n+1, l-d:l) = b;

    % right end
% WRONG (RMS ERRORS ARE OFF!)
%     for i = 1:dim
%         coeff = zeros(1, dim);
%         coeff(i) = 1;
%         [yd, cd] = derspl(d, y, coeff);
%         [yd2, cd2] = derspl(d-1, yd, cd);
%         l = findinterval(dim-2, yd, t(end));
%         M(n+2, i) = sval(d-2, yd2, cd2, t(end), l);
%     end
    l = findinterval(dim, y, t(end));
    % get d2/dt2 of all b-splines with non-zero value at 'b'
    b = bsplder(d, y, l, t(end), 2);
    % assemble row #(n+2) of [M]
    for j = 1:dim
        if j < (l-d)
            M(n+2, j) = 0;
        end
        if j > l
            M(n+2, j) = 0;
        end
    end
    M(n+2, l-d:l) = b;

        
    % setup z-vector
    zvector = zeros(dim, 1);
    zvector(1:n, 1) = z;
    % natural boundary conditions
    zvector(n+1) = 0;
    zvector(n+2) = 0;
    
    % compute coefficient vector
    c = M\zvector;
end