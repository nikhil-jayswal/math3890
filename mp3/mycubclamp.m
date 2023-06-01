function [y, c] = mycubclamp(t,z, alpha, beta)
    % Nikhil Jayswal
    % Clamped Spline Interpolation
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
    % iterate over all columns
    for i = 1:dim
        % since s(t) = \Sigma c_i N_i(t)
        % N_i(t) = s(t) when c_j = 1 iff j = i
        coeff = zeros(1, dim);
        coeff(i) = 1;
        M(1:n, i) = sval2(d, y, coeff, t);  
    end  
    % boundary conditions
    % clamped boundary condition
    % left end
    for i = 1:dim
        coeff = zeros(1, dim);
        coeff(i) = 1;
        [yd, cd] = derspl(d, y, coeff);
        l = findinterval(dim-1, yd, t(1));
        M(n+1, i) = sval(d-1, yd, cd, t(1), l);
    end
    % right end
    for i = 1:dim
        coeff = zeros(1, dim);
        coeff(i) = 1;
        [yd, cd] = derspl(d, y, coeff);
        l = findinterval(dim-1, yd, t(end));
        M(n+2, i) = sval(d-1, yd, cd, t(end), l);
   end
     
    % setup z-vector
    zvector = zeros(dim, 1);
    zvector(1:n, 1) = z;
    % for debugging
    % clamped spline
    zvector(n+1) = alpha;
    zvector(n+2) = beta;
 
    
    % compute coefficient vector
    c = M\zvector;
end