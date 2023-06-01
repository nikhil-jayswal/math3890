function c = lsqsplo(d, y, t, z)
    % Nikhil Jayswal
    % MATH 3890
    % Machine Problem 4
    % 15/2/2021
    
    % z = [z1 z2 ... zN] = data values
    % t = [t1 t2 ... tN] = sample points
    % y = extended knot vector
    % c = coefficient vector 
    % d = degree of spline
    
    % observation equations
    % z1 = c1*N1(t1) + c2*N2(t2) + ... 
    % ...
    % zN = c1*N1(tN) + c2*N2(tN) + ...
    % [z] = [N][c]
    
    % number of sample points
    N = length(t);
    % dimension of spline space
    dim = length(y) - d - 1;
    
    % assemble the matrix N
    % # of rows = # of sample points = N
    % # of cols = # of B-splines = dim
    Nmatrix = zeros(N, dim);
    for i = 1:N
        % find interval containing t(i)
        l = findinterval(dim, y, t(i));
        % get all b-splines with non-zero value at the point t(i)
        b = bspl(d, l, t(i), y);
        % assemble row #i of [N]
        for j = 1:dim
            if j < (l-d)
                Nmatrix(i, j) = 0;
            end
            if j > l
                Nmatrix(i, j) = 0;
            end
            Nmatrix(i, l-d:l) = b;
        end
    end
    % compute coefficient vector    
    c = Nmatrix\z';
end