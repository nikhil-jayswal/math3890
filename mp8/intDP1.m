function c = intDP1(d,x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,f,E,lambda)
    % Nikhil Jayswal
    % MATH 3890
    % Machine Problem 8
    % 16 April 2021

    % number of triangles
    nt = length(v1);

    % # of coefficients = nv + (d - 1)ne + (d-1)*(d-2)*nt/2
    % # of vertices = length(x/y)
    nv = length(x);
    % # of edges = length(ie1/ie2)
    ne = length(ie1);
    % # of triangles = length(v1/v2/v3/e1/e2/e3)
    nc = nv + (d-1)*ne + nchoosek(d-1, 2)*nt;

    % length of z = total # of domain points + # of rows of E
    ncond = size(E, 1);
    ntotal = nc + ncond;
    z = zeros(ntotal, 1);

    % get barycentric coordinates
    barycoords = domT(d);
    numDP = size(barycoords, 1);

    % setup system of observation equations at each domain point
    % initialise containers
    % Bz = [B z]
    Bz = zeros(nt*numDP, nc+1);
    index = 1;
    % find values of basis functions at domain points in a triangle
    B = basisv(d, barycoords(:, 1), barycoords(:, 2), barycoords(:, 3));
    for k = 1:nt
        % get indices of all coefficients associated with T_k
        clist = getindex(d, k, nv, v1, v2, v3, e1, e2, e3, ie1);
        for i = 1:size(B, 1)
            Bz(index, clist) = B(i, :);
            % last column = z-value
            a = [x(v1(k)) y(v1(k))];
            b = [x(v2(k)) y(v2(k))];
            c = [x(v3(k)) y(v3(k))];
            xy = a*barycoords(i, 1) + b*barycoords(i, 2) + ...
                        c*barycoords(i, 3);
            Bz(index, end) = f(xy(1), xy(2));
            index = index + 1;
        end
    end
    % remove redundant rows
    % unique() won't work (floats!)
    B = log10(10.^Bz);
    Bz = uniquetol([Bz;B], 'Byrows', true);
    z(1:nc) = Bz(:, end);
    B = Bz(:, 1:end-1);

    % lambda E = 0 conditions
    z(nc+1:end) = 0;

    B = [B; lambda*E];

    % solve for [c]
    c = B\z;
end
