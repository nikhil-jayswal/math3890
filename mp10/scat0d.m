function c = scat0d(d,x,y,z,v1,v2,v3,e1,e2,e3,ie1,ie2)
    % (x, y) = sample points
    % z = f(x, y) = values to be interpolated
    
    % number of triangles
    nt = length(v1);
    
    % # of coefficients = nv + (d - 1)ne + (d-1)*(d-2)*nt/2
    % # of vertices = length(x/y)
    nv = length(x);
    % # of edges = length(ie1/ie2)
    ne = length(ie1);
    % # of triangles = length(v1/v2/v3/e1/e2/e3)
    nc = nv + (d-1)*ne + nchoosek(d-1, 2)*nt;
    c = zeros(nc, 1);
    
    % for each triangle T
    for i = 1:nt
        
        % get indices of all coefficients associated with T
        list = getindex(d, i, nv, v1, v2, v3, e1, e2, e3, ie1);
        
        % coefficient at each vertex =  z(vertex)
        % get coefficients of vertices = first 3 smallest values in list
        % since vertices are stored first according to convention
        cT = zeros(size(list));
        [clist, index] = mink(list, 3);
        cT(index) = z(clist);
        
        % find 'm' points closest to barycenter of T
        % barycenter
        bTx = (x(v1(i)) + x(v2(i)) + x(v3(i)))/3;
        bTy = (y(v1(i)) + y(v2(i)) + y(v3(i)))/3;
        bT = [bTx, bTy];
        % m = 3*n_d, n_d = dimension of degree-d polynomials
        nd = nchoosek(d+2, 2); m = 3*nd; 
        % find closest points
        [vlist, ~] = knnsearch([x, y], bT, 'K', m); 
        % drop first three points = vertices of T
        vlist = vlist';
        vlist = vlist(4:end, 1);
        
        % for each point in vlist, write observation equations
        BT = zeros(length(vlist), nd);
        zT = zeros(size(vlist));
        for j = 1:length(vlist)
            % get barycentric coordinates
            x1 = x(v1(i)); y1 = y(v1(i));
            x2 = x(v2(i)); y2 = y(v2(i));
            x3 = x(v3(i)); y3 = y(v3(i));
            tx = x(vlist(j)); ty = y(vlist(j));
            [b1,b2,b3] = bcoord(x1,y1,x2,y2,x3,y3,tx,ty);
            % get matrix of basis functions
            BT(j, :) = basisv(d, b1, b2, b3);
            % fill in entries of zT
            zT(j) = z(vlist(j));
        end
        % compute coefficients
        zT = zT - BT(:, index)*z(clist);
        BT(:, index) = [];
        cT(cT == 0) = BT\zT;
        c(list) = cT;
    end
        
end