function [x, y, TRI] = refinecti(xo, yo, TRIo)

    % create lists
    [v1,v2,v3,e1,e2,e3,ie1,ie2,~] = mylists(xo ,yo ,TRIo);
    
    % number of triangles
    nt = size(TRIo, 1);
    
    % find incenters of triangles in TRIo
    bx = zeros(nt, 1);
    by = zeros(nt, 1);
    for i = 1:nt
        a = sqrt((xo(v3(i)) - xo(v2(i)))^2 + (yo(v3(i)) - yo(v2(i)))^2);
        b = sqrt((xo(v1(i)) - xo(v3(i)))^2 + (yo(v1(i)) - yo(v3(i)))^2);
        c = sqrt((xo(v2(i)) - xo(v1(i)))^2 + (yo(v2(i)) - yo(v1(i)))^2);
        bx(i) = (a*xo(v1(i)) + b*xo(v2(i)) + c*xo(v3(i)))/(a+b+c);
        by(i) = (a*yo(v1(i)) + b*yo(v2(i)) + c*yo(v3(i)))/(a+b+c);
    end
    
    % create x and y vectors
    x = [xo; bx];
    y = [yo; by];
    
    % construct refined connection matrix
    % number of vertices
    nv = max(TRIo, [], 'all');
    nt = 3*nt;
    TRI = zeros(nt, 3);
    itr = 1;
    i = 1;
    while(i <= nt)
        % get vertices of triangle in unrefined TRI 
        v = TRIo(itr, :);
        % create 3 new triangles
        TRI(i, :) = [v(1) v(2) nv+1];
        TRI(i+1, :) = [v(2) v(3) nv+1];
        TRI(i+2, :) = [v(3) v(1) nv+1];
        nv = nv + 1;
        itr = itr + 1;
        i = i+3;
    end
end