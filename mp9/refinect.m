function [x, y, TRI] = refinect(xo, yo, TRIo)

    % create lists
    [v1,v2,v3,~,~,~,~,~,~] = mylists(xo ,yo ,TRIo);
    
    % number of triangles
    nt = size(TRIo, 1);
    
    % find barycenters of triangles in TRIo
    bx = (xo(v1) + xo(v2) + xo(v3))/3;
    by = (yo(v1) + yo(v2) + yo(v3))/3;
    
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