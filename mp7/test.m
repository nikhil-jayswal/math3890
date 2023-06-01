% list all domain points of a triangulation

clc; clear; close all

% read a triangulation
[nv, x, y, nt, TRI] = readtri;

% setup lists
[nb,ne,nt,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,bdy,vadj,eadj, ...
    adjstart,tadj,tstart,area,TRI] = trilists(x,y,TRI);

% degree of spline S(d, 0)
d = 3;

% expected number of unique domain points
nc = length(x) + (d-1)*length(ie1) + nchoosek(d-1, 2)*length(v1);

% get barycentric coordinates
barycoords = domT(d);
numDP = size(barycoords, 1);

% domain points of the triangulation
DP = zeros(nt*numDP, 2);
index = 1;

for i = 1:nt
    a = [x(v1(i)) y(v1(i))];
    b = [x(v2(i)) y(v2(i))];
    c = [x(v3(i)) y(v3(i))];
    for j = 1:size(barycoords, 1)
        DP(index, :) = a*barycoords(j, 1) + b*barycoords(j, 2) + ...
                    c*barycoords(j, 3);
        index = index + 1;
    end
end

% remove redundant points
% unique() won't work (floats!)
B = log10(10.^DP);
DPunique = uniquetol([DP;B], 'Byrows', true);

% number of domain points == expected number of domain points ?
size(DPunique, 1) == nc
