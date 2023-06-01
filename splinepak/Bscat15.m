% Bscat15 
% SplinePAK: Copyright Larry Schumaker 2014
% Interpolate scattered data with a C1 quintic Argyris spline
%  Uses derivatives estimated by local polynomial least-squares

% Read in a triangulation
[n,x,y,nt,TRI] = readtri;
xmin = min(x); xmax = max(x); ymin = min(y); ymax = max(y);

% Calculate the triangulation lists
[nb,ne,nt,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,bdy,...
   vadj,eadj,adjstart,tadj,tstart,area,TRI] = trilists(x,y,TRI);

% Sample a function at the vertices
f = @(x,y) franke2(x,y);
z = f(x,y);

% Compute the coefs of the interpolating spline
tic
[der,ze] = derest15(x,y,z,v1,v2,v3,e1,e2,e3,ie1,ie2,adjstart,vadj);
c =  arg15(x,y,v1,v2,v3,e1,e2,e3,ie1,der,ze);
toc

% Check the C1 smoothness
d = 5;
c1ck(d,x,y,v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir,c)

% Evaluate the spline on a grid
ng = 51;
[xg,yg,g] =  valspgrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c,ng,xmin,xmax,ymin,ymax);

% Plot the spline
figure; surfl(xg,yg,g');  colormap(copper);

% Compute the max and RMS errors
e = errg(xg,yg,g,@franke2);
fprintf('Error:   emax =%5.2e, RMS = %5.2e\n',norm(e,inf),erms(e));

% Evaluate the x-derivative on the grid
u = [1,0];
[xg,yg,g] =  valspdergrid(d,x,y,v1,v2,v3,e1,e2,e3,ie1,c,ng,u,...
    xmin,xmax,ymin,ymax);

% Plot the x-derivative
figure; surfl(xg,yg,g');  colormap(copper);

% Calculate the error of the x-derivative
eder = errgder(xg,yg,g,@franke2d,2);
fprintf('Derivative: emax =%5.2e, RMS = %5.2e\n',norm(eder,inf),erms(eder));
