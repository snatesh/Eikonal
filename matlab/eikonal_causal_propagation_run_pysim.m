clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

P = load('../python/tri_points.txt');    
T = load('../python/tri_faces.txt') + 1; 
f_vals = load('../python/speed_field.txt');
DT = triangulation(T,P);

meshT = DT.ConnectivityList;
meshP = DT.Points';
nTri = size(meshT,1);

Nrs = length(W);
XX = zeros(Nrs*nTri,1);
YY = zeros(Nrs*nTri,1);
UU = zeros(Nrs*nTri,1);
for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    XX((iTri-1)*Nrs+1:Nrs*iTri) = X;
    YY((iTri-1)*Nrs+1:Nrs*iTri) = Y;
end
figure()
triplot(DT)

figure()
TRq = delaunayTriangulation(XX, YY);  % build tri mesh from quadrature points
trisurf(TRq.ConnectivityList, XX, YY, f_vals);
shading interp
view(2)   % Top-down (2D) view
axis equal
colorbar
title('Speed function f(x,y)')
drawnow; 

% compute min-time surface
m = 12; n = m+1; M = n*(n+1)/2; nquad = length(W);
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% dirchlet and neumann data/rhs
udirch = @(x,y) zeros(size(x));
rhs = 1./f_vals;
edges = freeBoundary(DT);    % E × 2 list of boundary vertex indices
pts = DT.Points;             % N × 2 matrix of vertex coordinates

% Choose a random edge
i = randi(size(edges, 1));   % random row index
v1 = edges(i,1);
v2 = edges(i,2);

% Get the coordinates of the vertices
p1 = pts(v1, :);  % 1×2
p2 = pts(v2, :);  % 1×2

% Compute a random point along the edge (not necessarily midpoint)
alpha = rand();              % scalar in [0,1]
x_anchor = (1 - alpha) * p1 + alpha * p2;  % point on edge
goal = x_anchor'

[cu_glb,V_abc,seed_tri,proxy] = eikonal_causal_propagation(DT, goal', m, R, S, W, udirch, rhs, 2);

% get quad grid and sol
Nrs = length(W);
XX = zeros(Nrs*nTri,1);
YY = zeros(Nrs*nTri,1);
UU = zeros(Nrs*nTri,1);
for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    XX((iTri-1)*Nrs+1:Nrs*iTri) = X;
    YY((iTri-1)*Nrs+1:Nrs*iTri) = Y;
    cu_abc = cu_glb((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    UU((iTri-1)*Nrs+1:iTri*Nrs) = u;
end

figure()
trisurf(TRq.ConnectivityList, XX, YY, UU);
shading interp
view(2)   % Top-down (2D) view
%axis equal
colorbar
title('Time surface u(x,y)')
drawnow; 
%%
triplot(DT)

cu = cu_glb;
x0 = [179,192]';
%x0 = [0.89, 0.92]';
tol = 1e-3;
%plot_sol(meshP,meshT,cu,V_abc,R,S,M);

step_size = 0.05;
step_tol = 1e-1;
max_steps = 1000;




% Plot contour of U over scattered points
% along with inverse velocity field
figure; hold on;
[Xg, Yg] = meshgrid(linspace(min(DT.Points(:,1)),max(DT.Points(:,1))),...
  linspace(min(DT.Points(:,2)),max(DT.Points(:,2))),200);



% Overlay contour lines
contourData = scatteredInterpolant(XX, YY, UU, 'natural', 'none');
f = scatteredInterpolant(XX,YY,rhs,'natural','none');
Zg = contourData(Xg, Yg);
contour(Xg, Yg, Zg, 20, 'k');

bbox = [min(DT.Points(:,1)), max(DT.Points(:,1));
        min(DT.Points(:,2)), max(DT.Points(:,2))];

fprintf('Bounding box: [%.2f, %.2f] x [%.2f, %.2f]\n', bbox(1,1), bbox(1,2), bbox(2,1), bbox(2,2));

H_abc = structure_factors_tri(n+1,a,b,c);

H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);

% derivative matrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
% eik_path = trace_path_greedy(x0, goal, cu, meshT, meshP, DT, ...
%                         Dx_a1bc1, Dy_ab1c1, H_a1bc1, H_ab1c1, n, ...
%                         step_size, step_tol, max_steps, H_abc);

eik_path = trace_path_greedy_proxy_pysim(cu_glb, meshT, meshP, DT, Dx_a1bc1, Dy_ab1c1, H_a1bc1, H_ab1c1, ...
                                   V_abc, proxy, f, x0, goal, ...
                                   H_abc, R, S, n);



[gradU,iTri,tx] = evaluate_grad_u(eik_path(:,end),cu,meshT,meshP,...
                            DT,Dx_a1bc1,Dy_ab1c1,H_a1bc1,H_ab1c1,n,H_abc);
z_offset = max(UU) + 1e-3;  % slight lift above the max surface
plot3(eik_path(1,:), eik_path(2,:), z_offset * ones(1,size(eik_path,2)), 'r-', 'LineWidth', 2);
plot3(eik_path(1,1), eik_path(2,1), z_offset, 'ko', 'MarkerFaceColor', 'k');  % start
plot3(eik_path(1,end), eik_path(2,end), z_offset, 'ko', 'MarkerFaceColor', 'k');  % goal