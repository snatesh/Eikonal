clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");
DT = delaunay_unit_square(20);
meshT = DT.ConnectivityList;
meshP = DT.Points';
nTri = size(meshT,1);

[centers, r0s, goal, x0] = generate_maze_circles(50, [-1 1 -1 1], [0.05 0.1], 0.02);

eps_obst = 1e-6;  % minimum speed inside obstacle

%[f_handle, obstacles, goal, x0] = generate_rectcirc_maze_speed(11, eps_obst);
%f = f_handle;
%[f_handle, obstacles, goal, x0] = generate_maze_speed(11, eps_obst);
%f = f_handle;

f = @(x,y) multi_obstacle_speed_varying_radii(x,y,centers,r0s,eps_obst);
%[f_handle, obstacles, goal, x0] = generate_maze_squares(10, 0.05, 0.15, eps_obst);
%f = f_handle;

[xg, yg] = meshgrid(linspace(-1,1,200));
fg = f(xg, yg);
figure();
contourf(xg, yg, fg); hold on;
%contourf(xg,yg,reshape(fg,200,200)); hold on;
drawnow;
% compute min-time surface
m = 10; n = m+1; M = n*(n+1)/2; nquad = length(W);
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% dirchlet and neumann data/rhs
udirch = @(x,y) zeros(size(x));
rhs = @(x,y) 1./f(x,y);
goal = [1;0.05];
[cu_glb,V_abc,seed_tri,proxy] = eikonal_causal_propagation(DT, goal', m, R, S, W, udirch, rhs, 1);

% get quad grid and sol
Nrs = length(W);
XX = zeros(Nrs*nTri,1);
YY = zeros(Nrs*nTri,1);
UU = zeros(Nrs*nTri,1);
[Xg, Yg] = meshgrid(linspace(-1, 1, 500));
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
%%
% Compute triangle midpoints
tris = DT.ConnectivityList;     % Nx3 triangle vertex indices
pts  = DT.Points;               % Nx2 point coordinates

midpts = (pts(tris(:,1),:) + pts(tris(:,2),:) + pts(tris(:,3),:)) / 3;

% 3D scatter plot
figure;
scatter3(midpts(:,1), midpts(:,2), proxy, 50, proxy, 'filled');
colorbar;
title('3D Proxy Values at Triangle Midpoints');
xlabel('x'); ylabel('y'); zlabel('proxy');
view(45, 30);
grid on;
axis tight;

%%
figure(); hold on; tol = 1e-3;
fgrid = f(Xg,Yg);
contourf(Xg, Yg, reshape(fgrid,500,500), 20, 'LineColor', 'none','FaceAlpha',0.6);  % background: f(x,y)
colormap(turbo); colorbar;
caxis([tol 1]);
contourData = scatteredInterpolant(XX, YY, UU, 'natural', 'none');
Zg = contourData(Xg, Yg);
contour(Xg, Yg, Zg, 20, 'k'); 
 

%%
cu = cu_glb;
x0 = [-0.81,0.77]';
%x0 = [0.89, 0.92]';
tol = 1e-3;
%plot_sol(meshP,meshT,cu,V_abc,R,S,M);

step_size = 0.05;
step_tol = 1e-1;
max_steps = 1000;




% Plot contour of U over scattered points
% along with inverse velocity field
figure; hold on;

fgrid = f(Xg,Yg);
contourf(Xg, Yg, reshape(fgrid,500,500), 30, 'LineColor', 'none','FaceAlpha',0.6);  % background: f(x,y)
colormap(turbo); colorbar;
caxis([tol 1]);  % to emphasize obstacle contrast
title('Time-to-goal surface with obstacles and path');
axis equal tight;


% Overlay contour lines
contourData = scatteredInterpolant(XX, YY, UU, 'natural', 'none');
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

eik_path = trace_path_greedy_proxy(cu_glb, meshT, meshP, DT, Dx_a1bc1, Dy_ab1c1, H_a1bc1, H_ab1c1, ...
                                   V_abc, proxy, f, x0, goal, ...
                                   H_abc, R, S, n, eps_obst);



[gradU,iTri,tx] = evaluate_grad_u(eik_path(:,end),cu,meshT,meshP,...
                            DT,Dx_a1bc1,Dy_ab1c1,H_a1bc1,H_ab1c1,n,H_abc);
z_offset = max(UU) + 1e-3;  % slight lift above the max surface
plot3(eik_path(1,:), eik_path(2,:), z_offset * ones(1,size(eik_path,2)), 'r-', 'LineWidth', 2);
plot3(eik_path(1,1), eik_path(2,1), z_offset, 'ko', 'MarkerFaceColor', 'k');  % start
plot3(eik_path(1,end), eik_path(2,end), z_offset, 'ko', 'MarkerFaceColor', 'k');  % goal


function plot_sol(meshP,meshT,cu,V_abc,R,S,M)
figure()
nTri = size(meshT,1);
for iTri = 1:nTri

    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    %T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    %semilogy(abs(cu_abc)); hold on;
    u = V_abc*cu_abc;
    %trisurf(T,X,Y,u); hold on;
    scatter3(X,Y,u,'.'); hold on;

end
drawnow;
hold off;
end