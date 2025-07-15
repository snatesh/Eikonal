clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");
DT = delaunay_unit_square(5);
meshT = DT.ConnectivityList;
meshP = DT.Points';
nTri = size(meshT,1);

[centers, r0s, goal, x0] = generate_maze_circles(50, [-1 1 -1 1], [0.05 0.1], 0.02);

eps_obst = 1e-6;  % minimum speed inside obstacle

% Fixed: smooth dip for obstacle
f = @(x,y) multi_obstacle_speed_varying_radii(x,y,centers,r0s,eps_obst);
[xg, yg] = meshgrid(linspace(-1,1,200));
fg = f(xg, yg);
contourf(xg, yg, fg); hold on;
plot(goal(1),goal(2),'ko','MarkerFaceColor','k');
plot(x0(1),x0(2),'ro','MarkerFaceColor','r');

% compute min-time surface
m = 14; n = m+1; M = n*(n+1)/2; nquad = length(W);
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% dirchlet and neumann data/rhs
udirch = @(x,y) zeros(size(x));
rhs = @(x,y) 1./f(x,y);
goal = [1;0];
eikonal_causal_propagation(DT, goal', m, R, S, W, udirch, rhs)
