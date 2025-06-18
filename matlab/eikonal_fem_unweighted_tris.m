clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

%DT = delaunayTriangulation([0 0; 1 0; 1 1; 0 1; 0.5 0.25; 0.5 0.75]);
%DT = delaunayTriangulation([0 0; 1 0; 0 1; 0.5 0.25; 0.5 0.75]);

%DT = delaunayTriangulation([0 0; 1 0; 0 1]);
%DT = delaunayTriangulation([0 0; 1 0; 1 1]);
%DT = delaunayTriangulation([0 0; 1 0; 1 1; 0 1]);
DT = delaunay_unit_square(5);
meshT = DT.ConnectivityList;
meshP = DT.Points';

Edges = [meshT(:, [1,2]); meshT(:, [2,3]); meshT(:, [3,1])];

% Sort edges so that (i,j) and (j,i) are considered the same
Edges = sort(Edges, 2);

% Count occurrences of each edge
[uniqueEdges, ~, ic] = unique(Edges, 'rows');
counts = accumarray(ic, 1);

% classify edges 
% edges that appear twice are interior
% edges that appear once are on boundary
intEdges = uniqueEdges(counts == 2, :);
bndEdges = uniqueEdges(counts == 1,:);
nTri = size(meshT,1);
figure(1);
for j = 1:size(bndEdges,1)
    p1 = meshP(:,bndEdges(j,1));
    p2 = meshP(:,bndEdges(j,2));
    plot([p1(1),p2(1)],[p1(2),p2(2)],'r-'); hold on;
end

for j = 1:size(intEdges,1)
    p1 = meshP(:,intEdges(j,1));
    p2 = meshP(:,intEdges(j,2));
    plot([p1(1),p2(1)],[p1(2),p2(2)],'b-'); hold on;
end

% there are M total polys
m = 20; n = m+1; M = n*(n+1)/2;
speed = 1;
xi = 0.01;
% manufactured sol and corresponding dirchlet data/rhs
udirch = @(x,y) zeros(size(x));
rhs = @(x,y) RHS(x,y,speed);
% finite element discretization for poisson
% on triangulated domain
[V_abc,dxP,dyP,Vl,Vb,Vh,...
    Vl_flip,Vb_flip,Vh_flip,...
    Rl,Sl,Rb,Sb,Rh,Sh,...
    Rl_flip,Sl_flip,Rb_flip,Sb_flip,...
    Rh_flip,Sh_flip,wleg,...
    intVlVl,intVbVb,intVhVh,...
    intVlVl_flip,intVbVb_flip,intVhVh_flip] = preAssemble_poisson1(n,R,S);
nleg = size(wleg,1);

[K_glb,Bint_glb,Bdirch_glb,...
    F_glb,G_glb,meshP,meshT,...
    sharedEdge_tri_map,...
    nintEdge,nbndEdge,...
    bndEdge_tri_map,dPdP] = assemble_poisson_pwc(n,R,S,W,...
                                            Rl,Sl,Rb,Sb,Rh,Sh,...
                                            Rl_flip,Sl_flip, ...
                                            Rb_flip,Sb_flip, ...
                                            Rh_flip,Sh_flip,wleg, ...
                                            V_abc,dxP,dyP,Vl,Vb,Vh,...
                                            Vl_flip,Vb_flip,Vh_flip,...
                                            rhs,udirch,DT);
% stack the boundary and inter-element continuity constraints
BB = [Bint_glb;Bdirch_glb];
% identify the nLambda independent rows of BB
[~,~,II] = qr(BB','vector');
nLambda = rank(BB);
BB_id = BB(II(1:nLambda),:);
Amat = zeros(M*nTri+nLambda);
Amat(1:M*nTri,1:M*nTri) = K_glb;
Amat(1:M*nTri,M*nTri+1:end) = BB_id';
Amat(M*nTri+1:end,1:M*nTri) = BB_id;
Gbc = [zeros(nleg*nintEdge,1);G_glb];
Gbc = Gbc(II(1:nLambda));
Fvec = [F_glb;Gbc];
% solve for solution modes on each tri
cu = Amat \ Fvec;

figure(2);
for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    scatter3(X,Y,u,'.'); hold on;
    %trisurf(T,X,Y,u); hold on;
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
end




% xi0 = 0.1;          % Initial penalty
% xi_min = 1e-2;      % Do not let xi drop below this
% xi = xi0;
% xi_decay = 0.5;     % Reduce xi by this factor
% xi_trigger = 1e-8;  % When residual norm drops below this, reduce xi
% Initial nonlinear term
[N_glb, H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                       V_abc, dPdP, ...
                                       meshT, meshP, K_glb);
H_glb = H_glb';

max_iter = 1000;
rtol = 1e-6;
tol = norm(N_glb + xi*K_glb*cu(1:M*nTri) - BB_id'*cu(M*nTri+1:end) - F_glb) * rtol;

fprintf("It \t (g,du) \t\t ||g|| \t\t xi\n");
figure(3);
for it = 1:max_iter

    alph = 0.2;
    % Assemble global system
    HHmat = zeros(M*nTri + nLambda);
    HHmat(1:M*nTri, 1:M*nTri) = H_glb;
    HHmat(1:M*nTri, M*nTri+1:end) = BB_id';
    HHmat(M*nTri+1:end, 1:M*nTri) = BB_id;

    % Compute residual
    brhs = [N_glb + xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end) - F_glb;
            BB_id*cu(1:M*nTri)];

    norm_g = norm(brhs);

    if norm_g < tol
        fprintf("Converged in %d iterations\n", it);
        break;
    end

    % Newton update (no line search)
    du = -HHmat \ brhs;
    cu = cu + alph*du;

    % Update residual and Jacobian
    [N_glb,H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                             V_abc, dPdP, ...
                                             meshT, meshP, K_glb);
    
    %brhs = [N_glb + xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end) - F_glb;
    %        BB_id*cu(1:M*nTri)];
    %norm_g = norm(brhs);

    H_glb = H_glb';

    % % Dynamically reduce xi if residual is small
    % if norm_g < xi_trigger && xi > xi_min
    %     xi = max(xi * xi_decay, xi_min);
    %     [N_glb, H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
    %                                                     V_abc, dPdP, ...
    %                                                     meshT, meshP, K_glb);
    %     H_glb = H_glb';
    % end

    fprintf("%d \t %.4e \t %.4e \t %.3e\n", it, -brhs'*du, norm_g, xi);

    for iTri = 1:nTri

        Pts_Ti = meshP(:,meshT(iTri,:));
        J = IncidenceMatrix(Pts_Ti);
        detJ = det(J);
        XYe = (J * [R,S]' + Pts_Ti(:,1))';
        X = XYe(:,1);
        Y = XYe(:,2);
        T = delaunay(X,Y);
        cu_abc = cu((iTri-1)*M+1:M*iTri);
        u = V_abc*cu_abc;
        trisurf(T,X,Y,u); hold on;
        %scatter3(X,Y,u,'.'); drawnow; hold on; 
    end
    drawnow;
    hold off;


end

%%
figure(3);
for iTri = 1:nTri

    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    %scatter3(X,Y,abs(u),'.'); hold on;
     trisurf(T,X,Y,u); hold on;
end

%%
[XX,YY] = meshgrid(linspace(0,1,10));
scatter3(XX(:),YY(:),distanceToUnitSquareBoundary(XX(:),YY(:)));
%%
% now do an inf-dim newton stepper for eikonal
[N_glb,H_glb] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                       V_abc, dPdP,...
                                       meshT,meshP,K_glb);
H_glb = H_glb';


max_iter = 1000;
rtol = 1e-6;
tol = norm(N_glb+xi*K_glb*cu(1:M*nTri) - BB_id'*cu(M*nTri+1:end)-F_glb)*rtol;


fprintf("It \t (g,du) \t ||g||l2 \t cond(H)\n")
for it = 1:max_iter

    alph = 1;
    HHmat = zeros(M*nTri+nLambda);
    HHmat(1:M*nTri,1:M*nTri) = H_glb;
    HHmat(1:M*nTri,M*nTri+1:end) = BB_id';
    HHmat(M*nTri+1:end,1:M*nTri) = BB_id;


    brhs = [N_glb+xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end)-F_glb; BB_id*cu(1:M*nTri)];
    if norm(brhs) < tol
        fprintf('Converged in %d Newton iterations\n', it);
        break;
    end

    du = -HHmat\brhs; 
    cu1 = cu + alph*du;

    cu = cu1;
    
    fprintf("%d \t %1.10f \t %1.10f \t %f\n", it, -brhs'*du, norm(brhs),cond(HHmat))

    [N_glb] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                             V_abc, dPdP,...
                                             meshT,meshP,K_glb);
end
figure(3);
for iTri = 1:nTri

    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    scatter3(X,Y,u,'.'); hold on;
end
%figure(3);
[XX,YY] = meshgrid(linspace(0,1,10));
scatter3(XX(:),YY(:),distanceToUnitSquareBoundary(XX(:),YY(:)));

%%
%%%%%%%%%%%%%%%%%%%%
function rhs = RHS(R,S,speed)
rhs = 1./(speed*ones(size(R)));
end


function d = distanceToUnitSquareBoundary(x, y)
    % Assumes (x,y) ∈ [0,1]^2
    % Returns the shortest distance to the boundary of the unit square
    
    % Distance to each side
    dx = min(x, 1 - x);
    dy = min(y, 1 - y);
    
    % Minimum distance to any of the four edges
    d = min(dx, dy);
end

function dt = delaunay_unit_square(N)
    % Generate N x N grid in [0,1]^2 with (N+1)^2 points
    x = linspace(0, 1, N+1);
    [X, Y] = meshgrid(x, x);
    p = [X(:), Y(:)];

    % Generate square elements and split into two triangles
    t = [];
    for i = 1:N
        for j = 1:N
            % Indices of square corners
            n0 = (i-1)*(N+1) + j;
            n1 = n0 + 1;
            n2 = n0 + (N+1);
            n3 = n2 + 1;

            % Split square into two triangles
            t1 = [n0, n1, n3];
            t2 = [n0, n3, n2];
            t = [t; t1; t2];
        end
    end

    % Optional: create Delaunay triangulation object
    dt = delaunayTriangulation(p);
end