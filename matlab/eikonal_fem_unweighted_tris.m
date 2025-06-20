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
%DT = delaunay_unit_square(2);
%DT = annulus();
DT = delaunay_disk(20,3);
meshT = DT.ConnectivityList;
meshP = DT.Points';
nTri = size(meshT,1);


% there are M total polys
m = 10; n = m+1; M = n*(n+1)/2;
speed = 1;
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

% IMEX loop
xi = 0.0001;
dt = 0.01;
dt_factor = 2;
max_tstep = 10000;
tol = 1e-7;
LHS = [eye(M*nTri) + dt*xi*K_glb BB_id';...
    BB_id zeros(nLambda)];
% get current residual
u = cu;
[N_glb, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                            V_abc, dPdP, ...
                                            meshT, meshP, K_glb);
res_curr = norm([N_glb + xi*K_glb*u(1:M*nTri)-F_glb+BB_id'*u(M*nTri+1:end);...
                BB_id*u(1:M*nTri)]);
for tstep = 1:max_tstep

    % Assemble global system

    [N_glb, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                                V_abc, dPdP, ...
                                                meshT, meshP, K_glb);

    rhs = [u(1:M*nTri) - dt*(N_glb-F_glb); Gbc];


    % solve
    sol = LHS \ rhs;
    
    % get new residual
    res_new = norm([N_glb + xi*K_glb*sol(1:M*nTri)-F_glb+BB_id'*sol(M*nTri+1:end);...
                BB_id*sol(1:M*nTri)]);
    fprintf("tstep %d: res_new = %.3e\t res_old=%.3e\n", tstep, res_new, res_curr);

    if res_new < res_curr
        fprintf('accept step\n');
        u = sol;
        res_curr = res_new;
        for iTri = 1:nTri
            Pts_Ti = meshP(:,meshT(iTri,:));
            Ixe = IncidenceMatrix(Pts_Ti);
            detJ = det(Ixe);
            XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
            X = XYe(:,1);
            Y = XYe(:,2);
            %T = delaunay(X,Y);
            cu_abc = u((iTri-1)*M+1:M*iTri);
            u_sol = V_abc*cu_abc;
            scatter3(X,Y,u_sol); hold on;
        end
        drawnow;
        hold off;
    else
        fprintf('step rejected. (res_new-res_old)=%.3e\n',abs(res_new-res_curr)')
        break;
        %dt = dt*dt_factor;
        %if (dt < 1e-8)
        %    fprintf('dt no longer reducible. (res_new-res_old)=%.3e\n',abs(res_new-res_curr));
        %    break;
        %end
        % reassemble sys mat with new dt
        %LHS = [eye(M*nTri) + dt*xi*K_glb BB_id';...
        %       BB_id zeros(nLambda)];
            
    end
end

%% Backward euler+newton loop
xi = 0.01;
dt = 1;
alph = 0.1;
max_iter = 500;
max_tstep = 10000;
tol = 1e-10;
u = cu;
% Initialize residual history
res_curr = Inf;
for tstep = 1:max_tstep
    u_n = u;
    u_m = u_n;
    fprintf("== Time Step %d ==\n", tstep)
    converged = false;
    for iter = 1:max_iter
        % Assemble global system
        [N_glb, H_glb] = assemble_eik_nonlin_unweighted(u_m, M, W, xi, ...
                                                        V_abc, dPdP, ...
                                                        meshT, meshP, K_glb);
        H_glb = H_glb';


        % Compute residual
        Ru = [N_glb + xi*K_glb*u_m(1:M*nTri) + BB_id'*u_m(M*nTri+1:end) - F_glb;
            BB_id*u_m(1:M*nTri)];

        F = [(u_m(1:M*nTri) - u_n(1:M*nTri)) / dt + Ru(1:M*nTri,:);...
              Ru(M*nTri+1:end,:)];
        JF = zeros(M*nTri + nLambda);
        JF(1:M*nTri, 1:M*nTri) = eye(M*nTri)/dt+H_glb;
        JF(1:M*nTri, M*nTri+1:end) = BB_id';
        JF(M*nTri+1:end, 1:M*nTri) = BB_id;
        du = -JF \ F;
        u_m = u_m + alph*du;

        rel_err = norm(du(1:M*nTri)) / max(norm(u_m(1:M*nTri)), 1e-14);
        %fprintf("  Newton iter %d: ||du||/||u|| = %.2e\t dt = %.2e \t ||R|| = %.2e\n", iter, rel_err,dt,norm(Ru));
        if rel_err < tol
            converged = true;
            break;
        end
    end
    % Compute residual norm for this time step
    res_new = norm(Ru);
    
    if ~converged || res_new > res_curr
        % Reject time step and reduce dt
        dt = dt / 2;
        fprintf("  Rejected step %d: residual = %.3e ↑: dt → %.3e\n", tstep, res_new, dt);
        % retry this time step with smaller dt
        continue;
    else
        % Accept step
        u = u_m;
        res_curr = res_new;
        % Optionally increase dt for efficiency
        dt = min(dt * 1.25, 1); % optional upper bound
        fprintf("  Accepted step %d: residual = %.3e ↓: dt → %.3e\n", tstep, res_new,dt);

    end

    for iTri = 1:nTri

        Pts_Ti = meshP(:,meshT(iTri,:));
        Ixe = IncidenceMatrix(Pts_Ti);
        detJ = det(Ixe);
        XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
        X = XYe(:,1);
        Y = XYe(:,2);
        T = delaunay(X,Y);
        cu_abc = u_m((iTri-1)*M+1:M*iTri);
        u_sol = V_abc*cu_abc;
        trisurf(T,X,Y,u_sol); hold on;
    end
    drawnow;
    hold off;

    if norm(u(1:M*nTri)- u_n(1:M*nTri)) / max(norm(u(1:M*nTri)), 1e-14) < tol
        fprintf("Steady state reached at step %d.\n", tstep);
        break;
    end
end
%%
hold on;
[XX,YY] = meshgrid(linspace(0,1,10));
scatter3(XX(:),YY(:),distanceToUnitSquareBoundary(XX(:),YY(:)));

%%
xi0 = 0.1;          % Initial penalty
xi_min = 1e-3;      % Do not let xi drop below this
xi = xi0;
xi_decay = 0.5;     % Reduce xi by this factor
[N_glb, H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                       V_abc, dPdP, ...
                                       meshT, meshP, K_glb);
H_glb = H_glb';


max_iter = 5000;
rtol = 1e-8;
tol = norm(N_glb + xi*K_glb*cu(1:M*nTri) - BB_id'*cu(M*nTri+1:end) - F_glb) * rtol;

xi_trigger = 10*rtol;  % When residual norm drops below this, reduce xi


fprintf("It \t (g,du) \t\t ||g|| \t\t xi\t\t cond(H)\n");
figure(3);
for it = 1:max_iter

    alph = 0.25;
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
    
    brhs = [N_glb + xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end) - F_glb;
           BB_id*cu(1:M*nTri)];
    norm_g = norm(brhs);

    H_glb = H_glb';

    % Dynamically reduce xi if residual is small
    if norm_g < xi_trigger && xi > xi_min
        xi = max(xi * xi_decay, xi_min);
        [N_glb, H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                                        V_abc, dPdP, ...
                                                        meshT, meshP, K_glb);
        H_glb = H_glb';
    end

    fprintf("%d \t %.4e \t %.4e \t %.3e\t %.4e\n", it, -brhs'*du, norm_g, xi, cond(HHmat));


    for iTri = 1:nTri

        Pts_Ti = meshP(:,meshT(iTri,:));
        J = IncidenceMatrix(Pts_Ti);
        detJ = det(J);
        XYe = (J * [R,S]' + Pts_Ti(:,1))';
        X = XYe(:,1);
        Y = XYe(:,2);
        %T = delaunay(X,Y);
        cu_abc = cu((iTri-1)*M+1:M*iTri);
        semilogy(abs(cu_abc)); hold on;
        %u = V_abc*cu_abc;
        %trisurf(T,X,Y,u); hold on;
        %scatter3(X,Y,u,'.'); hold on; 
    end
    drawnow;
    %hold off;


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
hold on;
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

function DT_annulus = annulus()
% Parameters
r_inner = 0.5;
r_outer = 1.0;
numPoints = 50;

% Generate points uniformly in the annulus
points = [];
while size(points,1) < numPoints
    % Generate points in square bounding box
    pts = (rand(numPoints*2,1)*2 - 1) * r_outer;
    pts = [pts, (rand(numPoints*2,1)*2 - 1) * r_outer];
    
    % Keep points between inner and outer radius
    r = sqrt(sum(pts.^2, 2));
    mask = (r >= r_inner) & (r <= r_outer);
    
    points = [points; pts(mask,:)];
    points = unique(points,'rows','stable');
end
points = points(1:numPoints,:);

% Create delaunay triangulation
DT = delaunayTriangulation(points);

% Filter triangles: keep only those with centroid inside annulus
triangles = DT.ConnectivityList;
pts = DT.Points;

centroids = (pts(triangles(:,1),:) + pts(triangles(:,2),:) + pts(triangles(:,3),:)) / 3;
r_cent = sqrt(sum(centroids.^2, 2));
keep = (r_cent >= r_inner) & (r_cent <= r_outer);

% Create filtered triangulation
DT_annulus = triangulation(triangles(keep,:), pts);

% Plot
%figure;
%triplot(DT_annulus, 'Color', 'b');
%axis equal;
%title('Delaunay Triangulation of Annulus');
%xlabel('x'); ylabel('y');
end 

function DT = delaunay_disk(N_boundary, N_radial)
    % Generate boundary points on unit circle
    theta = linspace(0, 2*pi, N_boundary + 1);
    theta(end) = [];  % Remove duplicate point at 2*pi
    xb = cos(theta);
    yb = sin(theta);

    % Generate interior points in polar coordinates
    r = linspace(0, 1, N_radial + 1);
    r = r(2:end);  % skip r=0 (handled separately)
    p_interior = [];

    for i = 1:length(r)
        Ni = round(N_boundary * r(i));  % number of points on i-th ring
        thetai = linspace(0, 2*pi, Ni+1); thetai(end) = [];
        xi = r(i) * cos(thetai);
        yi = r(i) * sin(thetai);
        p_interior = [p_interior; xi(:), yi(:)];
    end

    % Combine all points (add center point)
    p = [0, 0; p_interior; xb(:), yb(:)];

    % Delaunay triangulation
    dt = delaunayTriangulation(p);
    t = dt.ConnectivityList;
    p = dt.Points;

    DT = delaunayTriangulation(p);
end