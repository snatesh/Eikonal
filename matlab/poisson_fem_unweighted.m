clear all; close all; clc;
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
  default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
  set(groot, default_name,'latex');
end
set(groot, 'defaultLegendFontSize',30)
set(groot, 'defaultAxesFontSize',30)
set(groot, 'defaultLineLineWidth',1)



% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

ms = 4:20;
errs = zeros(length(ms),1);
for iM = 1:length(ms)
    % there are M total polys
    m = ms(iM); n = m+1; M = n*(n+1)/2;
    tripts = [2 1; 3 5; 1.5 4];
    %tripts = [0 0; 1 0; 0 1];
    DT = delaunayTriangulation(tripts);
    meshT = DT.ConnectivityList;
    meshP = DT.Points';
    J = IncidenceMatrix(meshP);
    detJ = det(J);
    usol = @(x,y) exp(x+y);
    udirch = @(x,y) usol(x,y);

    % finite element discretization for poisson
    rhs  = @(x,y) -2*exp(x+y);
    [V_abc,dxP,dyP,Vl,Vb,Vh,VlVl,VbVb,VhVh,Rl,Sl,Rb,Sb,Rh,Sh,wleg] = preAssemble_poisson(n,R,S);
    [K,B,F] = assemble_poisson_bulk(n,R,S,W,wleg, ...
                                    V_abc,dxP,dyP,...
                                    VlVl,VbVb,VhVh,...
                                    rhs,meshP);
    G = assemble_poisson_bnd(M,Rl,Sl,Rb,Sb,Rh,Sh,wleg,...
                             Vl,Vb,Vh,udirch,meshP);

    [~,~,II] = qr(B,0);
    nLambda = rank(B);
    BB = B(II(1:nLambda),:);
    Amat = zeros(M+nLambda);
    Amat(1:M,1:M) = K;
    Amat(1:M,M+1:end) = BB';
    Amat(M+1:end,1:M) = BB;
    Fvec = [F;G(II(1:nLambda))];
    cu = Amat\Fvec; cu_abc = cu(1:M);
    XYe = (J * [R,S]' + meshP(:,1))';
    errs(iM) = sqrt((W/2)'*((V_abc*cu_abc-usol(XYe(:,1),XYe(:,2))).^2 * detJ))

end

semilogy(ms,errs,'ro-'); hold on; semilogy(ms, exp(-2*ms), 'k--')
%%


%%
clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

DT = delaunayTriangulation([0 0; 1 0; 1 1; 0 1; 0.5 0.25; 0.5 0.75]);
%DT = delaunayTriangulation([0 0; 1 0; 0 1; 0.5 0.25; 0.5 0.75]);

%DT = delaunayTriangulation([0 0; 1 0; 0 1]);
%DT = delaunayTriangulation([0 0; 1 0; 1 1]);
%DT = delaunayTriangulation([0 0; 1 0; 1 1; 0 1]);
%DT = annulus();

meshT = DT.ConnectivityList;
meshP = DT.Points';

% Extract all edges from triangles[~,~,II] = qr(Bint_glb','vector');
% nLambda = rank(Bint_glb);
% BBint = Bint_glb(II(1:nLambda),:);
% 
% 
% Amat = zeros(M*nTri+nLambda);
% Amat(1:M*nTri,1:M*nTri) = K_glb;
% Amat(1:M*nTri,M*nTri+1:end) = BBint';
% Amat(M*nTri+1:end,1:M*nTri) = BBint;
% Gbc = zeros(nLambda,1);
% Fvec = [F_glb;Gbc];

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
m = 10; n = m+1; M = n*(n+1)/2;

% manufactured sol and corresponding dirchlet data/rhs
% usol = @(x,y) exp(x+y);
% udirch = @(x,y) usol(x,y);
% rhs  = @(x,y) -2*exp(x+y);
usol = @(x,y) sin(x.*y).*exp(x+y);
udirch = @(x,y) usol(x,y);
rhs = @(x,y) ((x.^2+y.^2).*sin(x.*y) - 2*(x+y).*cos(x.*y) - 2*sin(x.*y)).*exp(x+y);

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
    F_glb,G_glb,nintEdge,...
    nbndEdge,nbndTri] = assemble_poisson_pwc(n,R,S,W,...
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
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    scatter3(X,Y,u,'.'); hold on;
    %scatter3(X,Y,usol(X,Y),'.');
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
    %scatter3(XYl(:,1),XYl(:,2),Vl*cu_abc,'.'); hold on; 
    %scatter3(XYb(:,1),XYb(:,2),Vb*cu_abc,'.'); 
    %scatter3(XYh(:,1),XYh(:,2),Vh*cu_abc,'.'); 
    %scatter3(XYl(:,1),XYl(:,2),usol(XYl(:,1),XYl(:,2)),'o'); hold on; 
    %scatter3(XYb(:,1),XYb(:,2),usol(XYb(:,1),XYb(:,2)),'o'); 
    %scatter3(XYh(:,1),XYh(:,2),usol(XYh(:,1),XYh(:,2)),'o'); 
    errs(iTri) = sqrt((W/2)'*((V_abc*cu_abc-usol(X,Y)).^2 * detJ));

end


%%
[K_glb,Bint_glb,Bdirch_glb,...
 F_glb,G_glb,...
 nintEdge,nbndEdge,nbndTri] = assemble_poisson(n,R,S,W,...
                                                Rl,Sl,Rb,Sb,Rh,Sh,...
                                                V_abc,dxP,dyP,Vl,Vb,Vh,...
                                                intVlVl,intVbVb,intVhVh,...
                                                intVlVl_flip,...
                                                intVbVb_flip,...
                                                intVhVh_flip,wleg,...
                                                rhs,udirch,DT);
BB = [Bint_glb;Bdirch_glb];
% identify the nLambda independent rows of BB
[~,~,II] = qr(BB','vector');
nLambda = rank(BB);
BB_id = BB(II(1:nLambda),:);
Amat = zeros(M*nTri+nLambda);
Amat(1:M*nTri,1:M*nTri) = K_glb;
Amat(1:M*nTri,M*nTri+1:end) = BB_id';
Amat(M*nTri+1:end,1:M*nTri) = BB_id;
Gbc = [zeros(M*nintEdge,1);G_glb];
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
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    scatter3(X,Y,u,'o'); hold on;
    %scatter3(X,Y,usol(X,Y),'.');
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
    scatter3(XYl(:,1),XYl(:,2),Vl*cu_abc,'.'); hold on; 
    scatter3(XYb(:,1),XYb(:,2),Vb*cu_abc,'.'); 
    scatter3(XYh(:,1),XYh(:,2),Vh*cu_abc,'.'); 
    %scatter3(XYl(:,1),XYl(:,2),usol(XYl(:,1),XYl(:,2)),'.'); hold on; 
    %scatter3(XYb(:,1),XYb(:,2),usol(XYb(:,1),XYb(:,2)),'.'); 
    %scatter3(XYh(:,1),XYh(:,2),usol(XYh(:,1),XYh(:,2)),'.'); 
    errs(iTri) = sqrt((W/2)'*((V_abc*cu_abc-usol(X,Y)).^2 * detJ));

end










%%%%%%%%%%%%%%%%%%%%%% FEM ASSEMBLY %%%%%%%%%%%%%%%%%%%









function [K,B,F] = assemble_poisson_bulk(n,R,S,W,wleg, ...
                                         V_abc,dxP,dyP,...
                                         VlVl,VbVb,VhVh,...
                                         rhs,Pts_Ti)
J = IncidenceMatrix(Pts_Ti);
detJ = det(J);
invJ = inv(J);
lenE1 = norm(Pts_Ti(:,2)-Pts_Ti(:,1));
lenE2 = norm(Pts_Ti(:,3)-Pts_Ti(:,2));
lenE3 = norm(Pts_Ti(:,1)-Pts_Ti(:,3));
XY = (J * [R,S]' + Pts_Ti(:,1))';
X = XY(:,1);
Y = XY(:,2);
M = n*(n+1)/2;

% interior stiffness
K = zeros(M,M); 
for i = 1:M
    dxP_i = dxP(:,i);
    dyP_i = dyP(:,i);
    gradP_i = invJ'*[dxP_i';dyP_i'];
    dxP_i = gradP_i(1,:)';
    dyP_i = gradP_i(2,:)';
    for j = 1:M
        dxP_j = dxP(:,j);
        dyP_j = dyP(:,j);   
        gradP_j = invJ'*[dxP_j';dyP_j'];
        dxP_j = gradP_j(1,:)';
        dyP_j = gradP_j(2,:)';
        dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
        K(i,j) = (W/2)'*dPidPj*detJ;
    end
end

% boundary stiffness
% TODO: have to pass in mesh to this function
%       so that we can evaluate the integrals
%       for weak enforcement of dirichlet conditions
%       on the boundary edges only - not all edges
%       we need to know which parametric edge 
%       corresponds to the current boundary edge as well
B = zeros(M,M); nleg = length(wleg);
for i = 1:M
    for j = 1:M
        VlVlij = reshape(VlVl(:,i,j),nleg,1);
        VbVbij = reshape(VbVb(:,i,j),nleg,1);
        VhVhij = reshape(VhVh(:,i,j),nleg,1);
        int1 = wleg'*VlVlij;
        int2 = wleg'*VbVbij;
        int3 = wleg'*VhVhij;
        % as is, this evaluates the line integral 
        % over the whole triangle
        B(i,j) = int1*lenE1 + int2*lenE2 + int3*lenE3;
    end
end

% load vector
F = zeros(M,1);
for i = 1:M
    P_if = V_abc(:,i).*rhs(X,Y);
    F(i) = (W/2)'*P_if*detJ;
end


end

function G = assemble_poisson_bnd(M,Rl,Sl,Rb,Sb,Rh,Sh,wleg,...
                                  Vl,Vb,Vh,udirch,Pts_Ti,varargin)

J = IncidenceMatrix(Pts_Ti);
lenE1 = norm(Pts_Ti(:,2)-Pts_Ti(:,1));
lenE2 = norm(Pts_Ti(:,3)-Pts_Ti(:,2));
lenE3 = norm(Pts_Ti(:,1)-Pts_Ti(:,3));
% boundary load
G = zeros(M,1);
XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
gl = udirch(XYl(:,1),XYl(:,2));
gb = udirch(XYb(:,1),XYb(:,2));
gh = udirch(XYh(:,1),XYh(:,2));
for i = 1:M
    int1 = wleg'*(Vl(:,i).*gl);
    int2 = wleg'*(Vb(:,i).*gb);
    int3 = wleg'*(Vh(:,i).*gh);
    G(i) = int1*lenE1 + int2*lenE2 + int3*lenE3;
end


end




function [V_abc,dxP,dyP,Vl,Vb,Vh,VlVl,VbVb,VhVh,...
          Rl,Sl,Rb,Sb,Rh,Sh,wleg] = preAssemble_poisson(n,R,S)
a = 1/2; b = a; c = a;
M = n*(n+1)/2;

% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);

% vandermonde under (a,b,c), (a+1,b,c+1), (a,b+1,c+1)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
V_a1bc1 = jPoly_tri(R,S,H_a1bc1,n-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(R,S,H_ab1c1,n-1,a,b+1,c+1);

% derivative amatrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);


Nrs = length(R);

% interior basis derivatives at abscissa
dxP = zeros(Nrs,M);
dyP = zeros(Nrs,M);
for i = 1:M
    dxP(:,i) = V_a1bc1 * Dx_a1bc1(:,i);
    dyP(:,i) = V_ab1c1 * Dy_ab1c1(:,i);
end

nleg = 50;
[xleg,wleg,~] = gjQuad(nleg,0,0);
wleg = wleg'/2; xleg = (xleg+1)/2;
Rl = 0*xleg; Sl = xleg;
Rb = xleg; Sb = 0*xleg;
Rh = xleg; Sh = 1-Rh;

Vl = jPoly_tri(Rl,Sl,H_abc,n-1,a,b,c);
Vb = jPoly_tri(Rb,Sb,H_abc,n-1,a,b,c);
Vh = jPoly_tri(Rh,Sh,H_abc,n-1,a,b,c);
VlVl = zeros(nleg,M,M);
VbVb = zeros(nleg,M,M);
VhVh = zeros(nleg,M,M);


for i = 1:M
    for j = 1:M
        VlVl(:,i,j) = Vl(:,i).*Vl(:,j);
        VbVb(:,i,j) = Vb(:,i).*Vb(:,j);
        VhVh(:,i,j) = Vh(:,i).*Vh(:,j);
    end
end

end






function DT = annulus()
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
figure;
triplot(DT_annulus, 'Color', 'b');
axis equal;
title('Delaunay Triangulation of Annulus');
xlabel('x'); ylabel('y');
end 
