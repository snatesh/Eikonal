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
clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

DT = delaunayTriangulation([0 0; 1 0; 1 1; 0 1; 0.5 0.25; 0.5 0.75]);
meshT = DT.ConnectivityList;
meshP = DT.Points';

% Extract all edges from triangles
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


Edges = [meshT(:, [1,2]), meshT(:, [2,3]), meshT(:, [3,1])];

nbndEdge = size(bndEdges,1);
nintEdge = size(intEdges,1);
refedge_map = zeros(nTri,3);
bnd_map = zeros(nTri,3);

iTri = 4;
vT = meshP(:,meshT(iTri,:));
edgesT = [Edges(iTri,1:2);Edges(iTri,3:4);Edges(iTri,5:6)];
J = IncidenceMatrix(vT);
ve1 = meshP(:,edgesT(1,:));
ve2 = meshP(:,edgesT(2,:));
ve3 = meshP(:,edgesT(3,:));
rve1 = J\(ve1-vT(:,1));
rve2 = J\(ve2-vT(:,1));
rve3 = J\(ve3-vT(:,1));


[meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map] = process_mesh(DT);


%%

% there are M total polys
m = 10; n = m+1; M = n*(n+1)/2;

% manufactured sol and corresponding dirchlet data/rhs
usol = @(x,y) exp(x+y);
udirch = @(x,y) usol(x,y);
rhs  = @(x,y) -2*exp(x+y);

% finite element discretization for poisson
% on triangulated domain
[V_abc,dxP,dyP,Vl,Vb,Vh,...
 Rl,Sl,Rb,Sb,Rh,Sh,wleg,...
 intVlVl,intVbVb,intVhVh] = preAssemble_poisson1(n,R,S);
[K_glb,Bint_glb,Bdirch_glb,...
          F_glb,G_glb] = assemble_poisson(n,R,S,W,...
                                          Rl,Sl,Rb,Sb,Rh,Sh,wleg, ...
                                          V_abc,dxP,dyP,Vl,Vb,Vh,...
                                          intVlVl,intVbVb,intVhVh,...
                                          rhs,udirch,DT);



% interpolatory decomposition of Bdirch_glb
[~,~,IIdirch] = qr(Bdirch_glb,0);
nLambda_dirch = rank(Bdirch_glb);
BBdirch = Bdirch_glb(IIdirch(1:nLambda_dirch),:);
% interpolator decomposition of Bint_glb
[~,~,IIint] = qr(Bint_glb,0);
nLambda_int = rank(Bint_glb);
BBint = Bint_glb(IIint(1:nLambda_int),:);
Amat = zeros(M*nTri+nLambda_int+nLambda_dirch);
Amat(1:M*nTri,1:M*nTri) = K_glb;
BB = [BBint;BBdirch];

Amat(1:M*nTri,1:M*nTri) = K_glb;
Amat(1:M*nTri,M*nTri+1:end) = BB';
Amat(M*nTri+1:end,1:M*nTri) = BB;
%Fvec = [F;G(II(1:nLambda))];



%%%%%%%%%%%%%%%%%%%%%% FEM ASSEMBLY %%%%%%%%%%%%%%%%%%%


% TODO: should write functions which evaluates
%       inner products of edge restricted basis (VlVl,VbVb,VhVh)
%       in quadrature. then, we can just use the
%       precomputed evaluations when assembling Bint, Bdirch
%       - this can be added to preAssemble_poisson
%       - but, they will not include the edge Jacobian factor
% DONE

% TODO: we also want a function which maps each interior edge
%       to the two triangles which share it
% DONE

% TODO: once we have this,
% 1) loop over triangles
% 2) for each triangle, assemble the stiffness matrix K
% 3) identify if triangle has boundary edges
% 4) assemble the boundary dirichlet matrix contribution
%    which will be the sum of at most 2 line integrals 
%    with edge Jacobian factor 
% DONE

% Then, we do another loop

% 1) loop over the shared edge map
% 2) on a given shared edge, look up the two triangles which
%    share it
% 3) determine which edges in the reference triangle corresponds
%    to the edge in T1 and T2
% 4) pick the correct boundary continuity matrix contributions
%    from the precomputed values of int VlVl, ..
% 5) assemble it correctly in the global matrix
%    like [Bint_T1 0 0 Bint_T4];

% this ensures that we hit all interior edges once,
% and have correct incorporation of the continuity constraints
% to the global inter-element continuity block of the problem matrix.


function [K_glb,Bint_glb,Bdirch_glb,...
          F_glb,G_glb] = assemble_poisson(n,R,S,W,...
                                          Rl,Sl,Rb,Sb,Rh,Sh,wleg, ...
                                          V_abc,dxP,dyP,Vl,Vb,Vh,...
                                          intVlVl,intVbVb,intVhVh,...
                                          rhs,udirch,DT)


% get mesh points, connectivity list
% as well as interior/bndry edges
% and map taking shared edge to the two triangles which share it
[meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map] = process_mesh(DT);

% number of triangles in mesh
nTri = size(meshT,1);
% total number of polynomials
M = n*(n+1)/2;
% number of boundary and interior edges
nintEdge = size(intEdges,1);
nbndEdge = size(bndEdges,1);

% global stifness, dirichlet and inter-element continuity matrices
K_glb = zeros(M*nTri,M*nTri);
Bdirch_glb = zeros(M*nTri,M*nTri);
Bint_glb = zeros(M*nTri,M*nTri);
F_glb = zeros(M*nTri,1);
G_glb = zeros(M*nTri,1);


% now we loop over triangles 
% and assemble stiffness and dirichlet matrices
for iTri = 1:nTri

    % get tri pts, incidence matrix and edge jacobians
    Pts_Ti = meshP(:,meshT(iTri,:));

    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    invJ = inv(J);
    lenE = zeros(3,1);
    lenE(1) = norm(Pts_Ti(:,2)-Pts_Ti(:,1));
    lenE(2) = norm(Pts_Ti(:,3)-Pts_Ti(:,2));
    lenE(3) = norm(Pts_Ti(:,1)-Pts_Ti(:,3));
    % quadrature points mapped to phyiscal tri
    XY = (J * [R,S]' + Pts_Ti(:,1))';
    X = XY(:,1);
    Y = XY(:,2);
    % boundary quadrature points mapped to physical tri
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
  
    % local interior stiffness
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
    % save to global stiffness
    K_glb((iTri-1)*M+1:M*iTri,(iTri-1)*M+1:M*iTri) = K;


    % local load (rhs)
    F = zeros(M,1);
    for i = 1:M
        P_if = V_abc(:,i).*rhs(X,Y);
        F(i) = (W/2)'*P_if*detJ;
    end
    % save to global load
    F_glb((iTri-1)*M+1:M*iTri) = F;

    % local boundary stiffness
    Bdirch = zeros(M,M);
    % boundary load (dirichlet condition)
    G = zeros(M,1);
    % get edges of triangle, sorted to match bndry/int Edge lookup
    edgesT = sort([meshT(iTri, [1,2]);...
                   meshT(iTri, [2,3]);...
                   meshT(iTri, [3,1])],2);
    % for each edge of triangle
    for iedge = 1:3
        % get endpts
        ve = meshP(:,edgesT(iedge,:));
        % get parametric coords of endpt
        rve = J\(ve-Pts_Ti(:,1));
        % label them as left, bottom or hypotenuse edge
        left = false; bot = false; hyp = false;
        if (all(rve(:,1) == [0;0]) && all(rve(:,2) == [1;0])) || ...
           (all(rve(:,1) == [1;0]) && all(rve(:,2) == [0;0]))           
            bot = true;
        elseif (all(rve(:,1) == [1;0]) && all(rve(:,2) == [0;1])) || ...
               (all(rve(:,1) == [0;1]) && all(rve(:,2) == [1;0]))
            hyp = true;
        else
            left = true;
        end
        % test if it is a boundary edge
        for ibnd = 1:nbndEdge
            % if it is, compute and store
            % local dirichlet contribution
            if edgesT(iedge,:) == bndEdges(ibnd,:)
                intVV = 0; % needed for boundary stiffness
                V = 0; g = 0; % dirichlet values on bndry
                if bot
                    intVV = intVbVb;
                    g = udirch(XYb(:,1),XYb(:,2));
                    V = Vb;
                elseif hyp
                    intVV = intVhVh;
                    g = udirch(XYh(:,1),XYh(:,2));
                    V = Vh;
                elseif left
                    intVV = intVlVl;
                    g = udirch(XYl(:,1),XYl(:,2));
                    V = Vl;
                end
                
                for i = 1:M
                    for j = 1:M     
                        Bdirch(i,j) = Bdirch(i,j) + intVV(i,j)*lenE(iedge);
                    end
                    G(i) = G(i) + wleg'*(V(:,i).*g)*lenE(iedge);
                end
            end
        end
    end
    % now save to global dirichlet matrix
    Bdirch_glb((iTri-1)*M+1:M*iTri,(iTri-1)*M+1:M*iTri) = Bdirch; 
    % save dirichlet condition to global Gdirch
    G_glb((iTri-1)*M+1:M*iTri) = G;
end

% now we loop over shared edges 
% to assemble inter-element continuity matrix
for iedge = 1:nintEdge
    % get triangles sharing the edge
    shared_tris = sharedEdge_tri_map(iedge,:);
    tri1 = shared_tris(1);
    tri2 = shared_tris(2);
    % get vertices of each tri
    Pts_T1 = meshP(:,meshT(tri1,:));
    Pts_T2 = meshP(:,meshT(tri2,:));
    % get parametric map for each tri
    J1 = IncidenceMatrix(Pts_T1);
    J2 = IncidenceMatrix(Pts_T2);
    % get shared edge endpts
    ve = meshP(:,intEdges(iedge,:));
    % get edge length
    lenE = norm(ve(:,2)-ve(:,1));
    % get parametric coords of endpt on each tri
    rve1 = J1\(ve-Pts_T1(:,1));
    rve2 = J2\(ve-Pts_T2(:,1));
    % label them as left, bottom or hypotenuse edge
    intVV1 = 0; 
    % bottom
    if (all(rve1(:,1) == [0;0]) && all(rve1(:,2) == [1;0])) || ...
       (all(rve1(:,1) == [1;0]) && all(rve1(:,2) == [0;0]))
        intVV1 = intVbVb;
    % hyp
    elseif (all(rve1(:,1) == [1;0]) && all(rve1(:,2) == [0;1])) || ...
           (all(rve1(:,1) == [0;1]) && all(rve1(:,2) == [1;0]))
        intVV1 = intVhVh;
    % left
    else
        intVV1 = intVlVl;
    end
    intVV2 = 0;
    % bottom
    if (all(rve2(:,1) == [0;0]) && all(rve2(:,2) == [1;0])) || ...
       (all(rve2(:,1) == [1;0]) && all(rve2(:,2) == [0;0]))
        intVV2 = intVbVb;
    % hyp
    elseif (all(rve2(:,1) == [1;0]) && all(rve2(:,2) == [0;1])) || ...
           (all(rve2(:,1) == [0;1]) && all(rve2(:,2) == [1;0]))
        intVV2 = intVhVh;
    % left
    else
        intVV2 = intVlVl;
    end

    Bint1 = zeros(M,M);
    Bint2 = zeros(M,M);

    for i = 1:M
        for j = 1:M
            Bint1(i,j) = Bint1(i,j) + intVV1(i,j)*lenE;
            Bint2(i,j) = Bint2(i,j) + intVV2(i,j)*lenE;
        end
    end

    % now save to global inter-element continuity matrix
    Bint_glb((tri1-1)*M+1:M*tri1,(tri1-1)*M+1:M*tri1) = Bint1; 
    Bint_glb((tri2-1)*M+1:M*tri2,(tri2-1)*M+1:M*tri2) = Bint2; 

end

end



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

function [V_abc,dxP,dyP,Vl,Vb,Vh,...
          Rl,Sl,Rb,Sb,Rh,Sh,wleg,...
          intVlVl,intVbVb,intVhVh] = preAssemble_poisson1(n,R,S)
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


intVlVl = zeros(M,M);
intVbVb = zeros(M,M);
intVhVh = zeros(M,M);

for i = 1:M
    for j = 1:M
        VlVlij = reshape(Vl(:,i).*Vl(:,j),nleg,1);
        VbVbij = reshape(Vb(:,i).*Vb(:,j),nleg,1);
        VhVhij = reshape(Vh(:,i).*Vh(:,j),nleg,1);
        intVlVl(i,j) = wleg'*VlVlij;
        intVbVb(i,j) = wleg'*VbVbij;
        intVhVh(i,j) = wleg'*VhVhij;
    end
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

function Ixe = IncidenceMatrix(Xe)

Ixe = [Xe(:,2)-Xe(:,1), Xe(:,3)-Xe(:,1)];

end


function [meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map] = process_mesh(DT)

meshT = DT.ConnectivityList;
meshP = DT.Points';

% Extract all edges from triangles
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
nintEdge = size(intEdges,1);

sharedEdge_tri_map = zeros(nintEdge,2);


for iedge = 1:nintEdge

    sharedTris = edgeAttachments(DT,intEdges(iedge,:));
    sharedEdge_tri_map(iedge,:) = sharedTris{1};

end

% sort it to match vertex index ordering in intEdges
sharedEdge_tri_map = sort(sharedEdge_tri_map,2);

end