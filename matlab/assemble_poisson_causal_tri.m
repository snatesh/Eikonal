function [K,B,F,G,meshP,meshT] = assemble_poisson_causal_tri(iTri,DT,n,R,S,W,...
                                                             Rl,Sl,Rb,Sb,Rh,Sh, ...
                                                             Rl_flip,Sl_flip,...
                                                             Rb_flip,Sb_flip,...
                                                             Rh_flip,Sh_flip,wleg,...
                                                             V_abc,dxP,dyP,Vl,Vb,Vh,...
                                                             Vl_flip,Vb_flip,Vh_flip,...
                                                             rhs,udirch,mode,varargin)


% get mesh points, connectivity list
meshT = DT.ConnectivityList;
meshP = DT.Points';
% jacobi poly params
a = 1/2; b = a; c = a;
% total number of polynomials
M = n*(n+1)/2;
% size of edge quadrature
nleg = size(wleg,1);

% local stifness
K = zeros(M,M);
% local dirichlet 
% TODO: (need to add evaluator for non-seed)
B = zeros(nleg*1,M);
% local load
F = zeros(M,1);
% local dirichlet
G = zeros(nleg,1);


% get tri pts, incidence matrix and edge jacobians
Pts_Ti = meshP(:,meshT(iTri,:));

J = IncidenceMatrix(Pts_Ti);
detJ = det(J);
invJ = inv(J);
% quadrature points mapped to physical tri
XY = (J * [R,S]' + Pts_Ti(:,1))';
X = XY(:,1);
Y = XY(:,2);

% local interior stiffness
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

% local load (rhs)
for i = 1:M
    P_if = V_abc(:,i).*rhs(X,Y);
    F(i) = (W/2)'*P_if*detJ;
end


% TODO: Entire mode=0 block needs adapaptation for
%       causal/sequential assembly
if mode == 0
  cu_known_list = udirch;
  shared_edges = varargin{1};
  known_tris = varargin{2};
  %T_cur = varargin{2};
  nshared = size(shared_edges,2);
  B = zeros(nleg*nshared,M);
  G = zeros(nleg*nshared,1);

  for ishared = 1:nshared
    shared_edge = shared_edges{ishared};
    cu_known = cu_known_list{ishared};
    T_cur = known_tris(ishared);
    % get vertices of each tri
    Pts_T1 = meshP(:,meshT(T_cur,:));
    Pts_T2 = meshP(:,meshT(iTri,:));
    % get parametric map for each tri
    J1 = IncidenceMatrix(Pts_T1);
    J2 = IncidenceMatrix(Pts_T2);
    % get shared edge endpts
    ve = meshP(:, shared_edge);    % 2×2
    % get parametric coords of endpt on each tri
    rve1 = J1\(ve-Pts_T1(:,1));
    rve2 = J2\(ve-Pts_T2(:,1));
    % get the correct precomputed basis functions corresponding
    % to the parametric shared edge in each tri
    tol = 1e-13;
    edge1 = 0; edge2 = 0;
    rleg1 = 0; sleg1 = 0;
    rleg2 = 0; sleg2 = 0;
    % bottom
    if ((norm(rve1(:,1) - [0;0]) < tol) && (norm(rve1(:,2) - [1;0]) < tol)) || ...
        ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) -[0;0]) < tol))
      V1 = Vb;
      edge1 = 2;
      rleg1 = Rb; sleg1 = Sb;
      % hyp
    elseif ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) - [0;1]) < tol)) || ...
        ((norm(rve1(:,1) - [0;1])< tol) && (norm(rve1(:,2) - [1;0]) < tol))
      V1 = Vh;
      edge1 = 3;
      rleg1 = Rh; sleg1 = Sh;
      % left
    else
      V1 = Vl;
      edge1 = 1;
      rleg1 = Rl; sleg1 = Sl;
    end
    % bottom
    if ((norm(rve2(:,1) - [0;0]) < tol) && (norm(rve2(:,2) - [1;0]) < tol)) || ...
        ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;0]) < tol))
      V2 = Vb;
      edge2 = 2;
      rleg2 = Rb; sleg2 = Sb;
      % hyp
    elseif ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;1]) < tol)) || ...
        ((norm(rve2(:,1) - [0;1]) < tol) && (norm(rve2(:,2) - [1;0]) < tol))
      V2 = Vh;
      edge2 = 3;
      rleg2 = Rh; sleg2 = Sh;
      % left
    else
      V2 = Vl;
      edge2 = 1;
      rleg2 = Rl; sleg2 = Sl;
    end

    % Align the bases
    [Vknown, Vtarget] = getAlignedEdgeBases(V1, rleg1, sleg1, V2, rleg2, sleg2, Pts_T1, Pts_T2, edge1, edge2);
    B((ishared-1)*nleg+1:ishared*nleg,:) = Vtarget;
    G((ishared-1)*nleg+1:ishared*nleg) = Vknown * cu_known;

  end


elseif mode == 1 % seed triangle
  % get the goal point
  goal = varargin{1};
  goal_pts = [goal(1);goal(2)];
  ngoal = 1;
  % structure factors
  H_abc = varargin{2};
  % global dirichlet
  B = zeros(ngoal,M);
  % global dirichlet load
  G = zeros(ngoal,1);

  for ipt = 1:ngoal
  Pts_Ti = meshP(:,meshT(iTri,:));
  J = IncidenceMatrix(Pts_Ti);
  % get parametric coords of goal region
  goal_pt = goal_pts(:,ipt);
  rsgoal = J\(goal_pt-Pts_Ti(:,1));
  Rg = rsgoal(1); Sg = rsgoal(2);
  % evaluate basis at goal point 
  V = jPoly_tri(Rg,Sg,H_abc,n-1,a,b,c);
  g = udirch(goal_pt(1),goal_pt(2));
  B(ipt,:) = V; G(ipt) = g;
  end

end


% now we loop over shared edges 
% to assemble inter-element continuity matrix
%for iedge = 1:nintEdge
%
%    % get triangles sharing the edge
%    shared_tris = sharedEdge_tri_map(iedge,:);
%    tri1 = shared_tris(1);
%    tri2 = shared_tris(2);
%
%    % get vertices of each tri
%    Pts_T1 = meshP(:,meshT(tri1,:));
%    Pts_T2 = meshP(:,meshT(tri2,:));
%    % get parametric map for each tri
%    J1 = IncidenceMatrix(Pts_T1);
%    J2 = IncidenceMatrix(Pts_T2);
%    % get shared edge endpts
%    ve = meshP(:,intEdges(iedge,:));
%    % get parametric coords of endpt on each tri
%    rve1 = J1\(ve-Pts_T1(:,1));
%    rve2 = J2\(ve-Pts_T2(:,1));
%    % get the correct precomputed basis functions corresponding
%    % to the parametric shared edge in each tri
%    tol = 1e-13;
%    edge1 = 0; edge2 = 0;
%    rleg1 = 0; sleg1 = 0;
%    rleg2 = 0; sleg2 = 0;
%    % bottom
%    if ((norm(rve1(:,1) - [0;0]) < tol) && (norm(rve1(:,2) - [1;0]) < tol)) || ... 
%       ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) -[0;0]) < tol))
%        V1 = Vb;
%        edge1 = 2;
%        rleg1 = Rb; sleg1 = Sb;
%    % hyp
%    elseif ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) - [0;1]) < tol)) || ...
%           ((norm(rve1(:,1) - [0;1])< tol) && (norm(rve1(:,2) - [1;0]) < tol))
%        V1 = Vh;
%        edge1 = 3;
%        rleg1 = Rh; sleg1 = Sh;
%    % left
%    else
%        V1 = Vl;
%        edge1 = 1;
%        rleg1 = Rl; sleg1 = Sl;
%    end
%    % bottom
%    if ((norm(rve2(:,1) - [0;0]) < tol) && (norm(rve2(:,2) - [1;0]) < tol)) || ...
%       ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;0]) < tol))
%        V2 = Vb;
%        edge2 = 2;
%        rleg2 = Rb; sleg2 = Sb;
%    % hyp
%    elseif ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;1]) < tol)) || ...
%           ((norm(rve2(:,1) - [0;1]) < tol) && (norm(rve2(:,2) - [1;0]) < tol))
%        V2 = Vh;
%        edge2 = 3;
%        rleg2 = Rh; sleg2 = Sh;
%    % left
%    else
%        V2 = Vl;
%        edge2 = 1;
%        rleg2 = Rl; sleg2 = Sl;
%    end
%
%    % Align the bases
%    [V1a, V2a] = getAlignedEdgeBases(V1, rleg1, sleg1, V2, rleg2, sleg2, Pts_T1, Pts_T2, edge1, edge2);
%
%    % now save to global inter-element continuity matrix
%    Bint_glb((iedge-1)*nleg+1:nleg*iedge,(tri1-1)*M+1:M*tri1) = V1a; 
%    Bint_glb((iedge-1)*nleg+1:nleg*iedge,(tri2-1)*M+1:M*tri2) = -V2a; 
%end


end
