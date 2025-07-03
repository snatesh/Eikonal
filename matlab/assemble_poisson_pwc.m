function [K_glb,Bint_glb,Bdirch_glb,...
          F_glb,G_glb,...
          meshP, meshT,...
          sharedEdge_tri_map,...
          nintEdge,nbndEdge,...
          bndEdge_tri_map,bndEdges,normals,varargout] = assemble_poisson_pwc(n,R,S,W,...
                                                                             Rl,Sl,Rb,Sb,Rh,Sh, ...
                                                                             Rl_flip,Sl_flip,...
                                                                             Rb_flip,Sb_flip,...
                                                                             Rh_flip,Sh_flip,wleg,...
                                                                             V_abc,dxP,dyP,Vl,Vb,Vh,...
                                                                             Vl_flip,Vb_flip,Vh_flip,...
                                                                             rhs,udirch,DT, mode, varargin)


% get mesh points, connectivity list
% as well as interior/bndry edges
% and map taking shared edge to the two triangles which share it
[meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map,bndEdge_tri_map,normals] = process_mesh(DT);
% jacobi poly params
a = 1/2; b = a; c = a;
% number of triangles in mesh
nTri = size(meshT,1);
% total number of polynomials
M = n*(n+1)/2;
% number of boundary and interior edges
nintEdge = size(intEdges,1);
nbndEdge = size(bndEdges,1);
nleg = size(wleg,1);

% global stifness
K_glb = zeros(M*nTri,M*nTri);
% global dirichlet
Bdirch_glb = zeros(nleg*nbndEdge,M*nTri);
% global inter-element continuity
% there are nintEdge row blocks, each with M rows and M*nTri cols
Bint_glb = zeros(nleg*nintEdge,M*nTri);
F_glb = zeros(M*nTri,1);
G_glb = zeros(nleg*nbndEdge,1);


% now we loop over triangles 
% and assemble the stiffness matrix
hasdpdp = false; dPdP = 0;
if nargout >= 14
  dPdP = zeros(M,M,length(R));
end
for iTri = 1:nTri
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
            if nargout >= 14 && ~hasdpdp
              dPdP(i,j,:) = dPidPj;
            end 
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
end

% now we loop over boundary edges
% to assemble dirichlet bc matrix
% handling one boundary edge at a time
if mode == 0
  for ibndedge = 1:nbndEdge
    iTri = bndEdge_tri_map(ibndedge);
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    % boundary quadrature points mapped to physical tri
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
    XYl_flip = (J * [Rl_flip,Sl_flip]' + Pts_Ti(:,1))';
    XYb_flip = (J * [Rb_flip,Sb_flip]' + Pts_Ti(:,1))';
    XYh_flip = (J * [Rh_flip,Sh_flip]' + Pts_Ti(:,1))';
    % get edges of triangle, sorted to match bndry/int Edge lookup
    edgesT = sort([meshT(iTri, [1,2]);...
                   meshT(iTri, [2,3]);...
                   meshT(iTri, [3,1])],2);
    tol = 1e-13;
    V = 0; g = 0; % dirichlet values on bndry
    for iedge = 1:3
      if all(edgesT(iedge,:) == bndEdges(ibndedge,:))
        % get endpts
        ve = meshP(:,edgesT(iedge,:));
        % get parametric coords of endpt
        rve = J\(ve-Pts_Ti(:,1));
        % bottom and flip
        if (norm(rve(:,1) - [0;0]) < tol && norm(rve(:,2) - [1;0]) < tol)
            g = udirch(XYb(:,1),XYb(:,2));
            V = Vb;
        elseif (norm(rve(:,1) - [1;0]) < tol && norm(rve(:,2) - [0;0]) < tol)
            g = udirch(XYb_flip(:,1),XYb_flip(:,2));
            V = Vb_flip;
        % hyp and flip
        elseif (norm(rve(:,1) - [1;0]) < tol && norm(rve(:,2) - [0;1]) < tol)
            g = udirch(XYh(:,1),XYh(:,2));
            V = Vh;
        elseif (all(rve(:,1) - [0;1]) < tol && norm(rve(:,2) - [1;0]) < tol)
            g = udirch(XYh_flip(:,1),XYh_flip(:,2));
            V = Vh_flip;
        % left and flip
        elseif (norm(rve(:,1) - [0;1]) < tol && norm(rve(:,2) - [0;0]) < tol)
            left = true;
            g = udirch(XYl(:,1),XYl(:,2));
            V = Vl;
        elseif (norm(rve(:,1) - [0;0]) < tol && norm(rve(:,2) - [0;1]) < tol)
            left = true;
            g = udirch(XYl_flip(:,1),XYl_flip(:,2));
            V = Vl_flip;
        end
        break;
      else
        continue;
      end
    end
    % now save to global dirichlet matrix
    Bdirch_glb((ibndedge-1)*nleg+1:ibndedge*nleg,(iTri-1)*M+1:M*iTri) = V;
    % save dirichlet condition to global Gdirch
    G_glb((ibndedge-1)*nleg+1:ibndedge*nleg) = g;
  
  end
elseif mode == 1 %single goal point
  % get the goal point
  goal = varargin{1};
  ngoal = 0
  goal_pts = generate_goal_region(goal, 0.01, ngoal);

  % structure factors
  H_abc = varargin{2};
  % global dirichlet
  %Bdirch_glb = zeros(1,M*nTri);
  Bdirch_glb = zeros(ngoal+1,M*nTri);
  % global dirichlet load
  %G_glb = zeros(1,1);
  G_glb = zeros(ngoal+1,1);
  % find triangle containing goal
  %iTri = pointLocation(DT, goal);  
  %if isnan(iTri)
  %    error('Goal point is outside the mesh.');
  %end
  % find triangles containing goal region
  iTris = pointLocation(DT,goal_pts');
  if any(isnan(iTris))
    error('Goal region is outside the mesh.');
  end
  for igoal = 1:(ngoal+1) % goal_pts includes goal
    iTri = iTris(igoal);
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    goalpt = goal_pts(:,igoal);
    % get parametric coords of goal region
    rsgoal = J\(goalpt-Pts_Ti(:,1));
    Rg = rsgoal(1); Sg = rsgoal(2);
    % evaluate basis at goal point 
    V = jPoly_tri(Rg,Sg,H_abc,n-1,a,b,c);
    g = udirch(goalpt(1),goalpt(2));
    Bdirch_glb(igoal,(iTri-1)*M+1:M*iTri) = V;
    G_glb(igoal) = g;
  end

  %Pts_Ti = meshP(:,meshT(iTri,:));
  %J = IncidenceMatrix(Pts_Ti);
  %% get parametric coords of goal region
  %rsgoal = J\(goal'-Pts_Ti(:,1));
  %Rg = rsgoal(1); Sg = rsgoal(2);
  %% evaluate basis at goal point 
  %V = jPoly_tri(Rg,Sg,H_abc,n-1,a,b,c);
  %g = udirch(goal(1),goal(2));
  %Bdirch_glb(1,(iTri-1)*M+1:M*iTri) = V;
  %G_glb = g;

  %% weak neumann 
  % now handle Neumann enforcement
  Bneumm_glb = zeros(M*nbndEdge,M*nTri);
  Gneumm_glb = zeros(M*nbndEdge,1);
  % get deriv ops on edges
  dxPl = varargin{3};
  dyPl = varargin{4};
  dxPb = varargin{5};
  dyPb = varargin{6};
  dxPh = varargin{7};
  dyPh = varargin{8};
  dxPl_flip = varargin{9};
  dyPl_flip = varargin{10};
  dxPb_flip = varargin{11};
  dyPb_flip = varargin{12};
  dxPh_flip = varargin{13};
  dyPh_flip = varargin{14};
  uneumm = varargin{15};
  intgradPl_Pl = zeros(M,M);
  intgradPb_Pb = zeros(M,M);
  intgradPh_Ph = zeros(M,M);
  for ibndedge = 1:nbndEdge
    iTri = bndEdge_tri_map(ibndedge);
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    ne = normals(:,ibndedge);
    % compute \int_e \nabla \phi_i \cdot n \phi_j ds
    % without edge jacobian or parameterization sign (flip/unflipped)
    for i = 1:M

        dxPl_i = dxPl(:,i);
        dyPl_i = dyPl(:,i);
        gradPl_i = invJ'*[dxPl_i';dyPl_i'];
        gradPl_i_ne = gradPl_i' * ne;
        
        dxPb_i = dxPb(:,i);
        dyPb_i = dyPb(:,i);
        gradPb_i = invJ'*[dxPb_i';dyPb_i'];
        gradPb_i_ne = gradPb_i' * ne;
        
        dxPh_i = dxPh(:,i);
        dyPh_i = dyPh(:,i);
        gradPh_i = invJ'*[dxPh_i';dyPh_i'];
        gradPh_i_ne = gradPh_i' * ne;
      for j = 1:M
        intgradPl_Pl(i,j) = wleg'*(gradPl_i_ne .* Vl(:,j));
        intgradPb_Pb(i,j) = wleg'*(gradPb_i_ne .* Vb(:,j));
        intgradPh_Ph(i,j) = wleg'*(gradPh_i_ne .* Vh(:,j));
      end

    end 


    % boundary quadrature points mapped to physical tri
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
    % get edges of triangle, sorted to match bndry/int Edge lookup
    edgesT = sort([meshT(iTri, [1,2]);...
                   meshT(iTri, [2,3]);...
                   meshT(iTri, [3,1])],2);
    tol = 1e-13;
    intgradVV = 0; intgV = 0; % integral of grad basis * basis and weak neumann values on bndry
    for iedge = 1:3
      if all(edgesT(iedge,:) == bndEdges(ibndedge,:))
        % get endpts
        ve = meshP(:,edgesT(iedge,:));
        lenE = norm(ve(:,2)-ve(:,1));
        % get parametric coords of endpt
        rve = J\(ve-Pts_Ti(:,1));
        % bottom and flip
        if (norm(rve(:,1) - [0;0]) < tol && norm(rve(:,2) - [1;0]) < tol)
            intgV = (wleg'*(uneumm(XYb(:,1),XYb(:,2)) .* Vb))' * lenE;
            intgradVV = intgradPb_Pb * lenE;
        elseif (norm(rve(:,1) - [1;0]) < tol && norm(rve(:,2) - [0;0]) < tol)
            intgV = -(wleg'*(uneumm(XYb(:,1),XYb(:,2)) .* Vb))' * lenE;
            intgradVV = -intgradPb_Pb * lenE;
        % hyp and flip
        elseif (norm(rve(:,1) - [1;0]) < tol && norm(rve(:,2) - [0;1]) < tol)
            intgV = (wleg'*(uneumm(XYh(:,1),XYh(:,2)) .* Vh))' * lenE;
            intgradVV = intgradPh_Ph * lenE;
        elseif (all(rve(:,1) - [0;1]) < tol && norm(rve(:,2) - [1;0]) < tol)
            intgV = -(wleg'*(uneumm(XYh(:,1),XYh(:,2)) .* Vh))' * lenE;
            intgradVV = -intgradPh_Ph * lenE;
        % left and flip
        elseif (norm(rve(:,1) - [0;1]) < tol && norm(rve(:,2) - [0;0]) < tol)
            intgV = (wleg'*(uneumm(XYl(:,1),XYl(:,2)) .* Vl))' * lenE;
            intgradVV = intgradPl_Pl * lenE;
        elseif (norm(rve(:,1) - [0;0]) < tol && norm(rve(:,2) - [0;1]) < tol)
            intgV = -(wleg'*(uneumm(XYl(:,1),XYl(:,2)) .* Vl))' * lenE;
            intgradVV = -intgradPl_Pl * lenE;
        end
        break;
      else
        continue;
      end
    end
    % now save to global dirichlet matrix
    Bneumm_glb((ibndedge-1)*M+1:ibndedge*M,(iTri-1)*M+1:M*iTri) = intgradVV';
    % save dirichlet condition to global Gdirch
    Gneumm_glb((ibndedge-1)*M+1:ibndedge*M) = intgV;
  
  end
  %% pointwise neumann 
  % now handle Neumann enforcement
  %Bneum_glb = zeros(nleg*nbndEdge,M*nTri);
  %Gneumm_glb = zeros(nleg*nbndEdge,1);
  %% get deriv ops on edges
  %dxPl = varargin{3};
  %dyPl = varargin{4};
  %dxPb = varargin{5};
  %dyPb = varargin{6};
  %dxPh = varargin{7};
  %dyPh = varargin{8};
  %dxPl_flip = varargin{9};
  %dyPl_flip = varargin{10};
  %dxPb_flip = varargin{11};
  %dyPb_flip = varargin{12};
  %dxPh_flip = varargin{13};
  %dyPh_flip = varargin{14};
  %uneumm = varargin{15};
  %gradPl_ne = zeros(nleg,M);
  %gradPl_ne_flip = zeros(nleg,M);
  %gradPb_ne = zeros(nleg,M);
  %gradPb_ne_flip = zeros(nleg,M);
  %gradPh_ne = zeros(nleg,M);
  %gradPh_ne_flip = zeros(nleg,M);
  %for ibndedge = 1:nbndEdge
  %  iTri = bndEdge_tri_map(ibndedge);
  %  Pts_Ti = meshP(:,meshT(iTri,:));
  %  J = IncidenceMatrix(Pts_Ti);
  %  ne = normals(:,ibndedge);
  %  % compute gradient of each basis function
  %  % on regular and flipped edges
  %  for i = 1:M

  %      dxPl_i = dxPl(:,i);
  %      dyPl_i = dyPl(:,i);
  %      gradPl_i = invJ'*[dxPl_i';dyPl_i'];
  %      gradPl_ne(:,i) = gradPl_i' * ne;
  %      
  %      dxPl_flip_i = dxPl_flip(:,i);
  %      dyPl_flip_i = dyPl_flip(:,i);
  %      gradPl_flip_i = invJ'*[dxPl_flip_i';dyPl_flip_i'];
  %      gradPl_ne_flip(:,i) = gradPl_flip_i' * ne;
  %      
  %      dxPb_i = dxPb(:,i);
  %      dyPb_i = dyPb(:,i);
  %      gradPb_i = invJ'*[dxPb_i';dyPb_i'];
  %      gradPb_ne(:,i) = gradPb_i' * ne;
  %      
  %      dxPb_flip_i = dxPb_flip(:,i);
  %      dyPb_flip_i = dyPb_flip(:,i);
  %      gradPb_flip_i = invJ'*[dxPb_flip_i';dyPb_flip_i'];
  %      gradPb_ne_flip(:,i) = gradPb_flip_i' * ne;
  %      
  %      dxPh_i = dxPh(:,i);
  %      dyPh_i = dyPh(:,i);
  %      gradPh_i = invJ'*[dxPh_i';dyPh_i'];
  %      gradPh_ne(:,i) = gradPh_i' * ne;
  %      
  %      dxPh_flip_i = dxPh_flip(:,i);
  %      dyPh_flip_i = dyPh_flip(:,i);
  %      gradPh_flip_i = invJ'*[dxPh_flip_i';dyPh_flip_i'];
  %      gradPh_ne_flip(:,i) = gradPh_flip_i' * ne;

  %  end 


  %  % boundary quadrature points mapped to physical tri
  %  XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
  %  XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
  %  XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
  %  XYl_flip = (J * [Rl_flip,Sl_flip]' + Pts_Ti(:,1))';
  %  XYb_flip = (J * [Rb_flip,Sb_flip]' + Pts_Ti(:,1))';
  %  XYh_flip = (J * [Rh_flip,Sh_flip]' + Pts_Ti(:,1))';
  %  % get edges of triangle, sorted to match bndry/int Edge lookup
  %  edgesT = sort([meshT(iTri, [1,2]);...
  %                 meshT(iTri, [2,3]);...
  %                 meshT(iTri, [3,1])],2);
  %  tol = 1e-13;
  %  V = 0; g = 0; % grad basis and neumann values on bndry
  %  for iedge = 1:3
  %    if all(edgesT(iedge,:) == bndEdges(ibndedge,:))
  %      % get endpts
  %      ve = meshP(:,edgesT(iedge,:));
  %      % get parametric coords of endpt
  %      rve = J\(ve-Pts_Ti(:,1));
  %      % bottom and flip
  %      if (norm(rve(:,1) - [0;0]) < tol && norm(rve(:,2) - [1;0]) < tol)
  %          g = uneumm(XYb(:,1),XYb(:,2)) * normals(:,ibndedge);
  %          V = gradPb_ne;
  %      elseif (norm(rve(:,1) - [1;0]) < tol && norm(rve(:,2) - [0;0]) < tol)
  %          g = uneumm(XYb_flip(:,1),XYb_flip(:,2)) * normals(:,ibndedge);
  %          V = gradPb_ne_flip;
  %      % hyp and flip
  %      elseif (norm(rve(:,1) - [1;0]) < tol && norm(rve(:,2) - [0;1]) < tol)
  %          g = uneumm(XYh(:,1),XYh(:,2)) * normals(:,ibndedge);
  %          V = gradPh_ne;
  %      elseif (all(rve(:,1) - [0;1]) < tol && norm(rve(:,2) - [1;0]) < tol)
  %          g = uneumm(XYh_flip(:,1),XYh_flip(:,2)) * normals(:,ibndedge);
  %          V = gradPh_ne_flip;
  %      % left and flip
  %      elseif (norm(rve(:,1) - [0;1]) < tol && norm(rve(:,2) - [0;0]) < tol)
  %          g = uneumm(XYl(:,1),XYl(:,2)) * normals(:,ibndedge);
  %          V = gradPl_ne;
  %      elseif (norm(rve(:,1) - [0;0]) < tol && norm(rve(:,2) - [0;1]) < tol)
  %          g = uneumm(XYl_flip(:,1),XYl_flip(:,2)) * normals(:,ibndedge);
  %          V = gradPl_ne_flip;
  %      end
  %      break;
  %    else
  %      continue;
  %    end
  %  end
  %  % now save to global dirichlet matrix
  %  Bneumm_glb((ibndedge-1)*nleg+1:ibndedge*nleg,(iTri-1)*M+1:M*iTri) = V;
  %  % save dirichlet condition to global Gdirch
  %  Gneumm_glb((ibndedge-1)*nleg+1:ibndedge*nleg) = g;
  %
  %end


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
    [V1a, V2a] = getAlignedEdgeBases(V1, rleg1, sleg1, V2, rleg2, sleg2, Pts_T1, Pts_T2, edge1, edge2);

    % now save to global inter-element continuity matrix
    Bint_glb((iedge-1)*nleg+1:nleg*iedge,(tri1-1)*M+1:M*tri1) = V1a; 
    Bint_glb((iedge-1)*nleg+1:nleg*iedge,(tri2-1)*M+1:M*tri2) = -V2a; 
end

  if nargout >= 14
    varargout{1} = dPdP;
  end
  if nargout == 16
    varargout{2} = Bneumm_glb;
    varargout{3} = Gneumm_glb;
  end

end
