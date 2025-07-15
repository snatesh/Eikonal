function eikonal_causal_propagation(DT, x_anchor, m, R, S, W, udirch, rhs)

  % TODO: For every instance of proxy simulation for u,
  %       we need to add a solve
  %       - for each triangle sharing edge
  %         1) Get dirichlet data from edge
  %            and assemble pointwise Bdirch
  %         2) assemble stiffness, load and 
  %           nonlinear contribution as usual
  %           NOTE: there is no longer a need for local
  %                 to global mappings. we solve on each
  %                 triangle independently for M coeffs     
  %         3) Pass all to newton iteration after initializing
  %            with poisson solve
  %           NOTE: This is a newton iteration per triangle
  %                 - first try it without homotopy 
  %
  %
  %       REQUIRE: x_anchor should belong to only one triangle
  %                  i.e. it should not be a mesh vertex
  %                - otherwise, the initial solve will involve
  %                  multiple triangles, with no inter-element
  %                  continuity constraints. These tris will
  %                  never be revisited, which means the domain-level
  %                  solution will have pointwise discontinuities there.
  %                - if there is only one seed triangle, then all visited
  %                  triangles will have to match dirichlet data on edges
  %                  shared with previous triangles

  %       COMPLETED: Seed triangle initialization for poisson
  %       PLAN: For each neighbor of seed, identify edge which shares it
  %             1) Evaluate sol on that edge, which is now passed in as 
  %                udirch to the assemble_poisson_causal_tri, along with  
  %                the iTri
  %             2) Identify which reference edge in neighbor tri we need,
  %                and determine whether any flipping of basis evals is needed
  %             3) If flip, then we flip the dirichlet data as well. 

  % Inputs:
  %   DT       - delaunayTriangulation object
  %   x_anchor - 1x2 anchor point [x, y], must lie on the boundary between mesh vertices
  %   m        - polynomial total degree
  %   (R,S,W)  - quadrature rule on simplex
  %   udirch   - dirichlet data at x_anchor (typically == 0)
  %   rhs      - inverse speed function


  % Step 0: solver precomputation
  n = m+1; M = n*(n+1)/2;   
  [V_abc,dxP,dyP,Vl,Vb,Vh,...
   Vl_flip,Vb_flip,Vh_flip,...
   Rl,Sl,Rb,Sb,Rh,Sh,...
   Rl_flip,Sl_flip,Rb_flip,Sb_flip,...
   Rh_flip,Sh_flip,wleg,...
   intVlVl,intVbVb,intVhVh,...
   intVlVl_flip,intVbVb_flip,intVhVh_flip,...
   dxPl,dyPl,dxPb,dyPb,dxPh,dyPh,...
   dxPl_flip,dyPl_flip,dxPb_flip,dyPb_flip,...
   dxPh_flip,dyPh_flip] = preAssemble_poisson1(n,R,S);

  % normalization under (a,b,c)
  a = 1/2; b = a; c = a;
  H_abc = structure_factors_tri(n+1,a,b,c);


 
  % Step 1: Find triangle(s) whose edge contains x_anchor
  edges = DT.freeBoundary();
  pts = DT.Points;
  tris = DT.ConnectivityList; 
  % Find closest edge midpoint to x_anchor
  midpoints = (pts(edges(:,1),:) + pts(edges(:,2),:)) / 2;
  [~, idx] = min(vecnorm(midpoints - x_anchor, 2, 2));
  anchor_edge = edges(idx, :);
  
  % Find triangle(s) sharing this edge
  attached_tris = edgeAttachments(DT, anchor_edge);
  disp(attached_tris)
  queue = attached_tris{1};  % Initial "known" triangle(s)
  iTri = attached_tris{1}(1);
  [K,B,F,G,meshP,meshT] = assemble_poisson_causal_tri(iTri,DT,n,R,S,W,...
                                          Rl,Sl,Rb,Sb,Rh,Sh, ...
                                          Rl_flip,Sl_flip,...
                                          Rb_flip,Sb_flip,...
                                          Rh_flip,Sh_flip,wleg,...
                                          V_abc,dxP,dyP,Vl,Vb,Vh,...
                                          Vl_flip,Vb_flip,Vh_flip,...
                                          rhs,udirch,1,x_anchor,H_abc);
  % solve poisson on seed
  nLambda = 1;
  Amat = sparse(zeros(M+nLambda));
  Amat(1:M,1:M) = K;
  Amat(1:M,M+1:end) = B';
  Amat(M+1:end,1:M) = B;
  Gbc = G;
  Fvec = [F;Gbc];
  % solve for solution modes on each tri
  cu0 = Amat \ Fvec;
  cu_abc = cu0(1:M);
  u_sol = V_abc*cu_abc;
  % compute integral of u on tri
  intu_seed = (W/2)'*u_sol;
  fprintf("integral of u on seed %.3e\n", intu_seed);

  nTri = size(tris,1);

  cu_glb = zeros(M*nTri,1);
  cu_glb((iTri-1)*M+1:iTri*M) = cu_abc;
  

  
  % Initialize sets
  numTris = size(tris, 1);
  status = strings(numTris, 1);     % "Unknown", "Known", "Trial"
  status(:) = "Unknown";
  proxy = inf(numTris, 1);          % ∫ u value (for priority)
  parent = zeros(numTris, 1);       % Store source triangle (for arrows)
  
  % Assign random proxy to seed triangles and mark as Known
  for it = 1:length(queue)
    t = queue(it);
    status(t) = "Known";
    proxy(t) = intu_seed; % ∫ u
  end
  
  % Frontier: list of triangles eligible to update (naive FIFO)
  frontier = queue;
  
  while ~isempty(frontier)
    % Select triangle in frontier with minimum proxy
    [~, idx_min] = min(proxy(frontier));
    T_cur = frontier(idx_min);
    frontier(idx_min) = [];  % remove from frontier
    
    % Get neighbors of current triangle
    nbrs = neighbors(DT, T_cur);
    nbrs = nbrs(~isnan(nbrs));
    
    %cu_known = cu_glb((T_cur-1)*M+1:T_cur*M);

    for inbr = 1:length(nbrs)
      t_nbr = nbrs(inbr);
      
      if status(t_nbr) == "Unknown"
        % === (1) Identify all known neighbors of t_nbr ===
        nbrs2 = neighbors(DT, t_nbr);
        nbrs2 = nbrs2(~isnan(nbrs2));
        known_nbrs = nbrs2(status(nbrs2) == "Known");

        % === (2) For each known neighbor, compute shared edge and get cu_known ===
        shared_edges = {};     % cell array of 1×2 integer row vectors
        cu_known_list = {};    % cell array of M×1 modal vectors

        for ikn = 1:length(known_nbrs)
          t_known = known_nbrs(ikn);

          verts_known = tris(t_known, :);
          verts_nbr   = tris(t_nbr, :);

          edges_known = sort([verts_known([1 2]); verts_known([2 3]); verts_known([3 1])], 2);
          edges_nbr   = sort([verts_nbr([1 2]);   verts_nbr([2 3]);   verts_nbr([3 1])],   2);

          shared_edge = intersect(edges_known, edges_nbr, 'rows');  % 1×2

          shared_edges{end+1} = shared_edge;
          cu_known_list{end+1} = cu_glb((t_known-1)*M + (1:M));  % M×1
        end
        disp(size(shared_edges,2));
        known_tris = known_nbrs;
        % % Vertex indices of each triangle
        % verts_cur = tris(T_cur, :);
        % verts_nbr = tris(t_nbr, :);
        % 
        % % All edges of each triangle (as sorted pairs)
        % edges_cur = sort([verts_cur([1 2]); verts_cur([2 3]); verts_cur([3 1])], 2);
        % edges_nbr = sort([verts_nbr([1 2]); verts_nbr([2 3]); verts_nbr([3 1])], 2);
        % 
        % % Find shared edge: intersection of edge lists
        % shared_edge = intersect(edges_cur, edges_nbr, 'rows');
        [K,B,F,G,meshP,meshT] = ...
          assemble_poisson_causal_tri(t_nbr,DT,n,R,S,W,...
                                      Rl,Sl,Rb,Sb,Rh,Sh, ...
                                      Rl_flip,Sl_flip,...
                                      Rb_flip,Sb_flip,...
                                      Rh_flip,Sh_flip,wleg,...
                                      V_abc,dxP,dyP,Vl,Vb,Vh,...
                                      Vl_flip,Vb_flip,Vh_flip,...
                                      rhs,cu_known_list,0,shared_edges,known_tris);

        BB = B;
        % identify the nLambda independent rows of BB
        [~,~,II] = qr(BB','vector');
        nLambda = rank(BB);
        BB_id = BB(II(1:nLambda),:);
        Amat = sparse(zeros(M+nLambda));
        Amat(1:M,1:M) = K;
        Amat(1:M,M+1:end) = BB_id';
        Amat(M+1:end,1:M) = BB_id;
        Gbc = G(II(1:nLambda));
        Fvec = [F;Gbc];
        cu0 = Amat \ Fvec;
        cu_abc = cu0(1:M);
        u_sol = V_abc*cu_abc;
        cu_glb((t_nbr-1)*M+1:t_nbr*M) = cu_abc;

        % compute integral of u on tri
        intu_nbr = (W/2)'*u_sol;
        
        
        % Simulate causal update
        status(t_nbr) = "Known";
        proxy(t_nbr) = proxy(T_cur) + intu_nbr;  % causally increasing
        parent(t_nbr) = T_cur;
        frontier(end+1) = t_nbr;
      end
    end
  end
  
  % === Visualization ===
  figure();
  for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    Ixe = IncidenceMatrix(Pts_Ti);
    XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    %T = delaunay(X,Y);
    cu_abc = cu_glb((iTri-1)*M+1:M*iTri);
    u_sol = V_abc*cu_abc;
    scatter3(X,Y,u_sol,'.'); hold on;
  end



  figure; hold on; axis equal;
  triplot(DT, 'Color', [0.8 0.8 0.8]);
  
  % Plot arrows for causal propagation
  centers = incenter(DT);  % triangle centers
  for t = 1:numTris
    if parent(t) > 0
      p = parent(t);
      c1 = centers(p, :);
      c2 = centers(t, :);
      quiver(c1(1), c1(2), c2(1)-c1(1), c2(2)-c1(2), 0, ...
          'MaxHeadSize', 0.5, 'Color', 'b', 'LineWidth', 1.2);
    end
  end
  
  % Plot anchor point
  plot(x_anchor(1), x_anchor(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
  title('Causal Propagation Graph');
end
