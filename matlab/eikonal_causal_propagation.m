function [cu_glb, V_abc, seed_tri, proxy] = eikonal_causal_propagation(DT, x_anchor, m, R, S, W, udirch, rhs)

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
  disp(x_anchor)
  % Find triangle(s) sharing this edge
  attached_tris = edgeAttachments(DT, anchor_edge);
  disp(attached_tris)
  queue = attached_tris{1};  % Initial "known" triangle(s)
  iTri = attached_tris{1}(1);
  seed_tri = iTri;
  [K,B,F,G,meshP,meshT,dPdP] = assemble_poisson_causal_tri(iTri,DT,n,R,S,W,...
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
  nTri = size(tris,1);
  cu_glb = zeros(M*nTri,1);
  cu_glb((iTri-1)*M+1:iTri*M) = cu0(1:M); 


  
  xi0 = 0.1;
  alph = 1;
  max_iter = 10000;
  rtol = 1e-10;

  % newton + homotopic continuation
  xi_min = 7e-2;      % Do not let xi drop below this
  xi = xi0;
  xi_decay = 0.8;     % Reduce xi by this factor
  % from init
  cu = cu0;
  [N, H] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb, M, W, xi, ...
                                              V_abc, dPdP, ...
                                              meshT, meshP, K);
  H = H';


  BB_id = B;
  tol = norm(N + xi*K*cu(1:M) - BB_id'*cu(M+1:end) - F) * rtol;

  xi_trigger = 10*tol;  % When residual norm drops below this, reduce xi


  fprintf("It \t (g,du) \t\t ||g|| \t\t xi\n");
  HHmat = zeros(M + nLambda);
  for it = 1:max_iter

    % Assemble global system
    HHmat(1:M, 1:M) = H;
    HHmat(1:M, M+1:end) = BB_id';
    HHmat(M+1:end, 1:M) = BB_id;

    % Compute residual
    brhs = [N + xi*K*cu(1:M) + BB_id'*cu(M+1:end) - F;
            BB_id*cu(1:M)];

    norm_g = norm(brhs);

    if norm_g < tol
      fprintf("Converged in %d iterations\n", it);
      break;
    end

    % Newton update (no line search)
    du = -HHmat \ brhs;
    cu = cu + alph*du;
    cu_abc = cu(1:M);
    cu_glb((iTri-1)*M+1:iTri*M) = cu_abc;
    % Update residual and Jacobian
    [N,H] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb, M, W, xi, ...
                                               V_abc, dPdP, ...
                                               meshT, meshP, K);

    brhs = [N + xi*K*cu(1:M) + BB_id'*cu(M+1:end) - F;
            BB_id*cu(1:M)];
    norm_g = norm(brhs);

    H = H';

    fprintf("%d \t %.4e \t %.4e \t %.3e\n", it, -brhs'*du, norm_g, xi);

    % Dynamically reduce xi if residual is small
    if norm_g < xi_trigger && xi > xi_min
      xi = max(xi * xi_decay, xi_min);
      alph = max(0.01,0.9*alph);
      [N, H] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb, M, W, xi, ...
                                              V_abc, dPdP, ...
                                              meshT, meshP, K);
      H = H';
    end

  end

  cu_abc = cu(1:M);
  u_sol = V_abc*cu_abc;
  % compute integral of u on tri
  minu_seed = max(0,min(u_sol));
  intu_seed = (W/2)'*u_sol;
  cu_mat = reshape(cu_abc, [], 1, 1);           % M×1×1
  cu_mat_T = permute(cu_mat, [2 1 3]);        % 1×M×1

  % Batch: (1×M×1) × (M×M×Nq) × (M×1×1) → (1×1×Nq)
  grad_mag_sq = pagemtimes(cu_mat_T, pagemtimes(dPdP, cu_mat));  % 1×1×Nq

  % Reshape to Nq×1 vector
  grad_mag_sq = reshape(grad_mag_sq, [], 1);  % Nq×1
  weights = (W/2) .* sqrt(grad_mag_sq);
  weights = weights / sum(weights);
  cw_fac = 0.8;
  proxy_val_seed = (1-cw_fac)*sum(weights .* u_sol) + cw_fac*minu_seed;


  fprintf("proxy on seed %.3e\n", proxy_val_seed);
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
    proxy(t) = proxy_val_seed;
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
          cu_known_list{end+1} = cu_glb((t_known-1)*M+1:t_known*M);  % M×1
        end
        known_tris = known_nbrs;
        [K,B,F,G,meshP,meshT,dPdP] = ...
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
        cu = cu0;
        iTri = t_nbr;

        cu_glb((iTri-1)*M+1:iTri*M) = cu0(1:M);

        xi0 = 0.1;
        alph = 1;
        max_iter = 10000;

        % newton + homotopic continuation
        xi = xi0;
        xi_decay = 0.8;     % Reduce xi by this factor
        % from init
        cu = cu0;
        [N, H] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb, M, W, xi, ...
                                                    V_abc, dPdP, ...
                                                    meshT, meshP, K);
        H = H';


        tol = norm(N + xi*K*cu(1:M) - BB_id'*cu(M+1:end) - F) * rtol;

        xi_trigger = 10*tol;  % When residual norm drops below this, reduce xi


        fprintf("It \t (g,du) \t\t ||g|| \t\t xi\n");
        HHmat = zeros(M + nLambda);
        for it = 1:max_iter

          % Assemble global system
          HHmat(1:M, 1:M) = H;
          HHmat(1:M, M+1:end) = BB_id';
          HHmat(M+1:end, 1:M) = BB_id;

          % Compute residual
          brhs = [N + xi*K*cu(1:M) + BB_id'*cu(M+1:end) - F;
            BB_id*cu(1:M) - Gbc];

          norm_g = norm(brhs);

          if norm_g < tol
            fprintf("Converged in %d iterations\n", it);
            break;
          end

          % Newton update (no line search)
          du = -HHmat \ brhs;
          cu = cu + alph*du;
          cu_abc = cu(1:M);
          cu_glb((iTri-1)*M+1:iTri*M) = cu_abc;
          % Update residual and Jacobian
          [N,H] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb, M, W, xi, ...
            V_abc, dPdP, ...
            meshT, meshP, K);

          brhs = [N + xi*K*cu(1:M) + BB_id'*cu(M+1:end) - F;
            BB_id*cu(1:M) - Gbc];
          norm_g = norm(brhs);

          H = H';

          fprintf("%d \t %.4e \t %.4e \t %.3e\n", it, -brhs'*du, norm_g, xi);

          % Dynamically reduce xi if residual is small
          if norm_g < xi_trigger && xi > xi_min
            xi = max(xi * xi_decay, xi_min);
            alph = max(0.01,0.9*alph);
            [N, H] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb, M, W, xi, ...
              V_abc, dPdP, ...
              meshT, meshP, K);
            H = H';
          end

        end

        cu_abc = cu(1:M);
        u_sol = V_abc*cu_abc;
        % compute integral of u on tri
        minu_nbr = max(0,min(u_sol));
        intu_nbr = (W/2)'*u_sol;
        cu_mat = reshape(cu_abc, [], 1, 1);           % M×1×1
        cu_mat_T = permute(cu_mat, [2 1 3]);        % 1×M×1

        % Batch: (1×M×1) × (M×M×Nq) × (M×1×1) → (1×1×Nq)
        grad_mag_sq = pagemtimes(cu_mat_T, pagemtimes(dPdP, cu_mat));  % 1×1×Nq

        % Reshape to Nq×1 vector
        grad_mag_sq = reshape(grad_mag_sq, [], 1);  % Nq×1
        weights = (W/2) .* sqrt(grad_mag_sq);
        weights = weights / sum(weights);
        proxy_val_nbr = (1-cw_fac)*sum(weights .* u_sol) + cw_fac*minu_nbr;


        fprintf("proxy on nbr %.3e\n", proxy_val_nbr);
        cu_glb((iTri-1)*M+1:iTri*M) = cu_abc;

               
        % Simulate causal update
        status(t_nbr) = "Known";
        proxy(t_nbr) = proxy(T_cur) + proxy_val_nbr;  % causally increasing
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
