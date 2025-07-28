
function path = trace_path_greedy_proxy_pysim(cu_glb, meshT, meshP, DT, Dx, Dy, Hx, Hy, ...
                                              V_abc, proxy, f_handle, x0, x_goal, ...
                                              H_abc, R, S, n)

% INPUTS:
% cu_glb     : Global modal coefficients (M * nTris, 1)
% meshT      : Triangle connectivity (nTris, 3)
% meshP      : Vertex coordinates (2, nPts)
% DT         : DelaunayTriangulation object
% Dx, Dy     : Derivative operators
% Hx, Hy     : Physical grad maps for triangle
% V_abc      : Basis evaluations at quadrature points (Nq x M)
% proxy      : Proxy value per triangle (1 x nTris)
% f_handle   : Function handle @(x,y) for speed function
% x0         : Starting point
% x_goal     : Goal location (termination condition)

% PARAMETERS
MAX_ITERS = 2000;
TOL = 1e-2;
step_size = 0.001;

path = x0;
x = x0;
M = n*(n+1)/2;
a = 0.5; b = 0.5; c = 0.5;
for k = 1:MAX_ITERS
    % Evaluate gradient and locate triangle
    [g, iTri, uval] = evaluate_grad_u(x, cu_glb, meshT, meshP, DT, Dx, Dy, Hx, Hy, n, H_abc);
    if norm(x - x_goal) < TOL  || uval < TOL
        break;
    end
    
    if uval < 0
      d = x_goal - x;
      d = d / norm(d);
      x_new = x + step_size * d;
      path(:,end+1) = x_new;
      x = x_new;
      continue;
    end

    % Get neighbors
    nbrs = neighbors(DT, iTri);
    nbrs = nbrs(~isnan(nbrs));

    % Choose lowest proxy neighbor
    proxy_vals = proxy(nbrs);
    [~, i_min] = min(proxy_vals);
    tid = nbrs(i_min);
    % Map quadrature points to physical triangle
    verts = meshP(:, meshT(tid,:));
    J = IncidenceMatrix(verts);
    xy_pts = (J * [R,S]' + verts(:,1))';
    invJ = inv(J);

    % Evaluate u and f in neighbor
    cu_nbr = cu_glb((tid-1)*M+1 : tid*M);
    u_vals = V_abc * cu_nbr;
    % Evaluate speed f at quad pts
    f_vals = f_handle(xy_pts(:,1),xy_pts(:,2));
    % Score: combine u and inverse f
    score = u_vals; %+ 0.5 ./ (f_vals + 1e-6);
    [~, best_idx] = min(score);
    x_cand = xy_pts(best_idx,:)';
    x_new = x_cand;

    % Add to path
    path(:,end+1) = x_new;
    x = x_new;
end
end
