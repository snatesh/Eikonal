function eik_path = trace_path_greedy_quad(XY_quad, U_quad, x0, goal, r_search, tol, max_steps)
% XY_quad: (N × 2) physical quadrature points
% U_quad:  (N × 1) solution values at those points

eik_path = x0(:)';  % initialize eik_path
x = x0(:);

for step = 1:max_steps
    plot(x(1),x(2),'r.'); hold on; drawnow;
    % Euclidean distances to all quad points
    dists = vecnorm(XY_quad - x', 2, 2);

    % Find candidate points within radius
    in_radius = find(dists < r_search);
    if isempty(in_radius)
        disp('empty');
        break;  % nothing to move toward
    end

    % Among these, pick point with lowest value of u
    [~, idx_min] = min(U_quad(in_radius));
    next_x = XY_quad(in_radius(idx_min), :);

    % Stop if not moving significantly
    if norm(next_x - x') < 1e-14
        disp('not moving');
        break;
    end

    % Update
    x = next_x(:);
    eik_path = [eik_path; x'];

    if norm(x - goal(:)) < tol
        break;
    end
end
