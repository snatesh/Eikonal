function path = trace_path_greedy_local(x0, goal, XY_quad_cell, U_quad_cell, DT, r_search, tol, max_steps)

x = x0(:);                         % Ensure column vector
path = x';                         % Initialize path
M = size(U_quad_cell{1}, 1);       % Quadrature points per triangle
visited = [];
for step = 1:max_steps
    % Find triangle containing current point
    plot(x(1), x(2), 'r.'); hold on; drawnow;
    iTri = pointLocation(DT, x');
    if isnan(iTri)
        disp('Point exited mesh.');
        break;
    end

    % Quadrature points and values in current triangle
    XY_tri = XY_quad_cell{iTri};    % (M × 2)
    U_tri  = U_quad_cell{iTri};     % (M × 1)

    % Compute distances to all quad points in triangle
    dists = sqrt((XY_tri(:,1) - x(1)).^2 + (XY_tri(:,2) - x(2)).^2);
    % Exclude previously visited points
    tol_same = 1e-12;
    not_visited = true(length(dists), 1);
    for j = 1:size(visited, 1)
        not_visited = not_visited & vecnorm(XY_tri - visited(j,:), 2, 2) > tol_same;
    end
    
    in_radius = find((dists < r_search) & not_visited);

    if isempty(in_radius)
        disp('No quad points in radius.');
        break;
    end

    % Among nearby points, pick the one with lowest u
    [~, local_idx_min] = min(U_tri(in_radius));
    disp(min(U_tri(in_radius)));
    next_x = XY_tri(in_radius(local_idx_min), :)';

    % Stop if step is too small
    if norm(next_x - x) < 1e-10
        disp('Step too small.');
        break;
    end

    % Update path
    x = next_x;
    visited = [visited; x'];
    path = [path; x'];

    % Check for convergence
    if norm(x - goal(:)) < tol
        disp('Reached goal.');
        break;
    end
end
