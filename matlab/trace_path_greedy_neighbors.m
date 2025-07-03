function path = trace_path_greedy_neighbors(x0, goal, XY_quad_cell, U_quad_cell, DT, tri_neighbors, r_search, tol, max_steps)

x = x0(:);               % current point
path = x';               % initialize path
visited = [];            % list of visited points

for step = 1:max_steps
    % Locate triangle
    plot(x(1), x(2), 'r.'); hold on; drawnow;
    iTri = pointLocation(DT, x');
    if isnan(iTri)
        disp('Exited mesh.');
        break;
    end

    candidate_pts = [];  % (N × 2)
    candidate_vals = []; % (N × 1)

    % Check current triangle and its neighbors
    tri_list = [iTri, tri_neighbors(iTri, ~isnan(tri_neighbors(iTri,:)))];

    for iT = tri_list
        XY_tri = XY_quad_cell{iT};   % (nquad × 2)
        U_tri  = U_quad_cell{iT};    % (nquad × 1)

        % Compute distances to x
        dists = sqrt((XY_tri(:,1) - x(1)).^2 + (XY_tri(:,2) - x(2)).^2);

        % Remove previously visited points
        not_visited = true(size(dists));
        for j = 1:size(visited, 1)
            not_visited = not_visited & vecnorm(XY_tri - visited(j,:), 2, 2) > 1e-12;
        end

        in_radius = find((dists < r_search) & not_visited);

        % Append valid candidates
        candidate_pts = [candidate_pts; XY_tri(in_radius,:)];     % (k × 2)
        candidate_vals = [candidate_vals; U_tri(in_radius)];      % (k × 1)
    end

    % Stop if nothing to move toward
    if isempty(candidate_pts)
        disp('No valid move found.');
        break;
    end

    % Choose point with lowest u among candidates
    [~, idx_min] = min(candidate_vals);
    next_x = candidate_pts(idx_min, :)';

    % Stop if step is too small
    if norm(next_x - x) < 1e-10
        disp('Step too small.');
        break;
    end

    % Update
    x = next_x;
    visited = [visited; x'];
    path = [path; x'];

    % Check if goal reached
    if norm(x - goal(:)) < tol
        disp('Reached goal.');
        break;
    end
end
