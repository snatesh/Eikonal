function eik_path = trace_path_greedy_quad(XY_quad, U_quad, x0, goal, r_search, tol, max_steps)
% XY_quad: (N × 2) physical quadrature points
% U_quad:  (N × 1) solution values at those points

eik_path = x0(:)';  % Initialize path with starting point
x = x0(:);

for step = 1:max_steps
    % Optional: visualize current step
    plot(x(1), x(2), 'r.'); hold on; drawnow;

    % Compute Euclidean distances from current point to all quad points
    dists = sqrt((XY_quad(:,1) - x(1)).^2 + (XY_quad(:,2) - x(2)).^2);

    % Find quadrature points within search radius
    in_radius = find(dists < r_search);
    if isempty(in_radius)
        disp('No neighbors within search radius.');
        break;
    end

    % Among these, choose the point with minimal value of u
    [~, local_idx_min] = min(U_quad(in_radius));
    global_idx = in_radius(local_idx_min);
    next_x = XY_quad(global_idx, :);  % 1×2 row vector

    % Check if movement is significant
    step_size = norm(next_x(:) - x);
    if step_size < 1e-10
        disp('Step size below threshold; stopping.');
        break;
    end

    % Update current position and path
    x = next_x(:);                   % Make x a column vector again
    eik_path = [eik_path; x'];       % Append row

    % Check for convergence to goal
    if norm(x - goal(:)) < tol
        disp('Reached goal.');
        break;
    end
end
