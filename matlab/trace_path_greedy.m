function path = trace_path_greedy(x0, goal, cu_glb, meshT, meshP, DT, ...
                                   Dx_a1bc1, Dy_ab1c1, H_a1bc1, H_ab1c1, n, ...
                                   step_size, tol, max_steps)

    x = x0(:);               % current point
    goal = goal(:);          % ensure column vector
    path = x';               % initialize path
    M = n*(n+1)/2;

    for step = 1:max_steps
        % Evaluate ∇u at current point
        [g,iTri] = evaluate_grad_u(x, cu_glb, meshT, meshP, DT, ...
                            Dx_a1bc1, Dy_ab1c1, H_a1bc1, H_ab1c1, n);
        disp(norm(g))

        % Check if ∇u is invalid or degenerate
        if any(~isfinite(g)) || norm(g) < 1e-10
            warning("Step %d: using goal direction fallback at x = [%f, %f]", ...
                    step, x(1), x(2));

            d_goal = goal - x;
            if norm(d_goal) > 1e-10
                v = d_goal / norm(d_goal);  % unit vector toward goal
            else
                v = [0; 0];  % already at goal
            end
        else
            v = -g / norm(g);  % normalized descent direction
        end

        % Take step
        x_new = x + step_size * v;
        % Get triangle barycenter
        vtx = meshP(:, meshT(iTri,:));
        xc = mean(vtx, 2);
        
        % Nudge toward barycenter to avoid edge
        dir = xc - x_new;
        x_new = x_new + 1e-4 * dir / norm(dir);


        % Clip to domain bounding box
        x_new = min(1, max(-1, x_new));
        disp(x_new)

        % Save to path
        path(end+1, :) = x_new';

        % Stop if close to goal
        if norm(x_new - goal) < tol
            break;
        end

        % Update position
        x = x_new;
    end
end
