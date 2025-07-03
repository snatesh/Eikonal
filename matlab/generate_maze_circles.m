function [centers, r0s, goal, x0] = generate_maze_circles(N, domain, r_range, sep)
% Generate N circular obstacles in a maze-like pattern
% Inputs:
%   N        : number of circles
%   domain   : [xmin xmax ymin ymax]
%   r_range  : [rmin rmax]
%   sep      : spacing between circle boundaries
%
% Outputs:
%   centers  : 2 x <=N array of circle centers
%   r0s      : 1 x <=N array of radii

    xmin = domain(1); xmax = domain(2);
    ymin = domain(3); ymax = domain(4);
    rmin = r_range(1); rmax = r_range(2);

    % Estimate grid size (as square as possible)
    ncols = ceil(sqrt(N));
    nrows = ceil(N / ncols);

    % Grid spacing
    xspace = (xmax - xmin) / ncols;
    yspace = (ymax - ymin) / nrows;
    base_r = min(xspace, yspace) / 2 - sep;

    % Safety cap on radius
    base_r = min(base_r, rmax);

    centers = [];
    r0s = [];

    for i = 1:nrows
        for j = 1:ncols
            if length(r0s) >= N
                break;
            end

            % Maze effect: skip some cells to make corridors
            if mod(i + j, 3) == 0  % tweak this rule for different maze styles
                continue;
            end

            cx = xmin + (j - 0.5)*xspace;
            cy = ymin + (i - 0.5)*yspace;

            % Jitter for natural feel
            jitter_x = 0.1 * xspace * (2*rand - 1);
            jitter_y = 0.1 * yspace * (2*rand - 1);
            c = [cx + jitter_x; cy + jitter_y];

            % Random radius in safe bounds
            r = rmin + (min(base_r, rmax) - rmin)*rand;

            centers = [centers, c];
            r0s = [r0s, r];
        end
    end

    % Find a random point not inside any circle (with buffer)
    max_attempts = 1000; found_goal = false;
    for k = 1:max_attempts
        goal_candidate = [xmin + (xmax - xmin) * rand;
                          ymin + (ymax - ymin) * rand];
        dists = vecnorm(centers - goal_candidate, 2, 1);
        if all(dists > r0s + sep)
            goal = goal_candidate;
            found_goal = true;
        end
    end
    if ~found_goal
      error('Could not find a valid goal point outside obstacles after %d attempts.', max_attempts);
    end
    
    if nargin < 5
        dmin = 0.5 * norm([domain(2)-domain(1); domain(4)-domain(3)]);
    end

    % Now generate x0 far from goal, outside all obstacles
    for k = 1:max_attempts
        x0_candidate = [xmin + (xmax - xmin) * rand;
                        ymin + (ymax - ymin) * rand];
        dists = vecnorm(centers - x0_candidate, 2, 1);
        if all(dists > r0s + sep) && norm(x0_candidate - goal) > dmin
            x0 = x0_candidate;
            return;
        end
    end
    error('Could not find valid x0 far from goal after %d attempts.', max_attempts);

end
