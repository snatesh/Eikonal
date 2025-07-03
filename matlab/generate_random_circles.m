function [centers, r0s] = generate_random_circles(N, domain, r_range, sep, max_attempts)
% Generate N non-overlapping circles with separation sep between boundaries
%
% Inputs:
%   N         : number of circles
%   domain    : [xmin xmax ymin ymax]
%   r_range   : [rmin rmax]
%   sep       : required separation between circle boundaries
%   max_attempts (optional) : max number of attempts (default: 1000)
%
% Outputs:
%   centers   : 2 x N array of circle centers
%   r0s       : 1 x N array of radii

    if nargin < 5
        max_attempts = 1000;
    end

    xmin = domain(1); xmax = domain(2);
    ymin = domain(3); ymax = domain(4);
    rmin = r_range(1); rmax = r_range(2);

    centers = zeros(2, N);
    r0s = zeros(1, N);
    count = 0; attempts = 0;

    while count < N && attempts < max_attempts
        attempts = attempts + 1;

        % Propose new center and radius
        c_new = [xmin + (xmax - xmin) * rand;
                 ymin + (ymax - ymin) * rand];
        r_new = rmin + (rmax - rmin) * rand;

        % Check for sufficient separation from all existing circles
        ok = true;
        for j = 1:count
            dist = norm(c_new - centers(:,j));
            if dist < (r_new + r0s(j) + sep)
                ok = false;
                break;
            end
        end

        if ok
            count = count + 1;
            centers(:,count) = c_new;
            r0s(count) = r_new;
        end
    end

    if count < N
        warning('Only placed %d non-overlapping circles out of %d after %d attempts.', ...
                 count, N, attempts);
        centers = centers(:,1:count);
        r0s = r0s(1:count);
    end
end
