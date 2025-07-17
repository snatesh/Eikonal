function fvals = multi_obstacle_speed_varying_radii(X, Y, centers, r0s, tol)
% Speed function with k circular obstacles of varying radii.
%
% Inputs:
%   X, Y     : evaluation grid (matrices of same size)
%   centers : 2-by-k matrix of circle centers
%   r0s     : 1-by-k vector of radii for each obstacle
%   tol     : scalar speed inside any obstacle
%
% Output:
%   fvals   : speed values of same size as X and Y

    [rows, cols] = size(X);
    fvals = ones(rows, cols);  % default speed = 1
    k = size(centers, 2);
    for j = 1:k
        cx = centers(1, j);
        cy = centers(2, j);
        r2 = r0s(j)^2;

        % Mask for points inside this obstacle
        inside = (X - cx).^2 + (Y - cy).^2 < r2;
        fvals(inside) = tol;
    end
end
