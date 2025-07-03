function dt = delaunay_unit_square(N)
    % Generate N x N grid in [0,1]^2 with (N+1)^2 points
    x = linspace(-1, 1, N+1);
    [X, Y] = meshgrid(x, x);
    p = [X(:), Y(:)];

    % Generate square elements and split into two triangles
    t = [];
    for i = 1:N
        for j = 1:N
            % Indices of square corners
            n0 = (i-1)*(N+1) + j;
            n1 = n0 + 1;
            n2 = n0 + (N+1);
            n3 = n2 + 1;

            % Split square into two triangles
            t1 = [n0, n1, n3];
            t2 = [n0, n3, n2];
            t = [t; t1; t2];
        end
    end

    % Optional: create Delaunay triangulation object
    dt = delaunayTriangulation(p);
end
