function visualize_causal_propagation(DT, x_anchor)
    % Inputs:
    %   DT       - delaunayTriangulation object
    %   x_anchor - 1x2 anchor point [x, y], must lie on the boundary

    % Step 1: Find triangle(s) whose edge contains x_anchor
    edges = DT.freeBoundary();
    pts = DT.Points;
    tris = DT.ConnectivityList;
    
    % Find closest edge midpoint to x_anchor
    midpoints = (pts(edges(:,1),:) + pts(edges(:,2),:)) / 2;
    [~, idx] = min(vecnorm(midpoints - x_anchor, 2, 2));
    anchor_edge = edges(idx, :);

    % Find triangle(s) sharing this edge
    attached_tris = edgeAttachments(DT, anchor_edge);
    queue = attached_tris{1};  % Initial "known" triangle(s)

    % Initialize sets
    numTris = size(tris, 1);
    status = strings(numTris, 1);     % "Unknown", "Known", "Trial"
    status(:) = "Unknown";
    proxy = inf(numTris, 1);          % ∫ u value (for priority)
    parent = zeros(numTris, 1);       % Store source triangle (for arrows)
    
    % Assign random proxy to seed triangles and mark as Known
    for it = 1:length(queue)
        t = queue(it);
        status(t) = "Known";
        proxy(t) = rand();  % simulate ∫ u
    end

    % Frontier: list of triangles eligible to update (naive FIFO)
    frontier = queue;

    while ~isempty(frontier)
        % Select triangle in frontier with minimum proxy
        [~, idx_min] = min(proxy(frontier));
        T_cur = frontier(idx_min);
        frontier(idx_min) = [];  % remove from frontier
    
        % Get neighbors of current triangle
        nbrs = neighbors(DT, T_cur);
        nbrs = nbrs(~isnan(nbrs));
    
        for inbr = 1:length(nbrs)
            t_nbr = nbrs(inbr);
    
            if status(t_nbr) == "Unknown"
                % Simulate causal update
                status(t_nbr) = "Known";
                proxy(t_nbr) = proxy(T_cur) + rand();  % causally increasing
                parent(t_nbr) = T_cur;
                frontier(end+1) = t_nbr;
            end
        end
    end

    % === Visualization ===
    figure; hold on; axis equal;
    triplot(DT, 'Color', [0.8 0.8 0.8]);

    % Plot arrows for causal propagation
    centers = incenter(DT);  % triangle centers
    for t = 1:numTris
        if parent(t) > 0
            p = parent(t);
            c1 = centers(p, :);
            c2 = centers(t, :);
            quiver(c1(1), c1(2), c2(1)-c1(1), c2(2)-c1(2), 0, ...
                'MaxHeadSize', 0.5, 'Color', 'b', 'LineWidth', 1.2);
        end
    end

    % Plot anchor point
    plot(x_anchor(1), x_anchor(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    title('Causal Propagation Graph');
end
