function animate_causal_propagation(DT, x_anchor)
    % Inputs:
    %   DT       - delaunayTriangulation object
    %   x_anchor - [x, y] anchor point on boundary

    % Initialization
    edges = DT.freeBoundary();
    pts = DT.Points;
    tris = DT.ConnectivityList;
    numTris = size(tris, 1);
    centers = incenter(DT);  % triangle centers

    % Step 1: Find seed triangle from boundary
    midpoints = (pts(edges(:,1),:) + pts(edges(:,2),:)) / 2;
    [~, idx] = min(vecnorm(midpoints - x_anchor, 2, 2));
    anchor_edge = edges(idx,:);
    seed_tris = edgeAttachments(DT, anchor_edge);
    frontier = seed_tris{1};

    figure; hold on;
    triplot(DT, 'Color', [0.8 0.8 0.8]);  % base mesh lines
    plot(x_anchor(1), x_anchor(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);  % anchor
    xmin = min(pts(:,1));
    xmax = max(pts(:,1));
    ymin = min(pts(:,2));
    ymax = max(pts(:,2));
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    axis equal manual;


    % Initialize data
    status = strings(numTris, 1);       % "Unknown", "Known"
    status(:) = "Unknown";
    proxy = inf(numTris, 1);            % simulated ∫ u values
    parent = zeros(numTris, 1);         % for arrows
    update_order = zeros(numTris, 1);   % for coloring
    t_update = 1;

    % Seed triangle(s)
    for iseed = 1:length(frontier)
        t = frontier(iseed);
        status(t) = "Known";
        proxy(t) = rand();
        update_order(t) = t_update;
        t_update = t_update + 1;
        % === Animate seed triangle ===
        patch('Faces', tris(t,:), 'Vertices', pts, ...
              'FaceColor', 'flat', ...
              'FaceVertexCData', update_order(t), ...
              'EdgeColor', 'k');
        drawnow;
        pause(0.1);
    end


    % === Causal Marching Loop ===
    while ~isempty(frontier)
        [~, idx_min] = min(proxy(frontier));
        T_cur = frontier(idx_min);
        frontier(idx_min) = [];

        nbrs = neighbors(DT, T_cur);
        nbrs = nbrs(~isnan(nbrs));

        for inbr = 1:length(nbrs)
            t_nbr = nbrs(inbr);

            if status(t_nbr) == "Unknown"
                status(t_nbr) = "Known";
                proxy(t_nbr) = proxy(T_cur) + rand();  % simulate ∫ u
                parent(t_nbr) = T_cur;
                frontier(end+1) = t_nbr;

                update_order(t_nbr) = t_update;
                t_update = t_update + 1;

                % === Animate update ===
                patch('Faces', tris(t_nbr,:), 'Vertices', pts, ...
                      'FaceColor', 'flat', ...
                      'FaceVertexCData', update_order(t_nbr), ...
                      'EdgeColor', 'k'); hold on;
                % Arrow from parent to current triangle
                c1 = centers(parent(t_nbr), :);
                c2 = centers(t_nbr, :);
                quiver(c1(1), c1(2), c2(1)-c1(1), c2(2)-c1(2), 0, ...
                       'MaxHeadSize', 0.5, 'Color', 'w', 'LineWidth', 1);
                drawnow;
                pause(0.1);
            end
        end
    end

    % === Final Color Map ===
    figure; hold on; axis equal; title('Update Order Heatmap');
    colormap(parula);
    patch('Faces', tris, 'Vertices', pts, ...
          'FaceVertexCData', accumarray(repelem((1:numTris)',3), update_order(tris(:))) ./ 3, ...
          'FaceColor', 'flat', 'EdgeColor', [0.7 0.7 0.7]);
    colorbar;
    plot(x_anchor(1), x_anchor(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
end
