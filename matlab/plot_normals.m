function plot_normals(DT, bndEdges, normals)
figure;
triplot(DT, 'k'); hold on; axis equal;

meshP = DT.Points';  % [2 x nPts]
nbnd = size(bndEdges,1);

for ibnd = 1:nbnd
    edge = bndEdges(ibnd,:);
    v1 = meshP(:,edge(1));
    v2 = meshP(:,edge(2));
    
    % Midpoint of edge
    mid = 0.5 * (v1 + v2);
    
    % Normal vector
    n = normals(:,ibnd);
    
    % Scale for visibility
    scale = 0.1 * norm(v2 - v1);
    
    % Plot normal
    quiver(mid(1), mid(2), scale * n(1), scale * n(2), ...
           0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2);
end

title('Boundary Normals Overlaid on Triangulation');
end
