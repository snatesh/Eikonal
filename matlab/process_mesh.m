function [meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map,bndTris] = process_mesh(DT)

meshT = DT.ConnectivityList;
meshP = DT.Points';

% Extract all edges from triangles
Edges = [meshT(:, [1,2]); meshT(:, [2,3]); meshT(:, [3,1])];

% Sort edges so that (i,j) and (j,i) are considered the same
Edges = sort(Edges, 2);

% Count occurrences of each edge
[uniqueEdges, ~, ic] = unique(Edges, 'rows');
counts = accumarray(ic, 1);

% classify edges 
% edges that appear twice are interior
% edges that appear once are on boundary
intEdges = uniqueEdges(counts == 2, :);
bndEdges = uniqueEdges(counts == 1,:);
nintEdge = size(intEdges,1);

sharedEdge_tri_map = zeros(nintEdge,2);


for iedge = 1:nintEdge

    sharedTris = edgeAttachments(DT,intEdges(iedge,:));
    sharedEdge_tri_map(iedge,:) = sharedTris{1};

end

% sort it to match vertex index ordering in intEdges
sharedEdge_tri_map = sort(sharedEdge_tri_map,2);

% find triangles with boundary edge
triangleIDs = repmat((1:size(meshT,1))', 3, 1);  
[~, loc] = ismember(Edges, bndEdges, 'rows');
isBoundaryEdge = loc > 0;
bndTris = unique(triangleIDs(isBoundaryEdge));

end
