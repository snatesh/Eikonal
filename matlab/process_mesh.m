function [meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map,bndEdge_tri_map,normals] = process_mesh(DT)

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
nbndEdge = size(bndEdges,1);
bndEdge_tri_map = zeros(nbndEdge,1);

for i = 1:nbndEdge
    edge = bndEdges(i,:);
    tris = edgeAttachments(DT, edge);  % Only one triangle for boundary edge
    bndEdge_tri_map(i) = tris{1};
end

sharedEdge_tri_map = zeros(nintEdge,2);


for iedge = 1:nintEdge

    sharedTris = edgeAttachments(DT,intEdges(iedge,:));
    sharedEdge_tri_map(iedge,:) = sharedTris{1};

end

% sort it to match vertex index ordering in intEdges
sharedEdge_tri_map = sort(sharedEdge_tri_map,2);

% get boundary normals
normals = compute_outward_normals(DT, bndEdges, bndEdge_tri_map);

end
