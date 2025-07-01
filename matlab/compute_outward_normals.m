function normals = compute_outward_normals(DT, bndEdges, bndEdge_tri_map)

meshP = DT.Points';
meshT = DT.ConnectivityList';
nbnd = size(bndEdges,1);
normals = zeros(2, nbnd);  % one 2D normal per edge

for ibnd = 1:nbnd
  % Get triangle and edge vertex indices
  iTri = bndEdge_tri_map(ibnd);
  triVerts = meshT(:,iTri);
  edge = bndEdges(ibnd,:);
  
  % Get physical coordinates of edge vertices
  v1 = meshP(:,edge(1));
  v2 = meshP(:,edge(2));
  
  % Edge tangent vector
  t = v2 - v1;
  t = t / norm(t);  % unit tangent
  
  % Rotate by 90° counter-clockwise to get normal
  n = [0 -1; 1 0] * t;  % ⊥ to t
  
  % Check direction: should point *out* of triangle
  % Compute triangle center and edge midpoint
  triCenter = mean(meshP(:,triVerts), 2);
  edgeMid = 0.5 * (v1 + v2);
  
  % If normal points into triangle, flip it
  if dot(n, triCenter - edgeMid) > 0
    n = -n;
  end
  
  normals(:,ibnd) = n;
end

end
