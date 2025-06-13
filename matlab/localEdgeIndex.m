function idx = localEdgeIndex(tri, edge)

% tri is a 1x3 vector of triangle vertex indices
edges = [tri([1,2]); tri([2,3]); tri([3,1])];
for i = 1:3
    if isequal(sort(edges(i,:)), sort(edge))
        idx = i;
        return
    end
end
error('Edge not found');

end

