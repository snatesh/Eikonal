function xy = mapToPhysicalEdge(edgeID, r, s, P, J)
% Maps reference edge points to physical triangle using affine map
    switch edgeID
        case 1  % left: (0,0) to (0,1)
            rl = r;
            sl = s;
        case 2  % bottom: (0,0) to (1,0)
            rl = r;
            sl = s;
        case 3  % hypotenuse: (1,0) to (0,1)
            rl = r;
            sl = s;
    end
    % Local reference → physical coordinates
    xy = P(:,1) + J * [rl'; sl'];
end
