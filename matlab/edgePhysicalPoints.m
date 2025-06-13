function xy_phys = edgePhysicalPoints(Pts_T, whichEdge, rleg)
    % Input:
    % Pts_T    - (2×3) array: triangle vertex positions [x1 x2 x3; y1 y2 y3]
    % whichEdge - integer (1: left, 2: bottom, 3: hypotenuse)
    % rleg     - (1×N) or (N×1) array of parametric edge points ∈ [0,1]
    
    rleg = rleg(:)';  % ensure row vector
    N = length(rleg);
    
    switch whichEdge
        case 1  % left: from v1 to v3
            lambdas = [1 - rleg; zeros(1,N); rleg];
        case 2  % bottom: from v1 to v2
            lambdas = [1 - rleg; rleg; zeros(1,N)];
        case 3  % hypotenuse: from v2 to v3
            lambdas = [zeros(1,N); 1 - rleg; rleg];
        otherwise
            error('Invalid edge index: must be 1, 2, or 3');
    end

    xy_phys = Pts_T * lambdas;  % (2×3) * (3×N) = (2×N)
end
