function [V1_aligned, V2_aligned] = getAlignedEdgeBases(V1, rleg1, sleg1, V2, rleg2, sleg2, P1, P2, edge1, edge2)
% Align basis functions V1, V2 from two triangles T1 and T2 on a shared edge
% Inputs:
%   - V1, V2       : Basis functions evaluated at edge quadrature points on T1 and T2
%   - rleg1,sleg1  : Reference edge coordinates on T1's edge (edge1)
%   - rleg2,sleg2  : Reference edge coordinates on T2's edge (edge2)
%   - P1, P2       : 2x3 vertex matrices of triangles T1 and T2
%   - edge1, edge2 : Which edge (1=left, 2=bottom, 3=hypotenuse)
% Outputs:
%   - V1_aligned, V2_aligned : basis evaluations at matching physical points

    % Triangle-to-physical maps
    J1 = IncidenceMatrix(P1);
    J2 = IncidenceMatrix(P2);

    % Get physical coordinates of quadrature points on edge
    xy1 = mapToPhysicalEdge(edge1, rleg1, sleg1, P1, J1);
    xy2 = mapToPhysicalEdge(edge2, rleg2, sleg2, P2, J2);
    xy2_flip = mapToPhysicalEdge(edge2, flip(rleg2), flip(sleg2), P2, J2);

    % Compute alignment error
    err_no_flip = max(vecnorm(xy1 - xy2, 2, 1));
    err_flip    = max(vecnorm(xy1 - xy2_flip, 2, 1));

    % Decide if V2 needs to be flipped to match V1
    V1_aligned = V1;
    if err_flip < err_no_flip
        V2_aligned = flip(V2, 1); % flip quadrature order
    else
        V2_aligned = V2;
    end
end
