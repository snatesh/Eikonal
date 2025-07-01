function poiss_scale = scale_poisson(cu_glb,rhs,xi,M,R,S,W,dxP,dyP,meshT,meshP)

nTri = size(meshT,1); Nrs = length(W);

intNumer = 0;
intDenom = 0;

for iTri = 1:nTri
    % get tri pts, incidence matrix and edge jacobians
    Pts_Ti = meshP(:,meshT(iTri,:));

    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    invJ = inv(J);
    % quadrature points mapped to phyiscal tri
    XY = (J * [R,S]' + Pts_Ti(:,1))';
    X = XY(:,1);
    Y = XY(:,2);
    cu_T = cu_glb((iTri-1)*M+1:iTri*M);
    ur_T = dxP * cu_T;
    us_T = dyP * cu_T;
    gradu_T = invJ'*[ur_T';us_T'];
    ux_T = gradu_T(1,:)';
    uy_T = gradu_T(2,:)';
    grad_mag = sqrt(ux_T.^2 + uy_T.^2);
    f_T = rhs(X,Y);
    
    intNumer = intNumer + (W/2)'*((grad_mag + xi*f_T).*f_T)*detJ;
    intDenom = intDenom + (W/2)'*((grad_mag + xi*f_T).^2)*detJ;   
end

poiss_scale = intNumer / intDenom;

end
