function [K_glb,Bint_glb,Bdirch_glb,...
          F_glb,G_glb,...
          nintEdge,nbndEdge,...
          nbndTri] = assemble_poisson_pwc(n,R,S,W,...
                                          Rl,Sl,Rb,Sb,Rh,Sh, ...
                                          Rl_flip,Sl_flip,...
                                          Rb_flip,Sb_flip,...
                                          Rh_flip,Sh_flip,wleg,...
                                          V_abc,dxP,dyP,Vl,Vb,Vh,...
                                          Vl_flip,Vb_flip,Vh_flip,...
                                          rhs,udirch,DT)


% get mesh points, connectivity list
% as well as interior/bndry edges
% and map taking shared edge to the two triangles which share it
[meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map,bndTris] = process_mesh(DT);

% number of triangles in mesh
nTri = size(meshT,1);
% total number of polynomials
M = n*(n+1)/2;
% number of boundary and interior edges
nintEdge = size(intEdges,1);
nbndEdge = size(bndEdges,1);
nbndTri = size(bndTris,1);
nleg = size(wleg,1);

% global stifness
K_glb = zeros(M*nTri,M*nTri);
% global dirichlet
Bdirch_glb = zeros(nleg*nbndEdge,M*nTri);
% global inter-element continuity
% there are nintEdge row blocks, each with M rows and M*nTri cols
Bint_glb = zeros(nleg*nintEdge,M*nTri);
F_glb = zeros(M*nTri,1);
G_glb = zeros(nleg*nbndEdge,1);


% now we loop over triangles 
% and assemble stiffness and dirichlet matrices
ibndtri = 0;
for iTri = 1:nTri

    % get tri pts, incidence matrix and edge jacobians
    Pts_Ti = meshP(:,meshT(iTri,:));

    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    invJ = inv(J);
    lenE = zeros(3,1);
    lenE(1) = norm(Pts_Ti(:,2)-Pts_Ti(:,1));
    lenE(2) = norm(Pts_Ti(:,3)-Pts_Ti(:,2));
    lenE(3) = norm(Pts_Ti(:,1)-Pts_Ti(:,3));
    % quadrature points mapped to phyiscal tri
    XY = (J * [R,S]' + Pts_Ti(:,1))';
    X = XY(:,1);
    Y = XY(:,2);
    % boundary quadrature points mapped to physical tri
    XYl = (J * [Rl,Sl]' + Pts_Ti(:,1))';
    XYb = (J * [Rb,Sb]' + Pts_Ti(:,1))';
    XYh = (J * [Rh,Sh]' + Pts_Ti(:,1))';
    XYl_flip = (J * [Rl_flip,Sl_flip]' + Pts_Ti(:,1))';
    XYb_flip = (J * [Rb_flip,Sb_flip]' + Pts_Ti(:,1))';
    XYh_flip = (J * [Rh_flip,Sh_flip]' + Pts_Ti(:,1))';


    % local interior stiffness
    K = zeros(M,M);
    for i = 1:M
        dxP_i = dxP(:,i);
        dyP_i = dyP(:,i);
        gradP_i = invJ'*[dxP_i';dyP_i'];
        dxP_i = gradP_i(1,:)';
        dyP_i = gradP_i(2,:)';
        for j = 1:M
            dxP_j = dxP(:,j);
            dyP_j = dyP(:,j);
            gradP_j = invJ'*[dxP_j';dyP_j'];
            dxP_j = gradP_j(1,:)';
            dyP_j = gradP_j(2,:)';
            dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
            K(i,j) = (W/2)'*dPidPj*detJ;
        end
    end
    % save to global stiffness
    K_glb((iTri-1)*M+1:M*iTri,(iTri-1)*M+1:M*iTri) = K;

    % local load (rhs)
    F = zeros(M,1);
    for i = 1:M
        P_if = V_abc(:,i).*rhs(X,Y);
        F(i) = (W/2)'*P_if*detJ;
    end
    % save to global load
    F_glb((iTri-1)*M+1:M*iTri) = F;

    % determine how many boundary edges in triangle
    nbndEdge_T = 0;
    if ismember(iTri,bndTris)
        % get edges of triangle, sorted to match bndry/int Edge lookup
        edgesT = sort([meshT(iTri, [1,2]);...
                       meshT(iTri, [2,3]);...
                       meshT(iTri, [3,1])],2);
        % for each edge of triangle
        for iedge = 1:3
            % test if it is a boundary edge
            for ibnd = 1:nbndEdge
                if all(edgesT(iedge,:) == bndEdges(ibnd,:))
                    nbndEdge_T = nbndEdge_T+1;
                end
            end
        end
    end

    % local boundary stiffness
    Bdirch = zeros(nleg*nbndEdge_T,M);
    disp(nbndEdge_T);
    % boundary load (dirichlet condition)
    G = zeros(nleg*nbndEdge_T,1);
    % if triangle is a boundary triangle, we compute
    % the needed entries of Bdirch and G as well 
    % as the global index for global matrix/vec
    if ismember(iTri,bndTris)
        ibndtri = ibndtri+1;

        % get edges of triangle, sorted to match bndry/int Edge lookup
        edgesT = sort([meshT(iTri, [1,2]);...
                       meshT(iTri, [2,3]);...
                       meshT(iTri, [3,1])],2);
        % for each edge of triangle
        iibnd = 0;
        flips = [false false false];
        for iedge = 1:3
            % get endpts
            ve = meshP(:,edgesT(iedge,:));
            % get parametric coords of endpt
            rve = J\(ve-Pts_Ti(:,1));
            % label them as left, bottom or hypotenuse edge
            left = false; bot = false; hyp = false;
            if (all(rve(:,1) == [0;0]) && all(rve(:,2) == [1;0]))
                flips(iedge) = false;
                bot = true;
            elseif (all(rve(:,1) == [1;0]) && all(rve(:,2) == [0;0]))
                flips(iedge) = true;
                bot = true;
            elseif (all(rve(:,1) == [1;0]) && all(rve(:,2) == [0;1]))
                flips(iedge) = false;
                hyp = true;
            elseif (all(rve(:,1) == [0;1]) && all(rve(:,2) == [1;0]))
                flips(iedge) = true;
                hyp = true;
            elseif (all(rve(:,1) == [0;1]) && all(rve(:,2) == [0;0]))
                flips(iedge) = false;
                left = true;
            elseif (all(rve(:,1) == [0;0]) && all(rve(:,2) == [0;1]))
                flips(iedge) = true;
                left = true;
            end

            % now test if it is a boundary edge
            
            for ibnd = 1:nbndEdge
                % if it is, compute and store
                % local dirichlet contribution
                if all(edgesT(iedge,:) == bndEdges(ibnd,:))
                    V = 0; g = 0; % dirichlet values on bndry
                    iibnd = iibnd+1;
                    if bot
                        if ~flips(iedge)
                            g = udirch(XYb(:,1),XYb(:,2));
                            V = Vb;
                        else
                            g = udirch(XYb_flip(:,1),XYb_flip(:,2));
                            V = Vb_flip;
                        end
                    elseif hyp
                        if ~flips(iedge)
                            g = udirch(XYh(:,1),XYh(:,2));
                            V = Vh;
                        else
                            g = udirch(XYh_flip(:,1),XYh_flip(:,2));
                            V = Vh_flip;
                        end
                    elseif left
                        if ~flips(iedge)
                            g = udirch(XYl(:,1),XYl(:,2));
                            V = Vl;
                        else
                            g = udirch(XYl_flip(:,1),XYl_flip(:,2));
                            V = Vl_flip;
                        end
                    end
                    
                    for i = 1:nleg
                        Bdirch(i+(iibnd-1)*nleg,:) = V(i,:);
                        G(i+(iibnd-1)*nleg) = g(i);
                    end
                end
            end
        end
        % now save to global dirichlet matrix
        Bdirch_glb((ibndtri-1)*nleg*nbndEdge_T+1:ibndtri*nleg*nbndEdge_T,...
                   (iTri-1)*M+1:M*iTri) = Bdirch;
        % save dirichlet condition to global Gdirch
        G_glb((ibndtri-1)*nleg*nbndEdge_T+1:ibndtri*nleg*nbndEdge_T) = G;        
    end
end

% now we loop over shared edges 
% to assemble inter-element continuity matrix
for iedge = 1:nintEdge

    % get triangles sharing the edge
    shared_tris = sharedEdge_tri_map(iedge,:);
    tri1 = shared_tris(1);
    tri2 = shared_tris(2);

    % get vertices of each tri
    Pts_T1 = meshP(:,meshT(tri1,:));
    Pts_T2 = meshP(:,meshT(tri2,:));
    % get parametric map for each tri
    J1 = IncidenceMatrix(Pts_T1);
    J2 = IncidenceMatrix(Pts_T2);
    % get shared edge endpts
    ve = meshP(:,intEdges(iedge,:));
    % get parametric coords of endpt on each tri
    rve1 = J1\(ve-Pts_T1(:,1));
    rve2 = J2\(ve-Pts_T2(:,1));
    % get the correct precomputed basis functions corresponding
    % to the parametric shared edge in each tri
    tol = 1e-13;
    edge1 = 0; edge2 = 0;
    rleg1 = 0; sleg1 = 0;
    rleg2 = 0; sleg2 = 0;
    % bottom
    if ((norm(rve1(:,1) - [0;0]) < tol) && (norm(rve1(:,2) - [1;0]) < tol)) || ... 
       ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) -[0;0]) < tol))
        V1 = Vb;
        edge1 = 2;
        rleg1 = Rb; sleg1 = Sb;
    % hyp
    elseif ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) - [0;1]) < tol)) || ...
           ((norm(rve1(:,1) - [0;1])< tol) && (norm(rve1(:,2) - [1;0]) < tol))
        V1 = Vh;
        edge1 = 3;
        rleg1 = Rh; sleg1 = Sh;
    % left
    else
        V1 = Vl;
        edge1 = 1;
        rleg1 = Rl; sleg1 = Sl;
    end
    % bottom
    if ((norm(rve2(:,1) - [0;0]) < tol) && (norm(rve2(:,2) - [1;0]) < tol)) || ...
       ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;0]) < tol))
        V2 = Vb;
        edge2 = 2;
        rleg2 = Rb; sleg2 = Sb;
    % hyp
    elseif ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;1]) < tol)) || ...
           ((norm(rve2(:,1) - [0;1]) < tol) && (norm(rve2(:,2) - [1;0]) < tol))
        V2 = Vh;
        edge2 = 3;
        rleg2 = Rh; sleg2 = Sh;
    % left
    else
        V2 = Vl;
        edge2 = 1;
        rleg2 = Rl; sleg2 = Sl;
    end

    % Align the bases
    [V1a, V2a] = getAlignedEdgeBases(V1, rleg1, sleg1, V2, rleg2, sleg2, Pts_T1, Pts_T2, edge1, edge2);

    % now save to global inter-element continuity matrix
    Bint_glb((iedge-1)*nleg+1:nleg*iedge,(tri1-1)*M+1:M*tri1) = V1a; 
    Bint_glb((iedge-1)*nleg+1:nleg*iedge,(tri2-1)*M+1:M*tri2) = -V2a; 


    %xy1 = edgePhysicalPoints(Pts_T1, edge1, Sl);
    %xy2 = edgePhysicalPoints(Pts_T2, edge2, Sl);
    %if (norm(xy1-xy2) > 1e-13)
    %  xy2 = edgePhysicalPoints(Pts_T2,edge2,1-Sl);
    %  disp(norm(xy1 - xy2))
    %else
    %  disp(norm(xy1-xy2))
    %end
    

end

end
