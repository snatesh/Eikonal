function [K_glb,Bint_glb,Bdirch_glb,...
          F_glb,G_glb,...
          nintEdge,nbndEdge,nbndTri] = assemble_poisson(n,R,S,W,...
                                                Rl,Sl,Rb,Sb,Rh,Sh,...
                                                Rl_flip,Sl_flip,...
                                                Rb_flip,Sb_flip,...
                                                Rh_flip,Sh_flip,...
                                                V_abc,dxP,dyP,Vl,Vb,Vh,...
                                                Vl_flip,Vb_flip,Vh_flip,...
                                                intVlVl,intVbVb,intVhVh,...
                                                intVlVl_flip,...
                                                intVbVb_flip,...
                                                intVhVh_flip,wleg,...
                                                rhs,udirch,DT)


% get mesh points, connectivity list
% as well as interior/bndry edges
% and map taking shared edge to the two triangles which share it
[meshP,meshT,intEdges,bndEdges,sharedEdge_tri_map,bndTris] = process_mesh(DT);

disp(bndTris)
% number of triangles in mesh
nTri = size(meshT,1);
% total number of polynomials
M = n*(n+1)/2;
% number of boundary and interior edges
nintEdge = size(intEdges,1);
nbndEdge = size(bndEdges,1);
nbndTri = size(bndTris,1);

% global stifness
K_glb = zeros(M*nTri,M*nTri);
% global dirichlet
Bdirch_glb = zeros(M*nbndTri,M*nTri);
% global inter-element continuity
% there are nintEdge row blocks, each with M rows and M*nTri cols
Bint_glb = zeros(M*nintEdge,M*nTri);
F_glb = zeros(M*nTri,1);
G_glb = zeros(M*nbndTri,1);


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

    % if triangle is a boundary triangle, we compute
    % the needed entries of Bdirch and G
    if ismember(iTri,bndTris)
        ibndtri = ibndtri+1;
        % local boundary stiffness
        Bdirch = zeros(M,M);
        % boundary load (dirichlet condition)
        G = zeros(M,1);
        % get edges of triangle, sorted to match bndry/int Edge lookup
        edgesT = sort([meshT(iTri, [1,2]);...
                       meshT(iTri, [2,3]);...
                       meshT(iTri, [3,1])],2);
        % for each edge of triangle
        for iedge = 1:3
            % get endpts
            ve = meshP(:,edgesT(iedge,:));
            % get parametric coords of endpt
            rve = J\(ve-Pts_Ti(:,1));
            % label them as left, bottom or hypotenuse edge
            left = false; bot = false; hyp = false; flip = 1;
            if (all(rve(:,1) == [0;0]) && all(rve(:,2) == [1;0]))
               bot = true; 
               flip = 1;
            elseif (all(rve(:,1) == [1;0]) && all(rve(:,2) == [0;0]))
                bot = true;
                flip = -1;
            elseif (all(rve(:,1) == [1;0]) && all(rve(:,2) == [0;1]))
                hyp = true;
                flip = 1;
            elseif (all(rve(:,1) == [0;1]) && all(rve(:,2) == [1;0]))
                hyp = true;
                flip = -1;
            elseif (all(rve(:,1) == [0;1]) && all(rve(:,2) == [0;0]))
                left = true;
                flip = 1;
            elseif (all(rve(:,1) == [0;0]) && all(rve(:,2) == [0;1]))
                left = true;
                flip = -1;
            end
            % test if it is a boundary edge
            for ibnd = 1:nbndEdge
                % if it is, compute and store
                % local dirichlet contribution
                if all(edgesT(iedge,:) == bndEdges(ibnd,:))
                    intVV = 0; % needed for boundary stiffness
                    V = 0; g = 0; % dirichlet values on bndry
                    if bot
                        intVV = intVbVb;
                        if flip == 1
                          g = udirch(XYb(:,1),XYb(:,2));
                          V = Vb;
                        elseif flip == -1
                          g = udirch(XYb_flip(:,1),XYb_flip(:,2));
                          V = Vb_flip;
                        end

                    elseif hyp
                        intVV = intVhVh;
                        if flip == 1;
                          g = udirch(XYh(:,1),XYh(:,2));
                          V = Vh;
                        elseif flip == -1;
                          g = udirch(XYh_flip(:,1),XYh_flip(:,2));
                          V = Vh_flip;
                        end    
                    elseif left
                        intVV = intVlVl;
                        if flip == 1;
                          g = udirch(XYl(:,1),XYl(:,2));
                          V = Vl;
                        elseif flip == -1;
                          g = udirch(XYl_flip(:,1),XYl_flip(:,2));
                          V = Vl_flip;
                        end

                    end

                    for i = 1:M
                        for j = 1:M
                            Bdirch(i,j) = Bdirch(i,j) + intVV(i,j)*lenE(iedge);
                        end
                        G(i) = G(i) + wleg'*(V(:,i).*g)*lenE(iedge);
                    end
                end
            end
        end
        % now save to global dirichlet matrix
        Bdirch_glb((ibndtri-1)*M+1:M*ibndtri,(iTri-1)*M+1:M*iTri) = Bdirch;
        % save dirichlet condition to global Gdirch
        G_glb((ibndtri-1)*M+1:M*ibndtri) = G;
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
    % get edge length
    lenE = norm(ve(:,2)-ve(:,1));
    % get parametric coords of endpt on each tri
    rve1 = J1\(ve-Pts_T1(:,1));
    rve2 = J2\(ve-Pts_T2(:,1));
    % get the correct precomputed line integral corresponding
    % to the parametric shared edge in each tri
    intVV1 = 0; intVV2 = 0;

    tol = 1e-13;
    edge1 = 0; edge2 = 0;
    rleg1 = 0; sleg1 = 0;
    rleg2 = 0; sleg2 = 0;
    % bottom
    if ((norm(rve1(:,1) - [0;0]) < tol) && (norm(rve1(:,2) - [1;0]) < tol)) || ... 
       ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) -[0;0]) < tol))
        intVV1 = intVbVb;
        edge1 = 2;
        rleg1 = Rb; sleg1 = Sb;
    % hyp
    elseif ((norm(rve1(:,1) - [1;0]) < tol) && (norm(rve1(:,2) - [0;1]) < tol)) || ...
           ((norm(rve1(:,1) - [0;1])< tol) && (norm(rve1(:,2) - [1;0]) < tol))
        intVV1 = intVhVh;
        rleg1 = Rh; sleg1 = Sh;
        edge1 = 3;
    % left
    else
        intVV1 = intVlVl;
        rleg1 = Rl; sleg1 = Sl;
        edge1 = 1;
    end
    % bottom
    if ((norm(rve2(:,1) - [0;0]) < tol) && (norm(rve2(:,2) - [1;0]) < tol)) || ...
       ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;0]) < tol))
        intVV2 = intVbVb;
        rleg2 = Rb; sleg2 = Sb;
        edge2 = 2;
    % hyp
    elseif ((norm(rve2(:,1) - [1;0]) < tol) && (norm(rve2(:,2) - [0;1]) < tol)) || ...
           ((norm(rve2(:,1) - [0;1]) < tol) && (norm(rve2(:,2) - [1;0]) < tol))
        intVV2 = intVhVh;
        rleg2 = Rh; sleg2 = Sh;
        edge2 = 3;
    % left
    else
        intVV2 = intVlVl;
        rleg2 = Rl; sleg2 = Sl;
        edge2 = 1;
    end
    
    % Align the bases

    [intVV1a,intVV2a] = getAlignedEdgeBasesInt(intVV1, rleg1, sleg1, intVV2, ...
                                               rleg2, sleg2, Pts_T1, Pts_T2, edge1, edge2);
    %% bottom
    %if (all(rve1(:,1) == [0;0]) && all(rve1(:,2) == [1;0]))
    %    intVV1 = intVbVb;
    %elseif (all(rve1(:,1) == [1;0]) && all(rve1(:,2) == [0;0]))
    %    intVV1 = -intVbVb;
    %% hyp
    %elseif (all(rve1(:,1) == [1;0]) && all(rve1(:,2) == [0;1]))
    %    intVV1 = intVhVh;
    %elseif (all(rve1(:,1) == [0;1]) && all(rve1(:,2) == [1;0]))
    %    intVV1 = -intVhVh;
    %% left
    %elseif (all(rve1(:,1) == [0;1]) && all(rve1(:,2) == [0;0]))
    %    intVV1 = intVlVl;
    %elseif (all(rve1(:,1) == [0;0]) && all(rve1(:,2) == [0;1]))
    %    intVV1 = -intVlVl;
    %end
    %intVV2 = 0;
    %% bottom
    %if (all(rve2(:,1) == [0;0]) && all(rve2(:,2) == [1;0]))
    %    intVV2 = intVbVb;
    %elseif (all(rve2(:,1) == [1;0]) && all(rve2(:,2) == [0;0]))
    %    intVV2 = -intVbVb;
    %% hyp
    %elseif (all(rve2(:,1) == [1;0]) && all(rve2(:,2) == [0;1]))
    %    intVV2 = intVhVh;
    %elseif (all(rve2(:,1) == [0;1]) && all(rve2(:,2) == [1;0]))
    %    intVV2 = -intVhVh;
    %% left
    %elseif (all(rve2(:,1) == [0;1]) && all(rve2(:,2) == [0;0]))
    %    intVV2 = intVlVl;
    %elseif (all(rve2(:,1) == [0;0]) && all(rve2(:,2) == [0;1]))
    %    intVV2 = -intVlVl;
    %end
    %Bint1 = zeros(M,M);
    %Bint2 = zeros(M,M);
    %for i = 1:M
    %    for j = 1:M
    %        Bint1(i,j) = intVV1(i,j)*lenE;
    %        Bint2(i,j) = intVV2(i,j)*lenE;
    %    end
    %end
    % now save to global inter-element continuity matrix
    Bint_glb((iedge-1)*M+1:M*iedge,(tri1-1)*M+1:M*tri1) = intVV1a; 
    Bint_glb((iedge-1)*M+1:M*iedge,(tri2-1)*M+1:M*tri2) = intVV2a; 

end

end
