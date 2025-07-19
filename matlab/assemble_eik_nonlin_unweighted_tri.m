function [N_glb,detJ,varargout] = assemble_eik_nonlin_unweighted_tri(iTri,cu_glb,M,W,xi,...
                                                                V_abc,dPdP,...
                                                                meshT,meshP,K_glb)

% Assume: cu_glb [M*nTri × 1], V_abc [Nrs × M], dPdP [M × M × Nrs], W [Nrs × 1], detJ [nTri × 1]
% but now we only have one triangle
nTri = 1; Nrs = length(W);
% Precompute Jacobians and detJ for triangle
Pts_Ti = meshP(:,meshT(iTri,:));
J = IncidenceMatrix(Pts_Ti);
detJ = det(J);

% Reshape global coefficients on iTri
cu = reshape(cu_glb((iTri-1)*M+1:iTri*M), [M, nTri]); % [M × nTri]

% Broadcast for outer product cu_k * cu_j
cu1 = reshape(cu, [M, 1, 1, nTri]);               % [M × 1 × 1 × nTri]
cu2 = reshape(cu, [1, M, 1, nTri]);               % [1 × M × 1 × nTri]
dPdP4 = reshape(dPdP, [M, M, Nrs, 1]);            % [M × M × Nrs × 1]

% Compute gradient dot product term
G = squeeze(sum(sum(cu1 .* cu2 .* dPdP4, 1), 2)); % [Nrs × nTri]
sqrtG = sqrt(G);                                  % [Nrs × nTri]

% Expand basis functions and weights
V = permute(V_abc, [2 1]);                        % [M × Nrs]
VW = reshape(W/2, [1, Nrs]) .* V;                 % [M × Nrs]
VW = repmat(VW, [1 1 nTri]);                      % [M × Nrs × nTri]

% Tile sqrtG to match
sqrtG = reshape(sqrtG, [1, Nrs, nTri]);           % [1 × Nrs × nTri]
F = VW .* sqrtG;                                  % [M × Nrs × nTri]

% Integrate and multiply by detJ
N = squeeze(sum(F, 2));                           % [M × nTri]
N = N .* reshape(detJ, [1, nTri]);                % apply detJ
N_glb = reshape(N, [M*nTri, 1]);                  % flatten to [M*nTri × 1]


if nargout == 3
    H_glb = zeros(M*nTri);
    CU = reshape(cu_glb((iTri-1)*M+1:iTri*M), [M, nTri]);
    
    % Reshape global stiffness into 3D: [M × M × nTri]
    K_blocks = zeros(M, M, nTri);
    for t = 1:nTri
        idx = (t-1)*M+1 : t*M;
        K_blocks(:,:,t) = K_glb(idx, idx);
    end   

 
    % Precompute Jacobians and detJ for iTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    
    % 1. U1, U2 reshaped for broadcasting: [1 × M × 1 × nTri], [M × 1 × 1 × nTri]
    U1 = reshape(CU, [1, M, 1, nTri]);
    U2 = reshape(CU, [M, 1, 1, nTri]);
    
    % 2. Expand dPdP to [M × M × N × 1]
    dPdP4 = reshape(dPdP, [M, M, Nrs, 1]);
    
    % 3. Contract u^T dPdP u → [1 × 1 × N × nTri]
    quad_numer = sum(U1 .* sum(dPdP4 .* U2, 1), 2);  % → [1 × 1 × N × nTri]
    D_all = sqrt(squeeze(quad_numer));              % [N × nTri]
    
    % 4. Niq: [M × N × nTri]
    U_broadcast = reshape(CU, [1, M, nTri]);         % [1 × M × nTri]
    dPdP_batched = reshape(dPdP, [M, M, Nrs]);         % [M × M × N]
    Niq = zeros(M, Nrs, nTri);
    for t = 1:nTri
        Niq(:,:,t) = squeeze(sum(dPdP .* CU(:,t)', 2));  % sum over k
    end
    
    % 5. VW = φ_j(x_q) * w_q [N × M]
    VW = V_abc .* (W/2.0);  % constant across triangles
    
    % 6. Compute F: [M × N × nTri]
    F_all = Niq ./ permute(D_all, [3,1,2]);  % D_all: [N × nTri] → [1 × N × nTri]
    
    % 7. Compute H(:,:,iTri) = F(:,:,iTri) * VW * detJ(iTri)
    H_all = zeros(M, M, nTri);
    for t = 1:nTri
        H_all(:,:,t) = F_all(:,:,t) * VW * detJ(t) + xi * K_blocks(:,:,t);
    end
    
    for t = 1:nTri
        idx = (t-1)*M+1 : t*M;
        H_glb(idx, idx) = H_all(:,:,t);
    end

    varargout{1} = H_glb;
end

end
