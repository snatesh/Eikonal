function [N_glb,varargout] = assemble_eik_nonlin_unweighted(cu_glb,M,W,xi,...
                                                            V_abc,dPdP,...
                                                            meshT,meshP,K_glb)
nTri = size(meshT,1);
N_glb = zeros(nTri*M,1); Nrs = length(W);
for iTri = 1:nTri
  cu = cu_glb((iTri-1)*M+1:iTri*M);
  N = zeros(M,1);
  % get tri pts, incidence matrix and area jacobian
  Pts_Ti = meshP(:,meshT(iTri,:));
  J = IncidenceMatrix(Pts_Ti);
  detJ = det(J);
  for i = 1:M
      P_i = V_abc(:,i);
      sum1 = zeros(Nrs,1);
      for j = 1:M
          for k = 1:M
              dPkdPj = reshape(dPdP(k,j,:),Nrs,1);
              sum1 = sum1 + cu(k)*cu(j).*dPkdPj;
          end
      end
      sum1 = sqrt(sum1);
      N(i) = (W/2)'*(P_i.*sum1)*detJ;
  end
  N_glb((iTri-1)*M+1:iTri*M) = N;
end

if nargout == 2
    H_glb = zeros(M*nTri);
    for iTri = 1:nTri
      cu = cu_glb((iTri-1)*M+1:iTri*M);
      K = K_glb((iTri-1)*M+1:M*iTri,(iTri-1)*M+1:M*iTri);
      H = zeros(M);
      Pts_Ti = meshP(:,meshT(iTri,:));
      J = IncidenceMatrix(Pts_Ti);
      detJ = det(J);
      for i = 1:M
          for j = 1:M
              P_j = V_abc(:,j);
              sum1 = zeros(Nrs,1);
              sum2 = zeros(Nrs,1);
              for k = 1:M
                  dPidPk = reshape(dPdP(i,k,:),Nrs,1);
                  sum1 = sum1 + cu(k)*dPidPk;
                  for l = 1:M
                      dPkdPl = reshape(dPdP(k,l,:),Nrs,1);
                      sum2 = sum2 + cu(k)*cu(l).*dPkdPl;
                  end
              end
              sum1 = P_j.*sum1;
              sum2 = sqrt(sum2);
              H(i,j) = (W/2)'*(sum1./sum2)*detJ + xi*K(i,j);
          end
      end
      H_glb((iTri-1)*M+1:M*iTri,(iTri-1)*M+1:M*iTri) = H;
    end
    varargout{1} = H_glb;
end

end
