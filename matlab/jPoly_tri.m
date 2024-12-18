function V = jPoly_tri(X,Y,H,n,a,b,c)
% Make the Koornwinder Vandermonde matrix
% Input:
% n -  max total degree of interpolant
% X,Y - (x,y) nodes in triangle
% a,b,c > -1/2 - Jacoby poly parameters
%
% Output:
% V - (length(X),nchoosek(n+d,n)) interp matrix
%   - (P_0,...,P_n) where P_k = (P_k0,...,P_kk)


% dimension
d = 2;
% number of polynomials up to total degree n in d dimensions
% TODO: Need to add limiting value correction for x=1
N = nchoosek(n+d,n); 

V = zeros(length(X),N);
Pk = jPoly(2*Y./(1-X)-1,n+1,c-1/2,b-1/2); Pnmk = cell(n); 
for kk = 0:n
  Pnmk{kk+1} = jPoly(2*X-1,n+1,2*kk+b+c,a-1/2);
end

oneinds = find(X==1);

ind = 1;
for nn = 0:n
  for kk = 0:nn
    V(:,ind+kk) = 1/H(kk+1,nn+1) .* ...
                  Pnmk{kk+1}(:,nn-kk+1) .* ... 
                  (1-X).^kk.*Pk(:,kk+1);
    if ~isempty(oneinds)
        V(oneinds,ind+kk) = (1/H(kk+1,nn+1) .* gamma(nn+kk+b+c+1) ./ ...
                             (gamma(nn-kk+1) .* gamma(nn+kk+b+c-(nn-kk)+1))) .* ...
                             (pochhammer(kk+c+b,kk) / (2.^kk * gamma(kk+1))) .* ...
                             (2 * Y(oneinds)).^kk;
    end
  end
  ind = ind+nn+1;
end

