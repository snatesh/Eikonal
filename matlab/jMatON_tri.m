function [Jn1,Jn2,A1,A2,B1,B2,H] = jMatON_tri(n,a,b,c)

% Make the Jacobi coefficient matrix for 
% orthonormal Jacobi polys
% Input:
% n - max poly deg
% a,b,c - Jacoby poly parameters
% Output:
% Sparse Jacobi matrices Jn1,Jn2, and blocks thereof

% dimension
d = 2;
% size of polynomial space
N = nchoosek(n-1+d,n-1);

% a = a+1/2; b = b+1/2; c = c+1/2;
kap = abs(a+b+c);

% temporary upper triangular matrices for 
% constructing coefficient matrices in recurrence
A = zeros(n+1,n+1);
B = zeros(n+1,n+1);
C = zeros(n+1,n+1);
D = zeros(n+1,n+1);
E = zeros(n+1,n+1);
F = zeros(n+1,n+1);
G = zeros(n+1,n+1);
% normalization constants (also upper tri)
H = structure_factors_tri(n+3,a,b,c); %zeros(n+3,n+3);
% coefficient matrices in 3-term recurrence
A1 = cell(n+1,1);
B1 = cell(n+1,1);
A2 = cell(n+1,1);
B2 = cell(n+1,1);
% jacobi matrices for each dim
Jn1 = zeros(N); Jn2 = zeros(N);

for nn = 0:n
  kk = (0:nn)'; 
  A(kk+1,nn+1) = (H(kk+1,nn+2)./H(kk+1,nn+1)).*(nn-kk+1).*(nn+kk+kap+1/2) ./ ...
    ((2*nn+kap+1/2).*(2*nn+kap+3/2));
  
  B(kk+1,nn+1) = 1/2-((b+c+2*kk).^2-(a-1/2).^2) ./ ...
    (2*(2*nn+kap-1/2).*(2*nn+kap+3/2));
  
  if nn > 0
    kk1 = (1:nn)';
    C(kk1,nn+1) = ((H(kk1,nn+2)./H(kk1+1,nn+1)).*(nn-kk1+1).*(nn-kk1+2).* ...
      (kk1+b-1/2).*(kk1+c-1/2)) ./ ...
      ((2*nn+kap+1/2).*(2*nn+kap+3/2).*(2*kk1+b+c).*(2*kk1+b+c-1));          
  end
  
  % note: fixes typo for e_kn on pg 81 (a_nk->a_kn)
  E(kk+1,nn+1) = -(1+((b-1/2)^2-(c-1/2)^2) ./ ...
    ((2*kk+b+c+1).*(2*kk+b+c-1))).*A(kk+1,nn+1)/2;

  % note fixes typo for f_kn on pg 81 (b_nk->b_kn)
  F(kk+1,nn+1) = (1+((b-1/2)^2-(c-1/2)^2) ./ ...
    ((2*kk+b+c+1).*(2*kk+b+c-1))).*(1-B(kk+1,nn+1))/2;  
  
  if a == 1/2 && b == 1/2 && c == 1/2
    E(1,nn+1) = -A(1,nn+1)/2;
    F(1,nn+1) = (1-B(1,nn+1))/2;
  end
  
  D(kk+1,nn+1) = (H(kk+2,nn+2)./H(kk+1,nn+1)).*(nn+kk+kap+1/2) .* ...
    (nn+kk+kap+3/2).*(kk+1).*(kk+b+c) ./ ...
    ((2*nn+kap+1/2).*(2*nn+kap+3/2).*(2*kk+b+c).*(2*kk+b+c+1));

  G(kk+1,nn+1) = (-2*H(kk+2,nn+1)./(H(kk+1,nn+1))) .* ...
    (nn-kk+a-1/2).*(nn+kk+kap+1/2).*(kk+1).*(kk+b+c) ./ ...
    ((2*nn+kap-1/2).*(2*nn+kap+3/2).*(2*kk+b+c).*(2*kk+b+c+1));
end

for nn = 0:n 
  A1{nn+1} = [diag(A(1:nn+1,nn+1)),zeros(nn+1,1)];
  B1{nn+1} =  diag(B(1:nn+1,nn+1));  
  A2{nn+1} = [diag(E(1:nn+1,nn+1)),zeros(nn+1,1)] + ...
             [zeros(nn+1,1),diag(D(1:nn+1,nn+1))] + ...
             [[zeros(1,nn);diag(C(1:nn,nn+1))],zeros(nn+1,2)];
  B2{nn+1} = diag(F(1:nn+1,nn+1)) + ...
             diag(G(1:nn,nn+1),1) + ... 
             diag(G(1:nn,nn+1),-1); 
end

inds   = [1,1]; 
for nn = 1:n-1
  Jn1(inds(1):inds(2),inds(1):inds(2)) = B1{nn}; 
  Jn2(inds(1):inds(2),inds(1):inds(2)) = B2{nn};
  inds1 = inds + [nn,nn+1];
  Jn1(inds(1):inds(2),inds1(1):inds1(2)) = A1{nn} ;
  Jn1(inds1(1):inds1(2),inds(1):inds(2)) = A1{nn}';
  Jn2(inds(1):inds(2),inds1(1):inds1(2)) = A2{nn} ;
  Jn2(inds1(1):inds1(2),inds(1):inds(2)) = A2{nn}';
  inds = inds + [nn,nn+1];
end
Jn1(inds(1):inds(2),inds(1):inds(2)) = B1{n}; 
Jn2(inds(1):inds(2),inds(1):inds(2)) = B2{n};
  
%% uncomment to check necessary conditions on recurrence matrices
% rank check (3.2.7)
% for nn = 0:n
%   rn = nchoosek(nn+d-1,nn);
%   disp([nn,rank(full(A2{nn+1})),rn]);
% end

% commutativity conditions (3.4.2)
% err = zeros(3,n);
% for nn = 1:n
%   err(1,nn) = norm(A1{nn}*A2{nn+1}-A2{nn}*A1{nn+1});
%   err(2,nn) = norm((A1{nn}*B2{nn+1}+B1{nn}*A2{nn}) - ...
%                   (B2{nn}*A1{nn}+A2{nn}*B1{nn+1}));
% end
% for nn = 2:(n+1)
%  err(3,nn) = norm((A1{nn-1}'*A2{nn-1}+B1{nn}*B2{nn}+A1{nn}*A2{nn}') - ...
%                    (A2{nn-1}'*A1{nn-1}+B2{nn}*B1{nn}+A2{nn}*A1{nn}'));
%              
% end

end











