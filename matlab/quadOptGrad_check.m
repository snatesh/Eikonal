% Check the analytical gradient 
% of the quadrature objective function 
% with second order finite differences
% i.e. Verify eq. 3.16 of eikonal.pdf

clear all; close all; clc;

% load abscissa and weights
load './triquadLeg_17_28.mat';
X = Zk(1:N); Y = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
a = 0.5; b = 0.5; c = 0.5;
% analytical gradient
dF = gradFobj(X,Y,W,17,28,a,b,c);
% series of h for finite difference
hs = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10];
errs = zeros(size(hs));
% compute FD approx to F for each h
% and find relative error in FD approx
for kk = 1:length(hs)
  dF_FD = gradFobjFD(X,Y,W,17,28,a,b,c,hs(kk));
  errs(kk) = norm(dF_FD-dF)/norm(dF);
end

% lines should be parallel until
% roundoff error dominates fd approx
% (typically around h ~ 1e-6 to 1e-7)
loglog(hs,errs,'b-','linewidth',15); hold on;
loglog(hs,hs.^2,'r--','linewidth',15)



% objective function eq. 3.13
function F = Fobj(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
F = V_abc'*W; F(1) = F(1) - 1;

end

% objective gradient eq. 3.16
function dF = gradFobj(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
H_a1bc1 = structure_factors_tri(m+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(m+1,a,b+1,c+1);
Dx = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy = D1_tri(a,b,c,H_abc,H_ab1c1,1);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
V_a1bc1 = jPoly_tri(X,Y,H_a1bc1,m-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(X,Y,H_ab1c1,m-1,a,b+1,c+1);
dF = [Dx'*(V_a1bc1'.*W') Dy'*(V_ab1c1'.*W') V_abc'];

end

% 2nd-order finite difference approx to eq. 3.16
function dF = gradFobjFD(X,Y,W,n,m,a,b,c,h)

N = n*(n+1)/2;
M = m*(m+1)/2;
dF = zeros(M,3*N);
Z = [X;Y;W];
for ii = 1:3*N
  % f(z+h)
  Z(ii) = Z(ii)+h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  fkph = Fobj(x,y,w,n,m,a,b,c);
  % f(z-h)
  Z(ii) = Z(ii)-2*h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  fkmh = Fobj(x,y,w,n,m,a,b,c);
  % FD approx to gradF
  dF(:,ii) = (fkph-fkmh)/(2.0*h);
  % restore Z
  Z(ii) = Z(ii)+h;
end

end
