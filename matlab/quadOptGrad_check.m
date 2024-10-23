% Check the analytical gradient 
% of the quadrature objective function 
% with second order finite differences
% i.e. Verify eq. 3.17 of eikonal.pdf

clear all; close all; clc;
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
  default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
  set(groot, default_name,'latex');
end
set(groot, 'defaultLegendFontSize',30)
set(groot, 'defaultAxesFontSize',30)

% first we define points inside of the triangle 

% jacobi poly params
a = 0.5; b = 0.5; c = 0.5;
% source and target poly total degree (+1)
n = 10; m = 12; d = 2;
% jacobi matrices
[Jn1,Jn2,A1,A2,B1,B2,Hn] = jMatON_tri(n,a,b,c);
% structure factors in target basis
Hm = structure_factors_tri(m,a,b,c);
% number of polynomials in source basis
N = nchoosek(n-1+d,n-1);
% number of polynomials in target basis
M = nchoosek(m-1+d,m-1);

% abscissa and weights
X0 = eig(Jn1+1j*Jn2);
Y = imag(X0); X = real(X0);
Vm = jPoly_tri(X,Y,Hm,m-1,a,b,c);

Vm_pinv = pinv(Vm');
W = Vm_pinv(:,1);
% analytical gradient 
dF = gradFobj(X,Y,W,n,m,a,b,c); % vector objective
df = gradfobj(X,Y,W,n,m,a,b,c); % scalar objective
% analytical Hessian
hF = hessFobj(X,Y,W,n,m,a,b,c); % vector objective
hf = hessfobj(X,Y,W,n,m,a,b,c);
 

% series of h for finite difference
hs = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10];
errsdF = zeros(size(hs));
errsdf = zeros(size(hs));
errshF = zeros(size(hs));
errshf = zeros(size(hs));
% compute FD approx to F for each h
% and find relative error in FD approx
for kk = 1:length(hs)
  dF_FD = gradFobjFD(X,Y,W,n,m,a,b,c,hs(kk));
  df_FD = gradfobjFD(X,Y,W,n,m,a,b,c,hs(kk));
  hF_FD = hessFobjFD(X,Y,W,n,m,a,b,c,hs(kk));
  hf_FD = hessfobjFD(X,Y,W,n,m,a,b,c,hs(kk));
  errsdF(kk) = norm(dF_FD-dF)/norm(dF);
  errsdf(kk) = norm(df_FD-df)/norm(df);
  errshF(kk) = norm(hF_FD(:)-hF(:))/norm(hF(:));
  errshf(kk) = norm(hf_FD(:)-hf(:))/norm(hf(:))
end

% lines should be parallel until
% roundoff error dominates fd approx
% (typically around h ~ 1e-6 to 1e-7)
% err ~ O(h^2) => log||df-df_FD|| ~ 2log(h)
loglog(hs,errsdF,'bs-','linewidth',5,'Displayname','$\mathcal{D}(\bf{F})$'); hold on;
loglog(hs,errsdf,'cs-','linewidth',5,'Displayname','$\mathcal{D}(f)$');
loglog(hs,errshF,'rs-','linewidth',5,'Displayname','$\mathcal{D}^2(\bf{F})$')
loglog(hs,errshf,'ms-','linewidth',5,'Displayname','$\mathcal{D}^2(f)$')
loglog(hs,hs.^2,'k--','linewidth',5,'Displayname','$\mathcal{O}(h^2)$');
legend('location','best')
xlabel('$h$');
ylabel('Relative 2-Norm Error')


% objective function eq. 3.13
function F = Fobj(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
F = V_abc'*W; F(1) = F(1) - 1;

end

function f = fobj(X,Y,W,n,m,a,b,c)

f = norm(Fobj(X,Y,W,n,m,a,b,c))^2;

end

% objective gradient eq. 3.17
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

function df = gradfobj(X,Y,W,n,m,a,b,c)

df = 2*Fobj(X,Y,W,n,m,a,b,c)'*gradFobj(X,Y,W,n,m,a,b,c);

end

function hF = hessFobj(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
H_a1bc1 = structure_factors_tri(m+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(m+1,a,b+1,c+1);
H_a2bc2 = structure_factors_tri(m+1,a+2,b,c+2);
H_ab2c2 = structure_factors_tri(m+1,a,b+2,c+2);
H_a1b1c2 = structure_factors_tri(m+1,a+1,b+1,c+2);

Dx = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Dxx = D1_tri(a+1,b,c+1,H_a1bc1,H_a2bc2,0) * Dx;
Dyy = D1_tri(a,b+1,c+1,H_ab1c1,H_ab2c2,1) * Dy;

Dxy = D1_tri(a,b+1,c+1,H_ab1c1,H_a1b1c2,0) * Dy;
Dyx = D1_tri(a+1,b,c+1,H_a1bc1,H_a1b1c2,1) * Dx;

Px = jPoly_tri(X,Y,H_a1bc1,m-1,a+1,b,c+1);
Py = jPoly_tri(X,Y,H_ab1c1,m-1,a,b+1,c+1);
Pxx = jPoly_tri(X,Y,H_a2bc2,m-1,a+2,b,c+2);
Pyy = jPoly_tri(X,Y,H_ab2c2,m-1,a,b+2,c+2);
Pxy = jPoly_tri(X,Y,H_a1b1c2,m-1,a+1,b+1,c+2);

N = n*(n+1)/2;
M = m*(m+1)/2;
hF = zeros(M,3*N,3*N);

dxxF = W.*(Pxx*Dxx);
dyxF = W.*(Pxy*Dyx);
dwxF = Px*Dx;
dxyF = W.*(Pxy*Dxy);
dyyF = W.*(Pyy*Dyy);
dwyF = Py*Dy;
dxwF = dwxF;
dywF = dwyF;
dwwFj = zeros(N,N); 
for j = 1:M
  dxxFj = diag(dxxF(:,j));
  dyxFj = diag(dyxF(:,j));
  dwxFj = diag(dwxF(:,j));
  dxyFj = diag(dxyF(:,j));
  dyyFj = diag(dyyF(:,j));
  dwyFj = diag(dwyF(:,j));
  dxwFj = diag(dxwF(:,j));
  dywFj = diag(dywF(:,j));
 
  hF(j,:,:) = [ dxxFj dyxFj dwxFj;...
                dxyFj dyyFj dwyFj;...
                dxwFj dywFj dwwFj ];
end

end


function hf = hessfobj(X,Y,W,n,m,a,b,c)
 
N = n*(n+1)/2;
M = m*(m+1)/2;
Z = [X;Y;W];
F = Fobj(X,Y,W,n,m,a,b,c);
DF = gradFobj(X,Y,W,n,m,a,b,c);
D2F = hessFobj(X,Y,W,n,m,a,b,c);
hf = zeros(3*N,3*N);
for j = 1:M
  D2Fj = reshape(D2F(j,:,:),3*N,3*N);
  hf = hf + F(j)*D2Fj;
end
hf = hf + DF'*DF;
hf = 2*hf;
end

% 2nd-order finite difference approx to eq. 3.17
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

function df = gradfobjFD(X,Y,W,n,m,a,b,c,h)

N = n*(n+1)/2;
M = m*(m+1)/2;
df = zeros(1,3*N);
Z = [X;Y;W];
for ii = 1:3*N
  % f(z+h)
  Z(ii) = Z(ii)+h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  fkph = fobj(x,y,w,n,m,a,b,c);
  % f(z-h)
  Z(ii) = Z(ii)-2*h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  fkmh = fobj(x,y,w,n,m,a,b,c);
  % FD approx to gradF
  df(ii) = (fkph-fkmh)/(2.0*h);
  % restore Z
  Z(ii) = Z(ii)+h;
end

end

function hF = hessFobjFD(X,Y,W,n,m,a,b,c,h)

N = n*(n+1)/2;
M = m*(m+1)/2;
hF = zeros(M,3*N,3*N);
Z = [X;Y;W];

for ii = 1:(3*N)
  % df(z+h)
  Z(ii) = Z(ii)+h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  dfkph = gradFobj(x,y,w,n,m,a,b,c);
  % df(z-h)
  Z(ii) = Z(ii)-2*h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  dfkmh = gradFobj(x,y,w,n,m,a,b,c);
  ddf = (dfkph-dfkmh)/(2.0*h);
  hF(:,:,ii) = ddf;
  % restore Z
  Z(ii) = Z(ii)+h;
end

end

function hf = hessfobjFD(X,Y,W,n,m,a,b,c,h)

N = n*(n+1)/2;
hf = zeros(3*N,3*N);
Z = [X;Y;W];
for ii = 1:3*N
  % f(z+h)
  Z(ii) = Z(ii)+h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  fkph = gradfobj(x,y,w,n,m,a,b,c)';
  % f(z-h)
  Z(ii) = Z(ii)-2*h;
  x = Z(1:N); y = Z(N+1:2*N); w = Z(2*N+1:3*N);
  fkmh = gradfobj(x,y,w,n,m,a,b,c)';
  % FD approx to gradF
  hf(:,ii) = (fkph-fkmh)/(2.0*h);
  % restore Z
  Z(ii) = Z(ii)+h;
end

end

