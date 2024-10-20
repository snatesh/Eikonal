clear all; close all; clc;

% jacobi poly params
a = 0.5; b = 0.5; c = 0.5;
% source and target poly total degree (+1)
n = 20; m = 22; d = 2;
fname = strcat('triquadLeg_',num2str(n-1),'_',num2str(m-1),'.mat');
% jacobi matrices
[Jn1,Jn2,A1,A2,B1,B2,Hn] = jMatON_tri(n,a,b,c);
% structure factors in target basis
Hm = structure_factors_tri(m,a,b,c);
% number of polynomials in source basis
N = nchoosek(n-1+d,n-1);
% number of polynomials in target basis
M = nchoosek(m-1+d,m-1);
% initial nodes and weights for Newton iter
X0 = eig(Jn1+1j*Jn2);
Yk = imag(X0); Xk = real(X0);
Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
Vm_pinv = pinv(Vm');
Wk = Vm_pinv(:,1);
usewolfe = false;
%gradFk = gradFobj(Xk,Yk,Wk,n,m,a,b,c);
%hessfk = hessfobj(Xk,Yk,Wk,n,m,a,b,c);
T = delaunay(Xk,Yk);
Fk = Fobj(Xk,Yk,Wk,n,m,a,b,c);
trisurf(T,Xk,Yk,Fk(1:N))

%%
tol = 1e-14; rho = 0.9; gam = 1e-4; h = 1e-8;
pk = tol+1; Zk = [Xk;Yk;Wk];
iter = 0; maxiter = 10000; go = true;
while pk>tol && iter < maxiter && go
  iter = iter + 1;
  Fk = Fobj(Xk,Yk,Wk,n,m,a,b,c);
  % compute gradient
  gradFk = gradFobj(Xk,Yk,Wk,n,m,a,b,c);
  dZk = -gradFk\Fk;
 % hessFk = hessFobj(Xk,Yk,Wk,n,m,a,b,c); 
  % step direction
%  dZk = zeros(3*N,1);
%   for j = 1:M
%     hessFkj = reshape(hessFk(j,:,:),3*N,3*N)';
%     issymmetric(hessFkj)
%     gradFkj = gradFk(j,:)';
%     dZk = dZk + hessFkj\gradFkj;
%   end
  % linesearch with wolfe conditions
  alph = 0.001; 
  Zk1 = Zk+alph*dZk;
  Xk1 = Zk1(1:N); 
  Yk1 = Zk1(N+1:2*N); 
  Wk1 = Zk1(2*N+1:3*N);
  Fk1 = Fobj(Xk1,Yk1,Wk1,n,m,a,b,c);
  pk1 = norm(Fk1); iter_i = 0; 
  if (usewolfe)
    while (pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_i < maxiter && alph>eps)
      alph = rho*alph;
      Zk1 = Zk+alph*dZk;
      Xk1 = Zk1(1:N);
      Yk1 = Zk1(N+1:2*N);
      Wk1 = Zk1(2*N+1:3*N);
      Fk1 = Fobj(Xk1,Yk1,Wk1,n,m,a,b,c);
      pk1 = norm(Fk1);
      iter_i = iter_i+1;
      disp([pk1,alph])
    end
  end

  for j = 1:N
    if (Xk1(j) < 0 || Xk1(j) > 1 || ...
        Yk1(j) < 0 || Yk1(j) > 1 || ...
        Yk1(j) > 1-Xk1(j))
      go = false; 
      break;
    end
  end

  if (norm(Zk-Zk1)/norm(Zk) < 1e-14)
    break;
  end
  if (go)
  Xk = Xk1; 
  Yk = Yk1; 
  Wk = Wk1; 
  Zk = [Xk;Yk;Wk];
  Fk = Fk1; 
  pk = pk1;
  disp([pk,alph])
  end
end


Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
ftest = @(X,Y) sin(X.^2+Y.^2);%.*exp(cos(Y));
disp([pk,Wk'*ftest(Xk,Yk),sum(abs(Wk)), cond(Vm)])

tri = [0 1 0 0 0 1];
figure(1);
X0 = eig(Jn1+1j*Jn2);
plot(X0,'b.');
plot(Xk,Yk,'ro');
save(fname,'Zk','n','m','N','M');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% objective function eq. 3.13
function F = Fobj(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
F = V_abc'*W; F(1) = F(1) - 1;

end

function f = fobj(X,Y,W,n,m,a,b,c)

f = norm(Fobj(X,Y,W,n,m,a,b,c))^2;

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


function df = gradfobj(X,Y,W,n,m,a,b,c)

df = 2*Fobj(X,Y,W,n,m,a,b,c)'*gradFobj(X,Y,W,n,m,a,b,c);

end

function Hf = hessFobj(X,Y,W,n,m,a,b,c)

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
Hf = zeros(M,3*N,3*N);

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
 
  Hf(j,:,:) = [ dxxFj dyxFj dwxFj;...
                dxyFj dyyFj dwyFj;...
                dxwFj dywFj dwwFj ];
end

end

function Hf = hessFobjFD(X,Y,W,n,m,a,b,c,h)

N = n*(n+1)/2;
M = m*(m+1)/2;
Hf = zeros(M,3*N,3*N);
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
  Hf(:,:,ii) = ddf;
  % restore Z
  Z(ii) = Z(ii)+h;
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
