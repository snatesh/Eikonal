clear all; close all; clc;

% jacobi poly params
a = 0.5; b = 0.5; c = 0.5;
% source and target poly total degree (+1)
n = 17; m = 20; d = 2;
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

iter = 0;
tol = 15*eps; maxiter = 100000;
%Fk = Fobj(Xk,Yk,Wk,n,m,a,b,c);
Fk = fobj(Xk,Yk,Wk,n,m,a,b,c);
pk = norm(Fk); 
rho = 0.9; gam = 1e-4;   
Zk = [Xk;Yk;Wk];

usewolfe = false;

while pk>tol && iter < maxiter
  iter = iter + 1;
  % compute gradient
  gradFk = gradfobj(Xk,Yk,Wk,n,m,a,b,c);  
  hessFk = hessfobjFD(Xk,Yk,Wk,n,m,a,b,c,1e-8);
  % step direction
  %dZk = -gradFk\Fk;  
  
  dZk = -hessFk\gradFk';
  % linesearch with wolfe conditions
  alph = 0.0001; 
  Zk1 = Zk+alph*dZk;
  Xk1 = Zk1(1:N); 
  Yk1 = Zk1(N+1:2*N); 
  Wk1 = Zk1(2*N+1:3*N);
  Fk1 = fobj(Xk1,Yk1,Wk1,n,m,a,b,c);
  pk1 = norm(Fk1); iter_i = 0; 
  if (usewolfe)
    while (pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_i < maxiter && alph>eps)
      alph = rho*alph;
      Zk1 = Zk+alph*dZk;
      Xk1 = Zk1(1:N);
      Yk1 = Zk1(N+1:2*N);
      Wk1 = Zk1(2*N+1:3*N);
      Fk1 = fobj(Xk1,Yk1,Wk1,n,m,a,b,c);
      pk1 = norm(Fk1);
      iter_i = iter_i+1;
      disp([iter_i,pk,alph])
    end
  end
  Xk = Xk1; 
  Yk = Yk1; 
  Wk = Wk1; 
  Zk = [Xk;Yk;Wk];
  Fk = Fk1; 
  pk = pk1;
  disp([pk,alph])
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

function Hf = hessfobjFD(X,Y,W,n,m,a,b,c,h)

N = n*(n+1)/2;
Hf = zeros(3*N,3*N);
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
  Hf(:,ii) = (fkph-fkmh)/(2.0*h);
  % restore Z
  Z(ii) = Z(ii)+h;
end

end
