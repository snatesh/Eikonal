clear all; close all; clc;

% primal dual algorithm for the full quadrature problem.

% jacobi poly params
a = 0.5; b = 0.5; c = 0.5; dim = 2;
% source and target poly total degree (+1)
n = 16; m = 27; d = 2;
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
Wk1 = Vm_pinv(:,1);
Q = Vm*Vm';
Wk = Q\ones(N,1);
fobj(Xk,Yk,Wk1,n,m,a,b,c)
ftobj(Xk,Yk,n,m,a,b,c)
norm(Fobj(Xk,Yk,Wk1,n,m,a,b,c))
norm(Ftobj(Xk,Yk,n,m,a,b,c))

fun = @(X) ftobj(X(1:N),X(N+1:2*N),n,m,a,b,c);
%fun = @(X) fobj(X(1:N),X(N+1:2*N),X(2*N+1:3*N),n,m,a,b,c);
X0 = [Xk;Yk];
%X0 = [Xk;Yk;Wk];
[ieq,G,bvec] = ftieqc(Xk,Yk,n);

options = ...
  optimoptions('fmincon',...
               'Display','iter',...
               'Algorithm','sqp',...
               'FiniteDifferenceType','central',...
               'FiniteDifferenceStepSize',1e-6,...
               'MaxFunctionEvaluations',1e10,...
               'MaxIterations',1e10,...
               'ConstraintTolerance',1e-14,...
               'OptimalityTolerance',1e-14,...
               'StepTolerance', 1e-14, ...
               'UseParallel', false)
[X,fval] = fmincon(fun,X0,G,bvec,[],[],[],[],[],options);
%%
Xk = X(1:N); Yk = X(N+1:2*N); 
Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
Vm_pinv = pinv(Vm');
Wk = Vm_pinv(:,1);
ftest = @(X,Y) sin(X.^2+Y.^2);%.*exp(cos(Y));
disp([Wk'*ftest(Xk,Yk),sum(Wk), cond(Vm)]);
tri = [0 1 0 0 0 1];
figure(1);
X0 = eig(Jn1+1j*Jn2);
plot(X0,'b.');
plot(Xk,Yk,'ro');
Zk = [Xk;Yk;Wk];
save(fname,'Zk','n','m','N','M');

%%
usewolfe = false;
tol = 1e-14; rho = 0.9; gam = 1e-8; h = 1e-3;
pk = tol+1; Zk = [Xk;Yk];
iter = 0; maxiter = 10000; go = true;
while pk>tol && iter < maxiter && go
  iter = iter + 1;
  Fk = Ftobj(Xk,Yk,n,m,a,b,c);
  % compute gradient
  gradFk = gradFtobjFD(Xk,Yk,n,m,a,b,c,h);
  dZk = -gradFk\Fk;
  % linesearch with wolfe conditions
  alph = 0.00001;
  Zk1 = Zk+alph*dZk;
  Xk1 = Zk1(1:N);
  Yk1 = Zk1(N+1:2*N);
  Fk1 = Ftobj(Xk1,Yk1,n,m,a,b,c);
  pk1 = norm(Fk1); iter_i = 0;
  if (usewolfe)
    while (pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_i < maxiter && alph>eps)
      alph = rho*alph;
      Zk1 = Zk+alph*dZk;
      Xk1 = Zk1(1:N);
      Yk1 = Zk1(N+1:2*N);
      Fk1 = Ftobj(Xk1,Yk1,n,m,a,b,c);
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

  %if (norm(Zk-Zk1)/norm(Zk) < 1e-14)
  %  break;
  %end
  if (go)
      Xk = Xk1;
      Yk = Yk1;
      Zk = [Xk;Yk];
      Fk = Fk1;
      pk = pk1;
      disp([pk,alph])
  end
end





%%
t0 = 0; t1 = 0;
for rep = 1:100
tic;
fobj(Xk,Yk,Wk,n,m,a,b,c);
t0 = t0 + toc;
tic;
ftobj(Xk,Yk,n,m,a,b,c);
t1 = t1 + toc;
end

%%
[ieqc,G]= fieqc(Xk,Yk,Wk,n);

lambdak = -5./ieqc;
nuk = 0

maxiter = 1000;
mu = 20; epsfeas = 1e-13; epsgap = 1e-13;
alpha = 0.01; beta = 0.5;
disp(fobj(Xk,Yk,abs(Wk),n,m,a,b,c))

for j = 1:maxiter
    % determine t
    [ieqc,~]= fieqc(Xk,Yk,Wk,n);
    etagap = -ieqc'*lambdak;
    t = mu*N/etagap;
    % compute primal-dual search direction
    kkt = rt(Xk,Yk,Wk,lambdak,nuk,t,n,m,a,b,c);
    rpri = kkt(end);
    rdual = kkt(1:(d+1)*N);
    pk = norm(fobj(Xk,Yk,Wk,n,m,a,b,c));
    disp([pk,norm(rpri),norm(rdual)])
    dkkt = drt(Xk,Yk,Wk,lambdak,nuk,t,n,m,a,b,c);
    dZLN = dkkt\-kkt;
    dx = dZLN(1:N);
    dy = dZLN(N+1:2*N);
    dw = dZLN(2*N+1:3*N);
    dlambda = dZLN(3*N+1:3*N+(d+2)*N); 
    dnu = dZLN(end);
    % backtracking line search
    if isempty(lambdak(dlambda<0))
        smax = 1;
    else
        smax = min(1,min(-lambdak(dlambda<0)./dlambda(dlambda<0)));
    end
    if (smax < 1e-6)
        smax = 1;
    end
    s = 0.99*smax;
    xp = Xk + s*dx;
    yp = Yk + s*dy;
    wp = Wk + s*dw;
    lambdap = lambdak + s*dlambda;Wk = abs(Wk);

    nup = nuk + s*dnu;
    for ibt = 1:maxiter

        kktp = rt(xp,yp,wp,lambdap,nup,t,n,m,a,b,c);

        if (norm(kktp) > (1-alpha*s)*norm(kkt))
            s = beta*s;
            xp = Xk + s*dx;
            yp = Yk + s*dy;
            wp = Wk + s*dw;
            lambdap = lambdak + s*dlambda;
            nup = nuk + s*dnu;
        else
            Xk = xp;
            Yk = yp;
            Wk = wp;
            lambdak = lambdap;
            nuk = nup;
            break;
        end
        
    end
    if (norm(rpri) <= epsfeas && norm(rdual) <= epsfeas && abs(etagap) <= epsgap)
        break
    end
    pk1 = norm(fobj(Xk,Yk,Wk,n,m,a,b,c));
    if (norm(dZLN) < 1e-15)
        break;
    end
end

Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
ftest = @(X,Y) sin(X.^2+Y.^2);%.*exp(cos(Y));
disp([Wk'*ftest(Xk,Yk),sum(abs(Wk)), cond(Vm)]);
disp(Fobj(Xk,Yk,Wk,n,m,a,b,c));
tri = [0 1 0 0 0 1];
figure(1);
X0 = eig(Jn1+1j*Jn2);
plot(X0,'b.');
plot(Xk,Yk,'ro');
Zk = [Xk;Yk;Wk];
save(fname,'Zk','n','m','N','M');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% objective function eq. 3.13
function F = Fobj(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
F = V_abc'*W; F(1) = F(1) - 1;

end

function F = Ftobj(X,Y,n,m,a,b,c)

N = n*(n+1)/2;
M = m*(m+1)/2;

H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
one = ones(1,N);
V_pinv = (V_abc*V_abc')\V_abc;


F = one*V_pinv; F(1) = F(1) - 1;
F = F';

end

% 2nd-order finite difference approx to eq. 3.17
function dF = gradFtobjFD(X,Y,n,m,a,b,c,h)

N = n*(n+1)/2;
M = m*(m+1)/2;
dF = zeros(M,2*N);
Z = [X;Y];
for ii = 1:2*N
  % f(z+h)
  Z(ii) = Z(ii)+h;
  x = Z(1:N); y = Z(N+1:2*N);
  fkph = Ftobj(x,y,n,m,a,b,c);
  % f(z-h)
  Z(ii) = Z(ii)-2*h;
  x = Z(1:N); y = Z(N+1:2*N);
  fkmh = Ftobj(x,y,n,m,a,b,c);
  % FD approx to gradF
  dF(:,ii) = (fkph-fkmh)/(2.0*h);
  % restore Z
  Z(ii) = Z(ii)+h;
end

end



function f = fobj(X,Y,W,n,m,a,b,c)

f = 0.5*norm(Fobj(X,Y,W,n,m,a,b,c))^2;

end

function f = ftobj(X,Y,n,m,a,b,c)

N = n*(n+1)/2;
H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
one = ones(N,1);
Q = V_abc*V_abc';
f = -0.5*(one'*(Q\one) - 1);

end

function [ieqc,G,b] = ftieqc(X,Y,n)

dim = 2;
N = n*(n+1)/2;
G = zeros((dim+1)*N,dim*N);
one = ones(dim*N,1);
G(1:dim*N,1:dim*N) = -diag(one);
zeroOne = [zeros(dim*N,1); ones(N,1)];
I = eye(N);
for j = 1:dim
    G(dim*N+1:(dim+1)*N,N*(j-1)+1:N*j) = I;
end
b = zeroOne;
ieqc = G*[X;Y] - b;

end

function [eqc, A] = feqc(X,Y,W,n)

dim = 2;
N = n*(n+1)/2;
A = zeros(1,(dim+1)*N);
A(dim*N+1:(dim+1)*N) = 1;
eqc = A*[X;Y;W] - 1;

end

function [ieqc,G] = fieqc(X,Y,W,n)

dim = 2;
N = n*(n+1)/2;
G = zeros((dim+2)*N,(dim+1)*N);
one = ones((dim+1)*N,1);
G(1:(dim+1)*N,1:(dim+1)*N) = -diag(one);
zeroOne = [zeros((dim+1)*N,1); ones(N,1)];
I = eye(N);
for j = 1:dim
    G((dim+1)*N+1:(dim+2)*N,N*(j-1)+1:N*j) = I;
end
ieqc = G*[X;Y;W] - zeroOne;

end

function kkt = rt(X,Y,W,lambda,nu,t,n,m,a,b,c)

dim = 2;
N = n*(n+1)/2;
[eqc, A] = feqc(X,Y,W,n);
[ieqc, G] = fieqc(X,Y,W,n);
df = dfobj(X,Y,W,n,m,a,b,c);
kkt = [df' + G'*lambda + A'*nu;...
       -diag(lambda)*ieqc - (1/t)*ones((dim+2)*N,1);...
      eqc];

end

function dkkt = drt(X,Y,W,lambda,nu,t,n,m,a,b,c)

dim = 2;
N = n*(n+1)/2;
[~, A] = feqc(X,Y,W,n);
[ieqc, G] = fieqc(X,Y,W,n);

d2f = hessfobj(X,Y,W,n,m,a,b,c);

dkkt = [d2f', G', A';...
        -diag(lambda)*G, -diag(ieqc), zeros((dim+2)*N,1);...
        A, zeros(1,(dim+2)*N), 0];

end

function dF = dFobj(X,Y,W,n,m,a,b,c)

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

function df = dfobj(X,Y,W,n,m,a,b,c)

df = 2*Fobj(X,Y,W,n,m,a,b,c)'*dFobj(X,Y,W,n,m,a,b,c);

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


function hf = hessfobj(X,Y,W,n,m,a,b,c)
 
N = n*(n+1)/2;
M = m*(m+1)/2;
Z = [X;Y;W];
F = Fobj(X,Y,W,n,m,a,b,c);
DF = dFobj(X,Y,W,n,m,a,b,c);
D2F = hessFobj(X,Y,W,n,m,a,b,c);
hf = zeros(3*N,3*N);
for j = 1:M
  D2Fj = reshape(D2F(j,:,:),3*N,3*N);
  hf = hf + F(j)*D2Fj;
end
hf = hf + DF'*DF;
hf = 2*hf;

end


