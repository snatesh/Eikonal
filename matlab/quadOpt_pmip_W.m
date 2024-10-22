clear all; close all; clc;

% primal dual algorithm for the quadrature problem with only
% weights as free variables.

% jacobi poly params
a = 0.5; b = 0.5; c = 0.5;
% source and target poly total degree (+1)
n = 25; m = 30; d = 2;
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

lambdak = 1./Wk;
nuk = 0;

maxiter = 1000;
mu = 10; epsfeas = 1e-13; epsgap = 1e-13;
alpha = 0.01; beta = 0.5;
disp(fobj(Xk,Yk,abs(Wk),n,m,a,b,c))
for j = 1:maxiter
    % determine t
    etagap = -Wk'*lambdak;
    t = mu*N/etagap;
    % compute primal-dual search direction
    rdual = 2*Vm*Fobj(Xk,Yk,Wk,n,m,a,b,c)-lambdak+nuk*ones(N,1);
    rcent = diag(lambdak)*Wk-(1/t)*ones(N,1);
    rpri = ones(1,N)*Wk-1;
    kkt = [rdual;rcent;rpri];
    disp([j,fobj(Xk,Yk,Wk,n,m,a,b,c),norm(rpri),norm(rcent),norm(rdual)])

    Dkkt = [2*(Vm)*(Vm'), -eye(N), ones(N,1);...
            diag(lambdak), diag(Wk), zeros(N,1);...
            ones(1,N), zeros(1,N), 0];
    dY = Dkkt\-[rdual;rcent;rpri];
    dw = dY(1:N); 
    dlambda = dY(N+1:2*N); 
    dnu = dY(2*N+1:end);
    % backtracking line search
    if isempty(lambdak(dlambda<0))
        smax = 1;
    else
        smax = min(1,min(-lambdak(dlambda<0)./dlambda(dlambda<0)));
    end
    if (smax < 1e-16)
        smax = 1;
    end
    s = 0.99*smax;
    wp = Wk + s*dw;
    lambdap = lambdak + s*dlambda;
    nup = nuk + s*dnu;
    for ibt = 1:maxiter
        rdualp = 2*Vm*Fobj(Xk,Yk,wp,n,m,a,b,c)-lambdap + nup*ones(N,1);
        rcentp = diag(lambdap)*wp-(1/t)*ones(N,1);
        rprip = ones(1,N)*wp-1;
        kktp = [rdualp;rcentp;rprip];   
        if (norm(kktp) >= (1-alpha*s)*norm(kkt))
            s = beta*s;
            wp = Wk + s*dw;
            lambdap = lambdak + s*dlambda;
            nup = nuk + s*dnu;
        else
            Wk = wp;
            lambdak = lambdap;
            nuk = nup;
            break;
        end
    end
    if (norm(rpri) <= epsfeas && norm(rdual) <= epsfeas && abs(etagap) <= epsgap)
        break
    end
end

%%
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

function f = fobj(X,Y,W,n,m,a,b,c)

f = norm(Fobj(X,Y,W,n,m,a,b,c))^2;

end

% objective gradient eq. 3.16
function dF = gradFobj_W(X,Y,W,n,m,a,b,c)

H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
dF = V_abc';

end


function df = gradfobj_W(X,Y,W,n,m,a,b,c)

df = 2*Fobj(X,Y,W,n,m,a,b,c)'*gradFobj_W(X,Y,W,n,m,a,b,c);

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
