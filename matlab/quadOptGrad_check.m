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
n = 12; m = 14; d = 2;
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

% constraints
Aieq = zeros(5*N,3*N); bieq = zeros(5*N,1);
Aeq = zeros(1,3*N); Aeq(2*N+1:end) = 1; beq = 1;
ones = ones(3*N,1);
Aieq(1:3*N,1:3*N) = -diag(ones);
Aieq(3*N+1:4*N,1:N) = eye(N);
Aieq(4*N+1:5*N,1:N) = eye(N);
Aieq(4*N+1:5*N,N+1:2*N) = eye(N);
bieq(3*N+1:end) = 1;
lb = zeros(3*N,1); ub = lb + 1;

fun = @(Z) fobj_Z(Z,n,m,a,b,c);
Zk = [X;Y;W];

Fk = Fobj(Zk(1:N),Zk(N+1:2*N),Zk(2*N+1:end),n,m,a,b,c);
pk = norm(Fk)
iter_i = 0;
rho = 0.9; gam = 1e-4;
tol = 1e-13; tol_up = 1e3; maxiter = 1000;
while pk>tol && pk<tol_up && iter_i<maxiter
    iter_i = iter_i + 1;
    gradFk = gradFobj(Zk(1:N),Zk(N+1:2*N),Zk(2*N+1:end),n,m,a,b,c);
    % step direction
    dZk = -gradFk\Fk;
    % linesearch with wolf conditions
    alph = 1; Zk1 = Zk+alph*dZk;
    Fk1 = Fobj(Zk1(1:N),Zk1(N+1:2*N),Zk1(2*N+1:end),n,m,a,b,c);
    pk1 = norm(Fk1); iter_ii = 0;
    while pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_ii<maxiter
        alph = rho*alph;
        Zk1 = Zk+alph*dZk;
        Fk1 = Fobj(Zk1(1:N),Zk1(N+1:2*N),Zk1(2*N+1:end),n,m,a,b,c);
        pk1 = norm(Fk1); iter_ii = iter_ii+1;
    end
    Zk = Zk1; Fk = Fk1; pk = pk1;
    disp([pk,alph]);
end
X = Zk(1:N); Y = Zk(N+1:2*N); W = Zk(2*N+1:3*N);



%%

F = Fobj(X,Y,W,n,m,a,b,c);
% test function
f = @(x,y) (x+y).^10;
% eval test function on quadrature nodes
Fref = f(X,Y); FrefW = Fref.*W; normFref = norm(Fref);
nn = 11; NN = nn*(nn+1)/2;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(nn+1,a,b,c);
% vandermonde under (a,b,c), (a+1,b,c) etc.
V_abc = jPoly_tri(X,Y,H_abc,nn-1,a,b,c);
% make coeffs of fref under (a,b,c)
cfref_abc = V_abc'*FrefW;
disp(norm(V_abc*cfref_abc-Fref)/normFref);

%%

options = ...
    optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','sqp',...
    'SpecifyObjectiveGradient',true,...
    'MaxFunctionEvaluations',1e10,...
    'MaxIterations',1e10,...
    'ConstraintTolerance',1e-16,...
    'OptimalityTolerance',1e-16,...
    'StepTolerance', 1e-16, ...
    'UseParallel', true);
[Zk,fval] = fmincon(fun,Z0,Aieq,bieq,Aeq,beq,lb,ub,[],options);
Fk = Fobj(Zk(1:N),Zk(N+1:2*N),Zk(2*N+1:end),n,m,a,b,c);
pk = norm(Fk)
iter_i = 0;
rho = 0.9; gam = 1e-4;
tol = 15*eps; tol_up = 1e3; maxiter = 1000;
while pk>tol && pk<tol_up && iter_i<maxiter
    iter_i = iter_i + 1;
    gradFk = gradFobj(Zk(1:N),Zk(N+1:2*N),Zk(2*N+1:end),n,m,a,b,c);
    % step direction
    dZk = -gradFk\Fk;
    % linesearch with wolf conditions
    alph = 1; Zk1 = Zk+alph*dZk;
    Fk1 = Fobj(Zk1(1:N),Zk1(N+1:2*N),Zk1(2*N+1:end),n,m,a,b,c);
    pk1 = norm(Fk1); iter_ii = 0;
    while pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_ii<maxiter
        alph = rho*alph;
        Zk1 = Zk+alph*dZk;
        Fk1 = Fobj(Zk1(1:N),Zk1(N+1:2*N),Zk1(2*N+1:end),n,m,a,b,c);
        pk1 = norm(Fk1); iter_ii = iter_ii+1;
    end
    Zk = Zk1; Fk = Fk1; pk = pk1;
    disp([pk,alph]);
end


%%

%%
X = Zk(1:N); Y = Zk(N+1:2*N); W = Zk(2*N+1:end);
F = Fobj(X,Y,W,n,m,a,b,c);
% test function
f = @(x,y) (x+y).^8;
% eval test function on quadrature nodes
Fref = f(X,Y); FrefW = Fref.*W; normFref = norm(Fref);
nn = 7;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(nn+1,a,b,c);
% vandermonde under (a,b,c), (a+1,b,c) etc.
V_abc = jPoly_tri(X,Y,H_abc,nn-1,a,b,c);
% make coeffs of fref under (a,b,c)
cfref_abc = V_abc'*FrefW;
disp(norm(V_abc*cfref_abc-Fref)/normFref);

%% Gradient checks
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

function [f,df] = fobj_Z(Z,n,m,a,b,c)
d = 2;
N = nchoosek(n-1+d,n-1);
X = Z(1:N); 
Y = Z(N+1:2*N); 
W = Z(2*N+1:3*N);
f = fobj(X,Y,W,n,m,a,b,c);
df = gradfobj(X,Y,W,n,m,a,b,c);

end

function [eqc, Aeq, beq] = feqc(X,Y,W,n)

dim = 2;
N = n*(n+1)/2;
Aeq = zeros(1,(dim+1)*N);
Aeq(dim*N+1:(dim+1)*N) = 1;
eqc = Aeq*[X;Y;W] - 1;
beq = 1;
end

function [ieqc,G,zeroOne] = fieqc(X,Y,W,n)

dim = 2;
N = n*(n+1)/2;
G = zeros((dim+2)*N,(dim+1)*N);
one = ones((dim+1)*N,1);
G(1:(dim+1)*N,1:(dim+1)*N) = -diag(one);

%zeroOne = [zeros((dim+1)*N,1); ones(N,1)];
%I = eye(N);
%for j = 1:dim
%    G((dim+1)*N+1:(dim+2)*N,N*(j-1)+1:N*j) = I;
%end
%ieqc = G*[X;Y;W] - zeroOne;

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

    % for j = 1:N
    %     if (Zk1(j) < 0 || Zk1(j) > 1)
    %         Zk1(j) = Zk(j);
    %     end
    %     if (Zk1(j+N) < 0 || Zk1(j+N) > 1)
    %         Zk1(j+N) = Zk(j+N);
    %     end
    %     Fk1 = Fobj(Zk1(1:N),Zk1(N+1:2*N),Zk1(2*N+1:end),n,m,a,b,c);
    %     pk1 = norm(Fk1);  
    % end