clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
load './triquadLeg_17_28.mat';
X = Zk(1:N); Y = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
% weight function for (a+1,b+1,c+1)
w_a1b1c1 = gamma(a+1+b+1+c+1+3/2)/(gamma(a+1+1/2)*gamma(b+1+1/2)*gamma(c+1+1/2));
wa1b1c1 = @(x,y) x.^(a+1-1/2).*y.^(b+1-1/2).*(1-x-y).^(c+1-1/2)*w_a1b1c1/2;
% test function
f = @(x,y) (x+y).^10; 
dxf = @(x,y) 10*(x+y).^9;
dyf = @(x,y) 10*(x+y).^9;
lapf = @(x,y) 180*(x+y).^8;
% eval test function on quadrature nodes
F = f(X,Y); FW = F.*W; normF = norm(F);
dxF = dxf(X,Y); normdxF = norm(dxF);
dyF = dyf(X,Y); normdyF = norm(dyF);
lapF = lapf(X,Y); normlapF = norm(lapF);
syms x y
U(x,y) = x.*y.*(1-x-y).*(x+y).^8;
Us{1} = U;
U(x,y) = x.*y.*(1-x-y).*exp(y.*sin(x+y));
Us{2} = U;
U(x,y) = x.*y.*(1-x-y).*exp(-(x.^2+y.^2));
Us{3} = U;
% highest poly degree is m-1
for iF = 1:length(Us)
  
ms = 2:13;
%ms = 13;
for j = 1:length(ms)
m = ms(j); n = m+1;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc = structure_factors_tri(n+1,a+1,b,c);
H_a1b1c = structure_factors_tri(n+1,a+1,b+1,c);
H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
H_ab1c = structure_factors_tri(n+1,a,b+1,c);
H_abc1 = structure_factors_tri(n+1,a,b,c+1);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);
H_a2bc2 = structure_factors_tri(n+1,a+2,b,c+2);
H_a2b1c2 = structure_factors_tri(n+1,a+2,b+1,c+2);
H_ab2c2 = structure_factors_tri(n+1,a,b+2,c+2);
H_a1b2c2 = structure_factors_tri(n+1,a+1,b+2,c+2);
H_a2b2c2 = structure_factors_tri(n+1,a+2,b+2,c+2);

% vandermonde under (a,b,c), (a+1,b,c) etc.
V_abc = jPoly_tri(X,Y,H_abc,n-1,a,b,c);
V_a1bc = jPoly_tri(X,Y,H_a1bc,n-1,a+1,b,c);
V_ab1c = jPoly_tri(X,Y,H_ab1c,n-1,a,b+1,c);
V_abc1 = jPoly_tri(X,Y,H_abc1,n-1,a,b,c+1);
V_a1bc1 = jPoly_tri(X,Y,H_a1bc1,n-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(X,Y,H_ab1c1,n-1,a,b+1,c+1);
V_a2b2c2 = jPoly_tri(X,Y,H_a2b2c2,n-1,a+2,b+2,c+2);
V_a1b1c1 = jPoly_tri(X,Y,H_a1b1c1,n-1,a+1,b+1,c+1);
V_a1b1c1w = jPoly_tri_weighted(X,Y,H_a1b1c1,n-1,a+1,b+1,c+1);
% compute normalizations for weighted polys and normalize
%H_a1b1c1w = sqrt(sum(V_a1b1c1w.^2.*wa1b1c1(X,Y).*W));
%V_a1b1c1w = V_a1b1c1w./H_a1b1c1w;

% promotion matrices for (a,b,c) -> (a+1,b,c) etc.
K_a1bc = promotion_mat_tri(a,b,c,H_abc,H_a1bc,0);
K_a1b1c = promotion_mat_tri(a+1,b,c,H_a1bc,H_a1b1c,1);
K_a1b1c1 = promotion_mat_tri(a+1,b+1,c,H_a1b1c,H_a1b1c1,2);
K_ab1c = promotion_mat_tri(a,b,c,H_abc,H_ab1c,1);
K_abc1 = promotion_mat_tri(a,b,c,H_abc,H_abc1,2);
K_a2b1c2 = promotion_mat_tri(a+2,b,c+2,H_a2bc2,H_a2b1c2,1);
K_a2b11c2 = promotion_mat_tri(a+2,b+1,c+2,H_a2b1c2,H_a2b2c2,1);
K_a1b2c2 = promotion_mat_tri(a,b+2,c+2,H_ab2c2,H_a1b2c2,0);
K_a11b2c2 = promotion_mat_tri(a+1,b+2,c+2,H_a1b2c2,H_a2b2c2,0);


% derivative matrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dx_a2bc2 = D1_tri(a+1,b,c+1,H_a1bc1,H_a2bc2,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Dy_ab2c2 = D1_tri(a,b+1,c+1,H_ab1c1,H_ab2c2,1);
% laplacian
Lap_a2b2c2 = K_a2b11c2*K_a2b1c2*Dx_a2bc2*Dx_a1bc1 + ...\
             K_a11b2c2*K_a1b2c2*Dy_ab2c2*Dy_ab1c1;

% make coeffs of f under (a,b,c), (a+1,b,c) etc.
cf_abc = V_abc'*FW;
cf_a1bc = K_a1bc*cf_abc;
cf_ab1c = K_ab1c*cf_abc;
cf_abc1 = K_abc1*cf_abc;
cf_a1b1c1 = K_a1b1c1*(K_a1b1c*cf_a1bc);

% make coeffs of df/dx and df/dy
cdxf_a1bc1 = Dx_a1bc1*cf_abc;
cdyf_ab1c1 = Dy_ab1c1*cf_abc;
clapf_a2b2c2 = Lap_a2b2c2*cf_abc;

% check expansions, promotions and derivatives
% disp(norm(V_abc*cf_abc-F)/normF);
% disp(norm(V_a1bc*cf_a1bc-F)/normF);
% disp(norm(V_ab1c*cf_ab1c-F)/normF);
% disp(norm(V_abc1*cf_abc1-F)/normF);
% disp(norm(V_a1b1c1*cf_a1b1c1-F)/normF);
% 
% disp(norm(V_a1bc1*cdxf_a1bc1-dxF)/normdxF);
% disp(norm(V_ab1c1*cdyf_ab1c1-dyF)/normdyF);
% disp(norm(V_a2b2c2*clapf_a2b2c2-lapF)/normlapF);

% construct weigted laplacian
Wx = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
Lx = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
Dx = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Kx = promotion_mat_tri(a+1,b,c+1,H_a1bc1,H_a1b1c1,1);

Wy = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);
Ly = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);
Dy = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Ky = promotion_mat_tri(a,b+1,c+1,H_ab1c1,H_a1b1c1,0);

Dxx = Kx*Dx*Lx*Wx;
Dyy = Ky*Dy*Ly*Wy;
M = m*(m+1)/2;
A = Dxx(1:M,1:M) + Dyy(1:M,1:M);
% specify problem
%syms x y
%U(x,y) = x.*y.*(1-x-y).*sin(x+y);
U(x,y) = Us{iF};
dxU(x,y) = diff(U,'x',1);
dyU(x,y) = diff(U,'y',1);
lapU(x,y) = diff(dxU,'x',1)+diff(dyU,'y',1);
% convert to matlab funcs
f = matlabFunction(lapU);
u = matlabFunction(U);
dxu = matlabFunction(dxU);

H_abc = structure_factors_tri(n,a,b,c);
H_a1bc = structure_factors_tri(n,a+1,b,c);
H_a1b1c = structure_factors_tri(n,a+1,b+1,c);
H_a1b1c1 = structure_factors_tri(n,a+1,b+1,c+1);

V_abc = jPoly_tri(X,Y,H_abc,n-2,a,b,c);
V_a1b1c1 = jPoly_tri(X,Y,H_a1b1c1,n-2,a+1,b+1,c+1);
K_a1bc = promotion_mat_tri(a,b,c,H_abc,H_a1bc,0);
K_a1b1c = promotion_mat_tri(a+1,b,c,H_a1bc,H_a1b1c,1);
K_a1b1c1 = promotion_mat_tri(a+1,b+1,c,H_a1b1c,H_a1b1c1,2);
cf_abc = V_abc'*(f(X,Y).*W);
cf_a1b1c1 = K_a1b1c1*(K_a1b1c*(K_a1bc*cf_abc));
V_a1b1c1w = jPoly_tri_weighted(X,Y,H_a1b1c1,n-2,a+1,b+1,c+1);


%norm(V_abc*cf_abc-f(X,Y))
%norm(V_a1b1c1*cf_a1b1c1-f(X,Y))
cu_a1b1c1w = A\cf_a1b1c1;
errs(j,iF) = norm(V_a1b1c1w*cu_a1b1c1w-u(X,Y))/norm(u(X,Y))
Ns(j) = m*(m+1)/2;
end
end
%%
figure(1)
spy(A);
figure(2)
semilogy(Ns,errs(:,1),'o--','displayname','$xy(1-x-y)(x+y)^8$'); hold on;
semilogy(Ns,errs(:,2),'o--','displayname','$xy(1-x-y)\exp(y\sin(x+y))$');
semilogy(Ns,errs(:,3),'o--','displayname','$xy(1-x-y)\exp(-(x^2+y^2))$');

xlabs = {'(3,1)','(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)'};
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 14.5;
xlabel('$(N,n)$');
ax.XLabel.FontSize = 25;
ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
legend show;
%ylim([1e-17,1]);
%% now let's add boundary conditions to lap

Vx0 = jPoly_tri(0*X,Y,H_abc,n-1,a,b,c);
Vy0 = jPoly_tri(X,0*Y,H_abc,n-1,a,b,c);
V1mx = jPoly_tri(X,1-X,H_abc,n-1,a,b,c);
V_bnd = [Vx0 Vy0 V1mx];
nbnd_eq = rank(V_bnd); 
ntot_eq = size(Lap_a2b2c2,1); 
nlap_eq = ntot_eq - nbnd_eq;
[~, J] = id_decomp_hack(V_bnd',nbnd_eq);
J = J(1:nbnd_eq);
[~,J1] = id_decomp_hack(V_bnd(J,:),nbnd_eq);
J1 = J1(1:nbnd_eq);




