clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
load './triquadLeg_17_28.mat';
X = Zk(1:N); Y = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
[xx,yy] = meshgrid(linspace(0,1,30)); 
xx = xx(:); yy = yy(:);
indtri = yy<=1-xx;
xx = xx(indtri);
yy = yy(indtri);
dist = distToTri(xx,yy);
T = delaunay(xx,yy);
figure(1);
trisurf(T,xx,yy,dist);


% weight function for (a+1,b+1,c+1)
w_a1b1c1 = gamma(a+1+b+1+c+1+3/2)/(gamma(a+1+1/2)*gamma(b+1+1/2)*gamma(c+1+1/2));
wa1b1c1 = @(x,y) x.^(a+1-1/2).*y.^(b+1-1/2).*(1-x-y).^(c+1-1/2)*w_a1b1c1/2;
% rhs
f = @(x,y) ones(size(X));
u = @(x,y) distToTri(x,y);

% eval rhs and sol on quadrature nodes
% W is defined without the factor of 1/2
F = f(X,Y); FW = F.*W; 
U = u(X,Y); UW = U.*W;

% highest poly degree is m-1
ms = 12; j = 1; iF = 1;
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

K_a1bc1_a1b1c1 = promotion_mat_tri(a+1,b,c+1,H_a1bc1,H_a1b1c1,1);
K_ab1c1_a1b1c1 = promotion_mat_tri(a,b+1,c+1,H_ab1c1,H_a1b1c1,0);

% derivative matrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dx_a2bc2 = D1_tri(a+1,b,c+1,H_a1bc1,H_a2bc2,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Dy_ab2c2 = D1_tri(a,b+1,c+1,H_ab1c1,H_ab2c2,1);

Dx = K_a1bc1_a1b1c1*Dx_a1bc1;
Dy = K_ab1c1_a1b1c1*Dy_ab1c1;

% Eikonal operator
Eik_abc_a1b1c1 = Dx'*Dx + Dy'*Dy;

% make coeffs of f under (a,b,c), (a+1,b,c) etc.
cf_abc = V_abc'*FW;
cf_a1b1c1 = K_a1b1c1*(K_a1b1c*(K_a1bc * cf_abc));
cu_abc = V_abc'*UW;
cu_a1b1c1 = K_a1b1c1*(K_a1b1c*(K_a1bc * cu_abc))
T = delaunay(X,Y);


Np = n*(n+1)/2;
PPT = zeros(Np);
PPT1 = zeros(Np);
Wa1b1c1 = wa1b1c1(X,Y);
for j = 1:Np
    for i = 1:Np
        PPT(i,j) = W'*(V_a1b1c1(:,i).*V_a1b1c1(:,j).*Wa1b1c1);
        PPT1(i,j) = W'*(V_abc(:,i).*V_abc(:,j));
    end
end
PPT(PPT < 1e-10) = 0;
PPT1(PPT1 < 1e-10) = 0;



spy(PPT)

resTri = ((V_a1b1c1*(Dx*cu_abc)).^2 + (V_a1b1c1*(Dy*cu_abc)).^2 ...
        -(V_a1b1c1*cf_a1b1c1).^2).*Wa1b1c1;

resTriNorm = W'*resTri
figure(2)
%trisurf(T,X,Y,resTri);





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




