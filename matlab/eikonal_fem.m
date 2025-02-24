clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

% highest poly degree is m-1, and there are M total polys
m = 8; n = m+1; M = m*(m+1)/2;
speed = 1; 
xi = 0.001;
u = @(x,y) distToTri(x,y)/speed;

%u = @(x,y) x.*y.*(1-x-y).*sin(x+y);
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc = structure_factors_tri(n+1,a+1,b,c);
H_a1b1c = structure_factors_tri(n+1,a+1,b+1,c);
H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
H_ab1c = structure_factors_tri(n+1,a,b+1,c);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);
H_a2bc2 = structure_factors_tri(n+1,a+2,b,c+2);
H_a2b1c2 = structure_factors_tri(n+1,a+2,b+1,c+2);
H_ab2c2 = structure_factors_tri(n+1,a,b+2,c+2);
H_a1b1c2 = structure_factors_tri(n+1,a+1,b+1,c+2);
H_a1b2c2 = structure_factors_tri(n+1,a+1,b+2,c+2);
H_a2b2c2 = structure_factors_tri(n+1,a+2,b+2,c+2);
H_a2b1c1 = structure_factors_tri(n+1,a+2,b+1,c+1);
H_a2b2c1 = structure_factors_tri(n+1,a+2,b+2,c+1);

% vandermonde under (a,b,c)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
% vandermonde under(a+1,b+1,c+1)
V_a1b1c1 = jPoly_tri(R,S,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);
% weighted vandermonde under (a+1,b+1,c+1)
V_a1b1c1w = jPoly_tri_weighted(R,S,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);
% weight function for weighted basis (a+1,b+1,c+1)
w_a1b1c1 = @(R,S) R.*S.*(1-R-S);

% promotion matrices (a+2,b,c+2) -> (a+2,b+1,c+2) etc.
K_abc_a1bc = promotion_mat_tri(a,b,c,H_abc,H_a1bc,0);
K_a1bc_a1b1c = promotion_mat_tri(a+1,b,c,H_a1bc,H_a1b1c,1);
K_a1b1c_a1b1c1 = promotion_mat_tri(a+1,b+1,c,H_a1b1c,H_a1b1c1,2);
K_a1bc1_a1b1c1 = promotion_mat_tri(a+1,b,c+1,H_a1bc1,H_a1b1c1,1);
K_ab1c1_a1b1c1 = promotion_mat_tri(a,b+1,c+1,H_ab1c1,H_a1b1c1,0);
% unweighted derivative matrices
Dx_abc_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy_abc_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
% weighted derivative amatrices
Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);
% lowering matrices
Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);
% laplacian in weighted basis (3/2,/3/2,3/2)
Dxx_w = K_a1bc1_a1b1c1*Dx_abc_a1bc1*Lx_ab1c_abc*Wx_a1b1c1_ab1c;
Dyy_w = K_ab1c1_a1b1c1*Dy_abc_ab1c1*Ly_a1bc_abc*Wy_a1b1c1_a1bc;
Dxy_w = K_a1bc1_a1b1c1*Dx_abc_a1bc1*Ly_a1bc_abc*Wy_a1b1c1_a1bc;
Dxx_w = Dxx_w(1:M,1:M);
Dyy_w = Dyy_w(1:M,1:M);
Dxy_w = Dxy_w(1:M,1:M);
Lap_a1b1c1w_a1b1c1 = Dxx_w + Dyy_w;
% RHS for poisson
F = @(R,S) RHS(R,S,speed);
% coefficients of RHS in (1/2,1/2,1/2)
cF_abc = V_abc(:,1:M)'*((-F(R,S)).*W);
% promote coefficients to (3/2,3/2,3/2)
cF_a1b1c1 = K_a1b1c_a1b1c1(1:M,1:M)...
    *(K_a1bc_a1b1c(1:M,1:M)...
    *(K_abc_a1bc(1:M,1:M)*cF_abc));
% coeffs of sol in weighted basis (3/2,3/2,3/2)
cu_a1b1c1w = Lap_a1b1c1w_a1b1c1\cF_a1b1c1;

T = delaunay(R,S);
[xx,yy] = meshgrid(linspace(1e-10,1-1e-10,50)); 
xx = xx(:); yy = yy(:);
indtri = yy<=1-xx;
xx = xx(indtri);
yy = yy(indtri);
Tgrid = delaunay(xx,yy);
Vw_grid = jPoly_tri_weighted(xx,yy,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);


Dx_w = Lx_ab1c_abc*Wx_a1b1c1_ab1c; Dx_w_tr = Dx_w(1:M,1:M);
Dy_w = Ly_a1bc_abc*Wy_a1b1c1_a1bc; Dy_w_tr = Dy_w(1:M,1:M);

% 1/speed = cf(1)
rhs = @(R,S) RHS(R,S,speed);
cf_abc = V_abc(:,1:M)'*(rhs(R,S).*W);
Kap = K_a1b1c_a1b1c1(1:M,1:M) * ...
      K_a1bc_a1b1c(1:M,1:M) * ...
      K_abc_a1bc(1:M,1:M);
cf_a1b1c1 = Kap * cf_abc;




% 
% val = eik_nonlin_coeffs(cu_a1b1c1w,cf_a1b1c1,...
%                         Lap_a1b1c1w_a1b1c1,Kap,...
%                         Dx_w_tr,Dy_w_tr,V_abc(:,1:M),,W,xi);
% 
% val1 = eik_nonlin_coeffs1(cu_a1b1c1w,cf_abc,...
%                           Lap_a1b1c1w_a1b1c1,Dx_w_tr,Dy_w_tr,...
%                           V_abc(:,1:M),V_a1b1c1,W,xi);

% fun = @(cu) abs(eik_nonlin_coeffs(cu,cf_a1b1c1,...
%                               Lap_a1b1c1w_a1b1c1,Kap,...
%                               Dx_w_tr,Dy_w_tr,V_abc(:,1:M),W,xi));


fun = @(cu) eik_nonlin_coeffs1(cu,cf_abc,...
                               Lap_a1b1c1w_a1b1c1,Dx_w_tr,Dy_w_tr,...
                               V_abc(:,1:M),V_a1b1c1,W,xi);

%A = []; bvec = [];
A = -V_a1b1c1w; bvec = zeros(length(R),1);
%A = -(W/2)'*V_a1b1c1w; bvec = 0;
%A = -Vw_grid; bvec = zeros(length(xx),1);
x0 = cu_a1b1c1w;
Aeq = [];
beq = [];

options = ...
    optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','sqp',...
    'FiniteDifferenceType','central',...
    'FiniteDifferenceStepSize',1e-6,...
    'MaxFunctionEvaluations',1e10,...
    'MaxIterations',1e10,...
    'ConstraintTolerance',1e-8,...
    'OptimalityTolerance',1e-8,...
    'StepTolerance', 1e-8, ...
    'UseParallel', false);


cu_sol = fmincon(fun,x0,A,bvec,Aeq,beq,[],[],[],options);

figure()
tiledlayout(1,2);
%ax1 = nexttile;
%trisurf(Tgrid,xx,yy,Vw_grid*cu_a1b1c1w);
ax2 = nexttile;
trisurf(Tgrid,xx,yy,u(xx,yy))

ax3 = nexttile;
trisurf(Tgrid,xx,yy,Vw_grid*cu_sol)
linkaxes([ax2,ax3],'z')
%% now try constructing a finite element discretization for poisson
[K,B,F,dPdP,V_abc,V_a1bc1,V_ab1c1,Dx,Dy] = assemble_poisson_unweighted(n,R,S,W,rhs);
M = n*(n+1)/2;
[~,~,II] = qr(B,0);
nLambda = rank(B);
BB = B(II(1:nLambda),:);
Amat = zeros(M+nLambda);
Amat(1:M,1:M) = K;
Amat(1:M,M+1:end) = BB';
Amat(M+1:end,1:M) = BB;
Fvec = [F;zeros(nLambda,1)];
cu = Amat\Fvec; cu_abc = cu(1:M);
tiledlayout(1,2);
nexttile;
trisurf(T,R,S,V_abc*cu_abc);
nexttile
scatter3(R,S,V_abc*cu_abc);

%% strong form discretization with weak bc enforcement for poisson
M = n*(n+1)/2;
% promotion matrices
K_a2bc2_a2b1c2 = promotion_mat_tri(a+2,b,c+2,H_a2bc2,H_a2b1c2,1);
K_a2b1c2_a2b2c2 = promotion_mat_tri(a+2,b+1,c+2,H_a2b1c2,H_a2b2c2,1);
K_ab2c2_a1b2c2 = promotion_mat_tri(a,b+2,c+2,H_ab2c2,H_a1b2c2,0);
K_a1b2c2_a2b2c2 = promotion_mat_tri(a+1,b+2,c+2,H_a1b2c2,H_a2b2c2,0);
% unweighted derivative matrices
Dx_abc_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dx_a1bc1_a2bc2 = D1_tri(a+1,b,c+1,H_a1bc1,H_a2bc2,0);
Dy_abc_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Dy_ab1c1_ab2c2 = D1_tri(a,b+1,c+1,H_ab1c1,H_ab2c2,1);

% total number of equations
% laplacian in unweighted basis
Lap_abc_a2b2c2 = K_a2b1c2_a2b2c2*K_a2bc2_a2b1c2*Dx_a1bc1_a2bc2*Dx_abc_a1bc1 + ...\
                 K_a1b2c2_a2b2c2*K_ab2c2_a1b2c2*Dy_ab1c1_ab2c2*Dy_abc_ab1c1;  

% boundary stiffness
[xleg,wleg,~] = gjQuad(100,0,0);
wleg = wleg'/2;
xleg = (xleg+1)/2;
Xl = 0*xleg;
Yl = xleg;
Xb = xleg;
Yb = 0*xleg;
Xh = xleg;
Yh = 1-Xh;

Vl = jPoly_tri(Xl,Yl,H_abc,n-1,a,b,c);
Vb = jPoly_tri(Xb,Yb,H_abc,n-1,a,b,c);
Vh = jPoly_tri(Xh,Yh,H_abc,n-1,a,b,c);
BCrows = zeros(M,M);
for i = 1:M
    for j = 1:M
        int1 = wleg'*(Vl(:,i).*Vl(:,j));
        int2 = wleg'*(Vb(:,i).*Vb(:,j));
        int3 = wleg'*(Vh(:,i).*Vh(:,j));
        BCrows(i,j) = int1+int2+int3;
    end
end

Fbc = zeros(M,1);
dirch = @(R,S) zeros(size(R));
for i = 1:M
    int1 = wleg'*(Vl(:,i).*dirch(Xl,Yl));
    int2 = wleg'*(Vb(:,i).*dirch(Xb,Yb));
    int3 = wleg'*(Vh(:,i).*dirch(Xh,Yh));
    Fbc(i) = int1+int2+int3;
end
[~,~,II] = qr(BCrows,0);
nLambda = rank(BCrows);
BB = BCrows(II(1:nLambda),:);
AAmat = zeros(M+nLambda);
AAmat(1:M,1:M) = Lap_abc_a2b2c2;
AAmat(1:M,M+1:end) = BB';
AAmat(M+1:end,1:M) = BB;
cf_a2b2c2 = zeros(M,1); cf_a2b2c2(1) = -1./speed;
cF = [cf_a2b2c2;Fbc(II(1:nLambda))];
cu1_abc = AAmat \ cF;
tiledlayout(1,2);
nexttile;
trisurf(T,R,S,V_abc*cu1_abc(1:M));
nexttile
scatter3(R,S,V_abc*cu1_abc(1:M));


%%
% % check implementation of weak form by verifying that it 
% % is indeed the first variation of poisson energy functional
% % cf_abc here is actually cf_a1b1c1w.
% [err_grad,hs] = check_variations_poiss(n,cf_abc,M,R,S,W,rhs);
% figure()
% loglog(err_grad,hs,'r.-'); hold on; loglog(hs,hs,'k--');
% solve poisson problem for weighted coeffs cu_femw
M = (n-1)*n/2;
[K,F,V_abc,V_a1b1c1w,Dx_w,Dy_w,dPdP] = assemble_poisson(n,M,R,S,W,rhs);

cu_femw = K\F;
%%
% %  now check weak form for nonlinear eikonal problem
% [err_grad,err_hess,hs] = check_variations(n,cu_femw,M,R,S,W,xi,rhs);
% figure()
% loglog(hs,err_grad,'r.-'); hold on;
% loglog(hs,err_hess,'b.-');
% loglog(hs,hs,'k--');


%%
% now do an inf-dim newton stepper for eikonal
[N,K,F,H] = assemble_eik_nonlin(n,cu_femw,M,R,S,W,xi,rhs,...
                               V_abc, V_a1b1c1w,Dx_w,Dy_w,K,F,...
                               dPdP);
H = H';
% 
% [N,K,F] = assemble_eik_nonlin(n,cu_femw,M,R,S,W,xi,rhs,...
%                                 V_abc, V_a1b1c1w,Dx_w,Dy_w,K,F,...
%                                 dPdP);

max_iter = 100;
cu = cu_femw;
rtol = 1e-6;
tol = norm(N+xi*K*cu-F)*rtol;

fprintf("It \t Energy \t (g,du) \t ||g||l2\n")
for it = 1:max_iter
    cu_prev = cu;
    %if ~(rem(it,5))
        [N,K,F,H] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,...
                                        V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
        H = H'; alph = 1/2;
    %else
    %    [N,K,F] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,...
    %                                  V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
    %    alph = 0.001;
    %end
    
    brhs = N+xi*K*cu-F;
    if norm(brhs) < tol
        fprintf('Converged in %d Newton iterations\n', it);
        break;
    end

    du = -H\brhs; 
    cu1 = cu + alph*du;

    % bfgs update for hessian
    % alph = 0.5; 
    % sk = alph*du;
    % cu1 = cu+sk;
    % [N1,~,~] = assemble_eik_nonlin(n,cu1,M,R,S,W,xi,rhs,...
    %                               V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
    % brhs1 = N1+xi*K*cu1-F;
    % yk = brhs1 - brhs;
    % H = H + ((sk'*yk + yk'*(H*yk))*(sk*sk'))/((sk'*yk)^2) ...
    %       - ((H*yk)*sk' + sk*(yk'*H))/(sk'*yk);

    % line search
    % [N1,~,~] = assemble_eik_nonlin(n,cu1,M,R,S,W,xi,rhs,...
    %                               V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
    % brhs1 = N1+xi*K*cu1-F;
    % pk1 = norm(brhs1); iter_ii = 0;
    % rho = 0.9; gam = 1e-4;
    % while pk1>norm(brhs+gam*alph*H*du) && iter_ii<max_iter && alph > 1e-6
    %   alph = rho*alph;
    %   cu1 = cu+alph*du;
    %   [N1,~,~] = assemble_eik_nonlin(n,cu1,M,R,S,W,xi,rhs,...
    %                               V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
    %   brhs1 = N1+xi*K*cu1-F;
    %   pk1 = norm(brhs1); iter_ii = iter_ii+1;
    % end

    % fk = assemble_energy(cu,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w);
    % gradfk = brhs;
    % du = -pinv(gradfk)*fk; 
    % alph = 0.005; cu1 = cu+alph*du';
    % % fk1 = assemble_energy(cu1,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w);
    % % pk1 = fk1; iter_ii = 0;
    % % rho = 0.9; gam = 1e-4;
    % % while pk1>norm(fk+gam*alph*gradfk'*du') && iter_ii<max_iter && alph > 1e-3
    % %   alph = rho*alph;
    % %   cu1 = cu+alph*du';
    % %   fk1 = assemble_energy(cu1,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w);
    % %   pk1 = fk1; iter_ii = iter_ii+1;
    % %   %fprintf("%d \t %1.10f \n", it+iter_ii, fk1)
    % % end

    cu = cu1;
    
    fk = assemble_energy(cu,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w);
    fprintf("%d \t %1.10f \t %1.10f \t %1.10f \t %1.10f\n", it, fk, -brhs'*du, norm(brhs), norm(cu_prev-cu))
    %fprintf("%d \t %1.10f \t %1.10f \n", it, fk, norm(brhs))

end

figure()
tiledlayout(1,2);
ax2 = nexttile;
trisurf(Tgrid,xx,yy,u(xx,yy))

ax3 = nexttile;
trisurf(Tgrid,xx,yy,Vw_grid*cu)
linkaxes([ax2,ax3],'z')


%%%%%%%%%%%%%%%%%%%%%%%%%% FEM assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,B,F,dPdP,V_abc,V_a1bc1,V_ab1c1,Dx,Dy] = assemble_poisson_unweighted(n,R,S,W,rhs)
a = 1/2; b = a; c = a;
M = n*(n+1)/2;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);

% vandermonde under (a,b,c), (a+1,b,c+1), (a,b+1,c+1)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
V_a1bc1 = jPoly_tri(R,S,H_a1bc1,n-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(R,S,H_ab1c1,n-1,a,b+1,c+1);

% derivative amatrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);

% interior stiffness
K = zeros(M,M); 
dPdP = zeros(M,M,length(R));
for i = 1:M
    dxP_i = V_a1bc1 * Dx_a1bc1(:,i);
    dyP_i = V_ab1c1 * Dy_ab1c1(:,i);
    for j = 1:M
        dxP_j = V_a1bc1 * Dx_a1bc1(:,j);
        dyP_j = V_ab1c1 * Dy_ab1c1(:,j);   
        dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
        K(i,j) = (W/2)'*dPidPj;
        dPdP(i,j,:) = dPidPj;
    end
end

% boundary stiffness
[xleg,wleg,~] = gjQuad(50,0,0);
wleg = wleg'/2;
xleg = (xleg+1)/2;
Xl = 0*xleg;
Yl = xleg;
Xb = xleg;
Yb = 0*xleg;
Xh = xleg;
Yh = 1-Xh;

Vl = jPoly_tri(Xl,Yl,H_abc,n-1,a,b,c);
Vb = jPoly_tri(Xb,Yb,H_abc,n-1,a,b,c);
Vh = jPoly_tri(Xh,Yh,H_abc,n-1,a,b,c);
B = zeros(M,M);
for i = 1:M
    for j = 1:M
        int1 = wleg'*(Vl(:,i).*Vl(:,j));
        int2 = wleg'*(Vb(:,i).*Vb(:,j));
        int3 = wleg'*(Vh(:,i).*Vh(:,j));
        B(i,j) = int1+int2+int3;
    end
end

F = zeros(M,1);
for i = 1:M
    P_if = V_abc(:,i).*rhs(R,S);
    F(i) = (W/2)'*P_if;
end

Dx = Dx_a1bc1;
Dy = Dy_ab1c1;

end



function [K,F,varargout] = assemble_poisson(n,M,R,S,W,rhs)
a = 1/2; b = a; c = a;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc = structure_factors_tri(n+1,a+1,b,c);
H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
H_ab1c = structure_factors_tri(n+1,a,b+1,c);

% vandermonde under (a,b,c)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
% weighted vandermonde under (a+1,b+1,c+1)
V_a1b1c1w = jPoly_tri_weighted(R,S,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);



% weighted derivative amatrices
Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);
% lowering matrices
Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);
% derivative ops in weighted basis (output in (a,b,c))
Dx_w = Lx_ab1c_abc*Wx_a1b1c1_ab1c;
Dy_w = Ly_a1bc_abc*Wy_a1b1c1_a1bc;


K = zeros(M,M); 
dPdP = zeros(M,M,length(R));
for i = 1:M
    dxP_i = V_abc * Dx_w(:,i);
    dyP_i = V_abc * Dy_w(:,i);
    for j = 1:M
        dxP_j = V_abc * Dx_w(:,j);
        dyP_j = V_abc * Dy_w(:,j);   
        dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
        K(i,j) = (W/2)'*dPidPj;
        dPdP(i,j,:) = dPidPj;
    end
end

F = zeros(M,1);
for i = 1:M
    P_if = V_a1b1c1w(:,i).*rhs(R,S);
    F(i) = (W/2)'*P_if;
end

if nargout > 2
    varargout{1} = V_abc;
    varargout{2} = V_a1b1c1w;
    varargout{3} = Dx_w;
    varargout{4} = Dy_w;
    varargout{5} = dPdP;
end

end

function [N,K,F,varargout] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,varargin)

a = 1/2; b = a; c = a;

if nargin == 8
    % normalization under (a,b,c), (a+1,b,c) etc.
    H_abc = structure_factors_tri(n+1,a,b,c);
    H_a1bc = structure_factors_tri(n+1,a+1,b,c);
    H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
    H_ab1c = structure_factors_tri(n+1,a,b+1,c);



    % vandermonde under (a,b,c)
    V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
    % weighted vandermonde under (a+1,b+1,c+1)
    V_a1b1c1w = jPoly_tri_weighted(R,S,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);

    % weighted derivative amatrices
    Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
    Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);
    % lowering matrices
    Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
    Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);
    % derivative ops in weighted basis (output in (a,b,c))
    Dx_w = Lx_ab1c_abc*Wx_a1b1c1_ab1c;
    Dy_w = Ly_a1bc_abc*Wy_a1b1c1_a1bc;

    K = zeros(M,M);
    for i = 1:M
        for j = 1:M
            dxP_i = V_abc * Dx_w(:,i);
            dxP_j = V_abc * Dx_w(:,j);
            dyP_i = V_abc * Dy_w(:,i);
            dyP_j = V_abc * Dy_w(:,j);
            dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
            K(i,j) = (W/2)'*dPidPj;
        end
    end

    F = zeros(M,1);
    for i = 1:M
        Pf = V_a1b1c1w(:,i).*rhs(R,S);
        F(i) = (W/2)'*Pf;
    end

    N = zeros(M,1); Nrs = length(R);
    for i = 1:M
        P_i = V_a1b1c1w(:,i);
        sum = zeros(Nrs,1);
        for j = 1:M
            for k = 1:M
                dxP_k = V_abc * Dx_w(:,k);
                dxP_j = V_abc * Dx_w(:,j);
                dyP_k = V_abc * Dy_w(:,k);
                dyP_j = V_abc * Dy_w(:,j);
                dPkdPj = dxP_k.*dxP_j + dyP_k.*dyP_j;
                sum = sum + cu(k)*cu(j).*dPkdPj;
            end
        end
        sum = sqrt(sum);
        N(i) = (W/2)'*(P_i.*sum);
    end

    H = zeros(M);
    for i = 1:M
        dxP_i = V_abc * Dx_w(:,i);
        dyP_i = V_abc * Dy_w(:,i);
        for j = 1:M
            P_j = V_a1b1c1w(:,j);
            sum1 = zeros(Nrs,1);
            sum2 = zeros(Nrs,1);
            for k = 1:M
                dxP_k = V_abc * Dx_w(:,k);
                dyP_k = V_abc * Dy_w(:,k);
                dPidPk = dxP_k.*dxP_i + dyP_k.*dyP_i;
                sum1 = sum1 + cu(k)*dPidPk;
                for l = 1:M
                    dxP_l = V_abc * Dx_w(:,l);
                    dyP_l = V_abc * Dy_w(:,l);
                    dPkdPl = dxP_k.*dxP_l + dyP_k.*dyP_l;
                    sum2 = sum2 + cu(k)*cu(l).*dPkdPl;
                end
            end
            sum1 = P_j.*sum1;
            sum2 = sqrt(sum2);
            H(i,j) = (W/2)'*(sum1./sum2) + xi*K(i,j);
        end
    end
else
    V_abc = varargin{1};
    V_a1b1c1w = varargin{2};
    Dx_w = varargin{3};
    Dy_w = varargin{4};
    K = varargin{5};
    F = varargin{6};
    dPdP = varargin{7};
    
    N = zeros(M,1); Nrs = length(R);
    for i = 1:M
        P_i = V_a1b1c1w(:,i);
        sum = zeros(Nrs,1);
        for j = 1:M
            for k = 1:M
                dPkdPj = reshape(dPdP(k,j,:),Nrs,1);
                sum = sum + cu(k)*cu(j).*dPkdPj;
            end
        end
        sum = sqrt(sum);
        N(i) = (W/2)'*(P_i.*sum);
    end
    
    if nargout == 4
    tic
    H = zeros(M);
    for i = 1:M
        for j = 1:M
            P_j = V_a1b1c1w(:,j);
            sum1 = zeros(Nrs,1);
            sum2 = zeros(Nrs,1);
            for k = 1:M
                dPidPk = reshape(dPdP(i,k,:),Nrs,1);
                sum1 = sum1 + cu(k)*dPidPk;
                for l = 1:M
                    dPkdPl = reshape(dPdP(k,l,:),Nrs,1);
                    sum2 = sum2 + cu(k)*cu(l).*dPkdPl;
                end
            end
            sum1 = P_j.*sum1;
            sum2 = sqrt(sum2);
            H(i,j) = (W/2)'*(sum1./sum2) + xi*K(i,j);
        end
    end
    toc
    varargout{1} = H;
    end
end


end
  
function energy = assemble_energy_poiss(n,cu,M,R,S,W,rhs)
a = 1/2; b = a; c = a;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc = structure_factors_tri(n+1,a+1,b,c);
H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
H_ab1c = structure_factors_tri(n+1,a,b+1,c);



% vandermonde under (a,b,c)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
% weighted vandermonde under (a+1,b+1,c+1)
V_a1b1c1w = jPoly_tri_weighted(R,S,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);

% weighted derivative amatrices
Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);
% lowering matrices
Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);
% derivative ops in weighted basis (output in (a,b,c))
Dx_w = Lx_ab1c_abc*Wx_a1b1c1_ab1c;
Dy_w = Ly_a1bc_abc*Wy_a1b1c1_a1bc;
sum1 = 0;
sum3 = 0;
for i = 1:M
    P_i = V_a1b1c1w(:,i);
    for j = 1:M
        dxP_i = V_abc * Dx_w(:,i);
        dxP_j = V_abc * Dx_w(:,j);
        dyP_i = V_abc * Dy_w(:,i);
        dyP_j = V_abc * Dy_w(:,j);
        dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
        sum1 = sum1 + cu(i)*cu(j)*dPidPj;
    end
    sum3 = sum3 + cu(i)*P_i;
end
sum3 = sum3.*rhs(R,S);

energy = (W/2)'*(sum1/2.0 - sum3);



end

function energy = assemble_energy(cu,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w,varargin)

if nargin == 9
    sum1 = 0;
    sum3 = 0;
    Nrs = length(R);
    for i = 1:M
        P_i = V_a1b1c1w(:,i);
        for j = 1:M
            dPidPj = reshape(dPdP(i,j,:),Nrs,1);
            sum1 = sum1 + cu(i)*cu(j)*dPidPj;
        end
        sum3 = sum3 + cu(i)*P_i;
    end
    sum2 = sqrt(sum1);
    sum3 = sum3.*rhs(R,S);
    energy = (W/2)'*((xi/2.0)*sum1 + sum2 - sum3);
else
    V_a1b1c1 = varargin{1};
    Lap = varargin{2};
    sum1 = 0;
    sum3 = 0;
    Nrs = length(R);
    for i = 1:M
        P_i = V_a1b1c1w(:,i);
        for j = 1:M
            dPidPj = reshape(dPdP(i,j,:),Nrs,1);
            sum1 = sum1 + cu(i)*cu(j)*dPidPj;
        end
        sum3 = sum3 + cu(i)*P_i;
    end
    sum2 = sqrt(sum1);
    sum3 = rhs(R,S);
    Lapu = V_a1b1c1*(Lap*cu);
    energy = (W/2)'*((sum2 - Lapu - sum3).^2);

end

end

function [err_grad,err_hess,hs] = check_variations(n,cu,M,R,S,W,xi,rhs,varargin)

if nargin == 8
% compute gradient at cu
[K,F,V_abc,V_a1b1c1w,Dx_w,Dy_w,dPdP] = assemble_poisson(n,M,R,S,W,rhs);
% compute energy at cu
pi0 = assemble_energy(cu,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w);
[N,K,F,H0] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,...
                                 V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
grad0 = N + xi*K*cu - F;
% compute random direction uh within (T^2,W_a1b1c1w)
dir = randn(M,1);
% project gradient onto this direction
grad0_dir = grad0'*dir;
% project hessian onto this direction
H0_dir = H0'*dir;
nh = 32; hs = 1e-1*2.^(-(0:nh));
err_grad = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    piph = assemble_energy(cuph,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w);
    err_grad(i) = abs((piph-pi0)/hs(i) - grad0_dir);
end

err_hess = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    [N,K,F] = assemble_eik_nonlin(n,cuph,M,R,S,W,xi,rhs,...
                                  V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
    gradph = N + xi*K*cuph - F;
    gradFD = (gradph-grad0)/hs(i);
    err_hess(i) = norm(gradFD-H0_dir);
end

else
V_a1b1c1 = varargin{1};
Lap = varargin{2};
% compute gradient at cu
[K,F,V_abc,V_a1b1c1w,Dx_w,Dy_w,dPdP] = assemble_poisson(n,M,R,S,W,rhs);
% compute energy at cu
pi0 = assemble_energy(cu,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w,V_a1b1c1,Lap);
[N,K,F,H0] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,...
                                 V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
grad0 = N + xi*K*cu - F;
% compute random direction uh within (T^2,W_a1b1c1w)
dir = randn(M,1);
% project gradient onto this direction
grad0_dir = grad0'*dir;
% project hessian onto this direction
H0_dir = H0'*dir;
nh = 32; hs = 1e-1*2.^(-(0:nh));
err_grad = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    piph = assemble_energy(cuph,M,R,S,W,xi,rhs,dPdP,V_a1b1c1w,V_a1b1c1,Lap);
    err_grad(i) = abs((piph-pi0)/hs(i) - grad0_dir);
end

err_hess = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    [N,K,F] = assemble_eik_nonlin(n,cuph,M,R,S,W,xi,rhs,...
                                  V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F,dPdP);
    gradph = N + xi*K*cuph - F;
    gradFD = (gradph-grad0)/hs(i);
    err_hess(i) = norm(gradFD-H0_dir);
end
end

end

function [err_grad,hs] = check_variations_poiss(n,cu,M,R,S,W,rhs)
% compute energy at cu
pi0 = assemble_energy_poiss(n,cu,M,R,S,W,rhs);
% compute gradient at cu
[K,F] = assemble_poisson(n,M,R,S,W,rhs);
grad0 = K*cu - F;
% compute random direction uh within (T^2,W_a1b1c1w)
dir = randn(M,1);
% project gradient onto this direction
grad0_dir = grad0'*dir;
% project hessian onto this direction
nh = 32; hs = 1e-1*2.^(-(0:nh));
err_grad = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    piph = assemble_energy_poiss(n,cuph,M,R,S,W,rhs);
    err_grad(i) = abs((piph-pi0)/hs(i) - grad0_dir);
end



end

function rhs = RHS(R,S,speed)
rhs = 1./(speed*ones(size(R)));
%rhs(abs(R-1/2)<1e-2) = 0.001;
% rhs = sqrt(S.^2 .* (-R .* (-1 + R + S) .* ...
%             cos(R + S) - (-1 + 2*R + S).* ...
%             sin(R + S)).^2 + R.^2 .* ...
%             (-S .* (-1 + R + S) .* cos(R + S) - ...
%             (-1 + R + 2*S).*sin(R + S)).^2);
end

function val = eik_nonlin_coeffs(cu,cf,Lap,Kap,Dx_w,Dy_w,V_abc,W,xi)

normSq_cf = norm(cf)^2;
normSq_Lapcu = norm(Lap*cu)^2;

dxcu = Dx_w * cu;
dycu = Dy_w * cu;
gradSq_uW = ( (V_abc * dxcu).^2 + (V_abc * dycu).^2 ).*W;
du = Kap * (V_abc' * gradSq_uW);
normSq_du = norm(du)^2;

val = normSq_du - normSq_cf - xi^2*normSq_Lapcu;

end


function val = eik_nonlin_coeffs1(cu,cf,Lap,Dx_w,Dy_w,V_abc,V_a1b1c1,W,xi)
dxcu = Dx_w*cu;
dycu = Dy_w*cu;
gradNorm_u = ((V_abc * dxcu).^2 + (V_abc * dycu).^2).^(0.5);
Lap_u = V_a1b1c1 * (Lap * cu);
f = V_abc * cf;
val = (W/2)'*((gradNorm_u - f - xi*Lap_u).^2);
end