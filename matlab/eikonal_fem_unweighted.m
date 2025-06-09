clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");
T = delaunay(R,S);

% highest poly degree is m-1, and there are M total polys
m = 10; n = m+1; M = n*(n+1)/2;
speed = 2; 
xi = 0.01;
u = @(x,y) distToTri(x,y)/speed;
tripts = [2 1; 3 5; 1.5 4];
%tripts = [0 0; 1 0; 0 1];
DT = delaunayTriangulation(tripts);
meshT = DT.ConnectivityList;
meshP = DT.Points';
Ixe = IncidenceMatrix(meshP);

%%

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

%% now try constructing a finite element discretization for poisson
rhs = @(R,S) RHS(R,S,speed);
[K,B,F,dPdP,V_abc,Dx,Dy] = assemble_poisson_unweighted(n,R,S,W,rhs,meshT,meshP);
[~,~,II] = qr(B,0);
nLambda = rank(B);
BB = B(II(1:nLambda),:);
Amat = zeros(M+nLambda);
Amat(1:M,1:M) = K;
Amat(1:M,M+1:end) = BB';
Amat(M+1:end,1:M) = BB;
Fvec = [F;zeros(nLambda,1)];
cu = Amat\Fvec; cu_abc = cu(1:M);

tiledlayout(1,1);
nexttile;
XYe = (Ixe * [R,S]' + meshP(:,1))';
scatter3(XYe(:,1),XYe(:,2),V_abc*cu_abc);
%trisurf(T,R,S,V_abc*cu_abc);
%%
%[err_hess,hs] = check_variations(n,cu_abc,M,R,S,W,xi,rhs)

%% check poisson residual
V_a2b2c2 = jPoly_tri(R,S,H_a2b2c2,n-1,a+2,b+2,c+2);
wa2b2c2 = @(R,S) R.^(a+2-0.5).*S.^(b+2-0.5).*(1-R-S).^(c+2-0.5);
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
cf_a2b2c2 = zeros(M,1); cf_a2b2c2(1) = -1/speed;
err1 = (W/2)'*((V_a2b2c2*(Lap_abc_a2b2c2*cu_abc - cf_a2b2c2)).*wa2b2c2(R,S))
err2 = (W/2)'*((V_a2b2c2*(Lap_abc_a2b2c2*cu_abc - cf_a2b2c2)))


%% 
% now do an inf-dim newton stepper for eikonal
[N,K,F,H] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                            V_abc, Dx, Dy,...
                                            K, F, dPdP);
H = H';


max_iter = 1000;
rtol = 1e-6;
tol = norm(N+xi*K*cu(1:M) - BB'*cu(M+1:end)-F)*rtol;


fprintf("It \t (g,du) \t ||g||l2\n")
for it = 1:max_iter

    alph = 1;
    HHmat = zeros(M+nLambda);
    HHmat(1:M,1:M) = H;
    HHmat(1:M,M+1:end) = BB';
    HHmat(M+1:end,1:M) = BB;


    brhs = [N+xi*K*cu(1:M) + BB'*cu(M+1:end)-F; BB*cu(1:M)];
    if norm(brhs) < tol
        fprintf('Converged in %d Newton iterations\n', it);
        break;
    end

    du = -HHmat\brhs; 
    cu1 = cu + alph*du;

    cu = cu1;
    
    fprintf("%d \t %1.10f \t %1.10f \n", it, -brhs'*du, norm(brhs))

    %[N,K,F,H] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
    %                                           V_abc, Dx, Dy,...
    %                                           K, F, dPdP);
    [N,K,F] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                              V_abc, Dx, Dy,...
                                              K, F, dPdP);
end
%%
[xx,yy] = meshgrid(linspace(1e-10,1-1e-10,50)); 
xx = xx(:); yy = yy(:);
indtri = yy<=1-xx;
xx = xx(indtri);
yy = yy(indtri);
V_abc_grid = jPoly_tri(xx,yy,H_abc,n-1,a,b,c);
Tgrid = delaunay(xx,yy);
figure()
tiledlayout(1,2);
ax2 = nexttile;
trisurf(Tgrid,xx,yy,u(xx,yy))

ax3 = nexttile;
trisurf(Tgrid,xx,yy,V_abc_grid*cu(1:M))
linkaxes([ax2,ax3],'z')


%%
% check weak form for nonlinear eikonal problem
% [err_grad,err_hess,hs] = check_variations(n,cu_femw,M,R,S,W,xi,rhs);
% figure()
% loglog(hs,err_grad,'r.-'); hold on;
% loglog(hs,err_hess,'b.-');
% loglog(hs,hs,'k--');


%%



%%%%%%%%%%%%%%%%%%%%%%%%%% FEM assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,B,F,dPdP,V_abc,Dx,Dy] = assemble_poisson_unweighted(n,R,S,W,rhs,meshT,meshP)
a = 1/2; b = a; c = a;
Ixe = IncidenceMatrix(meshP);
detJ = det(Ixe);
invIxe = inv(Ixe);
lenE1 = norm(meshP(:,2)-meshP(:,1));
lenE2 = norm(meshP(:,3)-meshP(:,2));
lenE3 = norm(meshP(:,1)-meshP(:,3));
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
    gradP_i = invIxe'*[dxP_i';dyP_i'];
    dxP_i = gradP_i(1,:)';
    dyP_i = gradP_i(2,:)';
    for j = 1:M
        dxP_j = V_a1bc1 * Dx_a1bc1(:,j);
        dyP_j = V_ab1c1 * Dy_ab1c1(:,j);   
        gradP_j = invIxe'*[dxP_j';dyP_j'];
        dxP_j = gradP_j(1,:)';
        dyP_j = gradP_j(2,:)';
        dPidPj = dxP_i.*dxP_j + dyP_i.*dyP_j;
        K(i,j) = (W/2)'*dPidPj*detJ;
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
        int3 = wleg'*(Vh(:,i).*Vh(:,j))*sqrt(2);
        B(i,j) = int1*lenE1+int2*lenE2+int3*lenE3;
    end
end

F = zeros(M,1);
for i = 1:M
    P_if = V_abc(:,i).*rhs(R,S);
    F(i) = (W/2)'*P_if*detJ;
end

Dx = Dx_a1bc1;
Dy = Dy_ab1c1;

end

function [N,K,F,varargout] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                                            V_abc, Dx, Dy,...
                                                            K, F, dPdP)
N = zeros(M,1); Nrs = length(W);
for i = 1:M
    P_i = V_abc(:,i);
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
    H = zeros(M);
    for i = 1:M
        for j = 1:M
            P_j = V_abc(:,j);
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
    varargout{1} = H;
end

end
  
function [err_hess,hs] = check_variations(n,cu,M,R,S,W,xi,rhs)

% compute gradient at cu
[K,B,F,dPdP,V_abc,Dx,Dy] = assemble_poisson_unweighted(n,R,S,W,rhs);

[N,K,F,H0] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                            V_abc, Dx, Dy,...
                                            K, F, dPdP);

grad0 = N + xi*K*cu - F;
% compute random direction uh within (T^2,W_abc)
dir = randn(M,1);

% project hessian onto this direction
H0_dir = H0'*dir;
nh = 32; hs = 1e-1*2.^(-(0:nh));

err_hess = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    [N,K,F] = assemble_eik_nonlin_unweighted(cuph,M,W,xi,...
                                             V_abc, Dx, Dy,...
                                             K, F, dPdP);
    gradph = N + xi*K*cuph - F;
    gradFD = (gradph-grad0)/hs(i);
    err_hess(i) = norm(gradFD-H0_dir);
end

end


function rhs = RHS(R,S,speed)
rhs = 1./(speed*ones(size(R)));
end

function Ixe = IncidenceMatrix(Xe)

Ixe = [Xe(:,2)-Xe(:,1), Xe(:,3)-Xe(:,1)];

end

function A = triarea(t, p)
% A = TRIAREA(t, p) area of triangles in triangulation
Xt = reshape(p(t, 1), size(t)); % X coordinates of vertices in triangulation
Yt = reshape(p(t, 2), size(t)); % Y coordinates of vertices in triangulation
A = 0.5 * abs((Xt(:, 2) - Xt(:, 1)) .* (Yt(:, 3) - Yt(:, 1)) - ...
    (Xt(:, 3) - Xt(:, 1)) .* (Yt(:, 2) - Yt(:, 1)));
end