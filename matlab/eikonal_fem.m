clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

% highest poly degree is m-1, and there are M total polys
m = 10; n = m+1; M = m*(m+1)/2;
speed = 1;
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
F = @(R,S) -speed*ones(size(R));
% coefficients of RHS in (1/2,1/2,1/2)
cF_abc = V_abc(:,1:M)'*((1./F(R,S)).*W);
% promote coefficients to (3/2,3/2,3/2)
cF_a1b1c1 = K_a1b1c_a1b1c1(1:M,1:M)...
    *(K_a1bc_a1b1c(1:M,1:M)...
    *(K_abc_a1bc(1:M,1:M)*cF_abc));
% coeffs of sol in weighted basis (3/2,3/2,3/2)
cu_a1b1c1w = Lap_a1b1c1w_a1b1c1\cF_a1b1c1;
cu_abc = V_abc(:,1:M)'*((V_a1b1c1w*cu_a1b1c1w).*W);

T = delaunay(R,S);
% figure(1);
% subplot(2,1,1);
% trisurf(T,R,S,V_a1b1c1w*cu_a1b1c1w)
% subplot(2,1,2);
% trisurf(T,R,S,V_abc(:,1:M)*cu_abc);
% norm((W/2)'*(F(R,S)-V_a1b1c1*(Lap_a1b1c1w_a1b1c1*cu_a1b1c1w)))

% now do some tests for first derivs in weighted basis
% ops act on coefficients in the weighted basis
f_uw = @(R,S) 1;
f_w = @(R,S) w_a1b1c1(R,S).*f_uw(R,S);
Dx_w = Lx_ab1c_abc*Wx_a1b1c1_ab1c; Dx_w_tr = Dx_w(1:M,1:M);
Dy_w = Ly_a1bc_abc*Wy_a1b1c1_a1bc; Dy_w_tr = Dy_w(1:M,1:M);

% 1/speed = cf(1)
cf_abc = zeros(M,1); cf_abc(1) = 1/speed;

Eik = Dx_w'*Dx_w + Dy_w'*Dy_w; 
Eik = Eik(1:M,1:M);
Lap2 = Lap_a1b1c1w_a1b1c1'*Lap_a1b1c1w_a1b1c1;
xi = 0.01;

fun = @(cu) (cu'*Eik*cu - cf_abc'*cf_abc ...
            - xi^2*cu'*Lap2*cu ...
            - 2*xi*cf_abc'*(Lap_a1b1c1w_a1b1c1*cu))^2;
A = []; b = [];
x0 = cu_a1b1c1w;
Aeq = [];
beq = [];


options = ...
    optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FiniteDifferenceType','central',...
    'FiniteDifferenceStepSize',1e-3,...
    'MaxFunctionEvaluations',1e10,...
    'MaxIterations',1e10,...
    'ConstraintTolerance',1e-6,...
    'OptimalityTolerance',1e-6,...
    'StepTolerance', 1e-6, ...
    'UseParallel', false);


cu_sol = fmincon(fun,x0,A,b,Aeq,beq,[],[],[],options);
u = @(x,y) distToTri(x,y);
U = u(R,S); UW = U.*W;
% figure(2);
% tiledlayout(1,3);
% ax1 = nexttile;
% trisurf(T,R,S,V_a1b1c1w*cu_sol);
% ax2 = nexttile;
% trisurf(T,R,S,U*cf_abc(1))
% ax3 = nexttile;
% trisurf(T,R,S,V_a1b1c1w*cu_a1b1c1w)
% linkaxes([ax1,ax2,ax3],'z')

%%%% now try constructing a finite element discretization
rhs = @(R,S) 1./(speed*ones(size(R)));
% check implementation of weak form by verifying that it 
% is indeed the first variation of poisson energy functional
%cf_abc here is actually cf_a1b1c1w.
[err_grad,hs] = check_variations_poiss(n,cf_abc,M,R,S,W,rhs);
figure()
loglog(err_grad,hs,'r.-'); hold on; loglog(hs,hs,'k--');
% solve poisson problem for weighted coeffs cu_femw
[K,F,V_abc,V_a1b1c1w,Dx_w,Dy_w] = assemble_poisson(n,M,R,S,W,rhs);
cu_femw = K\F;
%energy = assemble_energy(n,cu_femw,M,R,S,W,xi,rhs)
% now check weak form for nonlinear eikonal problem
[err_grad, err_hess, hs] = check_variations(n,cu_femw,M,R,S,W,xi,rhs);
figure()
loglog(hs,err_grad,'r.-'); hold on;
loglog(hs,err_hess,'b.-');
loglog(hs,hs,'k--');


% now do an inf-dim newton stepper for eikonal
[N,K,F,H] = assemble_eik_nonlin(n,cu_femw,M,R,S,W,xi,rhs,...
                                V_abc, V_a1b1c1w,Dx_w,Dy_w,K,F);

max_iter = 1000;
cu = cu_femw;
rtol = 1e-6;
tol = norm(N+xi*K*cu-F)*rtol;

fprintf("It \t Energy \t (g,du) \t ||g||l2\n")
for it = 1:max_iter
    [N,K,F,H] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,...
                                    V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F);
    brhs = N+xi*K*cu-F;
    if norm(brhs) < tol
        fprintf('Converged in %d Newton iterations\n', it);
        break;
    end
    
    du = H'\brhs;
    cu = cu - du;
    energy = assemble_energy(n,cu,M,R,S,W,xi,rhs);
    fprintf("%d \t %1.10f \t %1.10f \t %1.10f\n", it, energy, -brhs'*du, norm(brhs))
end


%%%%%%%%%%%%%%%%%%%%%%%%%% FEM assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    P_if = V_a1b1c1w(:,i).*rhs(R,S);
    F(i) = (W/2)'*P_if;
end

if nargout > 2
    varargout{1} = V_abc;
    varargout{2} = V_a1b1c1w;
    varargout{3} = Dx_w;
    varargout{4} = Dy_w;
end

end

function [N,K,F,H] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,varargin)

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

function energy = assemble_energy(n,cu,M,R,S,W,xi,rhs)

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
sum2 = sqrt(sum1);
sum3 = sum3.*rhs(R,S);

energy = (W/2)'*((xi/2.0)*sum1 + sum2 - sum3);

end

function [err_grad,err_hess,hs] = check_variations(n,cu,M,R,S,W,xi,rhs)

% compute energy at cu
pi0 = assemble_energy(n,cu,M,R,S,W,xi,rhs);
% compute gradient at cu
[K,F,V_abc,V_a1b1c1w,Dx_w,Dy_w] = assemble_poisson(n,M,R,S,W,rhs);
[N,K,F,H0] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs,...
                                 V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F);
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
    piph = assemble_energy(n,cuph,M,R,S,W,xi,rhs);
    err_grad(i) = abs((piph-pi0)/hs(i) - grad0_dir);
end

err_hess = zeros(length(hs),1);
for i = 1:nh
    cuph = cu+hs(i)*dir;
    [N,K,F,~] = assemble_eik_nonlin(n,cuph,M,R,S,W,xi,rhs,...
                                 V_abc,V_a1b1c1w,Dx_w,Dy_w,K,F);
    gradph = N + xi*K*cuph - F;
    gradFD = (gradph-grad0)/hs(i);
    err_hess(i) = norm(gradFD-H0_dir);
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