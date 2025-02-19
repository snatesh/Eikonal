clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

% highest poly degree is m-1, and there are M total polys
m = 10; n = m+1; M = m*(m+1)/2;
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
F = @(R,S) -2*ones(size(R));
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
cf_abc = zeros(M,1); cf_abc(1) = 0.5;

Eik = Dx_w'*Dx_w + Dy_w'*Dy_w; 
Eik = Eik(1:M,1:M);
Lap2 = Lap_a1b1c1w_a1b1c1'*Lap_a1b1c1w_a1b1c1;
xi = 0.001;

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
figure(2);
tiledlayout(1,3);
ax1 = nexttile;
trisurf(T,R,S,V_a1b1c1w*cu_sol);
ax2 = nexttile;
trisurf(T,R,S,U*cf_abc(1))
ax3 = nexttile;
trisurf(T,R,S,V_a1b1c1w*cu_a1b1c1w)
linkaxes([ax1,ax2,ax3],'z')

%%%% now try constructing a finite element discretization
rhs = @(R,S) 1./(2*ones(size(R)));
[K,F] = assemble_poisson(n,M,R,S,W,rhs);
% solve poisson problem for weighted coeffs cu_femw
cu_femw = K\F;

%energy = assemble_energy(n,cu_femw,M,R,S,W,xi,rhs)
err = grad_check(n,cu_femw,M,R,S,W,xi,rhs);


[N,K,F,H] = assemble_eik_nonlin(n,cu_femw,M,R,S,W,xi,rhs);
%%


%err = hess_check(n,cu_femw,M,R,S,W,xi,rhs);

max_iter = 1000;
cu = cu_femw;
rtol = 1e-10;
tol = norm(N+xi*K*cu-F)*rtol;


for it = 1:max_iter
    [N,K,F,H] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs);
    if norm(N+xi*K*cu-F) < tol
        disp("Converged in ", it, "Newton iterations");
        break;
    end
    du = H\(N+xi*K*cu-F);
    cu = cu - 0.001*du;
    disp(norm(du))
end


%%%%%%%%%%%%%%%%%%%%%%%%%% FEM assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,F] = assemble_poisson(n,M,R,S,W,rhs)
a = 1/2; b = a; c = a;
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


end

function [N,K,F,H] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs)

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

function err = grad_check(n,cu,M,R,S,W,xi,rhs)
a = 1/2; b = a; c = a;


[N,K,F,~] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs);
grad0 = N + xi*K*cu - F;
disp(assemble_energy(n,cu,M,R,S,W,xi,rhs));
disp(norm(grad0))
hs = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
err = zeros(length(hs),1);

for i = 1:length(hs)
    gradFD = zeros(M,1);
    for j = 1:M
        cuph = cu; 
        cuph(j) = cuph(j) + hs(i);
        cumh = cu;
        cumh(j) = cumh(j) - hs(i);
        piph = assemble_energy(n,cuph,M,R,S,W,xi,rhs);
        pimh = assemble_energy(n,cumh,M,R,S,W,xi,rhs);
        gradFD(j) = (piph - pimh)/(2*hs(i));
    end
    err(i) = norm(gradFD-grad0);
end


end

function err = hess_check(n,cu,M,R,S,W,xi,rhs)

hs = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
err = zeros(length(hs),1);
[N,K,F,H] = assemble_eik_nonlin(n,cu,M,R,S,W,xi,rhs);
dPi_eval = N + xi*K*cu - F;
for i = 1:length(hs)
    cuph = cu+hs(i);
    [Nph,Kph,Fph] = assemble_eik_nonlin(n,cuph,M,R,S,W,xi,rhs);
    dPiph_eval = Nph + xi*Kph*cuph - Fph;
    diff_dPi = (dPiph_eval - dPi_eval)/hs(i);
    H_eval = H*(hs(i)*ones(length(F),1));
    err(i) = norm(diff_dPi - H_eval);
end

end