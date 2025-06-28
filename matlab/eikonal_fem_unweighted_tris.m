clear all; close all; clc;
% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
R = importdata("../bin/xtri_N496_n30_M1378_m51.txt");
S = importdata("../bin/ytri_N496_n30_M1378_m51.txt");
W = importdata("../bin/wtri_N496_n30_M1378_m51.txt");

%DT = delaunayTriangulation([0 0; 1 0; 1 1; 0 1; 0.5 0.25; 0.5 0.75]);
%DT = delaunayTriangulation([0 0; 1 0; 0 1; 0.5 0.25; 0.5 0.75]);

%DT = delaunayTriangulation([0 0; 1 0; 0 1]);
%DT = delaunayTriangulation([0 0; 1 0; 1 1]);
DT = delaunayTriangulation([-1 -1; 1 -1; 1 1; -1 1]);
%DT = delaunay_unit_square(4);
%DT = annulus();
%DT = delaunay_disk(10,3);
%DT = delaunay_square_with_hole();

meshT = DT.ConnectivityList;
meshP = DT.Points';
nTri = size(meshT,1);
triplot(DT);

% Parameters
c_hole = [0, 0];       % obstacle center
r0 = 0.3;         % radius
delta = 1e-2;     % transition width
eps_obst = 1e-3;  % minimum speed inside obstacle

% Fixed: smooth dip for obstacle
f = @(x,y) eps_obst + (1 - eps_obst) * 0.5 * (1 + tanh((sqrt((x - c_hole(1)).^2 + (y - c_hole(2)).^2) - r0) / delta));
[xg, yg] = meshgrid(linspace(-1,1,200));
fg = f(xg, yg);
surfl(xg, yg, fg); 

%%
m = 10; n = m+1; M = n*(n+1)/2; nquad = length(W);
goal = [-1,-0.4];

% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% vandermonde under (a,b,c), (a+1,b,c+1), (a,b+1,c+1)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
[u0,cu0] = euclidean_dist2_goal(goal,meshT,meshP,V_abc,R,S,W);
for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    Ixe = IncidenceMatrix(Pts_Ti);
    detJ = det(Ixe);
    XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    u0_T = u0((iTri-1)*nquad+1:(iTri)*nquad);
    cu0_T = cu0((iTri-1)*M+1:M*iTri);
    scatter3(X,Y,u0_T,'.'); hold on;
    scatter3(X,Y,V_abc*cu0_T,'o');
    errs(iTri) = (W'/2)*(V_abc*cu0_T-u0_T);
end
%%
%speed = 1;
m = 10; n = m+1; M = n*(n+1)/2; nquad = length(W);
goal = [-1,0];
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% manufactured sol and corresponding dirchlet data/rhs
udirch = @(x,y) zeros(size(x));
%rhs = @(x,y) RHS(x,y,speed);
rhs = @(x,y) 1./f(x,y);
% finite element discretization for poisson
% on triangulated domain
[V_abc,dxP,dyP,Vl,Vb,Vh,...
    Vl_flip,Vb_flip,Vh_flip,...
    Rl,Sl,Rb,Sb,Rh,Sh,...
    Rl_flip,Sl_flip,Rb_flip,Sb_flip,...
    Rh_flip,Sh_flip,wleg,...
    intVlVl,intVbVb,intVhVh,...
    intVlVl_flip,intVbVb_flip,intVhVh_flip] = preAssemble_poisson1(n,R,S);
nleg = size(wleg,1);

[K_glb,Bint_glb,Bdirch_glb,...
    F_glb,G_glb,meshP,meshT,...
    sharedEdge_tri_map,...
    nintEdge,nbndEdge,...
    bndEdge_tri_map,dPdP] = assemble_poisson_pwc(n,R,S,W,...
                                            Rl,Sl,Rb,Sb,Rh,Sh,...
                                            Rl_flip,Sl_flip, ...
                                            Rb_flip,Sb_flip, ...
                                            Rh_flip,Sh_flip,wleg, ...
                                            V_abc,dxP,dyP,Vl,Vb,Vh,...
                                            Vl_flip,Vb_flip,Vh_flip,...
                                            rhs,udirch,DT,1,goal,H_abc);
%%
% stack the boundary and inter-element continuity constraints
BB = [Bint_glb;Bdirch_glb];
% identify the nLambda independent rows of BB
[~,~,II] = qr(BB','vector');
nLambda = rank(BB);
BB_id = BB(II(1:nLambda),:);
Gbc = [zeros(nleg*nintEdge,1);G_glb];
Gbc = Gbc(II(1:nLambda));
Amat = sparse(zeros(M*nTri+nLambda));
Amat(1:M*nTri,1:M*nTri) = K_glb;
Amat(1:M*nTri,M*nTri+1:end) = BB_id';
Amat(M*nTri+1:end,1:M*nTri) = BB_id;
Gbc = [zeros(nleg*nintEdge,1);G_glb];
Gbc = Gbc(II(1:nLambda));
Fvec = [F_glb;Gbc];
% solve for solution modes on each tri
cu = Amat \ Fvec;
for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    Ixe = IncidenceMatrix(Pts_Ti);
    detJ = det(Ixe);
    XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    %T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u_sol = V_abc*cu_abc;
    scatter3(X,Y,u_sol,'.'); hold on;
end
cu0 = cu;
%%
% there are M total polys
m = 8; n = m+1; M = n*(n+1)/2;
speed = 1;
% manufactured sol and corresponding dirchlet data/rhs
udirch = @(x,y) zeros(size(x));
rhs = @(x,y) RHS(x,y,speed);
% finite element discretization for poisson
% on triangulated domain
[V_abc,dxP,dyP,Vl,Vb,Vh,...
    Vl_flip,Vb_flip,Vh_flip,...
    Rl,Sl,Rb,Sb,Rh,Sh,...
    Rl_flip,Sl_flip,Rb_flip,Sb_flip,...
    Rh_flip,Sh_flip,wleg,...
    intVlVl,intVbVb,intVhVh,...
    intVlVl_flip,intVbVb_flip,intVhVh_flip] = preAssemble_poisson1(n,R,S);
nleg = size(wleg,1);

[K_glb,Bint_glb,Bdirch_glb,...
    F_glb,G_glb,meshP,meshT,...
    sharedEdge_tri_map,...
    nintEdge,nbndEdge,...
    bndEdge_tri_map,dPdP] = assemble_poisson_pwc(n,R,S,W,...
                                            Rl,Sl,Rb,Sb,Rh,Sh,...
                                            Rl_flip,Sl_flip, ...
                                            Rb_flip,Sb_flip, ...
                                            Rh_flip,Sh_flip,wleg, ...
                                            V_abc,dxP,dyP,Vl,Vb,Vh,...
                                            Vl_flip,Vb_flip,Vh_flip,...
                                            rhs,udirch,DT);
% stack the boundary and inter-element continuity constraints
BB = [Bint_glb;Bdirch_glb];
% identify the nLambda independent rows of BB
[~,~,II] = qr(BB','vector');
nLambda = rank(BB);
BB_id = BB(II(1:nLambda),:);
Amat = sparse(zeros(M*nTri+nLambda));
Amat(1:M*nTri,1:M*nTri) = K_glb;
Amat(1:M*nTri,M*nTri+1:end) = BB_id';
Amat(M*nTri+1:end,1:M*nTri) = BB_id;
Gbc = [zeros(nleg*nintEdge,1);G_glb];
Gbc = Gbc(II(1:nLambda));
Fvec = [F_glb;Gbc];
% solve for solution modes on each tri
cu = Amat \ Fvec;
for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    Ixe = IncidenceMatrix(Pts_Ti);
    detJ = det(Ixe);
    XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    %T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u_sol = V_abc*cu_abc;
    scatter3(X,Y,u_sol,'.'); hold on;
end

%% ETDRK2
cu = cu0;
% sim params
xi = 0.1; dtfac = 0.5;
dt = dtfac/(xi*norm(K_glb));
dt_max = 10000;
% DAE system
Netd = @(F,N) F-N;
Letd = sparse(dt*(-xi*K_glb)); % use block sparse storage
% num taylor terms and tol for matrix exp, etdrk2 loop
Ntay_max = 30; 
tay_tol = 1e-15; 
step_tol = 1e-8;
res_tol = 1e-10;
cnstr_tol = 1e-13;
% constraint projector
Pv = eye(M*nTri) - BB_id'*((BB_id*BB_id')\BB_id);

% etdrk2 loop
sol_history = zeros(M*nTri,dt_max+1);
res_history = zeros(dt_max+1,1);
u = cu(1:M*nTri); sol_history(:,1) = u;
% initial residual
[Nu, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                         V_abc, dPdP, ...
                                         meshT, meshP, K_glb);
R_u = Netd(F_glb,Nu) + (1/dt)*Letd*u;
% solve for lagrange multipliers in least-squares
% to reconstruct the correct residual
Lambda = -(BB_id*BB_id')\(BB_id * R_u);
res_norm = norm(R_u + BB_id'*Lambda);
res_history(1) = res_norm;

for tstep = 1:dt_max
    % (i) exp step and projection
    phi0 = phij_etd(0,Letd,u(1:M*nTri),tay_tol,Ntay_max);
    Phi1 = Pv*phi0;
    % (ii) evaluate nonlinearity
    [N_glb1, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                            V_abc, dPdP, ...
                                            meshT, meshP, K_glb);
    N1 = Netd(F_glb,N_glb1);
    % (iii) midpoint update and projection
    phi2 = Phi1 + dt*phij_etd(1,Letd,N1,tay_tol,Ntay_max);
    Phi2 = Pv*phi2;
    % (iv) evaluate nonlinearity at midpoint
    [N_glb2, ~] = assemble_eik_nonlin_unweighted(Phi2, M, W, xi, ...
                                            V_abc, dPdP, ...
                                            meshT, meshP, K_glb);    
    N2 = Netd(F_glb,N_glb2);
    % (v) final update end projection
    u = Phi2 + Pv*(dt*phij_etd(2,Letd,N2-N1,tay_tol,Ntay_max));
    sol_history(:,tstep+1) = u;
    step_delta = norm(sol_history(:,tstep)-sol_history(:,tstep+1))/norm(sol_history(:,tstep));
    [Nu, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                             V_abc, dPdP, ...
                                             meshT, meshP, K_glb);
    R_u = Netd(F_glb,Nu) + (1/dt)*Letd*u;
    % solve for lagrange multipliers in least-squares
    % to reconstruct the correct residual
    Lambda = -(BB_id*BB_id')\(BB_id * R_u);
    res_norm = norm(R_u + BB_id'*Lambda);
    res_history(tstep+1) = res_norm;
    %if res_history(tstep+1) > res_history(tstep)
    %    fprintf(" global residual increase. terminating.\n")
    %    u = sol_history(:,tstep);
    %    break;
    %end
    cnstr_res = BB_id * u;
    cnstr_norm = norm(cnstr_res);
    fprintf("  Relative Δu: %.2e, Residual: %.2e, Constraint: %.2e\n", ...
            step_delta, res_norm, cnstr_norm);    
    % Check convergence
    if (step_delta < step_tol || res_norm < res_tol) && cnstr_norm < cnstr_tol
        fprintf(">>> Converged at step %d.\n", tstep);
        break;
    end

    if ~rem(tstep,10)
        for iTri = 1:nTri
            Pts_Ti = meshP(:,meshT(iTri,:));
            Ixe = IncidenceMatrix(Pts_Ti);
            detJ = det(Ixe);
            XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
            X = XYe(:,1);
            Y = XYe(:,2);
            %T = delaunay(X,Y);
            cu_abc = u((iTri-1)*M+1:M*iTri);
            u_sol = V_abc*cu_abc;
            %trisurf(T,X,Y,u_sol); hold on;
            scatter3(X,Y,u_sol); hold on;
        end
        drawnow;
        hold off;
    end
end



%% IMEX 
dt = 0.001; xi = 0.1;
%dt_factor = 2;
max_tstep = 100;
tol = 1e-7;
LHS = [eye(M*nTri) + dt*xi*K_glb BB_id';...
    BB_id zeros(nLambda)];
% get current residual
u = cu0;
%u = [u; zeros(nLambda,1)];
[N_glb, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                            V_abc, dPdP, ...
                                            meshT, meshP, K_glb);
res_curr = norm([N_glb + xi*K_glb*u(1:M*nTri)-F_glb+BB_id'*u(M*nTri+1:end);...
                BB_id*u(1:M*nTri)]);
for tstep = 1:max_tstep

    % Assemble global system

    [N_glb, ~] = assemble_eik_nonlin_unweighted(u, M, W, xi, ...
                                                V_abc, dPdP, ...
                                                meshT, meshP, K_glb);

    rhs = [u(1:M*nTri) - dt*(N_glb-F_glb); Gbc];


    % solve
    sol = LHS \ rhs;
    
    % get new residual
    res_new = norm([N_glb + xi*K_glb*sol(1:M*nTri)-F_glb+BB_id'*sol(M*nTri+1:end);...
                BB_id*sol(1:M*nTri)]);
    fprintf("tstep %d: res_new = %.3e\t res_old=%.3e\n", tstep, res_new, res_curr);
    % Optionally check convergence
    if norm(u - sol) / max(norm(u), 1e-14) < tol
        break;
    end
    if res_new < res_curr
        res_curr = res_new;
        u = sol;
    else
        fprintf('increase in global residual. terminating\n')
        break;
    end
    for iTri = 1:nTri
        Pts_Ti = meshP(:,meshT(iTri,:));
        Ixe = IncidenceMatrix(Pts_Ti);
        detJ = det(Ixe);
        XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
        X = XYe(:,1);
        Y = XYe(:,2);
        %T = delaunay(X,Y);
        cu_abc = u((iTri-1)*M+1:M*iTri);
        u_sol = V_abc*cu_abc;
        %trisurf(T,X,Y,u_sol); hold on;
        scatter3(X,Y,u_sol); hold on;
    end
    drawnow;
    hold off;
    

end

%% Backward euler (with newton inner loop for implicit solve)
xi = 0.01;
dt = 1;
alph = 0.1;
max_iter = 500;
max_tstep = 10000;
tol = 1e-10;
u = cu;
% Initialize residual history
res_curr = Inf;
for tstep = 1:max_tstep
    u_n = u;
    u_m = u_n;
    fprintf("== Time Step %d ==\n", tstep)
    converged = false;
    for iter = 1:max_iter
        % Assemble global system
        [N_glb, H_glb] = assemble_eik_nonlin_unweighted(u_m, M, W, xi, ...
                                                        V_abc, dPdP, ...
                                                        meshT, meshP, K_glb);
        H_glb = H_glb';


        % Compute residual
        Ru = [N_glb + xi*K_glb*u_m(1:M*nTri) + BB_id'*u_m(M*nTri+1:end) - F_glb;
            BB_id*u_m(1:M*nTri)];

        F = [(u_m(1:M*nTri) - u_n(1:M*nTri)) / dt + Ru(1:M*nTri,:);...
              Ru(M*nTri+1:end,:)];
        JF = zeros(M*nTri + nLambda);
        JF(1:M*nTri, 1:M*nTri) = eye(M*nTri)/dt+H_glb;
        JF(1:M*nTri, M*nTri+1:end) = BB_id';
        JF(M*nTri+1:end, 1:M*nTri) = BB_id;
        du = -JF \ F;
        u_m = u_m + alph*du;

        rel_err = norm(du(1:M*nTri)) / max(norm(u_m(1:M*nTri)), 1e-14);
        %fprintf("  Newton iter %d: ||du||/||u|| = %.2e\t dt = %.2e \t ||R|| = %.2e\n", iter, rel_err,dt,norm(Ru));
        if rel_err < tol
            converged = true;
            break;
        end
    end
    % Compute residual norm for this time step
    res_new = norm(Ru);
    
    if ~converged || res_new > res_curr
        % Reject time step and reduce dt
        dt = dt / 2;
        fprintf("  Rejected step %d: residual = %.3e ↑: dt → %.3e\n", tstep, res_new, dt);
        % retry this time step with smaller dt
        continue;
    else
        % Accept step
        u = u_m;
        res_curr = res_new;
        % Optionally increase dt for efficiency
        dt = min(dt * 1.25, 1); % optional upper bound
        fprintf("  Accepted step %d: residual = %.3e ↓: dt → %.3e\n", tstep, res_new,dt);

    end

    for iTri = 1:nTri

        Pts_Ti = meshP(:,meshT(iTri,:));
        Ixe = IncidenceMatrix(Pts_Ti);
        detJ = det(Ixe);
        XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
        X = XYe(:,1);
        Y = XYe(:,2);
        T = delaunay(X,Y);
        cu_abc = u_m((iTri-1)*M+1:M*iTri);
        u_sol = V_abc*cu_abc;
        trisurf(T,X,Y,u_sol); hold on;
    end
    drawnow;
    hold off;

    if norm(u(1:M*nTri)- u_n(1:M*nTri)) / max(norm(u(1:M*nTri)), 1e-14) < tol
        fprintf("Steady state reached at step %d.\n", tstep);
        break;
    end
end
%%
hold on;
[XX,YY] = meshgrid(linspace(0,1,10));
scatter3(XX(:),YY(:),distanceToUnitSquareBoundary(XX(:),YY(:)));

%% newton + homotopic continuation
xi0 = 0.5;          % Initial penalty
xi_min = 1e-3;      % Do not let xi drop below this
xi = xi0;
xi_decay = 0.5;     % Reduce xi by this factor
cu = cu0;
[N_glb, H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                       V_abc, dPdP, ...
                                       meshT, meshP, K_glb);
H_glb = H_glb';


max_iter = 10000;
rtol = 1e-8;
tol = norm(N_glb + xi*K_glb*cu(1:M*nTri) - BB_id'*cu(M*nTri+1:end) - F_glb) * rtol;

xi_trigger = 10*rtol;  % When residual norm drops below this, reduce xi


fprintf("It \t (g,du) \t\t ||g|| \t\t xi\t\t cond(H)\n");
figure(3);
for it = 1:max_iter

    alph = 0.01;
    % Assemble global system
    HHmat = zeros(M*nTri + nLambda);
    HHmat(1:M*nTri, 1:M*nTri) = H_glb;
    HHmat(1:M*nTri, M*nTri+1:end) = BB_id';
    HHmat(M*nTri+1:end, 1:M*nTri) = BB_id;

    % Compute residual
    brhs = [N_glb + xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end) - F_glb;
            BB_id*cu(1:M*nTri)];

    norm_g = norm(brhs);

    if norm_g < tol
        fprintf("Converged in %d iterations\n", it);
        break;
    end

    % Newton update (no line search)
    du = -HHmat \ brhs;
    cu = cu + alph*du;

    % Update residual and Jacobian
    [N_glb,H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                             V_abc, dPdP, ...
                                             meshT, meshP, K_glb);
    
    brhs = [N_glb + xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end) - F_glb;
           BB_id*cu(1:M*nTri)];
    norm_g = norm(brhs);

    H_glb = H_glb';

    % Dynamically reduce xi if residual is small
    if norm_g < xi_trigger && xi > xi_min
        xi = max(xi * xi_decay, xi_min);
        [N_glb, H_glb] = assemble_eik_nonlin_unweighted(cu, M, W, xi, ...
                                                        V_abc, dPdP, ...
                                                        meshT, meshP, K_glb);
        H_glb = H_glb';
    end

    fprintf("%d \t %.4e \t %.4e \t %.3e\t %.4e\n", it, -brhs'*du, norm_g, xi, cond(HHmat));


    for iTri = 1:nTri

        Pts_Ti = meshP(:,meshT(iTri,:));
        J = IncidenceMatrix(Pts_Ti);
        detJ = det(J);
        XYe = (J * [R,S]' + Pts_Ti(:,1))';
        X = XYe(:,1);
        Y = XYe(:,2);
        %T = delaunay(X,Y);
        cu_abc = cu((iTri-1)*M+1:M*iTri);
        %semilogy(abs(cu_abc)); hold on;
        u = V_abc*cu_abc;
        %trisurf(T,X,Y,u); hold on;
        scatter3(X,Y,u,'.'); hold on; 
    end
    drawnow;
    hold off;


end



%% straight forward newton iteration
cu = [u; zeros(nLambda,1)];
[N_glb,H_glb] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                       V_abc, dPdP,...
                                       meshT,meshP,K_glb);
H_glb = H_glb';


max_iter = 1000;
rtol = 1e-6;
tol = norm(N_glb+xi*K_glb*cu(1:M*nTri) - BB_id'*cu(M*nTri+1:end)-F_glb)*rtol;


fprintf("It \t (g,du) \t ||g||l2 \t cond(H)\n")
for it = 1:max_iter

    alph = 1;
    HHmat = zeros(M*nTri+nLambda);
    HHmat(1:M*nTri,1:M*nTri) = H_glb;
    HHmat(1:M*nTri,M*nTri+1:end) = BB_id';
    HHmat(M*nTri+1:end,1:M*nTri) = BB_id;


    brhs = [N_glb+xi*K_glb*cu(1:M*nTri) + BB_id'*cu(M*nTri+1:end)-F_glb; BB_id*cu(1:M*nTri)];
    if norm(brhs) < tol
        fprintf('Converged in %d Newton iterations\n', it);
        break;
    end

    du = -HHmat\brhs; 
    cu1 = cu + alph*du;

    cu = cu1;
    
    fprintf("%d \t %1.10f \t %1.10f \t %f\n", it, -brhs'*du, norm(brhs),cond(HHmat))

    [N_glb] = assemble_eik_nonlin_unweighted(cu,M,W,xi,...
                                             V_abc, dPdP,...
                                             meshT,meshP,K_glb);
end
figure(3);
for iTri = 1:nTri

    Pts_Ti = meshP(:,meshT(iTri,:));
    J = IncidenceMatrix(Pts_Ti);
    detJ = det(J);
    XYe = (J * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    T = delaunay(X,Y);
    cu_abc = cu((iTri-1)*M+1:M*iTri);
    u = V_abc*cu_abc;
    scatter3(X,Y,u,'.'); hold on;
end
%figure(3);
[XX,YY] = meshgrid(linspace(0,1,10));
scatter3(XX(:),YY(:),distanceToUnitSquareBoundary(XX(:),YY(:)));

%%
%%%%%%%%%%%%%%%%%%%%
function rhs = RHS(R,S,speed)
rhs = 1./(speed*ones(size(R)));
end


function d = distanceToUnitSquareBoundary(x, y)
    % Assumes (x,y) ∈ [0,1]^2
    % Returns the shortest distance to the boundary of the unit square
    
    % Distance to each side
    dx = min(x, 1 - x);
    dy = min(y, 1 - y);
    
    % Minimum distance to any of the four edges
    d = min(dx, dy);
end

function dt = delaunay_unit_square(N)
    % Generate N x N grid in [0,1]^2 with (N+1)^2 points
    x = linspace(-1, 1, N+1);
    [X, Y] = meshgrid(x, x);
    p = [X(:), Y(:)];

    % Generate square elements and split into two triangles
    t = [];
    for i = 1:N
        for j = 1:N
            % Indices of square corners
            n0 = (i-1)*(N+1) + j;
            n1 = n0 + 1;
            n2 = n0 + (N+1);
            n3 = n2 + 1;

            % Split square into two triangles
            t1 = [n0, n1, n3];
            t2 = [n0, n3, n2];
            t = [t; t1; t2];
        end
    end

    % Optional: create Delaunay triangulation object
    dt = delaunayTriangulation(p);
end

function DT_annulus = annulus()
% Parameters
r_inner = 0.5;
r_outer = 1.0;
numPoints = 100;

% Generate points uniformly in the annulus
points = [];
while size(points,1) < numPoints
    % Generate points in square bounding box
    pts = (rand(numPoints*2,1)*2 - 1) * r_outer;
    pts = [pts, (rand(numPoints*2,1)*2 - 1) * r_outer];
    
    % Keep points between inner and outer radius
    r = sqrt(sum(pts.^2, 2));
    mask = (r >= r_inner) & (r <= r_outer);
    
    points = [points; pts(mask,:)];
    points = unique(points,'rows','stable');
end
points = points(1:numPoints,:);

% Create delaunay triangulation
DT = delaunayTriangulation(points);

% Filter triangles: keep only those with centroid inside annulus
triangles = DT.ConnectivityList;
pts = DT.Points;

centroids = (pts(triangles(:,1),:) + pts(triangles(:,2),:) + pts(triangles(:,3),:)) / 3;
r_cent = sqrt(sum(centroids.^2, 2));
keep = (r_cent >= r_inner) & (r_cent <= r_outer);

% Create filtered triangulation
DT_annulus = triangulation(triangles(keep,:), pts);

% Plot
%figure;
%triplot(DT_annulus, 'Color', 'b');
%axis equal;
%title('Delaunay Triangulation of Annulus');
%xlabel('x'); ylabel('y');
end 

function DT = delaunay_disk(N_boundary, N_radial)
    % Generate boundary points on unit circle
    theta = linspace(0, 2*pi, N_boundary + 1);
    theta(end) = [];  % Remove duplicate point at 2*pi
    xb = cos(theta);
    yb = sin(theta);

    % Generate interior points in polar coordinates
    r = linspace(0, 1, N_radial + 1);
    r = r(2:end);  % skip r=0 (handled separately)
    p_interior = [];

    for i = 1:length(r)
        Ni = round(N_boundary * r(i));  % number of points on i-th ring
        thetai = linspace(0, 2*pi, Ni+1); thetai(end) = [];
        xi = r(i) * cos(thetai);
        yi = r(i) * sin(thetai);
        p_interior = [p_interior; xi(:), yi(:)];
    end

    % Combine all points (add center point)
    p = [0, 0; p_interior; xb(:), yb(:)];

    % Delaunay triangulation
    dt = delaunayTriangulation(p);
    t = dt.ConnectivityList;
    p = dt.Points;

    DT = delaunayTriangulation(p);
end

function DT = delaunay_square_with_hole()
% Parameters
r_hole = 0.3;
domain = [-1, 1];
n_circle = 10;      % Points around the hole
n_outer = 20;       % Outer points along square
n_radial = 2;       % Radial rings between hole and square

% --- 1. Points on circular hole boundary ---
theta = linspace(0, 2*pi, n_circle+1)';
theta(end) = [];  % Remove duplicate
x_hole = r_hole * cos(theta);
y_hole = r_hole * sin(theta);

% --- 2. Points on square boundary ---
x_square = [];
y_square = [];

% Bottom edge
x_square = [x_square; linspace(-1, 1, n_outer)'];
y_square = [y_square; -ones(n_outer,1)];

% Right edge
x_square = [x_square; ones(n_outer,1)];
y_square = [y_square; linspace(-1, 1, n_outer)'];

% Top edge
x_square = [x_square; linspace(1, -1, n_outer)'];
y_square = [y_square; ones(n_outer,1)];

% Left edge
x_square = [x_square; -ones(n_outer,1)];
y_square = [y_square; linspace(1, -1, n_outer)'];

% --- 3. Radial rings between hole and square ---
radial_pts = [];
for j = 1:n_radial
    r = r_hole + (1 - r_hole) * j / (n_radial + 1);
    n_ring = round(n_circle * (1 + j/2));  % increase points in outer rings
    t = linspace(0, 2*pi, n_ring + 1)';
    t(end) = [];
    radial_pts = [radial_pts; r * cos(t), r * sin(t)];
end

% --- 4. Combine all points ---
pts = unique([ ...
    [x_hole, y_hole]; ...
    [x_square, y_square]; ...
    radial_pts; ...
], 'rows');

% --- 5. Delaunay triangulation ---
DT_raw = delaunayTriangulation(pts);

% --- 6. Remove triangles inside the hole ---
tri_centers = incenter(DT_raw);
tri_mask = vecnorm(tri_centers, 2, 2) >= r_hole;
DT = triangulation(DT_raw.ConnectivityList(tri_mask, :), DT_raw.Points);


end


function [u0,cu0] = euclidean_dist2_goal(goal,meshT,meshP,V_abc,R,S,W)

nquad = length(W); 
M = size(V_abc,2); 
nTri = size(meshT,1);

u0 = zeros(nTri*nquad,1);
cu0 = zeros(M*nTri,1);
dist2_goal = @(X,Y,goal) sqrt((X-goal(1)).^2 + (X-goal(2)).^2);



for iTri = 1:nTri
    Pts_Ti = meshP(:,meshT(iTri,:));
    Ixe = IncidenceMatrix(Pts_Ti);
    XYe = (Ixe * [R,S]' + Pts_Ti(:,1))';
    X = XYe(:,1);
    Y = XYe(:,2);
    u0_T = dist2_goal(X,Y,goal);
    u0((iTri-1)*nquad+1:iTri*nquad) = u0_T;
    cu0((iTri-1)*M+1:M*iTri) = V_abc(:,1:M)'*(u0_T.*W);
end

end