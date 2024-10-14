clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2; kap = abs(a+b+c);
% normalization for weight
wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
% dimension
d = 2;
% max total poly degree
%n = 16; m = 27;
nm1 = 19; mm1 = 32;
n = nm1+1; m = mm1+1;
fname = strcat('triquadLeg_',num2str(n),'_',num2str(m),'.mat');
% jacobi matrices
[Jn1,Jn2,A1,A2,B1,B2,Hn] = jMatON_tri(n,a,b,c);
Hm = structure_factors_tri(m,a,b,c);
% number of polynomials in source basis
N = nchoosek(n-1+d,n-1)
% number of polynomials in target basis 
M = nchoosek(m-1+d,m-1)


% initialize the nodes for Newton iter
X0 = eig(Jn1+1j*Jn2); 
Yk = imag(X0); Xk = real(X0); 


%[tmp,D] = joint_diag([Jn1,Jn2],1e-8); %joint_diag([Jn1,Jn2], N,d,0.5,1e-8);
% X0 = diag(D(:,1:N)); Y0 = diag(D(:,N+1:end)); 
% Xk = X0; Yk = Y0;

% evaluate Vandermonde on initial nodes
Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
Vm_pinv = pinv(Vm');
% evaluate integrals in terms of numeric val of Vm
Ints = zeros(M,1); Ints(1) = 1;
% initialize the weights for Newton iter
Wk = Vm_pinv(:,1);
% minimize subject to inequality constraints preserving
% positivity of weights and keeping nodes in triangle
F_norm = @(Zk) norm((jPoly_tri(Zk(1:N),Zk(N+1:2*N),Hm,m-1,a,b,c))'*Zk(2*N+1:3*N)-Ints);
F = @(Zk) (jPoly_tri(Zk(1:N),Zk(N+1:2*N),Hm,m-1,a,b,c))'*Zk(2*N+1:3*N)-Ints;
G = zeros(5*N,3*N); 

% x,y,w > 0
G(1:N,1:N) = -eye(N); 
G(N+1:2*N,N+1:2*N) = -eye(N); 
G(2*N+1:3*N,2*N+1:3*N) = -eye(N);
% x < 1
G(3*N+1:4*N,1:N) = eye(N);
% x+y < 1
G(4*N+1:5*N,1:N) = eye(N); G(4*N+1:5*N,N+1:2*N) = eye(N);
% sum(w) = 1
Geq = zeros(1,3*N); Geq(2*N+1:end) = 1;
B = zeros(5*N,1); 
B(3*N+1:end) = 1-eps^(1/6); B(1:3*N) = -eps^(1/6);
Beq = 1;
Zk = [Xk;Yk;Wk];

options = ...
  optimoptions('fmincon',...
               'Display','iter',...
               'Algorithm','sqp',...
               'FiniteDifferenceType','central',...
               'FiniteDifferenceStepSize',1e-6,...
               'MaxFunctionEvaluations',1e10,...
               'MaxIterations',1e10,...
               'ConstraintTolerance',1e-14,...
               'OptimalityTolerance',1e-6,...
               'StepTolerance', 1e-6, ...
               'UseParallel', true)
Fk = F(Zk);
gradFk = zeros(M,3*N);
tol = 15*eps; tol_up = 1e3; maxiter_o = 10; maxiter = 1000;
h = 1e-7; pk = norm(Fk); iter_o = 0; lb = zeros(3*N,1); ub = ones(3*N,1);
while pk>tol && iter_o<maxiter_o
  iter_o = iter_o+1;
  % try solving with fmincon
  Zk = fmincon(F_norm,Zk,G,B,Geq,Beq,lb,ub,[],options); Z=Zk;
  Fk = F(Zk); pk = norm(Fk);
  iter_i = 0;
  rho = 0.9; gam = 1e-4;
  while pk>tol && pk<tol_up && iter_i<maxiter
    iter_i = iter_i + 1;
    Zkph = Zk; Zkmh = Zk;
    % compute gradient
    for jj = 1:3*N
      % evaluate above and below zk
      Zkph(jj) = Zkph(jj)+h;
      Zkmh(jj) = Zkmh(jj)-h;
      Fkph = (jPoly_tri(Zkph(1:N),Zkph(N+1:2*N),Hm,m-1,a,b,c))'*Zkph(2*N+1:3*N) - Ints;
      Fkmh = (jPoly_tri(Zkmh(1:N),Zkmh(N+1:2*N),Hm,m-1,a,b,c))'*Zkmh(2*N+1:3*N) - Ints;
      % compute dFk/dx_j
      gradFk(:,jj) = (Fkph-Fkmh)/(2*h);
      % revert to original zk
      Zkph(jj) = Zkph(jj)-h;
      Zkmh(jj) = Zkmh(jj)+h;
    end
    % step direction
    dZk = -gradFk\Fk;
    % linesearch with wolf conditions
    alph = 1; Zk1 = Zk+alph*dZk;
    Fk1 = (jPoly_tri(Zk1(1:N),Zk1(N+1:2*N),Hm,m-1,a,b,c))' * ...
      Zk1(2*N+1:3*N) - Ints;
    pk1 = norm(Fk1); iter_ii = 0;
    while pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_ii<maxiter
      alph = rho*alph;
      Zk1 = Zk+alph*dZk;
      Fk1 = (jPoly_tri(Zk1(1:N),Zk1(N+1:2*N),Hm,m-1,a,b,c))' * ...
        Zk1(2*N+1:3*N) - Ints;
      pk1 = norm(Fk1); iter_ii = iter_ii+1;
    end
    Zk = Zk1; Fk = Fk1; pk = pk1;
    disp([pk,alph]);  
  end
end


Xk = Zk(1:N); Yk = Zk(N+1:2*N); Wk = Zk(2*N+1:end)'; 
Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
ftest = @(X,Y) sin(X.^2+Y.^2);%.*exp(cos(Y));
disp([F_norm(Zk),Wk*ftest(Xk,Yk),sum(abs(Wk)), cond(Vm)])

tri = [0 1 0 0 0 1];
figure(1);
plot_tri(tri,'k-'); hold on;
X0 = eig(Jn1+1j*Jn2);
plot(X0,'b.');
plot(Xk,Yk,'ro');
save(fname,'Zk','n','m','N','M');


%%
Fk = (jPoly_tri(Zk(1:N),Zk(N+1:2*N),Hm,m-1,a,b,c))'*Zk(2*N+1:3*N)-Ints;



%%
gradFk = zeros(M,3*N);
tol = 15*eps; tol_up = 1e3;
step = 0; h = 1e-7; pk = norm(Fk); 
b0 = 3; bk = b0; q = 0.9; gamk = min(1,bk/pk); tolb = 1e-2; tolb_up = 1e-1;
iter_i = 0; maxiter = 100000; iter_i = 0;
rho = 0.9; gam = 1e-4;
while pk>tol && pk<tol_up && iter_i<maxiter && iter_i<maxiter
  % only compute grad outside of implicit inner loop
  %if iter_i == 0
    Zkph = Zk; Zkmh = Zk;
    % compute gradient
    for jj = 1:3*N
      % evaluate above and below zk
      Zkph(jj) = Zkph(jj)+h;
      Zkmh(jj) = Zkmh(jj)-h;
      Fkph = (jPoly_tri(Zkph(1:N),Zkph(N+1:2*N),Hm,m-1,a,b,c))'*Zkph(2*N+1:3*N) - Ints; 
      Fkmh = (jPoly_tri(Zkmh(1:N),Zkmh(N+1:2*N),Hm,m-1,a,b,c))'*Zkmh(2*N+1:3*N) - Ints; 
      % compute dFk/dx_j
      gradFk(:,jj) = (Fkph-Fkmh)/(2*h);
      % revert to original zk
      Zkph(jj) = Zkph(jj)-h;
      Zkmh(jj) = Zkmh(jj)+h;
    end
    % step direction
    dZk = -gradFk\Fk;
  %end

  % linesearch with wolf conditions
  alph = 1; 
  Zk1 = Zk+alph*dZk;
  Fk1 = (jPoly_tri(Zk1(1:N),Zk1(N+1:2*N),Hm,m-1,a,b,c))' * ...
          Zk1(2*N+1:3*N) - Ints;
  pk1 = norm(Fk1); 
  while(pk1 > norm(Fk+gam*alph*gradFk*dZk))
    %phi0p = gradFk*dZk;
    alph = rho*alph;%norm(alph*phi0p./(2*(Fk1-Fk-phi0p*alph)));
    Zk1 = Zk+alph*dZk;
    Fk1 = (jPoly_tri(Zk1(1:N),Zk1(N+1:2*N),Hm,m-1,a,b,c))' * ...
            Zk1(2*N+1:3*N) - Ints;
    pk1 = norm(Fk1); 
  end
%   % if inside tri 
%  if min(Zk1(1:N))>0 && max(Zk1(1:N))<1 && ...
%     min(Zk1(N+1:2*N))>0 && min(1-Zk1(1:N)-Zk1(N+1:2*N))>0     
    Zk = Zk1; Fk = Fk1; pk = pk1; iter_i = 0;
    disp([pk,alph]);
%   else % replace with prev coord
%     for jj = 1:N
%       if Zk1(jj) < 0
%         Zk1(jj) = Zk(jj);
%       elseif Zk1(jj) > 1
%         Zk1(jj) = Zk(jj);
%       elseif Zk1(N+jj) < 0 
%         Zk1(N+jj) = Zk(N+jj);
%       elseif (1-Zk1(jj)-Zk1(N+jj)) < 0
%         Zk1(N+jj) = Zk(N+jj);
%       end
%     end
%     Zk = Zk1; 
%   end
 
  % step length control
  %gamk = min(1,bk/pk);
  %Zk1 = Zk-gamk*dZk;
  
  % if inside tri 
%   if min(Zk1(1:N))>0 && max(Zk1(1:N))<1 && ...
%      min(Zk1(N+1:2*N))>0 && min(1-Zk1(1:N)-Zk1(N+1:2*N))>0    
%     Fk1 = (jPoly_tri(Zk1(1:N),Zk1(N+1:2*N),Hm,m-1,a,b,c))'*Zk1(2*N+1:3*N) - Ints;
%     %F_fun(Zk1,Ints,N,m,a,b,c);
%     pk1 = norm(Fk1); 
%     % and satisfy conditions
%     if (gamk < 1 && pk1 < pk - bk/2) || ...
%       (gamk == 1 && pk1 < pk^2/(2*bk))   
%       % update
%       Zk = Zk1; Fk = Fk1; pk = pk1;
%       iter_o = iter_o + 1; iter_i = 0;     
%       disp([pk,gamk]);
%       if gamk < tolb || gamk > tolb_up
%         bk = b0;
%       end
%     % otherwise search again
%     else 
%       bk = q*bk;
%       iter_i = iter_i+1;% cycles step length ctrlr
%     end
%   % otherwise search again
%   else  
%     bk = q*bk;
%     iter_i = iter_i+1; % cycles step length ctrlr
%   end
end


%%
Xk = Z(1:N); Yk = Z(N+1:2*N); Wk = Z(2*N+1:end)';

Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);

ftest = @(X,Y) sin(X.^2).*exp(cos(Y));

Wk*ftest(Xk,Yk)
sum(abs(Wk))
%W*(ftest(X,Y).*X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc)


tri = [0 1 0 0 0 1];
figure(1);
subplot(1,2,1);
plot_tri(tri,'r-'); X0 = eig(Jn1+1j*Jn2);
plot(X0,'b.');
subplot(1,2,2);
plot_tri(tri,'r-');
plot(Xk,Yk,'b.');
%%
save(fname,'Z','n','m','N','M');


% 3m^2 is number of initial interpolation points
% n_interp = 30;
% % for each quad in T, get maps to reference quad
% % and Jacobian of map
% [X_fun,Y_fun,Jxy_fun] = get_Param();
% % initial quadrature to evaluate integrals
% [X,Y,W] = ljT2Q_Quad(n_interp,X_fun,Y_fun,Jxy_fun);
% compute J_n using quadrature
% Vn = jPoly_tri(X,Y,n-1,a,b,c);
% kap = abs(a+b+c); Jvr = zeros(N);
% wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
% for ii = 1:N
%   for jj = 1:N
%     Jvr(ii,jj) = W*((X+1j*Y).*Vn(:,ii).*Vn(:,jj) .* ...
%                      X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc);
%   end
% end


% J = zeros(2*M);
% J(1:M,1:M) = Jx1; J(1:M,M+1:end) = -Jx2;
% J(M+1:end,1:M) = Jx2; J(M+1:end,M+1:end) = Jx1;
% 
% plot_quad(xya,'ks-');
% plot_quad(xyb,'ks-');
% plot_quad(xyc,'ks-');

% evaluate integrals
%Vm = jPoly_tri(X,Y,m,a,b,c); 
% for jj = 1:M 
%   Ints(jj) = W*(Vm(:,jj).*X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc);  
% end


%compute J_n using quadrature
% Vn = jPoly_tri(X,Y,n-1,a,b,c);
% kap = abs(a+b+c); Jvr = zeros(N);
% wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
% for ii = 1:N
%   for jj = 1:N
%     Jvr(ii,jj) = W*((X+1j*Y).*Vn(:,ii).*Vn(:,jj) .* ...
%                      X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc);
%   end
% end

%%

function [X,Y,W] = ljT2Q_Quad(m,X_fun,Y_fun,Jxy_fun)
  tri = [0 1 0 0 0 1];
  % barycenter coords
  x_cent = mean(tri(1:3));
  y_cent = mean(tri(4:end));
  % quad a coords
  xya1 = [tri(1),tri(4)];
  xya2 = (xya1 + [tri(2),tri(5)]) / 2;
  xya3 = [x_cent,y_cent];
  xya4 = (xya1 + [tri(3),tri(6)]) / 2;
  xya = [xya1;xya2;xya3;xya4];
  % quad b coords
  xyb1 = xya2;
  xyb2 = [tri(2),tri(5)];
  xyb3 = (xyb2 + [tri(3),tri(6)]) / 2;
  xyb4 = xya3;
  xyb = [xyb1;xyb2;xyb3;xyb4];
  % quad c coords
  xyc1 = xya3;
  xyc2 = xyb3;
  xyc3 = [tri(3),tri(6)];
  xyc4 = xya4;
  xyc = [xyc1;xyc2;xyc3;xyc4];
  Q2QRef = 0.25 * [1 1 1 1; -1 1 1 -1; -1 -1 1 1; 1 -1 1 -1];
  % quad ref coords
  alpha_a = Q2QRef * xya;
  alpha_b = Q2QRef * xyb;
  alpha_c = Q2QRef * xyc;

  [~,~,wa,Xa,Ya] = eval_Param(X_fun,Y_fun,Jxy_fun,alpha_a,m);
  [~,~,wb,Xb,Yb] = eval_Param(X_fun,Y_fun,Jxy_fun,alpha_b,m);
  [~,~,wc,Xc,Yc] = eval_Param(X_fun,Y_fun,Jxy_fun,alpha_c,m);
  X = [Xa;Xb;Xc]; Y = [Ya;Yb;Yc]; W = [wa,wb,wc];
end

function [rr,ss,w,Xa,Ya] = eval_Param(X_fun,Y_fun,Jxy_fun,quad,N)
  [X,Wx] = gjQuad(N,0,0); 
  % remove endpoints and take product
  [rr,ss] = meshgrid(X); rr = rr(:); ss = ss(:);
  Xa = X_fun(rr,ss,quad(:,1));
  Ya = Y_fun(rr,ss,quad(:,2));
  [wx,wy] = meshgrid(Wx); 
  detJxy = zeros(length(X));
  for ii = 1:length(rr)
    detJxy(ii) = det(Jxy_fun(rr(ii),ss(ii),quad(:,1),quad(:,2)));
  end
  w = wx(:).*wy(:).*detJxy(:); w = w';
end

function [X_fun, Y_fun, Jxy_fun] = get_Param()
  % get jacobian factors for pde
  syms x y r s alpha0 alpha1 alpha2 alpha3 beta0 beta1 beta2 beta3;
  X_eq(r,s,alpha0,alpha1,alpha2,alpha3,beta0,beta1,beta2,beta3) = alpha0 + alpha1 * r + alpha2 * s + alpha3 * r .* s;
  Y_eq(r,s,alpha0,alpha1,alpha2,alpha3,beta0,beta1,beta2,beta3) = beta0 + beta1 * r + beta2 * s + beta3 * r .* s;
  drX = diff(X_eq,r);
  dsX = diff(X_eq,s);
  drY = diff(Y_eq,r);
  dsY = diff(Y_eq,s);
  Jxy = [drX dsX; drY dsY];
  % detJxy = det(Jxy);
  % invdet = 1/detJxy;
  % dxR = invdet .* dsY;
  % dyR = -invdet .* dsX;
  % dxS = -invdet .* drY;
  % dyS = invdet .* drX;
  % Jrs = [dxR dyR; dxS dyS];
  X_fun_tmp = matlabFunction(X_eq); 
  X_fun = @(r,s,alpha) X_fun_tmp(r,s,alpha(1),alpha(2),alpha(3),alpha(4),0,0,0,0);
  Y_fun_tmp = matlabFunction(Y_eq); 
  Y_fun = @(r,s,beta) Y_fun_tmp(r,s,0,0,0,0,beta(1),beta(2),beta(3),beta(4));
  Jxy_fun_tmp = matlabFunction(Jxy);
  Jxy_fun = @(r,s,alpha,beta) Jxy_fun_tmp(r,s,alpha(1),alpha(2),alpha(3),alpha(4),...
                                                beta(1),beta(2),beta(3),beta(4));
end

function plot_quad(xa, spec)
xa1 = xa(1,:); xa2 = xa(2,:);
xa3 = xa(3,:); xa4 = xa(4,:);
plot([xa1(1),xa2(1)],[xa1(2),xa2(2)],spec,'markerfacecolor',spec(1)); hold on;
plot([xa2(1),xa3(1)],[xa2(2),xa3(2)],spec,'markerfacecolor',spec(1));
plot([xa3(1),xa4(1)],[xa3(2),xa4(2)],spec,'markerfacecolor',spec(1));
plot([xa4(1),xa1(1)],[xa4(2),xa1(2)],spec,'markerfacecolor',spec(1));
end

function plot_tri(tri,col)
plot([tri(1),tri(2)], [tri(4),tri(5)], col); hold on;
plot([tri(2),tri(3)], [tri(5),tri(6)], col);
plot([tri(3),tri(1)], [tri(6),tri(4)], col);
end

%check inner products
% kap = abs(a+b+c); inner_prods = zeros(M);
% wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
% for jj = 1:M
%   for ii = 1:M
%     inner_prods(ii,jj) = W*(V(:,jj).*V(:,ii).*X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc);  
%   end
% end


% this is how many cols we actually need from Vandermonde
% the remaining are for P_n contribution to Jn
% M = nchoosek(n-1+d,n-1);
% V = jPoly_tri(X,Y,n,a,b,c);
% Vnm1 = V(:,1:M);
% [Jn1,Jn2,A1,A2,~,~,H] = jMatON_tri(n,a,b,c);
% Jx1 = Jn1 + (Vnm1'*Vnm1)\(Vnm1'*[zeros(3*m^2,M-n), V(:,M+1:end)*A1{n}']);
% Jx2 = Jn2 + (Vnm1'*Vnm1)\(Vnm1'*[zeros(3*m^2,M-n), V(:,M+1:end)*A2{n}']);
% 

% compute J_n using quadrature
% kap = abs(a+b+c); Jvr = zeros(M);
% wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
% for ii = 1:M
%   for jj = 1:M
%     Jvr(ii,jj) = W*((X+1j*Y).*V(:,ii).*V(:,jj) .* ...
%                      X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc);
%   end
% end
