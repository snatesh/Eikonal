clear all; close all; clc;
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
  default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
  set(groot, default_name,'latex');
end
set(groot, 'defaultLegendFontSize',30)
set(groot, 'defaultAxesFontSize',30)

% primal dual algorithm for the full quadrature problem.

% jacobi poly params
a = 0.5; b = 0.5; c = 0.5; dim = 2;
% source and target poly total degree (+1)
n = 7; m = 9; d = 2;
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

% 
dU = dPTP([Xk,Yk],n,m,a,b,c);
U = PTP([Xk,Yk],n,m,a,b,c);
dg = dgU(U);
dgdz = dgUdz([Xk,Yk],n,m,a,b,c);
v = inv(U);
dh = dhV(v);
dV = dvZ([Xk,Yk],n,m,a,b,c);
dhdz = dhVdz([Xk,Yk],n,m,a,b,c);
d2U = d2PTP([Xk,Yk],n,m,a,b,c);
d2gdz2 = d2gUdz2([Xk,Yk],n,m,a,b,c);

Z = [Xk;Yk];
phi = phi_Z(Z,n);
dphi = dphi_Z(Z,n);
d2phi = d2phi_Z(Z,n);

hs = [5e-1 4e-1 3e-1 2e-1 1e-1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4 5e-5 1e-5 5e-6 1e-6 5e-7 1e-7];
errsdphi = zeros(size(hs));
errsd2phi = zeros(size(hs));
for kk = 1:length(hs)
    dphi_FD = dphi_Z_FD(Z,n,hs(kk));
    d2phi_FD = d2phi_Z_FD(Z,n,hs(kk));
    errsdphi(kk) = norm(dphi-dphi_FD)/norm(dphi)
    errsd2phi(kk) = norm(d2phi-d2phi_FD,'fro')/norm(d2phi,'fro')
end

loglog(hs,errsdphi,'r-','linewidth',5,'DisplayName','$\nabla\phi$'); hold on;
loglog(hs,errsd2phi,'b-','linewidth',5,'DisplayName','$\nabla^2\phi$');
loglog(hs,hs.^2,'k--','linewidth',5,'DisplayName', '$\mathcal{O}(h^2)$');
legend('location','best')
xlabel('$h$');
ylabel('Relative 2-Norm Error')


% 
% Xk = X(1:N); Yk = X(N+1:2*N); 
% Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
% Vm_pinv = pinv(Vm');
% Wk = Vm_pinv(:,1);
% ftest = @(X,Y) sin(X.^2+Y.^2);%.*exp(cos(Y));
% disp([Wk'*ftest(Xk,Yk),sum(Wk), cond(Vm)]);
% tri = [0 1 0 0 0 1];
% figure(1);
% X0 = eig(Jn1+1j*Jn2);
% plot(X0,'b.');
% plot(Xk,Yk,'ro');
% Zk = [Xk;Yk;Wk];
% %save(fname,'Zk','n','m','N','M');

%%
usewolfe = false;
tol = 1e-14; rho = 0.9; gam = 1e-8; %h = 1e-3;
pk = tol+1; Zk = [Xk;Yk];
iter = 0; maxiter = 10000; go = true;
while pk>tol && iter < maxiter && go
  iter = iter + 1;
  Fk = Ftobj(Xk,Yk,n,m,a,b,c);
  % compute gradient
  gradFk = gradFtobjFD(Xk,Yk,n,m,a,b,c,h);
  dZk = -gradFk\Fk;
  % linesearch with wolfe conditions
  alph = 0.00001;
  Zk1 = Zk+alph*dZk;
  Xk1 = Zk1(1:N);
  Yk1 = Zk1(N+1:2*N);
  Fk1 = Ftobj(Xk1,Yk1,n,m,a,b,c);
  pk1 = norm(Fk1); iter_i = 0;
  if (usewolfe)
    while (pk1>norm(Fk+gam*alph*gradFk*dZk) && iter_i < maxiter && alph>eps)
      alph = rho*alph;
      Zk1 = Zk+alph*dZk;
      Xk1 = Zk1(1:N);
      Yk1 = Zk1(N+1:2*N);
      Fk1 = Ftobj(Xk1,Yk1,n,m,a,b,c);
      pk1 = norm(Fk1);
      iter_i = iter_i+1;
      disp([pk1,alph])
    end
  end

  for j = 1:N
    if (Xk1(j) < 0 || Xk1(j) > 1 || ...
        Yk1(j) < 0 || Yk1(j) > 1 || ...
        Yk1(j) > 1-Xk1(j))
      go = false;
      break;
    end
  end

  %if (norm(Zk-Zk1)/norm(Zk) < 1e-14)
  %  break;
  %end
  if (go)
      Xk = Xk1;
      Yk = Yk1;
      Zk = [Xk;Yk];
      Fk = Fk1;
      pk = pk1;
      disp([pk,alph])
  end
end




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

function [f, df, d2f] = ftobj(X,Y,n,m,a,b,c)

N = n*(n+1)/2;
H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
one = ones(N,1);
Q = V_abc*V_abc';
f = -0.5*(one'*(Q\one) - 1);
df = 0.5*reshape(dgUdz([X,Y],n,m,a,b,c),N*2,1);
d2f = 0.5*reshape(d2gUdz2([X,Y],n,m,a,b,c),N*2,N*2); 
end

function [f, df, d2f] = ftobj_Z(Z,n,m,a,b,c)

N = n*(n+1)/2;
[f, df. d2f] = ftobj(Z(1:N),Z(N+1:2*N),n,m,a,b,c);

end

function [tfphi, dtfphi, d2tfphi] = tftobj_phi_Z(t,Z,n,m,a,b,c)

[f,df,d2f] = ftobj_Z(Z,n,m,a,b,c);
[phi,dphi,d2phi] = phi_Z(Z,n);
tfphi = t*f + phi;
dtfphi = t*df + dphi;
d2tfphi = t*d2f + d2phi;

end

function [ieqc,G,b] = ftieqc(X,Y,n)

dim = 2;
N = n*(n+1)/2;
G = zeros((dim+1)*N,dim*N);
one = ones(dim*N,1);
G(1:dim*N,1:dim*N) = -diag(one);
zeroOne = [zeros(dim*N,1); ones(N,1)];
I = eye(N);
for j = 1:dim
    G(dim*N+1:(dim+1)*N,N*(j-1)+1:N*j) = I;
end
b = zeroOne;
ieqc = G*[X;Y] - b;

end

function [ieqc,G,b] = ftieqc_Z(Z,n)

N = n*(n+1)/2;
[ieqc,G,b] = ftieqc(Z(1:N),Z(N+1:2*N),n);

end

function [phi,dphi,d2phi] = phi_Z(Z,n)

[ieqc,G,b] = ftieqc_Z(Z,n);
phi = -sum(log(-ieqc));
dphi = dphi_Z(Z,n);
d2phi = d2phi_Z(Z,n);

end

function dphi = dphi_Z(Z,n)

[ieqc,G,~] = ftieqc_Z(Z,n);
dphi = -G'*(1./ieqc); 
end

function d2phi = d2phi_Z(Z,n)

dim = 2;
N = n*(n+1)/2;
[ieqc,G,~] = ftieqc_Z(Z,n);
d2phi = zeros(dim*N,dim*N);
invftsq = (1./ieqc).^2;
d2phi = G'*(invftsq.*G);

end

function dphi = dphi_Z_FD(Z,n,h)

dim = 2;
N = n*(n+1)/2;
dphi = zeros(dim*N,1);
for j = 1:(dim*N)
  % phi(z+h)
  Z(j) = Z(j) + h;
  phiph = phi_Z(Z,n);
  % phi(z-h)
  Z(j) = Z(j) - 2*h;
  phimh = phi_Z(Z,n);
  % FD approx to dphi/dzj
  dphi(j) = (phiph-phimh)/(2.0*h);
  % restore z
  Z(j) = Z(j) + h;
end

end

function d2phi = d2phi_Z_FD(Z,n,h)

dim = 2;
N = n*(n+1)/2;
d2phi = zeros(dim*N,dim*N);
for j = 1:(dim*N)
  % dphi(z+h)
  Z(j) = Z(j) + h;
  dphiph = dphi_Z(Z,n);
  % dphi(z-h)
  Z(j) = Z(j) - 2*h;
  dphimh = dphi_Z(Z,n);
  % FD approx to d2phi/dzjdzk
  d2phi(:,j) = (dphiph-dphimh)/(2.0*h);
  % restore Z
  Z(j) = Z(j) + h;
end

end

function dU = dPTP(Z,n,m,a,b,c)

N = n*(n+1)/2;
M = m*(m+1)/2;
dim = 2;
X = Z(:,1); Y = Z(:,2);
H_abc = structure_factors_tri(m+1,a,b,c);
H_a1bc1 = structure_factors_tri(m+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(m+1,a,b+1,c+1);

Dx = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy = D1_tri(a,b,c,H_abc,H_ab1c1,1);
D = cell(2); D{1} = Dx; D{2} = Dy;

V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
V_a1bc1 = jPoly_tri(X,Y,H_a1bc1,m-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(X,Y,H_ab1c1,m-1,a,b+1,c+1);
V = cell(2); V{1} = V_a1bc1; V{2} = V_ab1c1;

dU = zeros(N,dim,N,N);
for i = 1:N
    for j = 1:dim
        for k = 1:N
            for l = 1:N
                if (i ~= l && i ~= k)
                    dukldzij = 0;
                elseif (i == l && i ~=k)
                    dukldzij = V_abc(k,:)*(V{j}(l,:)*D{j})';
                elseif (i ~= l && i == k)
                    dukldzij = V_abc(l,:)*(V{j}(k,:)*D{j})';
                elseif (i == l && i == k)
                    dukldzij = 2*V_abc(k,:)*(V{j}(k,:)*D{j})';
                end
                dU(i,j,k,l) = dukldzij;
            end
        end
    end
end
       
end

function d2U = d2PTP(Z,n,m,a,b,c)

N = n*(n+1)/2;
M = m*(m+1)/2;
dim = 2;
X = Z(:,1); Y = Z(:,2);

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

P = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
Px = jPoly_tri(X,Y,H_a1bc1,m-1,a+1,b,c+1);
Py = jPoly_tri(X,Y,H_ab1c1,m-1,a,b+1,c+1);
Pxx = jPoly_tri(X,Y,H_a2bc2,m-1,a+2,b,c+2);
Pyy = jPoly_tri(X,Y,H_ab2c2,m-1,a,b+2,c+2);
Pxy = jPoly_tri(X,Y,H_a1b1c2,m-1,a+1,b+1,c+2);
Pyx = Pxy;

D = cell(2); DD = cell(2,2);
V = cell(2); VV = cell(2,2);
D{1} = Dx; D{2} = Dy;
V{1} = Px; V{2} = Py; 
DD{1,1} = Dxx; DD{1,2} = Dyx;
DD{2,1} = Dxy; DD{2,2} = Dyy;
VV{1,1} = Pxx; VV{1,2} = Pyx;
VV{2,1} = Pxy; VV{2,2} = Pyy;

d2U = zeros(N,dim,N,dim,N,N);
for s = 1:N
    for t = 1:dim
        for i = 1:N
            for j = 1:dim
                for k = 1:N
                    for l = 1:N
                        if (i ~= l && i ~= k)
                            d2u = 0;
                        elseif (s ~= l && s ~= k)
                            d2u = 0;
                        elseif (((i == l) && (i == k)) && (s ~=k))
                            d2u = 0;
                        elseif (i == l && i ~= k && s == l)
                            d2u = P(k,:)*(DD{t,j}'*(VV{t,j}(l,:))');
                        elseif (i == l && i ~= k && s == k)
                            d2u = (D{j}'*(V{j}(l,:))')'*(D{t}'*(V{t}(k,:))');
                        elseif (i ~= l && i == k && s == l)
                            d2u = (D{j}'*(V{j}(k,:))')'*(D{t}'*(V{t}(l,:))');
                        elseif (i ~= l && i == k && s == k)
                            d2u = P(l,:)*(DD{t,j}'*(VV{t,j}(k,:))');
                        elseif (i == l && i == k && s == k)
                            d2u = 2*P(k,:)*(DD{t,j}'*(VV{t,j}(k,:))') +...
                                2*((D{j}'*(V{j}(k,:))')'*(D{t}'*(V{t}(k,:))'));     
                       end
                        d2U(s,t,i,j,k,l) = d2u;
                    end
                end
            end
        end
    end
end


end


function U = PTP(Z,n,m,a,b,c)

X = Z(:,1); Y = Z(:,2);
H_abc = structure_factors_tri(m+1,a,b,c);
V_abc = jPoly_tri(X,Y,H_abc,m-1,a,b,c);
U = V_abc*V_abc';

end

function gU = g(U)

[N,~] = size(U);
one = ones(N,1);
gU = -one'*(U\one);

end

function dg = dgU(U)

[N,~] = size(U);
one = ones(N,1);
UTinvOne = (U')\one;
dg = UTinvOne*UTinvOne';

end

function hV = h(V)

[N,~] = size(V);
one = ones(N,1);
Vone = V*one;
hV = Vone*Vone';

end

function dh = dhV(V)

[N,~] = size(V);
one = ones(N,1);
dh = zeros(N,N,N,N);
VoneoneT = (V*one)*one';
oneoneTV = VoneoneT';
for i = 1:N
    for j = 1:N
        for k = 1:N
            for l = 1:N
                if (l ~= j && k ~= j)
                    dhkldvij = 0;
                elseif (l == j && k ~= j)
                    dhkldvij = VoneoneT(k,i);
                elseif (l ~= j && k == j)
                    dhkldvij = oneoneTV(i,l);
                elseif (l == j && k == j)
                    dhkldvij = VoneoneT(k,i) + oneoneTV(i,l);
                end
                dh(i,j,k,l) = dhkldvij;
            end
        end
    end
end

end

function vZ = V(Z,n,m,a,b,c)

U = PTP(Z,n,m,a,b,c);
vZ = inv(U);

end

function dV = dvZ(Z,n,m,a,b,c)

dim = 2;
N = n*(n+1)/2;
dV = zeros(N,dim,N,N);
dU = dPTP(Z,n,m,a,b,c);
invU = V(Z,n,m,a,b,c);
for s = 1:N
    for t = 1:dim
        dudzst = reshape(dU(s,t,:,:),N,N);
        dV(s,t,:,:) = -invU*dudzst*invU;
    end
end

end

function dgdz = dgUdz(Z,n,m,a,b,c)

[N,dim] = size(Z);
dgdz = zeros(size(Z));
U = PTP(Z,n,m,a,b,c);
dg = dgU(U); 
dU = dPTP(Z,n,m,a,b,c);
for i = 1:N
    for j = 1:dim
        du = reshape(dU(i,j,:,:),N,N);
        dgdz(i,j) = trace(dg'*du);
    end
end

end

function d2gdz2 = d2gUdz2(Z,n,m,a,b,c)
[N,dim] = size(Z);
d2gdz2 = zeros(N,dim,N,dim);
dU = dPTP(Z,n,m,a,b,c);
d2U = d2PTP(Z,n,m,a,b,c);
dhdz = dhVdz(Z,n,m,a,b,c);
hV = h(V(Z,n,m,a,b,c));
for s = 1:N
    for t = 1:dim
        for i = 1:N
            for j = 1:dim
                dh = reshape(dhdz(s,t,:,:),N,N);
                du = reshape(dU(i,j,:,:),N,N);
                d2u = reshape(d2U(s,t,i,j,:,:),N,N);
                d2gdz2(s,t,i,j) = trace(dh*du + hV*d2u);
            end
        end
    end
end

end

function dhdz = dhVdz(Z,n,m,a,b,c)

[N,dim] = size(Z);
% z is N x dim
% dV/dzst is N x N
% dh/dV is N x N x N x N
% dh/dzst is N x dim x N x N
dhdz = zeros(N,dim,N,N);
dh = dhV(V(Z,n,m,a,b,c));
dV = dvZ(Z,n,m,a,b,c);
for s = 1:N
    for t = 1:dim
        dvdzst = reshape(dV(s,t,:,:),N,N);
        for k = 1:N
            for l = 1:N
                dhdvkl = reshape(dh(:,:,k,l),N,N);
                dhdz(s,t,k,l) = trace(dhdvkl'*dvdzst);
            end
        end
    end
end

end
