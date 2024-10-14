clear all; close all; clc;

% jacobi poly params
a = 1; b = 1; c = 1;
% dimension
d = 2;
% max total poly degree
n = 5;
% m^2 is number of initial interpolation points
m = n+1;

% maps for quad to tri
[X_fun, Y_fun, dXYdRS_fun] = Q2T();
tri = [0 1 0 0 0 1];
% these are the initial interpolation points
% all inside of T
[rr,ss,w,X,Y] = eval_Param(X_fun,Y_fun,dXYdRS_fun,tri,m);
% create Vandermonde matrix
V = jPoly_tri(X,Y,n,a,b,c); 
% this is how many cols we actually need from Vandermonde
% the remaining are for P_n contribution to Jn
M = nchoosek(n-1+d,n-1);

% check inner products
kap = abs(a+b+c);
wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
for jj = 1:M
  for ii = 1:M
    inner_prods(ii,jj) = w*(V(:,jj).*V(:,ii).*X.^(a-1/2).*Y.^(b-1/2).*(1-X-Y).^(c-1/2)*wabc);  
  end
end




Vnm1 = V(:,1:M);
[Jn1,Jn2,A1,A2,~,~,H] = jMatON_tri(n,a,b,c);
Jx1 = Jn1 + (Vnm1'*Vnm1)\(Vnm1'*[zeros(m^2,M-n), V(:,M+1:end)*A1{n}']);
Jx2 = Jn2 + (Vnm1'*Vnm1)\(Vnm1'*[zeros(m^2,M-n), V(:,M+1:end)*A2{n}']);

J_vr = zeros(M);
for ii = 1:M
  for jj = 1:M
    J_vr(ii,jj) = w*((X+1j*Y).*V(:,ii).*V(:,jj));
  end
end

X1 = eig(Jx1+1j*Jx2);
X2 = eig(Jn1+1j*Jn2);

V1n = jPoly_tri(real(X1),imag(X1),n,a,b,c);
V1nm1 = V1n(:,1:M);

figure(1);
plot_tri(tri,'r-');
plot(real(X1),imag(X1),'k.'); 
plot(real(X2),imag(X2),'b.');
% J = zeros(2*M);
% J(1:M,1:M) = Jx1; J(1:M,M+1:end) = -Jx2;
% J(M+1:end,1:M) = Jx2; J(M+1:end,M+1:end) = Jx1;


function [rr,ss,w,Xa,Ya] = eval_Param(X_fun,Y_fun,Jxy_fun,tri,N)
[X,Wx] = gjQuad(N,0,0); 
[rr,ss] = meshgrid(X); rr = rr(:); ss = ss(:);
Xa = X_fun(rr,ss,tri(1),tri(2),tri(3),tri(4),tri(5),tri(6));
Ya = Y_fun(rr,ss,tri(1),tri(2),tri(3),tri(4),tri(5),tri(6));
[wx,wy] = meshgrid(Wx); 
detJxy = zeros(N);
for ii = 1:N^2
  detJxy(ii) = det(Jxy_fun(rr(ii),ss(ii),tri(:,1),tri(:,2),...
                    tri(3),tri(4),tri(5),tri(6)));
end
w = wx(:).*wy(:).*detJxy(:); w = w';
return
end



function plot_tri(tri,col)
plot([tri(1),tri(2)], [tri(4),tri(5)], col); hold on;
plot([tri(2),tri(3)], [tri(5),tri(6)], col);
plot([tri(3),tri(1)], [tri(6),tri(4)], col);
end
