%% LOOK AT ULTRASPHERICAL 1D PAPER FOR MULT OP

clear all; close all; clc;

Ms = [3,4,15];
gs{1} = @(x) x.^2;
gs{2} = @(x) x.^3;
gs{3} = @(x) sin(x);

for j = 1:length(Ms)
N = 20; M = Ms(j); a = 3; b = 5;
f = @(x) exp(x.*cos(x));
g = gs{j};%@(x) x.^2;
gf = @(x) g(x) .* f(x);
% f lives on N pt grid
[Xf,Wf] = gjQuad(N,a,b);
[Vf,hf] = jPoly(Xf,N,a,b);
% g lives on M pt grid (M <= N)
[Xg,Wg] = gjQuad(M,a,b);
[Vg,hg] = jPoly(Xg,M,a,b);
% N coeffs of f
cf = (Vf.' * (f(Xf) .* Wf.'))./hf;
% M coeffs of g
cg = (Vg.' * (g(Xg) .* Wg.'))./hg;
% coeffs of cgf w/o using mult op
cgf = (Vf.' * (gf(Xf) .* Wf.'))./hf;

% multiplication by x operator
Jx = multByX(Xf,N,a,b);
% multiplication by "any func" operator
Mop = multBy(Jx,M,a,b);
% multiplication by g operator
Mg = multByFunc(cg,Mop);
cgf_new = Mg * cf;
norm(Vf*cgf - gf(Xf))
norm(Vf*cgf_new-gf(Xf))
subplot(1,3,j)
spy(Mg)
title(strcat('$g=x^{',num2str(M-1),'}$'));
end
function Jx = multByX(x,N,a,b)
% Generate the multipilcation operator for 
% coeffs of Jacobi poly expansion

% Input:
% N-1 - max poly deg
% a,b - Jacoby poly parameters
% x - grid from -1 to 1

% Output:
% Multiplication operator M s.t. if f is vector 
% of function values on grid X with f = V*c 
% for Vandermonde V and coeffs c, then
% X.*f = V*(M*c)
% i.e. M*c are the Jacoby coeffs of xf(x)

% Jacobi matrix with p_n contrib
[J, anm1] = jMatOG(N,a,b);
% vandermonde up to Nth deg poly
[V,~] = jPoly(x,N+1,a,b); 
% p_n contrib for mult mat
en = V(:,end)*anm1; 
Jx = J'; Jx(:,end) = Jx(:,end) + V(:,1:end-1) \ en;
end

function Mop = multBy(Jx,M,a,b)
% Generate the operator which can be used 
% to create the multiplication by f operator
% for any g represented by at most a degree M-1 poly
%
% Inputs:
% Jx - NxN multiplication by X op (for deg(f)<= N-1)
% M-1 - max poly degree for g (for f*g, with deg(g) <= M-1)
% a,b - Jacobi poly params
%
% Output:
% Mop - Multiplication by "any g" with deg(g) <= M-1
%       i.e. if c_g are M coeffs of g, 
%               c_f are N coeffs of f, 
%               and I is NxN identity,
%       then M = kron(c_g,I)' * Mop is such that 
%            M*g = c_{gf}, the N coeffs of g*f
N = size(Jx,1);
L = clenshawOM(Jx,M,a,b);
In = eye(N); e0 = (1:M==1)';
Mop = L \ kron(e0,In);
end

function Mg = multByFunc(gm,Mop)
N = size(Mop,2);
In = eye(N); 
Mg = kron(gm,In)' * Mop;
end

function L = clenshawOM(Jx,M,a,b)
% Generate the block lower triangular matrix
% for Clenshaw's algorithm for operator valued input

% Inputs:
% Jx - NxN matrix
% M-1 - max poly degree
% a,b - Jacobi params

% Output:
% L - MN x MN block lower triangular matrix
%     s.t. if c are the M coeffs of f, then
%     f(Jx) = kron(c,I)' * (L \ kron(e0,I))
%     with e0 1st basis vec in R^M 
%     and I the NxN identity    

N = size(Jx,1);
[J,~] = jMatOG(M,a,b); I = speye(N);
L = zeros(N*M,N*M);

cvec = diag(J,-1); 
bvec = diag(J,0);
avec = diag(J,1);
avec = [1;avec]; cvec = [cvec;0];

% fill main block diagonal (NxN blocks)
for j = 1:M
    L((j-1)*N+1:j*N,(j-1)*N+1:j*N) = avec(j) * I;
end

% fill sub1 block diagonal 
for j = 1:(M-1)
    L(j*N+1:(j+1)*N,(j-1)*N+1:j*N) = bvec(j)*I - Jx;
end

% fill sub2 block diagonal
for j = 1:(M-2)
   L((j+1)*N+1:(j+2)*N,(j-1)*N+1:j*N) = cvec(j) * I; 
end
L = sparse(L);
end

function L = clenshawSM(X,M,a,b)
% Generate the lower triangular matrix
% for Clenshaw's algorithm on entire x grid

% Inputs:
% x - grid from -1 to 1
% M-1 - max poly degree
% a,b - Jacobi params

% Output:
% L - length(x) cell array for MxM Clenshaw matrices 
%     s.t. if c are coeffs of f, and v = L{j}' \ c,
%          then f(x(j)) = v(1).


[J,~] = jMatOG(M,a,b);
L = {}; 
cvec = diag(J,-1); 
bvec = diag(J,0);
avec = diag(J,1);
avec = [1;avec]; cvec = [cvec;0];
for j = 1:length(X)
    dm1 = bvec-X(j);
    L{j} = spdiags([cvec dm1 avec],-2:0,M,M);
end
end

function [f,varargout] = clenshawEval(cf,L,op)
if ~op
N = size(L,2);
f = zeros(N,1);
for j = 1:N
   v = L{j}' \ cf; 
   f(j) = v(1);
end
else
    N = length(cf); M = int(size(L,1)/N);
    In = eye(N); e0 = (1:M==1)';
    Mop = L \ kron(e0,In);
    f = kron(cf,In)' * Mop;
    if nargout > 1
        varargout{1} = Mop;
    end
end
end
% %%
% % test multiplication by x
% xf = @(x) x.* f(x);
% 
% [X,W] = gjQuad(N,a,b);
% [V,n] = jPoly(X,N,a,b);
% cf = (V.' * (f(X) .* W.'))./n;
% cxf = (V.' * (xf(X) .* W.'))./n;
% MbyX = multByX(X,N,a,b);
% cxf_new = MbyX*cf;
% 
% norm(V*cxf - xf(X))
% norm(V*cxf_new - xf(X))
% 
% % test multiplication g*f
% gf = @(x) g(x) .* f(x);
% [X,W] = gjQuad(N,a,b);
% [V,n] = jPoly(X,N,a,b);
% cgf = (V.' * (gf(X) .* W.'))./n;
% 
% L = clenshawSM(X,N,a,b);
% gf_new = clenshawEval(cgf,L,false);
% norm(gf(X)-gf_new)
% 
% 
% Mop = multBy(N,M,a,b);
% Mg = multByFunc(g,Mop,N,M,a,b);
% cgf_new = Mg * cf;
% 
% norm(V*cgf - gf(X))
% norm(V*cgf_new-gf(X))
% 
