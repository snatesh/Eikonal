% test a second order differential equation
% right now everything must live in [-1,1]

function test_2ndorder

% N = order of the approximation
N = 10;

CC = cell(2, 2);
Sx = convertmat(N+1,2,1);
Sy = convertmat(N+1,0,1);
Dx = diffmat(N+1,2);
Dy = diffmat(N+1,0);
CC{1,1} = Sy * Dy;
CC{1,2} = Sx * Dx;

Sx = convertmat(N+1,0,1);
Sy = convertmat(N+1,2,1);
Dx = diffmat(N+1,0);
Dy = diffmat(N+1,2);
CC{2,1} = Sy * Dy;
CC{2,2} = Sx * Dx;

Bx = [(-1).^(0:N); ones(1,(N+1))];
By = [(-1).^(0:N); ones(1,(N+1))];
Gx = zeros(2, (N+1), 4*(N+1));
Gy = zeros(2, (N+1), 4*(N+1));
% Canonicalize:
[By, Gy, Py] = canonicalBC(By, Gy);
[Bx, Gx, Px] = canonicalBC(Bx, Gx);
BC = zeros((N+1), (N+1), 4*(N+1));

for k = 1:size(CC, 1)
    [CC{k,1}, BC] = zeroDOF(CC{k,1}, CC{k,2}, BC, By, Gy);
    [CC{k,2}, BC] = zeroDOF(CC{k,2}, CC{k,1}, BC, Bx, Gx, true);
end

% Remove degrees of freedom.
for k = 1:size(CC, 1)
    CC{k,1} = CC{k,1}(1:N-1, 3:(N+1));
    CC{k,2} = CC{k,2}(1:N-1, 3:(N+1));
end

BC = BC(1:N-1, 1:N-1, :);
BC = reshape(BC, (N-1)^2, 4*(N+1));
A = spalloc((N-1)^2, (N-1)^2, (N+1)*(N-1)^2 + (N-1)^2);
for k = 1:size(CC, 1)
    Ak = kron(CC{k,2}, CC{k,1});
    A = A + Ak;
end
rhs = zeros((N-1)^2,1);
rhs(1,:) = -1;

S22 = A \ [BC, rhs];
Gx(:,:,end+1) = zeros(2, (N+1));
Gy(:,:,end+1) = zeros(2, (N+1));
S = imposeBCs(S22, Px, Py, Bx, By, Gx, Gy, (N+1));

bc = zeros(4*(N+1),1);
u = S * [bc ; 1]; 
u = reshape(u, (N+1)*[1,1]);
[xk, yk] = meshgrid(linspace(-1, 1, 200));
nplotpts = 200;
pmax = 100;
Eval = zeros(nplotpts, pmax);
x = linspace(-1, 1, nplotpts).';
c = zeros(pmax, 1);
for k = 1:nplotpts
   c(k) = 1;
   Eval(:,k) = clenshaw_scl(x, c);
   c(k) = 0;
end
[px,py] = size(u);
V = Eval(:,1:py) * u * Eval(:,1:px)';
disp(norm(V));
surf(xk,yk,V);
shading interp;

%keyboard

return

function y = clenshaw_scl(x, c)
% Clenshaw scheme for scalar-valued functions.
bk1 = 0*x; 
bk2 = bk1;
x = 2*x;
n = size(c,1)-1;
for k = (n+1):-2:3
    bk2 = c(k) + x.*bk1 - bk2;
    bk1 = c(k-1) + x.*bk2 - bk1;
end
if ( mod(n, 2) )
    tmp = bk1;
    bk1 = c(2) + x.*bk1 - bk2;
    bk2 = tmp;
end
y = c(1) + .5*x.*bk1 - bk2;
return

function S = imposeBCs(S22, Px, Py, Bx, By, Gx, Gy, n)
%IMPOSEBCS   Impose the boundary conditions on the solution.

By = By * Py;
Bx = Bx * Px;

nc = size(Gx, 3);
S = zeros(n^2, nc);
for k = 1:nc
    % Recombine the boundary conditions.
    X22 = reshape(S22(:,k), n-2, n-2);
    X12 = By(:,1:2) \ (Gy(:,3:n,k) - By(:,3:n)*X22);
    X = [ X12; X22 ];
    X2 = Bx(:,1:2) \ (Gx(:,1:n,k) - Bx(:,3:n)*X.');
    X = [ X2.' X ];
    S(:,k) = X(:);
end

return

function P = nonsingularPermute(B)
%NONSINGULARPERMUTE   Permute the columns of B to ensure that the principal
%m*m submatrix of B is nonsingular, where m = size(B, 1).
%
% Note: This is needed for solving the matrix equations with linear
% constraints, see DPhil thesis of Alex Townsend (section 6.5).

m = size(B, 1);
k = 1;

% [TODO]: improve this check.
% Try each mxm block in a linear fashion:
while ( rank(B(:,k:m+k-1)) < m )
    k = k+1;
    if ( m+k > size(B, 2) )
        error('ULTRASEM:buildSolOp:nonsingularPermute:BCs', ...
            'Boundary conditions are linearly dependent.');
    end
end

P = speye(size(B, 2));
P = P(:,[k:m+k-1, 1:k-1, m+k:end]);

return

function [B, G, P] = canonicalBC(B, G)
%CANONICALBC   Form a linear combination of the boundary conditions
%so that they can be used for imposing on the PDE.

P = nonsingularPermute(B);
B = B*P;
[L, B] = lu(B);

% Scale so that B is unit upper triangular.
if ( min(size(B)) > 1 )
    D = diag(1./diag(B));
elseif ( ~isempty(B) )
    D = 1./B(1,1);
else
    D = []; % No boundary conditions.
end
B = D*B;

for k = 1:size(G,3)
    G(:,:,k) = L \ G(:,:,k);
    G(:,:,k) = D*G(:,:,k);
end

return

function [C1, E] = zeroDOF(C1, C2, E, B, G, trans)
%ZERODOF  Eliminate some degrees of freedom in the matrix equation can be removed.

if ( nargin < 6)
    trans = false;
end

tol = 100*eps;
n = size(C1,1);
perm = [2 1 3];
if ( trans )
    perm = [1 2 3];
end

for ii = 1:size(B, 1) % For each boundary condition, zero a column.
    C1ii = C1(:,ii); % Constant required to zero entry out.
    if ( ~any( abs(C1ii) > tol ) ), continue, end
    C1 = C1 - C1ii*sparse(B(ii,:));
    Gii = permute(G(ii,:,:), [2 3 1]);
    C2Gii = C2*Gii;
    R = repelem(full(C1ii), n, 1) .* repmat(C2Gii, n, 1);
    R = reshape(R, n, n, 4*n);
    R = permute(R, perm);
    E = E - R;
end

C1(abs(C1) < tol) = 0;

return

function b = build_Dir_rhs(n,f,ua,ub)
[~,x] = LOCAL_cheb(n);
% x = fliplr(x')';
fx = f(x);

%make chebychev coefficients
a=fcgltran(fx,1);
% convert to the next order
C = ConvertCheby2ultra1(n+1);
a2 = C*a;
S = convertmat(n+1,1, 1);
a2 = S*a2;

b = [ua; ub;a2(1:end-2)];


return

function A = build_Dir_system(N)
% A is the left hand side of the system
% ud is the left end point boundary condition

m = 2;

D = diffmat(N+1, m);

B = ones(2,N+1);
B(1,2:2:end) = -1;

A = [B;D(1:end-2,:)];

return

function u = eval_ultralam(x,coef,lam)
% x = places to evaluate at
% coef = coeficents of the expansion
% lam = superscript on ultraspherical expansion

n = length(coef)-1;

P = zeros(n+1);

P(:,1) = ones(size(x));

A0 = 2*lam;
P(:,2) = A0*x;


for j = 2:n
    Aj = 2*(j+lam-1)/j;
    Cj = (j+2*lam-2)/j;
    
    P(:,j+1) = Aj*x.*P(:,j)-Cj*P(:,j-1);
end

u = zeros(size(x));

for j = 1:n+1
    u = u+coef(j)*P(:,j);
end

return

function C = ConvertCheby2ultra1(n)
%convert chebychev to ultraspherical of 1st kind

C = diag([1 0.5*ones(1,n-1)])+diag(-1/2*ones(1,n-2),2);

return


function B=fcgltran(A,direction)

% Fast Chebyshev Transform 
%
% Performs the fast transform of data sampled at the
% Chebyshev-Gauss-Lobatto Nodes x=cos(k*pi/N);
%
% A - original data in columns
% B - transformed data in columns
% direction - set equal to 1 for nodal to spectral
%             anything else for spectral to nodal
%
% Written by Greg von Winckel 03/08/04  
% Contact: gregvw@chtm.unm.edu

[N,M]=size(A);

if direction==1 % Nodal-to-spectral
    F=ifft([A(1:N,:);A(N-1:-1:2,:)]);
    B=([F(1,:); 2*F(2:(N-1),:); F(N,:)]);
else            % Spectral-to-nodal
%     F=fft([A(1,:); [A(2:N,:);A(N-1:-1:2,:)]/2]);
F = fft([A(1,:); [A(2:N-1,:); A(N,:)*2; A(N-1:-1:2,:)]/2]);
B=(F(1:N,:));
end
return


function D = diffmat(n, m)
%DIFFMAT   Differentiation matrices for ultraspherical spectral method.
%   D = DIFFMAT(N, M) returns the differentiation matrix that takes N Chebyshev
%   coefficients and returns N C^{(M)} coefficients that represent the derivative
%   of the Chebyshev series. Here, C^{(K)} is the ultraspherical polynomial basis
%   with parameter K.
%  n= number of chebychev coefficients
% m is the order of the derivative.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    m = 1;
end

% Create the differentation matrix.
if ( m > 0 )
    D = spdiags((0 : n - 1)', 1, n, n);
    for s = 1:m-1
        D = spdiags(2*s*ones(n, 1), 1, n, n) * D;
    end
else
    D = speye(n);
end

return


function T = spconvert(n, lam)
% ONLY WORKS WITH ULTRA Spherical.
%SPCONVERT   Compute sparse representation for conversion operators. 
%   CONVERMAT(N, LAM) returns the truncation of the operator that transforms
%   C^{lam} (Ultraspherical polynomials) to C^{lam+1}.  The truncation gives
%   back a matrix of size n x n.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Relation is: C_n^(lam) = (lam/(n+lam))(C_n^(lam+1) - C_{n-2}^(lam+1))

if ( lam == 0 )
    dg = .5*ones(n - 2, 1);
    T = spdiags([1 0 ; .5 0 ; dg -dg], [0 2], n, n);
else
    dg = lam./(lam + (2 : n - 1))';
    T = spdiags([1 0 ; lam./(lam + 1) 0 ; dg -dg], [0 2], n, n);
end

return

function S = convertmat(n, K1, K2)
%CONVERTMAT  Conversion matrix used in the ultraspherical spectral method.
%   S = CONVERTMAT(N, K1, K2) computes the N-by-N matrix realization of the
%   conversion operator between two bases of ultrapherical polynomials.  The
%   matrix S maps N coefficients in a C^{(K1)} basis to N coefficients in a
%   C^{(K2 + 1)} basis, where, C^{(K)} denotes ultraspherical polynomial basis
%   with parameter K.  If K2 < K1, S is the N-by-N identity matrix.
%
%   This function is meant for internal use only and does not validate its
%   inputs.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Create the conversion matrix.
S = speye(n);
for s = K1:K2
    S = spconvert(n, s) * S;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Chebyshev nodes on [-1,1].
% It also computes a differentiation operator.
% page 54 of the text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D,x] = LOCAL_cheb(N)
if N==0
    D=0;
    x=1;
    return
end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
return




