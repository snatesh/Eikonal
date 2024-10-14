function [V,varargout] = jPoly(x,N,a,b)
% Make the Jacobi Vandermonde matrix
% Input:
% N-1 - degree of interpolant
% a,b - Jacoby poly parameters
% Output:
% V - (length(x),N) interp matrix
% scale - normalization for cols
apb = a + b; aa  = a * a; bb  = b * b;
N = N-1;
V = zeros(length(x),N+1); 
V(:,1) = 1;
if N > 0
V(:,2) = 0.5*(2*(a + 1) + (apb + 2)*(x - 1));
end
for k = 2:N
    k2 = 2*k;
    k2apb = k2 + apb;
    q1 =  k2*(k + apb)*(k2apb - 2);
    q2 = (k2apb - 1)*(aa - bb);
    q3 = (k2apb - 2)*(k2apb - 1)*k2apb;
    q4 =  2*(k + a - 1)*(k + b - 1)*k2apb;
    V(:,k+1) = ((q2 + q3*x).*V(:,k) - q4*V(:,k-1)) / q1;
end
% Scaling
if nargout == 2
varargout{1} = structure_factors(N+1,a,b);
end
end

