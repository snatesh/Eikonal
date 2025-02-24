function [J,varargout] = jMatON(N,a,b)
% Make the Jacobi coefficient matrix for 
% orthonormal Jacobi polys
% Input:
% N-1 - max poly deg
% a,b - Jacoby poly parameters
% Output:
% Sparse Jacobi matrix J s.t. xP(x) = JP(x) + a_{n-1}p_n(x)e_n,
% where P = (p_0,p_1,...,p_{n-1})
NN = (0:(N-1))';
Avec = (2*NN+a+b+1).*(2*NN+a+b+2)./(2*(NN+1).*(NN+a+b+1));
Bvec = (a^2-b^2)*(2*NN+a+b+1)./(2*(NN+1).*(NN+a+b+1).*(2*NN+a+b)); 
Avec(1) = 0.5*(a+b)+1; 
Bvec(1) = 0.5*(a-b); 
bvec = -Bvec./Avec;
avecON = 2./(a+b+2*NN+2).*((a+NN+1).*(b+NN+1).*(NN+1).*(a+b+NN+1) ...
            ./((a+b+2*NN+1).*(a+b+2*NN+3))).^(1/2);
avecON(1) = 2/(a+b+2)*((a+1)*(b+1)/(a+b+3))^(1/2);
J = spdiags([avecON,bvec,[0;avecON(1:N-1)]],-1:1,N,N); 
if nargout > 1
   varargout{1} = avecON(end); 
end
end
