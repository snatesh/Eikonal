% Given an input matrix A, this function determines an index 
% vector I and a "coefficient matrix" T such that
%
%   A(:,I) \approx A(:,I(1:k))*[I,T]
%
% The function has two modes of operating:
%
% (1) If the second input is smaller than 1, then the function
%     determines the rank k adaptively.
%
% (2) If the second input is larger than of equal to 1, then it
%     sets k = INPUT2 

function [T, I] = id_decomp_hack(A,INPUT2)

if (INPUT2 < 1)
  acc = INPUT2;
  [T,I] = skel_col2(A,acc);
else
  k = INPUT2;
  [T,I] = skel_col3(A,k);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, I] = skel_col2(A,acc)

ss = svd(A);
k  = sum(ss > acc);
[~, R, I] = qr(A,0);

[U,D,V] = svd(R(1:k,1:k));
q = sum(diag(D) > acc);
T = V(:,1:q)*(D(1:q,1:q)\(U(:,1:q)'))*R(1:k, (k+1):size(R,2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, I] = skel_col3(A,k)

[~, R, I] = qr(A,0);
[U,D,V] = svd(R(1:k,1:k));
q = sum(diag(D) > 1e-12);
T = V(:,1:q)*(D(1:q,1:q)\(U(:,1:q)'))*R(1:k, (k+1):size(R,2));

return
