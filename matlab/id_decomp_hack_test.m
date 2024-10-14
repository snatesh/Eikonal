A = rand(100,10);
B = rand(10,200);
C = A*B;
[T,J] = id_decomp_hack(C',10);
Cp = C(J(1:10),:)' * [eye(10),T];
%keyboard;