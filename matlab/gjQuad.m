function [X,w,J] = gjQuad(N,a,b)
% N point Gauss quad rule for jacoby poly with params a,b>-1
% \int_{-1}^1 f dmu ~ w * f(X), where dmu = (1-x)^a(1+x)^b dx.
%
% Input:
% N - num pts in quad rule
% a,b - Jacoby poly parameters
% Output:
% X - quadradture grid
% w - weights
% J - Jacobi matrix
%
% Example Usage: integrate e^xsin(x)dx from -1 to 1 
% N = 10; a = 0; b = 0;
% [X,w,~] = gjQuad(N,a,b)
% f = @(x) exp(x).*sin(x);
% If = (cos(1) + sin(1) + exp(1)^2 * (-cos(1) + sin(1)))/(2*exp(1));
% disp(w*f(X)-If)

% normalization for w(x) = (1-x)^a(1+x)^b on [-1,1]
h0 = 2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
% jacobi matrix
J = jMatON(N,a,b);
% Golub-Welsch with measure normalization
[V,D] = eigs(J,N);
[X,I] = sort(diag(D));
w = h0*V(1,I).^2;
end
