clear all; close all; clc;

N = 10; Mmax= 6;
A = randn(N,Mmax);
b = randn(N,1);

x = A\b;

[U,S,V] = svd(A);
invS = S;
invS(1:Mmax,:) = diag(1./diag(S));

invA = V*invS'*U';

xsvd = invA*b;