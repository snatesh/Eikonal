#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll[P,W,Wab,H,h,a,b];

(* first, we generate the 1D polynomials in closed form *)
P[a_,b_,n_,x_] := 
Assuming[

		Element[n,Integers] &&
		Element[a,Reals] &&
		Element[b,Reals] &&
		Element[x,Reals] &&
		n >= 0 &&
		a > -1 && 
		b > -1 &&
		x >= -1 &&
		x <= 1,
		
			((-1)^n)/(2^n * n!) * 
			(1-x)^(-a)*(1+x)^(-b) * 
			D[(1-x)^(a+n)*(1+x)^(b+n),{x,n}]
			
		];
(* Pi,Pj are orthogonal in L2([-1,1],W)= L2mu weighted inner product *)			
W[a_,b_,x_] := 
Assuming[

		Element[a,Reals] &&
		Element[b,Reals] &&
		Element[x,Reals] &&
		a > -1 && 
		b > -1 &&
		x >= -1 &&
		x <= 1,
		
			((1-x)^a) * ((1+x)^b)
			
		];
(* <1,1> = 1 only if we divide by below Wab *)
Wab[a_,b_] := Integrate[W[a,b,x],{x,-1,1}];
(* Use H to check orthogonality, 
and compute structural constants for Pk *)
H[a_,b_,i_,j_] := Integrate[P[a,b,i,x]*P[a,b,j,x]*W[a,b,x], {x,-1,1}]/Wab[a,b];
(* normalize the polynomials *)
p[a_,b_,n_,x_] := P[a,b,n,x] / Sqrt[H[a,b,n,n]];
(* check if we've done everything right so far *)
a = 1; b = 3;
h[i_,j_] := Integrate[p[a,b,i,x]*p[a,b,j,x]*W[a,b,x],{x,-1,1}]/Wab[a,b];
(* expect identity *)
MatrixForm[Table[h[i,j],{i,0,5},{j,0,5}]]
