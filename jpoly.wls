#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Generation of Jacobi Matrices on the Simplex*)
(**)


(* ::Subtitle:: *)
(*d=1*)


ClearAll[P,W,Wab,Wabc,H,a,b,c,h,Ax,Bx,Ay,By,An1,An2,Bn1,Bn2,x,y];

(* Orthogonal Jacobi polynomials in 1D *)
P[a_,b_,n_,x_] := 
Assuming[

	Element[n, Integers] &&
	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[x, Reals] &&
	n >= 0 &&
	a > -1 && 
	b > -1 &&
	x >= -1 &&
	x <= 1,
		
	((-1)^n)/(2^n * n!) * 
	(1-x)^(-a)*(1+x)^(-b) * 
	D[(1-x)^(a+n)*(1+x)^(b+n),{x,n}] 
	
];
		
(*  Pi,Pj are orthogonal in L2([-1,1],W) weighted inner product *)			
W[a_,b_,x_] := 
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[x, Reals] &&
	a > -1 && 
	b > -1 &&
	x >= -1 &&
	x <= 1,
		
	((1-x)^a) * ((1+x)^b)
			
];

(*  <1,1> = 1 only if we divide inner products by Wab *)
Wab[a_,b_] := 
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	a > -1 &&
	b > -1,
	
	Integrate[W[a,b,x],{x,-1,1}]
	
];

(* Use H to check orthogonality, and 
compute structural constants for Pk *)
H[a_,b_,i_,j_] := 
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[i, Integers] &&
	Element[j, Integers] &&
	a > -1 &&
	b > -1 &&
	i >= 0 &&
	j >= 0,

	Integrate[P[a,b,i,x]*P[a,b,j,x]*W[a,b,x], {x,-1,1}]/Wab[a,b]
	
];

(* normalize the polynomials *)
p[a_,b_,n_,x_] := 
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	a > -1 &&
	b > -1 &&
	n >= 0 &&
	x >= -1 &&
	x <= 1,
	
	P[a,b,n,x] / Sqrt[H[a,b,n,n]]
];





(* ::Text:: *)
(* *)


(* check if we've done everything right so far *)
a = 1; b = 3;
h[i_,j_] := Integrate[p[a,b,i,x]*p[a,b,j,x]*W[a,b,x],{x,-1,1}]/Wab[a,b];
(* expect identity *)
MatrixForm[Table[h[i,j],{i,0,5},{j,0,5}]]


(* ::Text:: *)
(* *)


(* ::Subtitle:: *)
(*d = 2*)


(* now move on to Tri *)
ClearAll[a,b,c,x,y,n,k,Sz]
(* Normalization for the measure  *)
Wabc[a_,b_,c_] := 
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&	
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2,
	
	Gamma[a+b+c+3/2]/(Gamma[a+1/2]*Gamma[b+1/2]*Gamma[c+1/2])
	
];
(* Normalized measure *)
W[a_,b_,c_,x_,y_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	x >= 0 &&
	y >= 0 &&
	1-x-y >= 0,
	
	x^(a-1/2)*y^(b-1/2)*(1-x-y)^(c-1/2)*Wabc[a,b,c]
	
];

(* structure constants *)

H[a_,b_,c_,k_,n_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	k >= 0 &&
	n >= 0,

	Piecewise[
	{{Sqrt[Wabc[a,b,c]/((2*n+a+b+c+1/2)*(2*k+b+c)) *
	Gamma[n+k+b+c+1]*Gamma[n-k+a+1/2]*Gamma[k+b+1/2]*Gamma[k+c+1/2] /
	((n-k)!*k!*Gamma[n+k+a+b+c+1/2]*Gamma[k+b+c])], k <= n},{0,k>n}}]
	
];

(* Orthonormal Koornwinder polynomials *)
P[a_,b_,c_,k_,n_,x_,y_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	k >= 0 &&
	n >= 0 &&
	x >= 0 &&
	y >= 0 &&
	1-x-y >= 0,
	
	H[a,b,c,k,n]^(-1) * 
	(P[2*k+b+c,a-1/2,n-k,z]/.z->2*x-1) *
	(1-x)^k *
	(P[c-1/2,b-1/2,k,z]/.z->2*y/(1-x)-1)
	
];






(* ::Text:: *)
(* *)


(* ::Text:: *)
(* *)


(* check if we've done everything right so far *)

(* expect upper trimat*)
MatrixForm[Table[H[1/2,1/2,1/2,k,n],{k,0,10},{n,0,10}]] 
a = 1/2; b = 1/2; c = 1/2 ;
h[i_,j_] := Integrate[P[a,b,c,i,i,x,y]*P[a,b,c,j,j,x,y]*W[a,b,c,x,y],{x,0,1},{y,0,1-x}];
(* expect identity *)
MatrixForm[Table[h[i,j],{i,0,5},{j,0,5}]]


(* our lex order is (x^n,x^n-1y,..,y^n) *)
ClearAll[a,b,c,x,y,n,k,Sz,Jn1,Jn2,B1,B2,inds,inds1,An1,An2,Bn1,Bn2,A1,A2]

P[a_,b_,c_,n_,x_,y_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0 &&
	x >= 0 &&
	y >= 0 &&
	1-x-y >= 0,
	
	Table[P[a,b,c,k,n,x,y],{k,0,n}]
					
]

(* Blocks of Jacobi matrix*)
Ax[a_,b_,c_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	x*P[a,b,c,k,x,y],
	P[a,b,c,k+1,x,y]]*W[a,b,c,x,y],
	{x,0,1},{y,0,1-x}],{k,0,n-2}]

];

Ay[a_,b_,c_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	y*P[a,b,c,k,x,y],
	P[a,b,c,k+1,x,y]]*W[a,b,c,x,y],
	{x,0,1},{y,0,1-x}],{k,0,n-2}]

];

Bx[a_,b_,c_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	x*P[a,b,c,k,x,y],
	P[a,b,c,k,x,y]]*W[a,b,c,x,y],
	{x,0,1},{y,0,1-x}],{k,0,n-1}]

];

By[a_,b_,c_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	y*P[a,b,c,k,x,y],
	P[a,b,c,k,x,y]]*W[a,b,c,x,y],
	{x,0,1},{y,0,1-x}],{k,0,n-1}]

];

MakeBlockDiag[An_,n_]:=
Assuming[

	Element[n,Integers] &&
	n >= 0,
	
	Sz = n *(n+1)/2;
	
	Antmp = BlockDiagonalMatrix[Table[Append[An[[nn]],
	ConstantArray[0,Dimensions[An[[nn,1]]]]],{nn,1,n-1}]] // Normal;
	Join[Table[Join[{0},Antmp[[i]]],{i,1,Sz-1}],
	ConstantArray[0,{1,Sz}]]
	
];

JacobiMat[a_,b_,c_,n_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0,
	
	B1 = BlockDiagonalMatrix[Bx[a, b, c, n]];
	A1 = MakeBlockDiag[Ax[a, b, c, n],n];
	A2 = MakeBlockDiag[Ay[a, b, c, n],n];
	B2 = BlockDiagonalMatrix[By[a, b, c, n]];

	{B1 + A1 + Transpose[A1],
	 B2 + A2 + Transpose[A2]}
];


a = 1/2; b = 1/2; c = 1/2; n =4;
{Jx,Jy} = JacobiMat[a,b,c,n];


Jx // MatrixForm
Jy // MatrixForm



