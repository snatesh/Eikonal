#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Generation of Jacobi Matrices on the Simplex*)


(* ::Subtitle:: *)
(*d=1*)


ClearAll[P,W,Wab,Wabc,H,a,b,c,h,Ax,Bx,Ay,By,A1,A2,B1,B2,x,y,Jx,Jy,Antmp];

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
ClearAll[a,b,c,x,y,n,nn,k,Sz,Jn1,Jn2,B1,B2,inds,inds1,An1,An2,Bn1,Bn2,A1,A2]

(* Vectorization for space of homogeneous polynomials of degree n *)
P[a_,b_,c_,n_,x_,y_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
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

(* Blocks of Jacobi matrix *)
Ax[a_,b_,c_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[n, Integers] &&
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
	Element[n, Integers] &&
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
	Element[n, Integers] &&
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
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	y*P[a,b,c,k,x,y],
	P[a,b,c,k,x,y]]*W[a,b,c,x,y],
	{x,0,1},{y,0,1-x}],{k,0,n-1}]

];




(* generate block recurrence matrices *)
a = 1/2; b = 1/2; c = 1/2; n = 3;
B1 = Bx[a,b,c,n];
A1 = Ax[a,b,c,n];
B2 = By[a,b,c,n];
A2 = Ay[a,b,c,n];



(* assemble into Jacobi matrices for each dim *)
n = 3;
Sz = n *(n+1)/2;
Jn1 = ConstantArray[0,{Sz,Sz}];
Jn2 = ConstantArray[0,{Sz,Sz}];
inds = {1,1}; 
For [nn = 1, nn <= n-1, nn++,
	Jn1[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B1[[nn]];
	Jn2[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B2[[nn]];
	inds1 = inds + {nn,nn+1};
	Jn1[[inds[[1]];;inds[[2]],inds1[[1]];;inds1[[2]]]] = A1[[nn]];
	Jn1[[inds1[[1]];;inds1[[2]],inds[[1]];;inds[[2]]]] = Transpose[A1[[nn]]];
	Jn2[[inds[[1]];;inds[[2]],inds1[[1]];;inds1[[2]]]] = A2[[nn]];
	Jn2[[inds1[[1]];;inds1[[2]],inds[[1]];;inds[[2]]]] = Transpose[A2[[nn]]];
	inds = inds + {nn,nn+1}
];
Jn1[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B1[[n]];
Jn2[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B2[[n]];


(* ::Input:: *)
(*MatrixPlot[Jn1]*)
(*MatrixPlot[Jn2]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Text:: *)
(**)


(* ::Subtitle:: *)
(*d = 3, the tetrahedron awaits!*)


(* Now we define orthogonal polynomials on the tetrahedron 

We proceed analogously to 2D in that we only need to multiply
the 2D polynomials by a specific form of the 1D polynomials
to generate the 3D polynomials. However, we just unroll this 
all in terms of the 1D polynomials.
*)


(* un-normalized measure *)
W[a_,b_,c_,d_,x_,y_,z_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	Element[z, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	x >= 0 &&
	y >= 0 &&
	z >= 0 &&
	1-x-y-z >= 0,
	
	x^(a-1/2)*y^(b-1/2)*z^(c-1/2)*(1-x-y-z)^(d-1/2)
	
];

(* Normalization for the measure  *)
Wabcd[a_,b_,c_,d_] := 
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals]
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2,
	
	Integrate[W[a,b,c,d,x,y,z],{x,0,1},{y,0,1-x},{z,0,1-x-y}]
	
];
(* normalized measure *)
w[a_,b_,c_,d_,x_,y_,z_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	Element[z, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	x >= 0 &&
	y >= 0 &&
	z >= 0 &&
	1-x-y-z >= 0,
	
	x^(a-1/2)*y^(b-1/2)*z^(c-1/2)*(1-x-y-z)^(d-1/2)/Wabcd[a,b,c,d]
	
	
];
a = 1/2; b = 1/2; c = 1/2; d = 1/2;




(* structure constants *)

H[a_,b_,c_,d_,k_,j_,n_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[k, Integers] &&
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	k >= 0 &&
	j >= 0 &&
	n >= 0,

	alphaj = {n-k,k-j,j};
	kj = {a,b,c,d};
	kjl1 = {a+b+c+d,b+c+d,c+d,d};
	alphjl1 = {n,k,j,0};
	aj = {2*k+b+c+d+1/2,2*j+c+d,d-1/2};
	Sqrt[Product[
	Pochhammer[kj[[l]]+1/2,alphaj[[l]]]*
	Pochhammer[kjl1[[l+1]]+(1/2)*(3-l+1),alphjl1[[l]]+alphjl1[[l+1]]]/
	(alphaj[[l]]!*Pochhammer[kjl1[[l]]+(1/2)*(3-l+2),alphjl1[[l]]+alphjl1[[l+1]]])*
	(2*(aj[[l]]+kj[[l]]+alphaj[[l]])+1)/
	(2*(aj[[l]]+kj[[l]]+2*alphaj[[l]])+1),
	{l,1,3}]]
	
];



 (* Orthonormal Jacobi polynomials on the tetrahedron *)
 
 P[a_,b_,c_,d_,n_,k_,j_,x_,y_,z_]:=
 Assuming[
 
    Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[n, Integers] &&
	Element[k, Integers] &&
	Element[j, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	Element[z, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	k >= 0 &&
	n >= 0 &&
	x >= 0 &&
	y >= 0 &&
	z >= 0 &&
	1-x-y-z >= 0,

	H[a,b,c,d,k,j,n]^(-1) *
	(P[2*k+b+c+d+1/2,a-1/2,n-k,q]/.q->2x-1) *
	(1-x)^(k-j)*(P[2*j+c+d,b-1/2,k-j,q]/.q->(2*y/(1-x)-1)) *
	(1-x-y)^j*(P[d-1/2,c-1/2,j,q]/.q->(2*z/(1-x-y)-1))
 ];


(* check that we've done everything right so far *)
htet[i_,j_] := Integrate[P[a,b,c,d,i,0,0,x,y,z]*P[a,b,c,d,j,0,0,x,y,z]*w[a,b,c,d,x,y,z],{x,0,1},{y,0,1-x},{z,0,1-x-y}];
(* expect identity *)
MatrixForm[Table[htet[i,j],{i,0,5},{j,0,5}]]



(* Vectorization for space of homogeneous polynomials of degree n *)

P[a_,b_,c_,d_,n_,x_,y_,z_]:=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[n, Integers] &&
	Element[x, Reals] &&
	Element[y, Reals] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	n >= 0 &&
	x >= 0 &&
	y >= 0 &&
	z >= 0
	1-x-y-z >= 0,
	
	ArrayFlatten[Table[P[a,b,c,d,n,k,j,x,y,z],{k,0,n},{j,0,k}],1]
					
];

(* Blocks of Jacobi matrix *)
Ax[a_,b_,c_,d_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	x*P[a,b,c,d,k,x,y,z],
	P[a,b,c,d,k+1,x,y,z]]*W[a,b,c,d,x,y,z],
	{x,0,1},{y,0,1-x},{z,0,1-x-y}],{k,0,n-2}]

];

Ay[a_,b_,c_,d_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	y*P[a,b,c,d,k,x,y,z],
	P[a,b,c,d,k+1,x,y,z]]*W[a,b,c,d,x,y,z],
	{x,0,1},{y,0,1-x},{z,0,1-x-y}],{k,0,n-2}]

];

Bx[a_,b_,c_,d_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	x*P[a,b,c,d,k,x,y,z],
	P[a,b,c,d,k,x,y,z]]*W[a,b,c,d,x,y,z],
	{x,0,1},{y,0,1-x},{z,0,1-x-y}],{k,0,n-1}]

];

By[a_,b_,c_,d_,n_] :=
Assuming[

	Element[a, Reals] &&
	Element[b, Reals] &&
	Element[c, Reals] &&
	Element[d, Reals] &&
	Element[n, Integers] &&
	a > -1/2 &&
	b > -1/2 &&
	c > -1/2 &&
	d > -1/2 &&
	n >= 0,
	
	Table[Integrate[Outer[Times,
	y*P[a,b,c,d,k,x,y,z],
	P[a,b,c,d,k,x,y,z]]*W[a,b,c,d,x,y,z],
	{x,0,1},{y,0,1-x},{z,0,1-x-y}],{k,0,n-1}]

];



(* generate block recurrence matrices *)
a = 1/2; b = 1/2; c = 1/2; d = 1/2; n = 3;
B1 = Bx[a,b,c,d,n];
A1 = Ax[a,b,c,d,n];
B2 = By[a,b,c,d,n];
A2 = Ay[a,b,c,d,n];



(* assemble into Jacobi matrices for each dim *)
n = 3; 
Sz = Binomial[n-1+3,n-1];
Jn1 = ConstantArray[0,{Sz,Sz}];
Jn2 = ConstantArray[0,{Sz,Sz}];
inds = {1,1}; 
For [nn = 1, nn <= n-1, nn++,
	Jn1[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B1[[nn]];
	Jn2[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B2[[nn]];
	inds1 = inds + {nn,nn+1};
	Jn1[[inds[[1]];;inds[[2]],inds1[[1]];;inds1[[2]]]] = A1[[nn]];
	Jn1[[inds1[[1]];;inds1[[2]],inds[[1]];;inds[[2]]]] = Transpose[A1[[nn]]];
	Jn2[[inds[[1]];;inds[[2]],inds1[[1]];;inds1[[2]]]] = A2[[nn]];
	Jn2[[inds1[[1]];;inds1[[2]],inds[[1]];;inds[[2]]]] = Transpose[A2[[nn]]];
	inds = inds + {nn,nn+1}
];
Jn1[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B1[[n]];
Jn2[[inds[[1]];;inds[[2]],inds[[1]];;inds[[2]]]] = B2[[n]];



