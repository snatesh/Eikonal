#!/usr/bin/env wolframscript
(* ::Package:: *)

Quit


u[x1_,x2_] := 
Assuming[

	Element[x1, Reals] &&
	Element[x2, Reals]  &&
	x1 >= 0 && 
	x2 >= 0 && 
	x1 <= 1 &&
	x2 <= 1-x1,

	Min[x1,x2,
		EuclideanDistance[
			{x1,x2},{(x1-x2+1)/2, (x2-x1+1)/2}
						]
		]
];


Plot3D[u[x1,x2],{x1,0,1},{x2,0,1-x1}]


Eik[x1_,x2_]:=
Assuming[
	Element[x1, Reals] &&
	Element[x2, Reals]  &&
	x1 >= 0 && 
	x2 >= 0 && 
	x1 <= 1 &&
	x2 <= 1-x1,
	
	(D[u[x,y],x]^2+D[u[x,y],y]^2) /. x ->x1 /. y->x2
	
]





Plot3D[u[x1,x2],{x1,0,1},{x2,0,1-x1}]
Plot3D[Eik[x1,x2],{x1,0,1},{x2,0,1-x1}]



