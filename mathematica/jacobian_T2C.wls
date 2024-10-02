#!/usr/bin/env wolframscript
(* ::Package:: *)

(* we define the mapping from standard Simplex to the 
   standard cube in the first quadrant *)
  << ToMatlab.m


mapT2C[x_,y_,z_]:= 
Assuming[
	Element[x, Reals] &&
	Element[y, Reals] &&
	Element[z, Reals],

	
		{x,y,z}*(Sqrt[x^2]+Sqrt[y^2]+Sqrt[z^2]) / Max[Sqrt[x^2],Sqrt[y^2],Sqrt[z^2]]
		

];


JacobianT2C[x_,y_,z_] :=
Assuming[
	Element[x, Reals] &&
	Element[y, Reals] &&
	Element[z, Reals],
	J = ConstantArray[0,{3,3}];
	J[[1]] = D[mapT2C[x,y,z],x];
	J[[2]] = D[mapT2C[x,y,z],y];
	J[[3]] = D[mapT2C[x,y,z],z];
	J
] // Evaluate


(* ::Input:: *)
(**)


JacobianT2C[x,y,z] 
