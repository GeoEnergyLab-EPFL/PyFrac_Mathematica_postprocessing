(* ::Package:: *)

(* --------------------------------------------- 

Wolfram Language Package created by Carlo Peruzzo on Tuesday 26th October 2020
EPFL ENAC IIC GEL
GC B1 385 (Batiment GC)
Station 18, 1015 Lausanne, Switzerland
Phone: +41 (0)21 693 76 47
carlo.peruzzo@epfl.ch
Geo-energy Laboratory EPFL 

--------------------------------------------- *)
BeginPackage["PostprocessPyFrac`"];

Needs[ "rectangularMesh`"]; (* when loaded via needs the context is not exposed *)
Needs[ "getDOUBLEfrac`"];
Needs[ "getSINGLEfrac`"];
Needs[ "getGENERALfrac`"];

postProcessPyfrac::usage = "postProcessPyfrac[ Simpath,Case] it returns an association";
meshPlan::usage =  "meshPlan[ Lx, Nx, Ly, Ny] It computes: {hx,hy,Nelts,vertexes,conn,centers,Neighbors}";

(*meshPlan::usage =  "meshPlan[ Lx, Nx, Ly, Ny] It computes: conn, centers, Nelts, ElemNei";
getNeighbors::usage = "getNeighbors[ element, Nx, Ny] It returns {left,right,bottom,up} of element.";*)

Begin["`Private`"];
postProcessPyfrac[Simpath_, Case_] := Which[Case == "double_fracture",postprocDOUBLEfrac[Simpath],
											Case == "single_fracture",postprocSINGLEfrac[Simpath],
	                                                             True,Print[" Nothing donee"]];

postprocDOUBLEfrac[Simpath_] := Module[{result, data},
	result = Association[{"Simpath"-> Simpath}];
	data = Import[Simpath, "RawJSON"];
	result = AssociateTo[result, getDOUBLEfrac`getFr[data]];
	
	result = AssociateTo[result, "NumRinFrLeft"     ->getDOUBLEfrac`getNumRinFrLeft[data,result]];
	result = AssociateTo[result, "NumRinFrRight"    ->getDOUBLEfrac`getNumRinFrRight[data,result]];
	result = AssociateTo[result, "NumRoutFrLeft"    ->getDOUBLEfrac`getNumRoutFrLeft[data,result]];
	result = AssociateTo[result, "NumRoutFrRight"   ->getDOUBLEfrac`getNumRoutFrRight[data,result]];
	result = AssociateTo[result, "NumDeltatilde"    ->getDOUBLEfrac`getNumDeltatilde[data,result]];
	
	result = AssociateTo[result, "LastMesh"         ->getGENERALfrac`getLastMesh[data]];
	
	varString="w_at_my_point"; timeString="time_list_W_at_my_point";
	result = AssociateTo[result, "w_A(t)"           ->getGENERALfrac`getVARofT[data,varString,timeString]];
	
	varString="pf_at_my_point_A"; timeString="time_list_pf_at_my_point_A";
	result = AssociateTo[result, "pf_A(t)"           ->getGENERALfrac`getVARofT[data,varString,timeString]];
	
	varString="pf_at_my_point_B"; timeString="time_list_pf_at_my_point_B";
	result = AssociateTo[result, "pf_B(t)"           ->getGENERALfrac`getVARofT[data,varString,timeString]];
	
	
	
	result
	]

postprocSINGLEfrac[Simpath_] := Module[{result, data,varString,timeString},
	result = Association[{"Simpath"-> Simpath}];
	data = Import[Simpath, "RawJSON"];
	result = AssociateTo[result, getSINGLEfrac`getFr[data]];
	result = AssociateTo[result, "LastMesh"         ->getGENERALfrac`getLastMesh[data]];
	
	varString="w_at_my_point"; timeString="time_list_W_at_my_point";
	result = AssociateTo[result, "w_A(t)"           ->getGENERALfrac`getVARofT[data,varString,timeString]];
	
	varString="pf_at_my_point_A"; timeString="time_list_pf_at_my_point_A";
	result = AssociateTo[result, "pf_A(t)"           ->getGENERALfrac`getVARofT[data,varString,timeString]];
	
	varString="pf_at_my_point_B"; timeString="time_list_pf_at_my_point_B";
	result = AssociateTo[result, "pf_B(t)"           ->getGENERALfrac`getVARofT[data,varString,timeString]];
	
	
	result
	]

End[]; (* End Private Context *)

EndPackage[];
