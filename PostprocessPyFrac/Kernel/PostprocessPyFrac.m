(* ::Package:: *)

(* --------------------------------------------- 

Wolfram Language Package created by Carlo Peruzzo on Tuesday 26th October 2020
EPFL ENAC IIC GEL
GC B1 385 (Batiment GC)
Station 18, 1015 Lausanne, Switzerland
Phone: +41 (0)21 693 76 47
carlo.peruzzo@epfl.ch
Geo-energy Laboratory EPFL 

Contributors: Andreas M\[ODoubleDot]ri (GEL, EPFL)

--------------------------------------------- *)
BeginPackage["PostprocessPyFrac`"];

Needs[ "rectangularMesh`"]; (* when loaded via needs the context is not exposed *)
Needs[ "getDOUBLEfrac`"];
Needs[ "getSINGLEfrac`"];
Needs[ "getGENERALfrac`"];
Needs[ "DescriptionUtilities`"];

postProcessPyfrac::usage = "postProcessPyfrac[ Simpath,Case] it returns an association";
meshPlan::usage =  "meshPlan[ Lx, Nx, Ly, Ny] It computes: {hx,hy,Nelts,vertexes,conn,centers,Neighbors}";
plotFractureFP::usage = "plotFractureFP[meshData_,prop_,opt_,elts_:True,dimless_:True,color_:White]";

(*meshPlan::usage =  "meshPlan[ Lx, Nx, Ly, Ny] It computes: conn, centers, Nelts, ElemNei";
getNeighbors::usage = "getNeighbors[ element, Nx, Ny] It returns {left,right,bottom,up} of element.";*)

Begin["`Private`"];
postProcessPyfrac[Simpath_, Case_] := Which[Case == "double_fracture",postprocDOUBLEfrac[Simpath],
											Case == "single_fracture",postprocSINGLEfrac[Simpath],
											Case == "single_fracture_oldFR",postprocSINGLEfracOF[Simpath],
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
	

postprocSINGLEfracOF[Simpath_] := Module[{result, data,varString,timeString},
	result = Association[{"Simpath"-> Simpath}];
	data = Import[Simpath, "RawJSON"];
	result = AssociateTo[result, getSINGLEfrac`getFrOF[data]];
	result = AssociateTo[result, "LastMesh"         ->getGENERALfrac`getLastMesh[data]];
		
	result
	]

plotFractureFP[meshData_,prop_,opt_,elts_:{},bck_:True,dimless_:True,color_:{Transparent,Black}] := Module[
{elPlt, fps, myMesh, bckground, iter, ind, vertices, indLast, colorPlot},
	If[Length[elts]==0,
	elPlt = Range[1,Length[meshData["time"]]];,
	elPlt = elts;];

	fps = ConstantArray[0,{Length[elPlt],1}];
	
	If[Length[color[[2]]]==1,
	Print["why?"];
	colorPlot = ConstantArray[color[[2]],{Length[elPlt],1}];,
	If[Length[color[[2]]] != Length[elPlt],
	colorPlot = ConstantArray[color[[2,1]],{Length[elPlt],1}];,
	colorPlot = color[[2]];];
	];

	For[iter = 1, iter<=Length[elPlt],iter++,
		ind = elPlt[[iter]];
		myMesh = rectangularMesh`meshPlan[{},{},{},{},meshData["mesh_info"][[ind]]];
		If[dimless,
			vertices=myMesh["vertexes"]/((KIc/\[CapitalDelta]\[Gamma])^(2/3)/.prop);
			fps[[iter]] = ListPlot[{meshData["fp_x"][[ind]],meshData["fp_y"][[ind]]}\[Transpose]/((KIc/\[CapitalDelta]\[Gamma])^(2/3)/.prop),Joined->True,
			PlotStyle->{Thickness[.02],colorPlot[[iter]]}];
			,
			vertices=myMesh["vertexes"];
			fps[[iter]] = ListPlot[{meshData["fp_x"][[ind]],meshData["fp_y"][[ind]]}\[Transpose],Joined->True,
			PlotStyle->{Thickness[.02],colorPlot[[iter]]}];
		];
		

	];
	If[bck,
	indLast = Max[{Length[meshData["time"]]+Max[Select[elPlt,#<0&]],Max[Select[elPlt,#>0&]]}];
	myMesh = rectangularMesh`meshPlan[{},{},{},{},meshData["mesh_info"][[indLast]]];
	bckground= Graphics[{EdgeForm[{Thin,Gray}],color[[1]],GraphicsComplex[vertices,
				Polygon[myMesh["conn"][[1;;]]]]}];
	Show[bckground,fps,opt]	
	,
	Show[fps,opt]	
	]
	
]


End[]; (* End Private Context *)

EndPackage[];
