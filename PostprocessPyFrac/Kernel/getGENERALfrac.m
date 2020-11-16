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

BeginPackage["getGENERALfrac`"]
Needs[ "rectangularMesh`"];

(* ----------- point evolution ------------ *)
getVARofT::usage =  "to be documented";
getLastMesh::usage =  "to be documented";

Begin["`Private`"] 

getLastMesh[data_]:= Module[{meshinfoLst,lastmesh,Lx,Ly,Nx,Ny},
	meshinfoLst = data["mesh_info"][[1]];
	lastmesh = meshinfoLst[[-1]];
	{Lx,Ly,Nx,Ny} = lastmesh;
	meshPlan[Lx, Nx, Ly, Ny] ](* end module *)
	
(* ----------- point evolution ------------ *)
getVARofT[data_,varString_,timeString_]:= Module[{times,var,i,result},
	var=data[varString][[1]];
	times=data[timeString][[1]];
	result=Table[{times[[i]], var[[i]]}, {i, 1, Length[times]}];
	result](* end module *)
     

     
End[] (* End Private Context *)

EndPackage[]