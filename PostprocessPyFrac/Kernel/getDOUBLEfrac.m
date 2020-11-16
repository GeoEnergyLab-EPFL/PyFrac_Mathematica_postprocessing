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

BeginPackage["getDOUBLEfrac`"]

getFr::usage =  "to be documented";

getNumDeltatilde::usage =  "to be documented";

getNumRinFrLeft::usage =  "to be documented";

getNumRinFrRight::usage =  "to be documented";

getNumRoutFrLeft::usage =  "to be documented";

getNumRoutFrRight::usage =  "to be documented";

Begin["`Private`"] 

(* ----------- Front functions ------------ *)
(*this Module gets a list "Frpoints" having two fronts inside. the function separates the two fronts *)
convertTOxyTWOFronts[Frpoints_] := Module[{i,j, 
										   x0 = {Frpoints[[1]][[1]], Frpoints[[1]][[2]]}, (* memorize the first point of the list *)
										   indexes = {}, 
										   Check = True, 
										   NofPoints = Length[Frpoints], 
										   FractureA,
										   FractureB,
										   return},
										   
   (* loop over the points and execute the first test that yelds true*)
   FractureA = Table[ 
   	Which[ Check && x0 == {Frpoints[[i]][[3]], Frpoints[[i]][[4]]},
   	       Check = False; {Frpoints[[i]][[1]], Frpoints[[i]][[2]]},
   	       (*(AA) if true it means that you closed the loop of fracture A by finding a recurrence of the first point of the list *)
           Check && ! x0 == {Frpoints[[i]][[3]],Frpoints[[i]][[4]]},
           {Frpoints[[i]][[1]], Frpoints[[i]][[2]]},
           (*(BB) if true it means that you did not yed closed the loop of fracture A *)
      		! Check && ! x0 == {Frpoints[[i]][[3]], Frpoints[[i]][[4]]}, 
      		indexes = Join[indexes, {i}]; Nothing
      		(*(CC) if true it means that you already have found (AA) True *)
                    ], {i, 1, NofPoints}];
   FractureA = Join[FractureA, {x0}];(*returning the first point and closing the fracture *)
   
   FractureB = Table[j = indexes[[i]]; {Frpoints[[j]][[1]], Frpoints[[j]][[2]]}, {i, 1, Length[indexes]}];
   FractureB = Join[FractureB, {FractureB[[1]]}]; (*returning the first point and closing the fracture *)
   return = {FractureA,FractureB};
   return
   ]
   
(* this function converts a front from pyfrac removing the repeated coordinates*)
convertTOxySingleFront[Frpoints_] := Table[{Frpoints[[i]][[1]], Frpoints[[i]][[2]]}, {i, 1, Length[Frpoints]}];
    
separateFronts[noffLst_,frLst_,nOfTstps_]:= Module[{FracturesLeft = {},  (* list of left fronts *)
													FracturesRight = {}, (* list of right fronts *)
													SingleFracture,
													FractureA,
										  			FractureB,
													front,i,return},
	(* take out the single front *)
	SingleFracture = Table[ If[noffLst[[i]] == 1, 
	                           front = convertTOxySingleFront[frLst[[i]]];
	                           Join[front , {front [[1]]}]
	                           ,
	                           Nothing], {i, 1, nOfTstps}]; 
    
    (* decide where to append each of the two fractures you got *)
    Do[ If[noffLst[[i]] == 2, 
    	{FractureA,FractureB} = convertTOxyTWOFronts[frLst[[i]]];
		If[Mean[FractureA][[1]] < 0 && Mean[FractureB][[1]] > 0,
   			FracturesLeft  = Join[FracturesLeft,  {FractureA}];
   			FracturesRight = Join[FracturesRight, {FractureB}];,
   			FracturesLeft  = Join[FracturesLeft,  {FractureB}];
   			FracturesRight = Join[FracturesRight, {FractureA}];],Nothing]
   			, {i, 1, nOfTstps}];
    return = <|"FracturesLeft"->FracturesLeft, "FracturesRight"->FracturesRight, "SingleFractures"->SingleFracture|>;
    return
	](* end module *)
getFr[data_]:= Module[{noffLst,frLst,nOfTstps,return},
    noffLst  = data["Number_of_fronts"][[1]];(* number of fracture fronts *)
    frLst    = data["Fr_list"][[1]]; (* fracture lists *)
    nOfTstps = Length[noffLst];  (* number of lists *)
    return   = separateFronts[noffLst,frLst,nOfTstps];
	return
	](* end module *)
getNumRinFrLeft[ data_, prevResult_]:= Module[{},{}
	];(* end module *)

getNumRinFrRight[ data_, prevResult_]:= Module[{},{}
	];(* end module *)

getNumRoutFrLeft[ data_, prevResult_]:= Module[{},{}
	];(* end module *)

getNumRoutFrRight[ data_, prevResult_]:= Module[{},{}
	];(* end module *)

getNumDeltatilde[ data_, prevResult_]:= Module[{},{}
	];(* end module *)
	

     
     
End[] (* End Private Context *)

EndPackage[]
