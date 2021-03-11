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

BeginPackage["getSINGLEfrac`"]

getFr::usage =  "to be documented";
getFrOF::usage =  "to be documented";

Begin["`Private`"] 

(* ----------- Front functions ------------ *)
(* this function converts a front from pyfrac removing the repeated coordinates*)

convertTOxySingleFront[Frpoints_] := Table[{Frpoints[[i]][[1]], Frpoints[[i]][[2]]}, {i, 1, Length[Frpoints]}]
    
separateFronts[frLst_,nOfTstps_]:= Module[{front,i,return},
	
	return = Table[front = convertTOxySingleFront[frLst[[i]]]; Join[front, {front[[1]]}], {i, 1, nOfTstps}];
    return](* end module *)

getFr[data_]:= Module[{ frLst,
						nOfTstps,
						singlefronts,
						xinj,
						yinj,
						times,
						xmax,xmin,ymax,ymin,lst,lst2,j,i,
						rRight,rLeft,rUp,rBottom,
						rRightVsT,
						rLeftVsT,
						rUpVsT,
						rBottomVsT,
						return},
    frLst    = data["Fr_list"][[1]]; (* fracture lists *)
    nOfTstps = Length[frLst];  (* number of lists *)
    singlefronts   = separateFronts[frLst,nOfTstps];
    xinj = data["simul_info"]["sources_coordinates_lastFR"][[1]][[1]];
	yinj = data["simul_info"]["sources_coordinates_lastFR"][[1]][[2]];
	times = data["time_srs_of_Fr_list"][[1]];
	
	xmax = Table[lst = singlefronts[[i]]; 
		   lst2=Table[lst[[j,1]],{j,1,Length[lst]}];
		   Max[lst2],{i,1,Length[singlefronts]}];(* x max *)
		   
	xmin = Table[lst = singlefronts[[i]]; 
		   lst2=Table[lst[[j,1]],{j,1,Length[lst]}];
		   Min[lst2],{i,1,Length[singlefronts]}];(* x max *)
		   	   
	ymax = Table[lst = singlefronts[[i]]; 
		   lst2=Table[lst[[j,2]],{j,1,Length[lst]}];
		   Max[lst2],{i,1,Length[singlefronts]}];(* x max *)
		   
	ymin = Table[lst = singlefronts[[i]]; 
		   lst2=Table[lst[[j,2]],{j,1,Length[lst]}];
		   Min[lst2],{i,1,Length[singlefronts]}];(* x max *)
	
	rRight = Abs[xmax - xinj];
	rLeft = Abs[xmin - xinj];
	rUp = Abs[ymax - yinj];
	rBottom = Abs[ymin - yinj];
	
	rRightVsT = Table[{times[[j]], rRight[[j]] },{j,1,Length[times]}];
	rLeftVsT = Table[{times[[j]], rLeft[[j]] },{j,1,Length[times]}];
	rUpVsT = Table[{times[[j]], rUp[[j]] },{j,1,Length[times]}];
	rBottomVsT = Table[{times[[j]], rBottom[[j]] },{j,1,Length[times]}];		
    return = <|"SingleFractures"->singlefronts,"times"->times,"xinj"->xinj,"yinj"->yinj,"rRightVsT"->rRightVsT,"rLeftVsT"->rLeftVsT,"rUpVsT"->rUpVsT,"rBottomVsT"->rBottomVsT|>;
	return](* end module *)   
    
getFrOF[data_]:= Module[{ frLst,
						nOfTstps,
						times,
						return},
    frLst    = data["Fr_list"][[1]]; (* fracture lists *)
	times = data["time_srs_of_Fr_list"][[1]];
			
    return = <|"SingleFractures"->frLst,"times"->times|>;
	return](* end module *)   
	
	
	
End[] (* End Private Context *)

EndPackage[]
