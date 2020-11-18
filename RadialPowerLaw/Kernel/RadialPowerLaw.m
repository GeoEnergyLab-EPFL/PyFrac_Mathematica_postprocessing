(* ::Package:: *)

BeginPackage["RadialPowerLaw`"]

Get[FileNameJoin[{"RadialPowerLaw","VertexSolutions.m"}]];
Get[FileNameJoin[{"RadialPowerLaw","NearVertexSolutions.m"}]];
Get[FileNameJoin[{"RadialPowerLaw","RadialScaling.m"}]];

Needs["DescriptionUtilities`"];

EarlySolution::usage = "Gives the early time solution for a radial HF driven by a power-law fluid"
MiddleSolution::usage = "Gives the middle time solution for a radial HF driven by a power-law fluid"
LaterSolution::usage = "Gives the large time solution for a radial HF driven by a power-law fluid"
CriticalOpening::usage = "Gives the Adimensionl critical opening"

Begin["`Private`"]
 
Get[FileNameJoin[{"RadialPowerLaw","InsidePowerLaw.txt"}]];
Get[FileNameJoin[{"RadialPowerLaw","MKEdgePowerLaw.txt"}]];


(* Early time MK EDGE *)
EarlySolution[n_]:= 
	Module[
			{
				times,Ns=40,
				ScalarSolution
			},
			
			Off[InterpolatingFunction::dmval];
			
			times = Exp[Table[Log[10^-8] + (i - 1)/Ns (Log[10^5] - Log[10^-8]), {i, 1, Ns + 1}] // N ];
			
			ScalarSolution =
			{
				Time -> times,
				FractureLength -> Transpose[{times, length["MK"][ToString[n]][times]}],
				FractureEfficiency -> Transpose[{times, efficiency["MK"][ToString[n]][times]}],
				InletPressure -> Transpose[{times, Pi1["MK"][ToString[n]][times]}],
				InletOpening -> Transpose[{times, Omega1["MK"][ToString[n]][times]}]
			};
			

	 		{
	 			(# -> (# /. ScalarSolution) ) 
	 			& /@ {Time, FractureLength,FractureEfficiency, InletOpening, InletPressure}
	 		}

		];

(* MEDIUM TIME simulation *)
MiddleSolution[n_,phi_]:=
	Module[
		{
			tmin,tmax,times,Ns=40,ro,rp,
			ScalarSolution,ProfileSolution
		},
		
		Off[InterpolatingFunction::dmval];
		
		ro = Table[(2*i - 1.)/(2.*30), {i, 1, 30}];
		rp = Table[ (2*i - 1.)/(2.*30), {i, 1, 29}];
		
		tmax = FinalTime["Inside"][ToString[n]];
		tmin = InitialTime["Inside"][ToString[n]];
		
		times = Exp[Table[Log[tmin] + (i - 1)/Ns (Log[tmax] - Log[tmin]), {i, 1, Ns + 1}] // N ];
		
		ScalarSolution = 
	     {
	     	Time -> times, 
   			FractureLength -> Transpose[{times,length["Inside"][ToString[n]][phi, times]}], 
   			FractureEfficiency -> Transpose[{times,efficiency["Inside"][ToString[n]][phi, times]}], 
   			InletPressure -> Transpose[{times,Pi1["Inside"][ToString[n]][phi, times]}], 
   			InletOpening -> Transpose[{times, Omega1["Inside"][ToString[n]][phi, times]}] 
   		};
   		
   		ProfileSolution = 
   		{
   			Time -> #, 
     		Opening -> (Transpose[{ro, Omegax["Inside"][ToString[n]][phi, #, ro] // Re}]  /. {Indeterminate -> 0}), 
     		NetPressure -> Transpose[{rp, Pix["Inside"][ToString[n]][phi, #, rp] // Re }]
     	} & /@ times;
		
		{
 			(# -> ( # /. ScalarSolution) ) & /@ {Time, FractureLength, FractureEfficiency,InletOpening, InletPressure},
  			ProfileSolution
  		}
		
	];
	
LaterSolution[n_,phi_]:=
	Module[
		{
			times,Ns=40,
			KTtv,KTtoN,MTtv,MTtoN,x,ro,rp,
			ScalarSolution
		},
		
		times = Exp[Table[Log[10^-3] + (i - 1)/Ns (Log[10^8] - Log[10^-3]), {i, 1, Ns + 1}] // N ];
		
		ro = Table[(2*i - 1.)/(2.*30), {i, 1, 30}];
		rp = Table[ (2*i - 1.)/(2.*30), {i, 1, 29}];
		
		If[
			n > 0.5,
			{
   				(* Large time  Kt vertex *)
  				KTtv = KTILDEVERTEXSOLUTION[x];
  				KTtoN = ToNumericalScaling["KT", n, times, phi],
  
  				ScalarSolution = 
  					{
  						Time -> times,
    					FractureLength -> Transpose[{times, (Lstar /. KTtoN) (FractureLength /. KTtv)}],
    					FractureEfficiency -> Transpose[{times, Table[10^(-10),{Length[times]}]}],
    					InletPressure ->  Transpose[{times, (pstar /. KTtoN) ( (NetPressure /. KTtv) /. {x -> 0.05}) }], 	(* correct x !!!!*)
    					InletOpening -> Transpose[{times, (wstar /. KTtoN) ( (Opening /. KTtv) /. {x -> 0.05}) }] 			(* correct x !!!!*)
    				};
  
			},
  			{  
		    	(* Large time  Mt vertex *)
		   		MTtv = MTILDEVERTEXSOLUTION[x, n ];
		   		MTtoN = ToNumericalScaling["MT", n, times, phi];
		    
		   		ScalarSolution = 
		   			{
		   				Time -> times,
		     			FractureLength -> Transpose[{times, (Lstar /. MTtoN) (FractureLength /. MTtv)}],
		     			FractureEfficiency -> Transpose[{times, Table[10^(-10),{Length[times]}]}],
		     			InletPressure ->  Transpose[{times, (pstar /. MTtoN) ( (NetPressure /. MTtv) /. {x -> 0.05}) }],
		     			InletOpening -> Transpose[{times, (wstar /. MTtoN) ( (Opening /. MTtv) /. {x -> 0.05}) }] 
		     		};
		   
		       }
     
	     ];
	     
	     {
  			(# -> (# /. ScalarSolution)) 
  			& /@ {Time, FractureLength,FractureEfficiency, InletOpening, InletPressure}
  		 }	
		
	];

CriticalOpening[ASolution_,OmegaCritic_,t_]:=
	Module[
		{
			opgdata,opg,xl,time,k
		},
		
		Off[FindRoot::reged];
		
		k=1;
		time=Time /. ASolution[[1]];
		While[time[[k]]<t,null;k++];
		
		xl = ConstantArray[ 0., k ];
		
		Do[
  			{
   				opgdata = Opening /. ASolution[[2]][[i]];
   
   				If[
    				opgdata[[1, 2]] > OmegaCritic,
    				{
     					opg = Interpolation[opgdata, InterpolationOrder -> 1];
     					xl[[i]] = x /. FindRoot[opg[x] == OmegaCritic, {x, 1., 0., 1.}]
     				},
    				xl[[i]] = 0.
    			]
   			},
  			{i, 1, k}
  		];
  		
  		{xl,k}
	];
	
End[]

EndPackage[]
