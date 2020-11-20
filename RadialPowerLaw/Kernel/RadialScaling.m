(* ::Package:: *)

BeginPackage["RadialScaling`"]

Needs["DescriptionUtilities`"];
Needs["UtilityForScalings`"];


transitionMKScales::usage =	"Gives the time-independent characteristic scales - call sequence is DimensionalScales[Eprime,KIc,CL,M,Q,n]"

timeParameters::usage =	"Gives phi and the dimensionless time along the edges - call sequence is TimeParameters[Eprime,KIc,CL,M,Q,n]"

toNumericalScaling::usage = "Scales transition to numerical scaling from vertex scaling ToNumericalScaling[V_?StringQ,n_,tau_,phi_1]"

vertexScaling::usage = "Vertex scalings Power-law (self-similar scaling) - call sequence is VertexScaling[V_?StringQ,Ep_,Kp_,Mp_,Qo_,n_,t_] 
with Kp=(32/pi)^0.5 , Mp= alpha M with alpha=2^(1 + n) n^-n (1 + 2 n)^n"


Begin["`Private`"]

vertexScaling[V_?StringQ,inpData_,t_,prime_:True] := Module[{data},

data = inputDataTransformation[inpData,prime];

Switch[V,
	"M",
			{Lstar->Ep^(1/(6 + 3 n)) Mp^(-(1/(6 + 3 n))) Qo^(1/3) t^((2 (1 + n))/( 3 (2 + n))),
			wstar->Ep^(-(2/(6 + 3 n))) Mp^(2/(6 + 3 n)) Qo^(1/3) t^((2 - n)/(6 + 3 n)),
			pstar->Ep^((1 + n)/(2 + n)) Mp^(1/(2 + n)) t^(-(n/(2 + n))),
			qstar->Ep^(-(1/(6 + 3 n))) Mp^(1/(6 + 3 n)) Qo^( 2/3) t^(-((2 (1 + n))/(3 (2 + n)))),
			\[ScriptCapitalK] -> (Kp t^(1/9))/(Ep^(13/18) Mp^(5/18) Qo^(1/6)),
			\[ScriptCapitalC] -> (Cp Ep^(2/9) t^(7/18))/(Mp^(2/9) Qo^(1/3))
			} /.data//N,
			"K",
		{Lstar->(Ep^(2/5) Qo^(2/5) t^(2/5))/Kp^(2/5),
		wstar->(Kp^(4/5) Qo^(1/5) t^(1/5))/Ep^(4/5),
		pstar->Kp^(6/5)/(Ep^(1/5) Qo^(1/5) t^(1/5)),
		qstar->(Kp^(2/5) Qo^(3/5))/(Ep^(2/5) t^(2/5)),
			\[ScriptCapitalM] -> (Ep^(13/5) Mp Qo^(3/5))/(Kp^(18/5) t^(2/5)),
			\[ScriptCapitalC] -> (Cp Ep^(4/5) t^(3/10))/(Kp^(4/5) Qo^(1/5))
			} /.data//N,
		"Mt",{Lstar->(Sqrt[Qo] t^(1/4))/Sqrt[Cp],
		wstar->Cp^((-2 + n)/(4 (1 + n))) Ep^(-(1/(2 + 2 n))) Mp^(1/(2 + 2 n)) Qo^(( 2 + n)/(4 + 4 n)) t^((2 - n)/(8 + 8 n)),
 		pstar->Cp^((3 n)/(4 + 4 n)) Ep^(1 - 1/(2 + 2 n)) Mp^(1/( 2 + 2 n)) Qo^(-(n/(4 + 4 n))) t^(-((3 n)/(8 + 8 n))),
 		qstar->(Sqrt[Cp] Sqrt[Qo])/t^(1/4),
			\[ScriptCapitalK] -> (Kp t^(1/16))/(Cp^(1/8) Ep^(3/4) Mp^(1/4) Qo^(1/8)),
			\[ScriptCapitalV] -> (Cp^(9/8) Ep^(1/4) t^(7/16))/(Mp^(1/4) Qo^(3/8))
			} /.data//N,
		"Kt",{Lstar->(Sqrt[Qo] t^(1/4))/Sqrt[Cp],
		wstar->(Kp Qo^(1/4) t^(1/8))/(Cp^(1/4) Ep),
		pstar->(Cp^(1/4) Kp)/(Qo^(1/4) t^(1/8)),
		qstar->(Sqrt[Cp] Sqrt[Qo])/t^(1/4),
			\[ScriptCapitalM] -> (Sqrt[Cp] Ep^3 Mp Sqrt[Qo])/(Kp^4 t^(1/4)),
			\[ScriptCapitalV] -> (Cp^(5/4) Ep t^(3/8))/(Kp Qo^(1/4))
			} /.data//N
]
];


transitionMKScales[inpData_,prime_:True] := 
Module[{data},
   
	data = inputDataTransformation[inpData,prime];
         
   		{
	   		Lstar-> Qo^(2/5 + (2 (2 + n))/(5 (-2 + 4 n)))* 
	    			Ep^(2/5 + (2 (7 + 6 n))/(5 (-2 + 4 n)))*
	      			Kp^(-(2/5) - (12 (2 + n))/(5 (-2 + 4 n))) Mp^(2/(-2 + 4 n)),
	      	wstar-> Qo^(1/5 + (2 + n)/(5 (-2 + 4 n)))*
	    			Ep^(-(4/5) + (7 + 6 n)/(5 (-2 + 4 n)))*
	    			Kp^(4/5 - (6 (2 + n))/(5 (-2 + 4 n)))*
	      			Mp^(1/(-2 + 4 n)),
	      	pstar-> Qo^(-(1/5) - (2 + n)/(5 (-2 + 4 n)))* 
	    			Ep^(-(1/5) - (7 + 6 n)/(5 (-2 + 4 n)))*
	    			Kp^(6/5 + (6 (2 + n))/(5 (-2 + 4 n)))*
	      			Mp^(-(1/(-2 + 4 n))),
	      	qstar -> Qo^(3/5 - (2 (2 + n))/(5 (-2 + 4 n)))*
	      			 Ep^(-(2/5) - (2 (7 + 6 n))/(5 (-2 + 4 n)))*
	      			 Kp^(2/5 + (12 (2 + n))/(5 (-2 + 4 n)))*
	      			 Mp^(-(2/(-2 + 4 n))),
	      	tstar -> ((Mp^5 Ep^(6 n + 7) Qo^(n + 2))/Kp^(6 (n + 2)))^(1/(4 n - 2))
      	} /.data//N
];


timeParameters[inpData_,prime_:True] := 
Module[
   		{data,phi,tmk,tmmt,tmtkt,tkkt},
   		data = inputDataTransformation[inpData,prime];
   
   		phi =
   		 If[
   			(Cp/.data)==0,
   			0,
   			((Ep^(10 n + 1) Mp^3 Cp^(4 (2 n - 1)) Qo^(2 - n))/Kp^(2 (2 + 5 n)))^(1/(2 n - 1))/.data
   		];
   
      	(* Times along the edges of the rectangular space MKK~M~ *)
   		
   		tmk = ((Mp^5 Ep^(6 n + 7) Qo^(n + 2))/Kp^(6 (n + 2)))^(1/(4 n - 2))/.data;
   		tmmt = If[(Cp/.data)== 0, Infinity, ((Mp^4 Qo^(2 (n + 2)))/(Ep^4 Cp^(6 (n + 2))))^(1/(5 n + 2))/.data];
   		tmtkt = If[(Cp/.data)==0,0,((Mp^4 Ep^(4 (2 n + 1)) Cp^(2 (2 n - 1)) Qo^2)/Kp^(8 (n + 1)))^(1/(2 n - 1))/.data];	
		tkkt = If[(Cp/.data) == 0, Infinity, ((Kp^8 Qo^2)/(Ep^8 Cp^10))^(1/3)/.data];

   		{phi, tmk, tmmt, tmtkt, tkkt}
   
];


toNumericalScaling[V_?StringQ,n_,tau_,phi_:1] :=
Switch[
		V,
		"M",
			{
				Lstar-> tau ^ (2(n+1)/(3(n+2))),
				wstar-> tau ^ ((2-n)/(3(n+2))),
				pstar-> tau ^ -(n/(n+2)),
				qstar-> tau ^ -(2(n+1)/(3(n+2))) 
			},
		"K",
			{
				Lstar-> tau ^ (2/5),
				wstar-> tau ^ (1/5),
				pstar-> tau ^ -(1/5),
				qstar-> tau ^ -(2/5)
			},
		"Mt",
			{
				Lstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ (1/4),
				wstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ ((2-n)/(8(1+n))),
				pstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ -(3n/(8(1+n))),
				qstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ -(1/4)
			},
		"Kt",
			{
				Lstar-> (tau/(phi^(1/2))) ^ (1/4),
				wstar-> (tau/(phi^(1/2))) ^ (1/8),
				pstar-> (tau/(phi^(1/2))) ^ -(1/8),
				qstar-> (tau/(phi^(1/2))) ^ -(1/4) 
			}
];


End[]

EndPackage[]
