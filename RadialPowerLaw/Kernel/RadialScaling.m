(* ::Package:: *)

BeginPackage["RadialScaling`"]

Needs["DescriptionUtilities`"];
Needs["UtilityForScalings`"];


transitionMKScales::usage =	"Gives the time-independent characteristic scales - call sequence is transitionMKScales[inpData,prime] with inpData an Array containg the material properties
and prime a boolean (true if Kp is used [default], false to use KIc)"

timeParameters::usage =	"Gives phi and the dimensionless time along the edges - call sequence is TimeParameters[Eprime,KIc,CL,M,Q,n] output is 
{phi, tmk, tmmt, tmtkt, tkkt, tmb, tkb, tb} where phi is the trajectory parameter, tmk to tkkt are the timescales from the parametric triangle,
tmb and tkb timescales to buoyancy and tb the timescale to use accounting for the dimensionless viscosity in the head \[ScriptCapitalM]_b"

toNumericalScaling::usage = "Scales transition to numerical scaling from vertex scaling ToNumericalScaling[Regime,n,tau,phi]. Only works for the
Vertexes and does not give dimensionless parameters."

vertexScaling::usage = "Vertex scalings Power-law (self-similar scaling) - call sequence is VertexScaling[Regime,inpData,t,prime]
with Regime the string indicating the scaling (M,Mt,K,Kt,Bh,Bt}, inpData an Array containg the material properties, t the time when the solution is calculated
and prime a boolean (true if Kp is used [default], false to use KIc)"


Begin["`Private`"]

vertexScaling[V_?StringQ,inpData_,t_,prime_:True] := Module[{data},

data = inputDataTransformation[inpData,prime];

Switch[V,
	"M",If[(n/.data)!=1.,
			{Lstar->Ep^(1/(6 + 3 n)) Mp^(-(1/(6 + 3 n))) Qo^(1/3) t^((2 (1 + n))/( 3 (2 + n))),
			wstar->Ep^(-(2/(6 + 3 n))) Mp^(2/(6 + 3 n)) Qo^(1/3) t^((2 - n)/(6 + 3 n)),
			pstar->Ep^((1 + n)/(2 + n)) Mp^(1/(2 + n)) t^(-(n/(2 + n))),
			qstar->Ep^(-(1/(6 + 3 n))) Mp^(1/(6 + 3 n)) Qo^( 2/3) t^(-((2 (1 + n))/(3 (2 + n)))),
			\[ScriptCapitalK] -> (Kp t^(1/9))/(Ep^(13/18) Mp^(5/18) Qo^(1/6)),
			\[ScriptCapitalC] -> (Cp Ep^(2/9) t^(7/18))/(Mp^(2/9) Qo^(1/3))
			} /.data//N,
			{Lstar->Ep^(1/(6 + 3 n)) Mp^(-(1/(6 + 3 n))) Qo^(1/3) t^((2 (1 + n))/( 3 (2 + n))),
			wstar->Ep^(-(2/(6 + 3 n))) Mp^(2/(6 + 3 n)) Qo^(1/3) t^((2 - n)/(6 + 3 n)),
			pstar->Ep^((1 + n)/(2 + n)) Mp^(1/(2 + n)) t^(-(n/(2 + n))),
			qstar->Ep^(-(1/(6 + 3 n))) Mp^(1/(6 + 3 n)) Qo^( 2/3) t^(-((2 (1 + n))/(3 (2 + n)))),
			\[ScriptCapitalK] -> (Kp t^(1/9))/(Ep^(13/18) Mp^(5/18) Qo^(1/6)),
			\[ScriptCapitalC] -> (Cp Ep^(2/9) t^(7/18))/(Mp^(2/9) Qo^(1/3)),
			\[ScriptCapitalB] -> (Qo^(1/3) t^(7/9) \[CapitalDelta]\[Gamma])/(Ep^(5/9) Mp^(4/9))
			} /.data//N],
			"K",If[(n/.data)!=1.,
		{Lstar->(Ep^(2/5) Qo^(2/5) t^(2/5))/Kp^(2/5),
		wstar->(Kp^(4/5) Qo^(1/5) t^(1/5))/Ep^(4/5),
		pstar->Kp^(6/5)/(Ep^(1/5) Qo^(1/5) t^(1/5)),
		qstar->(Kp^(2/5) Qo^(3/5))/(Ep^(2/5) t^(2/5)),
			\[ScriptCapitalM] -> (Ep^(13/5) Mp Qo^(3/5))/(Kp^(18/5) t^(2/5)),
			\[ScriptCapitalC] -> (Cp Ep^(4/5) t^(3/10))/(Kp^(4/5) Qo^(1/5))
			} /.data//N,
			{Lstar->(Ep^(2/5) Qo^(2/5) t^(2/5))/Kp^(2/5),
		wstar->(Kp^(4/5) Qo^(1/5) t^(1/5))/Ep^(4/5),
		pstar->Kp^(6/5)/(Ep^(1/5) Qo^(1/5) t^(1/5)),
		qstar->(Kp^(2/5) Qo^(3/5))/(Ep^(2/5) t^(2/5)),
			\[ScriptCapitalM] -> (Ep^(13/5) Mp Qo^(3/5))/(Kp^(18/5) t^(2/5)),
			\[ScriptCapitalC] -> (Cp Ep^(4/5) t^(3/10))/(Kp^(4/5) Qo^(1/5)),
			\[ScriptCapitalB] -> (Ep^(3/5) Qo^(3/5) t^(3/5) \[CapitalDelta]\[Gamma])/Kp^(8/5)
			} /.data//N
			],
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
		,
		"Bh",{Lstarv->Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
		Lstarh-> Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
		wstar->Kp^(4/3)/(Ep \[CapitalDelta]\[Gamma]^(1/3)),
		pstar-> Kp^(2/3) \[CapitalDelta]\[Gamma]^(1/3),
		qstar-> Qo,
			\[ScriptCapitalM] -> (Ep^3 Qo \[CapitalDelta]\[Gamma]^(2/3) Mp)/Kp^(14/3),
			Vh -> Kp^(8/3)/(Ep \[CapitalDelta]\[Gamma]^(5/3))
			} /.data//N,
		"Bt",{Lstarv->(Qo^(2/3) t \[CapitalDelta]\[Gamma]^(7/9))/(Kp^(4/9) Mp^(1/3)),
		Lstarh-> Kp^(2/3)/\[CapitalDelta]\[Gamma]^(2/3),
		wstar->(Qo^(1/3) Mp^(1/3))/(Kp^(2/9) \[CapitalDelta]\[Gamma]^(1/9)),
		pstar-> (Qo^(2/3) t \[CapitalDelta]\[Gamma]^(16/9))/(Kp^(4/9) Mp^(1/3)),
		qstar-> Qo,
			\[ScriptCapitalK] -> Kp^(14/9)/(Ep Qo^(1/3) \[CapitalDelta]\[Gamma]^(2/9) Mp^(1/3)),
			Vt -> Qo t-Kp^(8/3)/(Ep \[CapitalDelta]\[Gamma]^(5/3))
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
   		{data,phi,tmk,tmmt,tmtkt,tkkt,tmb,tkb,tb,Mb,preFac},
   		
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
		   tmb = If[(\[CapitalDelta]\[Gamma]/.data) == 0, Infinity, (Ep^(5/7) Mp^(4/7))/(Qo^(3/7) \[CapitalDelta]\[Gamma]^(9/7))/.data];
		   tkb = If[(\[CapitalDelta]\[Gamma]/.data) == 0, Infinity, Kp^(8/3)/(Ep Qo \[CapitalDelta]\[Gamma]^(5/3))/.data];
		   
			Mb =\[ScriptCapitalM]/.vertexScaling["Bh",inpData,0.1,False];
			preFac = If[Mb<10^-3,5 10^-3 Mb^(-4/7),
			If[Mb>500,0.5,1/(2 10^(2/7))+(Mb-10^-3)*(0.5-1/(2 10^(2/7)))/(500-10^-3)]];
			tb = preFac*tmb;
			
   		{phi, tmk, tmmt, tmtkt, tkkt, tmb, tkb, tb}
   
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
