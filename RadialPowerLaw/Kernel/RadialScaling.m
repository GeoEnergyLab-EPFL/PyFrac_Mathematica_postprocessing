(* ::Package:: *)

BeginPackage["RadialPowerLaw`RadialScaling`"]

Needs["DescriptionUtilities`"];

DimensionalScales::usage =	"Gives the time-independent characteristic scales - call sequence is DimensionalScales[Eprime,KIc,CL,M,Q,n]"

TimeParameters::usage =	"Gives phi and the dimensionless time along the edges - call sequence is TimeParameters[Eprime,KIc,CL,M,Q,n] "

ToNumericalScaling::usage = "Scales transition to numerical scaling from vertex scaling ToNumericalScaling[V_?StringQ,n_,tau_,phi_:1]"

VertexScaling::usage = "Vertex scalings Power-law (self-similar scaling) - call sequence is VertexScaling[V_?StringQ,Ep_,Kp_,Mp_,Qo_,n_,t_] 
with Kp=(32/pi)^0.5 , Mp= alpha M with alpha=2^(1 + n) n^-n (1 + 2 n)^n" 

Begin["`Private`"]

VertexScaling[V_?StringQ,Ep_,Kp_,Mp_,Qo_,n_,t_] =
Switch[V,
	"M",
			{Lstar->Ep^(1/(6 + 3 n)) Mp^(-(1/(6 + 3 n))) Qo^(1/3) t^((2 (1 + n))/( 3 (2 + n))),
			wstar->Ep^(-(2/(6 + 3 n))) Mp^(2/(6 + 3 n)) Qo^(1/3) t^((2 - n)/(6 + 3 n)),
			pstar->Ep^((1 + n)/(2 + n)) Mp^(1/(2 + n)) t^(-(n/(2 + n))),
			qstar->Ep^(-(1/(6 + 3 n))) Mp^(1/(6 + 3 n)) Qo^( 2/3) t^(-((2 (1 + n))/(3 (2 + n))))
			},
			"K",
		{Lstar->(Ep^(2/5) Qo^(2/5) t^(2/5))/Kp^(2/5),
		wstar->(Kp^(4/5) Qo^(1/5) t^(1/5))/Ep^(4/5),
		pstar->Kp^(6/5)/(Ep^(1/5) Qo^(1/5) t^(1/5)),
		qstar->(Kp^(2/5) Qo^(3/5))/(Ep^(2/5) t^(2/5))
		},
		"Mt",{Lstar->(Sqrt[Qo] t^(1/4))/Sqrt[Cp],
		wstar->Cp^((-2 + n)/(4 (1 + n))) Ep^(-(1/(2 + 2 n))) Mp^(1/(2 + 2 n)) Qo^(( 2 + n)/(4 + 4 n)) t^((2 - n)/(8 + 8 n)),
 		pstar->Cp^((3 n)/(4 + 4 n)) Ep^(1 - 1/(2 + 2 n)) Mp^(1/( 2 + 2 n)) Qo^(-(n/(4 + 4 n))) t^(-((3 n)/(8 + 8 n))),
 		qstar->(Sqrt[Cp] Sqrt[Qo])/t^(1/4)
		},
		"Kt",{Lstar->(Sqrt[Qo] t^(1/4))/Sqrt[Cp],
		wstar->(Kp Qo^(1/4) t^(1/8))/(Cp^(1/4) Ep),
		pstar->(Cp^(1/4) Kp)/(Qo^(1/4) t^(1/8)),
		qstar->(Sqrt[Cp] Sqrt[Qo])/t^(1/4)
		}
];



DimensionalScales[Eprime_, KIc_, CL_, M_, Q_, n_] = 
Module[
   		{
   			Kprime, Cprime, Mprime, alpha    		
    	},
   
   		alpha = (2^(n + 1) (2 n + 1)^n)/n^n;
   		Kprime = 4 (2/\[Pi])^(1/2) KIc;
   		Cprime = 2 CL;
   		Mprime = alpha M;
         
   		{
	   		Lstar-> Q^(2/5 + (2 (2 + n))/(5 (-2 + 4 n)))* 
	    			Eprime^(2/5 + (2 (7 + 6 n))/(5 (-2 + 4 n)))*
	      			Kprime^(-(2/5) - (12 (2 + n))/(5 (-2 + 4 n))) Mprime^(2/(-2 + 4 n)),
	      	wstar-> Q^(1/5 + (2 + n)/(5 (-2 + 4 n)))*
	    			Eprime^(-(4/5) + (7 + 6 n)/(5 (-2 + 4 n)))*
	    			Kprime^(4/5 - (6 (2 + n))/(5 (-2 + 4 n)))*
	      			Mprime^(1/(-2 + 4 n)),
	      	pstar-> Q^(-(1/5) - (2 + n)/(5 (-2 + 4 n)))* 
	    			Eprime^(-(1/5) - (7 + 6 n)/(5 (-2 + 4 n)))*
	    			Kprime^(6/5 + (6 (2 + n))/(5 (-2 + 4 n)))*
	      			Mprime^(-(1/(-2 + 4 n))),
	      	qstar -> Q^(3/5 - (2 (2 + n))/(5 (-2 + 4 n)))*
	      			 Eprime^(-(2/5) - (2 (7 + 6 n))/(5 (-2 + 4 n)))*
	      			 Kprime^(2/5 + (12 (2 + n))/(5 (-2 + 4 n)))*
	      			 Mprime^(-(2/(-2 + 4 n))),
	      	tstar -> ((Mprime^5 Eprime^(6 n + 7) Q^(n + 2))/Kprime^(6 (n + 2)))^(1/(4 n - 2))
      	}
];


TimeParameters[Eprime_, KIc_, CL_, M_, Q_, n_] := 
Module[
   		{
   			Kprime, Cprime, Mprime, alpha,
    		phi, tmk, tmmt, tmtkt, tkkt
    	},
   		
   		alpha = (2^(n + 1) (2 n + 1)^n)/n^n;
   		Kprime = 4 (2/\[Pi])^(1/2) KIc;
   		Cprime = 2 CL; 
   		Mprime = alpha M;
   
   		phi =
   		 If[
   			Cprime==0,
   			0,
   			((Eprime^(10 n + 1) Mprime^3 Cprime^(4 (2 n - 1)) Q^(2 - n))/Kprime^(2 (2 + 5 n)))^(1/(2 n - 1))
   		];
   
      	(* Times along the edges of the rectangular space MKK~M~ *)
   		
   		tmk = ((Mprime^5 Eprime^(6 n + 7) Q^(n + 2))/Kprime^(6 (n + 2)))^(1/(4 n - 2));
   		tmmt = If[Cprime == 0, Infinity, ((Mprime^4 Q^(2 (n + 2)))/(Eprime^4 Cprime^(6 (n + 2))))^(1/(5 n + 2))];
   		tmtkt = If[Cprime==0,0,((Mprime^4 Eprime^(4 (2 n + 1)) Cprime^(2 (2 n - 1)) Q^2)/Kprime^(8 (n + 1)))^(1/(2 n - 1))];	
		tkkt = If[Cprime == 0, Infinity, ((Kprime^8 Q^2)/(Eprime^8 Cprime^10))^(1/3)];

   		{phi, tmk, tmmt, tmtkt, tkkt}
   
];


ToNumericalScaling[V_?StringQ,n_,tau_,phi_:1] =
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
		"MT",
			{
				Lstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ (1/4),
				wstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ ((2-n)/(8(1+n))),
				pstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ -(3n/(8(1+n))),
				qstar-> (tau/(phi^(3(n+2)/(2(5n+2))))) ^ -(1/4)
			},
		"KT",
			{
				Lstar-> (tau/(phi^(1/2))) ^ (1/4),
				wstar-> (tau/(phi^(1/2))) ^ (1/8),
				pstar-> (tau/(phi^(1/2))) ^ -(1/8),
				qstar-> (tau/(phi^(1/2))) ^ -(1/4) 
			}
];

End[]

EndPackage[]
