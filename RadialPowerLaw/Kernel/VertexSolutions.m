(* ::Package:: *)

BeginPackage["VertexSolutions`"]

Needs["DescriptionUtilities`"];

MVERTEXSOLUTION::usage = "The M vertex solution correspond to the asymptotic case 
of a penny-shaped fracture in an impermeable elastic rock with zero toughness
MVERTEXSOLUTION[rho,n]"

KVERTEXSOLUTION::usage = "The K vertex solution correspond to the problem of a 
fracture driven by an inviscid fluid in an impermeable medium KVERTEXSOLUTION[rho]"

MTILDEVERTEXSOLUTION::usage = "The M~ vertex solution correspond to the asymptotic 
solution of a penny shaped fracture in a permeable medium with finite toughness
MTILDEVERTEXSOLUTION[rho,n]"

KTILDEVERTEXSOLUTION::usage = "The K~ vertex solution correspond to the problem of a 
fracture propagating in the toughness-leak off dominated KTILDEVERTEXSOLUTION[rho]"

Begin["`Private`"]

(* Auxiliary functions for M and M~ vertex*)
G[i_, p_, q_, rho_] :=
    (i!*Gamma[p + i])/Gamma[p + 2*i]*JacobiP[i, p - q, q - 1, 2*rho - 1]
h[i_, p_, q_] :=
    (i!*Gamma[i + q]*Gamma[i + p]*Gamma[i + p - q + 1])/((2 i + p)*Gamma[2 i + p]^2)

(*-------------------------------------M VERTEX SOLUTION-------------------------------------*)

CoefficientsM = (*Coefficients for M vertex solution*)
{{b -> 0.505448, a[1] -> 4.14615, a[2] -> -0.0512996, a[3] -> 0.00190012, a[4] -> 0.00057974, 
  c[1] -> 0.496277, c[2] -> 0.0940619, c[3] -> 0.012608, c[4] -> 0.00350545, c[5] -> 0.00137, 
  w[1] -> 1.49158, w[2] -> 0.719954, w[3] -> 0.597918, w[4] -> 0.508402}, 
 {b -> 0.412098, a[1] -> 2.0444, a[2] -> -0.0403466, a[3] -> 0.00312301, a[4] -> 0.00056157,
  c[1] -> 0.52214, c[2] -> 0.0907525, c[3] -> 0.0103306, c[4] -> 0.00293439, c[5] -> 0.00113799, 
  w[1] -> 1.57264, w[2] -> 0.726771, w[3] -> 0.602933, w[4] -> 0.511255}, 
 {b -> 0.340435, a[1] -> 1.33963, a[2] -> -0.0311094, a[3] -> 0.00391093, a[4] -> 0.00049646,
  c[1] -> 0.544558, c[2] -> 0.0875195, c[3] -> 0.00839435, c[4] -> 0.00247058, c[5] -> 0.00094958, 
  w[1] -> 1.65809, w[2] -> 0.73291, w[3] -> 0.607599, w[4] -> 0.513862}, 
 {b ->  0.282418, a[1] -> 0.985152, a[2] -> -0.0228307, a[3] -> 0.00431844, a[4] -> 0.0004034, 
  c[1] -> 0.565351, c[2] -> 0.0843429, c[3] -> 0.00669772, c[4] -> 0.002092, c[5] -> 0.00079279,
  w[1] -> 1.74876, w[2] -> 0.73844, w[3] -> 0.611944, w[4] -> 0.516275}, 
 {b -> 0.234717, a[1] -> 0.771342, a[2] -> -0.0153711, a[3] -> 0.00440752, a[4] -> 0.0003025, 
  c[1] -> 0.585091, c[2] -> 0.0812624, c[3] -> 0.00520872, c[4] -> 0.00178367, c[5] -> 0.00066089, 
  w[1] -> 1.84562, w[2] -> 0.743409, w[3] -> 0.616002, w[4] -> 0.518497},
 {b -> 0.195237, a[1] -> 0.628173, a[2] -> -0.00869054, a[3] -> 0.00424317, a[4] -> 0.00021074, 
  c[1] -> 0.604109, c[2] -> 0.0783232, c[3] -> 0.00390739, c[4] -> 0.00153293, c[5] -> 0.00054886,
  w[1] -> 1.9499, w[2] -> 0.747871, w[3] -> 0.619804, w[4] -> 0.520557}, 
 {b -> 0.162433, a[1] -> 0.525514, a[2] -> -0.00275709, a[3] -> 0.00388914, a[4] -> 0.00014137, 
  c[1] -> 0.622687, c[2] -> 0.0755662, c[3] -> 0.00277627, c[4] -> 0.00132884, c[5] -> 0.0004527, 
  w[1] -> 2.06317, w[2] -> 0.751892, w[3] -> 0.623393, w[4] -> 0.52249}, 
 {b -> 0.1351, a[1] -> 0.448275, a[2] -> 0.00245555, a[3] -> 0.00340339, a[4] -> 0.00010331, 
  c[1] -> 0.641086, c[2] -> 0.0730238, c[3] -> 0.00179633, c[4] -> 0.00116064, c[5] -> 0.00036837, 
  w[1] -> 2.1873, w[2] -> 0.755445, w[3] -> 0.626727, w[4] -> 0.524247}, 
 {b -> 0.112268, a[1] -> 0.387937, a[2] -> 0.00699536, a[3] -> 0.00284055, a[4] -> 0.00010226, 
  c[1] -> 0.659608, c[2] -> 0.070731, c[3] -> 0.00095222, c[4] -> 0.001021, c[5] -> 0.00029381, 
  w[1] -> 2.32491, w[2] -> 0.758633, w[3] -> 0.629889, w[4] -> 0.525906}, 
 {b -> 0.0931746, a[1] -> 0.338851, a[2] -> 0.0112089, a[3] -> 0.00232964, a[4] -> 0.00016, 
  c[1] -> 0.678415, c[2] -> 0.0686695, c[3] -> 0.00020376, c[4] -> 0.00089093, c[5] -> 0.00022263,
  w[1] -> 2.47937, w[2] -> 0.761466, w[3] -> 0.632875, w[4] -> 0.527456}};
  
gammaM = { 0.730912, 0.728744, 0.726022, 0.722864, 0.719351, 0.715535, 0.711441, 0.70708, 0.702436, 0.6976};  (* from 0.1 to 1. by steps of .1 *)

(* OPENING FOR M VERTEX SOLUTION *)
(* Asymptotic behavior at the Fracture Tip *)
OpeningTipM[rho_,n_,i_?IntegerQ] = (1 - rho)^(2/(2 + n))/Sqrt[h[i - 1, (2 (4 + n))/(2 + n), 2]] G[i - 1, (2 (4 + n))/(2 + n), 2, rho];
 
(* Particular soltuion having the proper singularity at the source *)
OpeningSourceM[rho_,n_] = If[ n == 1,
                              8/\[Pi] (Sqrt[1 - rho^2] - rho ArcCos[rho]),
                              4/(Sqrt[\[Pi]] (1 - n)) Gamma[(3 - n)/2]/ Gamma[(4 - n)/2](Sqrt[1 - rho^2] -
                              (Sqrt[\[Pi]] rho^(2 - n) Gamma[(n - 2)/2])/(2 Gamma[(n - 1)/2]) -
                               Hypergeometric2F1[1/2, (n - 2)/2, n/2, rho^2]/(2 - n))
                          ];

(* M Vertex opening solution *)
OpeningM[rho_,n_] :=
    gammaM[[IntegerPart[10 n]]] (b OpeningSourceM[rho,n] + Sum[c[i] OpeningTipM[rho,n,i],{i,1,5}]) /.(CoefficientsM[[IntegerPart[10 n]]]);

(* PRESSURE FOR M VERTEX SOLUTION *)
(* Asymptotic behavior  at the Fracture Tip *)
PressureTipM[rho_,n_,i_?IntegerQ] :=
    -((1 - rho)^(-(n/(2 + n)))/Sqrt[h[i - 1, 4/(2 + n), 2]]) G[i - 1, 4/(2 + n), 2, rho] + w[i] /.(CoefficientsM[[IntegerPart[10 n]]]);

(* Particular solution having the proper singularity at the source *)
PressureSourceM[rho_,n_] = If[ n == 1,
                               Log[2/rho] - 1,
                               -(rho^(1 - n)/(1 - n)) + Sqrt[\[Pi]]/(2 (1 - n)) Gamma[(3 - n)/2]/Gamma[(4 - n)/2]
                           ];

(* M Vertex pressure solution *)
PressureM[rho_,n_] :=
    (b PressureSourceM[rho,n] + Sum[a[i] PressureTipM[rho,n,i],{i,1,4}]) /. (CoefficientsM[[IntegerPart[10 n]]]);

MVERTEXSOLUTION[rho_,n_] :={ 
								FractureLength	-> gammaM[[IntegerPart[10 n]]],
      							Opening     	-> OpeningM[rho, n], 
      							NetPressure     -> PressureM[rho, n]
      						}

(*-------------------------------------K VERTEX SOLUTION-------------------------------------*)

KVERTEXSOLUTION[rho_] = {    
							FractureLength	-> (3/(\[Pi] Sqrt[2]))^(2/5), 
                            Opening     	-> (3/(8 \[Pi]))^(1/5) (1 - rho^2)^(1/2), 
                            NetPressure     -> \[Pi]/8 (\[Pi]/12)^(1/5),
                            FluidFlux       -> 1/(5 (2^4 3^2 \[Pi]^3)^(1/5)) ((5 + rho^2) Sqrt[1 - rho^2])/rho
                            }

(*-------------------------------------M~VERTEX SOLUTION-------------------------------------*)

CoefficientsMT = (*Coefficients for M~ vertex solution*)
{{d -> 0.527462, e[1] -> 3.06401, e[2] -> -0.0533871, e[3] -> 0.0019528, e[4] -> 0.00056779, 
  f[1] -> 0.546224, f[2] -> 0.100386, f[3] -> 0.0129076, f[4] -> 0.00362441, f[5] -> 0.0014039, 
  w[1] -> 1.52865, w[2] -> 0.723201, w[3] -> 0.600289, w[4] -> 0.509753}, 
 {d -> 0.45366, e[1] -> 1.62233, e[2] -> -0.0422339, e[3] -> 0.00383398, e[4] -> 0.00062295, 
  f[1] -> 0.612047, f[2] -> 0.101213, f[3] -> 0.01052, f[4] -> 0.00301115, f[5] -> 0.0011566, 
  w[1] -> 1.64559, w[2] -> 0.732073, w[3] -> 0.606952, w[4] -> 0.513505}, 
 {d -> 0.391675, e[1] -> 1.12346, e[2] -> -0.0312839, e[3] -> 0.00508731, e[4] -> 0.00057867, 
  f[1] -> 0.66944, f[2] -> 0.10069, f[3] -> 0.00828782, f[4] -> 0.00249066, f[5] -> 0.00094629,
  w[1] -> 1.76627, w[2] -> 0.739404, w[3] -> 0.612718, w[4] -> 0.5167}, 
 {d -> 0.337804, e[1] -> 0.86302, e[2] -> -0.0205496, e[3] -> 0.00576674, e[4] -> 0.0004715, 
  f[1] -> 0.720327, f[2] -> 0.0992223, f[3] -> 0.00622569, f[4] -> 0.00205813, f[5] -> 0.00076775, 
  w[1] -> 1.89202, w[2] -> 0.745497, w[3] -> 0.61776, w[4] -> 0.519453}, 
 {d -> 0.290804, e[1] -> 0.699762, e[2] -> -0.0103062, e[3] -> 0.0059519, e[4] -> 0.00033974, 
  f[1] -> 0.765807, f[2] -> 0.0971481, f[3] -> 0.00436178, f[4] -> 0.00170569, f[5] -> 0.00061643, 
  w[1] -> 2.02429, w[2] -> 0.750586, w[3] -> 0.622207, w[4] -> 0.521846}, 
 {d -> 0.249854, e[1] -> 0.58627, e[2] -> -0.00077367, e[3] -> 0.00573595, e[4] -> 0.00021694, 
  f[1] -> 0.806748, f[2] -> 0.0947268, f[3] -> 0.00270953, f[4] -> 0.00142254, f[5] -> 0.00048749, 
  w[1] -> 2.16482, w[2] -> 0.754852, w[3] -> 0.626158, w[4] -> 0.523947}, 
 {d -> 0.214236, e[1] -> 0.501966, e[2] -> 0.00791363, e[3] -> 0.00521268, e[4] -> 0.0001302, 
  f[1] -> 0.843911, f[2] -> 0.0921519, f[3] -> 0.0012696, f[4] -> 0.00119707, f[5] -> 0.00037626, 
  w[1] -> 2.31566, w[2] -> 0.758441, w[3] -> 0.629693, w[4] -> 0.525803}, 
 {d -> 0.183289, e[1] -> 0.43637, e[2] -> 0.0156862, e[3] -> 0.00446932, e[4] -> 0.00009971, 
  f[1] -> 0.877984, f[2] -> 0.0895681, f[3] -> 0.00003461, f[4] -> 0.00101828, f[5] -> 0.00027862, 
  w[1] -> 2.47937, w[2] -> 0.761466, w[3] -> 0.632875, w[4] -> 0.527456}, 
  {d -> 0.156415, e[1] -> 0.38353, e[2] -> 0.0225133, e[3] -> 0.00358335, e[4] -> 0.00013853, 
   f[1] -> 0.909598, f[2] -> 0.087087, f[3] -> -0.00100725, f[4] -> 0.00087648, f[5] -> 0.00019123, 
   w[1] -> 2.65922, w[2] -> 0.76402, w[3] -> 0.635754, w[4] -> 0.528936}};

(* OPENING FOR M~ VERTEX SOLUTION *)
(* Asymptotic behavior  at the Fracture Tip *)
OpeningTipMT[rho_,n_,i_?IntegerQ] = (1 - rho)^((n + 4)/(4 (n + 1)))/Sqrt[h[i - 1, (5 n + 8)/(2 (n + 1)), 2]]G[i - 1, (5 n + 8)/(2 (n + 1)), 2, rho];

(* Particular soltuion having the proper singularity at the source *)
OpeningSourceMT[rho_,n_] = If[ n == 1,
                               8./\[Pi] (Sqrt[1 - rho^2] - rho ArcCos[rho]),
                               4/(Sqrt[\[Pi]] (1 - n)) Gamma[(3 - n)/2]/Gamma[(4 - n)/2] (Sqrt[1 - rho^2] - 
                               (Sqrt[\[Pi]] rho^(2 - n) Gamma[(n - 2)/2])/(2 Gamma[(n - 1)/2]) - 
                               Hypergeometric2F1[1/2, (n - 2)/2, n/2, rho^2]/(2 - n))
                           ];
 
 (* M~ Vertex opening solution *)
OpeningMT[rho_,n_] :=
    2^(1/2)/\[Pi] (d OpeningSourceMT[rho,n] + Sum[f[i] OpeningTipMT[rho,n,i],{i,1,5}])/. (CoefficientsMT[[IntegerPart[10 n]]]);

(* PRESSURE FOR M~ VERTEX SOLUTION *)
(* Asymptotic behavior  at the Fracture Tip *)
PressureTipMT[rho_,n_,i_?IntegerQ] :=
    -((1 - rho)^(-((3 n)/(4 (n + 1))))/Sqrt[h[i - 1, (n + 4)/(2 (n + 1)), 2]]) G[i - 1, (n + 4)/(2 (n + 1)),  2, rho] + w[i] /. (CoefficientsMT[[IntegerPart[10 n]]]);

(* Particular soltuion having the proper singularity at the source *)
PressureSourceMT[rho_,n_] = If[ n == 1,
                                Log[2/rho] - 1,
                                -(rho^(1 - n)/(1 - n)) + Sqrt[\[Pi]]/(2 (1 - n)) Gamma[(3 - n)/2]/Gamma[(4 - n)/2]
                            ];

(* M~ Vertex pressure solution *)
PressureMT[rho_,n_] :=
    (d PressureSourceMT[rho,n] + Sum[e[i] PressureTipMT[rho,n,i],{i,1,4}])/. (CoefficientsMT[[IntegerPart[10 n]]]);

MTILDEVERTEXSOLUTION[rho_,n_] := {    
									FractureLength  -> 2^(1/2)/\[Pi], 
     								Opening    		-> OpeningMT[rho, n], 
     								NetPressure     -> PressureMT[rho, n]
     								}

(*-------------------------------------K~VERTEX SOLUTION-------------------------------------*)

KTILDEVERTEXSOLUTION[rho_] = {  
								FractureLength	-> 2^(1/2)/\[Pi], 
                                Opening    		-> 1/(2^(1/4) \[Pi]^(1/2)) Sqrt[1 - rho^2],
                                NetPressure     -> \[Pi]^(3/2)/2^(15/4)
                                }
                             
End[]

EndPackage[]




