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

BeginPackage["rectangularMesh`"]

meshPlan::usage =  "meshPlan[ Lx, Nx, Ly, Ny] It computes: conn, centers, Nelts, ElemNei";
getNeighbors::usage = "getNeighbors[ element, Nx, Ny] It returns {left,right,bottom,up} of element.";

Begin["`Private`"] 

(* ----------- Mesh functions ------------ *)
(* get the cell size *)
geth[L_, N_] := 2 L/ N

(* get the Neighbors of one element *)
getElemNeighbors[elem_, Nx_, Ny_] := 
 Module[{j, i, left, right, up, bottom},
  j = Quotient[elem, Nx];
  i =  Mod[elem, Nx];
  If[i == 0, left = elem, left = j*Nx + i - 1];
  If[i == Nx - 1, right = elem, right = j*Nx + i + 1];
  If[j == 0, bottom = elem, bottom = (j - 1)*Nx + i];
  If[j == Ny - 1, up = elem, up = (j + 1)*Nx + i];
  {left, right, bottom, up}]

(* get the Neighbors of all the elements *)
getNeighbors[Nx_, Ny_] := 
  getElemNeighbors[#, Nx, Ny] & /@ Range[Nx Ny];

(* get the coordinates of all the vertexes of the mesh *)
getVertexesCoords[Lx_, Nx_, Ly_, Ny_] := 
 With[{hx = geth[Lx, Nx], hy = geth[Ly, Ny]}, 
  Flatten[Transpose[
    Table[  {-Lx + i hx, -Ly + j hy}, {i, 0, Nx}, {j, 0, Ny}] ] , 1]]

(* get the connectivity matrix between one element and the \
corresponding vertexes *)
getElemVertexesConnectivity[i_, j_, Nx_, 
  Ny_] := {i + (j - 1 ) (Nx + 1), (i + 1) + (j - 1 ) (Nx + 1), 
  i + 1 + j (Nx + 1), (i) + j (Nx + 1)}

(* get the connectivity matrix between all the elements and the \
corresponding vertexes *)
getElemtsVertexesConnectivity[Nx_, Ny_] := 
 Flatten[Transpose[
   Table[getElemVertexesConnectivity[i, j, Nx, Ny], {i, 1, Nx}, {j, 1,
      Ny}] ], 1] 

(* get the coordinates of the center of the elements *)
getElemtsCenters[Lx_, Nx_, Ly_, Ny_] := Module[{conn, vertexcoords},
  conn = getElemtsVertexesConnectivity[Nx, Ny]; 
  vertexcoords = getVertexesCoords[Lx, Nx, Ly, Ny];
  Mean[vertexcoords[[conn[[#]]]] ] & /@ Range[Nx Ny]]
  
  
meshPlan[Lx_, Nx_, Ly_, Ny_] := Module[{hx,hy,vertexes,conn,centers,nelts,neighbors,mesh},

   hx = geth[Lx, Nx]; 
   hy = geth[Lx, Nx]; (* hx, hy element sizes in x & y directions*)
  
    vertexes =  getVertexesCoords[Lx, Nx, Ly, Ny];
   conn =   getElemtsVertexesConnectivity[Nx, Ny];
   centers = getElemtsCenters[Lx, Nx, Ly, Ny];
   nelts = Nx Ny;
   neighbors = getNeighbors[Nx, Ny];
   mesh = <|"hx" -> hx, "hy" -> hy, "Nelts" -> nelts, 
     "vertexes" -> vertexes, "conn" -> conn, "centers" -> centers, 
     "Neighbors" -> neighbors |>;
     
   mesh
   ];
     
     
End[] (* End Private Context *)

EndPackage[]