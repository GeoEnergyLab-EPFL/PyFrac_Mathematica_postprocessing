(* ::Package:: *)

(* Paclet Info File *)

(* created 2020/11/02*)

Paclet[
    "Name" -> "PyFracPost", 
    "Version" -> "0.0.2",
    "MathematicaVersion" -> "10+",
    "SystemID" -> {"MacOSX-x86-64","Linux-x86-64","Windows-x-86-64"},
    "Description" -> "Some functions to import easily PyFrac results into Mathematica",
    "Creator" -> "Carlo Peruzzo",
	"Support" -> "carlo.peruzzo@epfl.ch",
	"Organization" -> "https://www.epfl.ch/labs/gel/",
    "Loading" -> "Automatic",
    "Extensions" -> { 
            {"Kernel", Root -> "Kernel", "Context" ->{"PostprocessPyFrac`","getDOUBLEfrac`","getGENERALfrac`","getSINGLEfrac`","rectangularMesh`" }}
        }
]


