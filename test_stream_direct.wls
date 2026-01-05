#!/usr/bin/env wolframscript

Get["EpidCRN.wl"];

(* Exact same RHS as SI2Co *)
rhN = {0.1 - 0.1*s - 0.27*i12*s - 0.15*i2*s,
       -0.1*i2 - 0.1*i12*i2 + 0.15*i12*s + 0.15*i2*s,
       -0.1*i12 + 0.1*i12*i2 + 0.12*i12*s};

Print["=== Testing StreamPlot3D directly ==="];

(* Test with explicit ranges *)
Print["\nTest 1: Direct StreamPlot3D"];
p1 = StreamPlot3D[rhN,
  {s, 0, 1.3},
  {i2, 0, 0.6},
  {i12, 0, 0.4},
  StreamPoints -> 30,
  StreamMarkers -> "Arrow3D"
];

Print["Head: ", Head[p1]];
Print["Graphics length: ", Length[p1[[1]]]];
str1 = ToString[InputForm[p1]];
Print["Has Arrow: ", StringContainsQ[str1, "Arrow"]];
Print["Has Tube: ", StringContainsQ[str1, "Tube"]];
Print["String length: ", StringLength[str1]];

(* Check PlotRange *)
pr = PlotRange /. AbsoluteOptions[p1, PlotRange];
Print["PlotRange: ", pr];

(* Evaluate vector field at a test point *)
testPt = {s -> 0.5, i2 -> 0.3, i12 -> 0.1};
vecVal = rhN /. testPt;
Print["\nVector at {0.5, 0.3, 0.1}: ", vecVal];
Print["Vector magnitude: ", Norm[vecVal]];

(* Try with MaxRecursion *)
Print["\nTest 2: With MaxRecursion"];
p2 = StreamPlot3D[rhN,
  {s, 0, 1.3},
  {i2, 0, 0.6},
  {i12, 0, 0.4},
  StreamPoints -> 30,
  StreamMarkers -> "Arrow3D",
  MaxRecursion -> 3
];
str2 = ToString[InputForm[p2]];
Print["String length: ", StringLength[str2]];

(* Export to see what's in there *)
If[StringLength[str1] > 1000,
  Print["\nGraphics object seems substantial"],
  Print["\nWARNING: Graphics object seems too small"]
];
