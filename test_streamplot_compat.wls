#!/usr/bin/env wolframscript

(* Test StreamPlot3D compatibility for Mathematica 13.3 *)

Print["Mathematica version: ", $Version];

(* Simple test system *)
testRHS = {0.1 - 0.1*x - 0.15*y*x,
           -0.1*y + 0.15*y*x,
           -0.1*z + 0.1*y*z};

Print["\n=== Test 1: Basic StreamPlot3D ==="];
p1 = StreamPlot3D[testRHS, {x, 0, 2}, {y, 0, 2}, {z, 0, 2}];
Print["Head: ", Head[p1]];
Print["Length: ", Length[p1[[1]]]];

Print["\n=== Test 2: With Arrow3D markers ==="];
p2 = StreamPlot3D[testRHS, {x, 0, 2}, {y, 0, 2}, {z, 0, 2},
  StreamMarkers -> "Arrow3D"];
Print["Head: ", Head[p2]];
Print["Length: ", Length[p2[[1]]]];

Print["\n=== Test 3: With Automatic markers ==="];
p3 = StreamPlot3D[testRHS, {x, 0, 2}, {y, 0, 2}, {z, 0, 2},
  StreamMarkers -> Automatic];
Print["Head: ", Head[p3]];
Print["Length: ", Length[p3[[1]]]];

Print["\n=== Test 4: Checking for Arrow primitives ==="];
str1 = ToString[InputForm[p1]];
str2 = ToString[InputForm[p2]];
str3 = ToString[InputForm[p3]];

Print["Test 1 has Arrow: ", StringContainsQ[str1, "Arrow"]];
Print["Test 2 has Arrow: ", StringContainsQ[str2, "Arrow"]];
Print["Test 3 has Arrow: ", StringContainsQ[str3, "Arrow"]];

Print["Test 1 has Tube: ", StringContainsQ[str1, "Tube"]];
Print["Test 2 has Tube: ", StringContainsQ[str2, "Tube"]];
Print["Test 3 has Tube: ", StringContainsQ[str3, "Tube"]];
