#!/usr/bin/env wolframscript

(* Test StreamPlot3D directly without any wrapping *)

Print["=== Direct StreamPlot3D test ==="];

(* Simple 3D system *)
rhN = {0.1 - 0.1*s - 0.27*i12*s - 0.15*i2*s,
       -0.1*i2 - 0.1*i12*i2 + 0.15*i12*s + 0.15*i2*s,
       -0.1*i12 + 0.1*i12*i2 + 0.12*i12*s};

Print["Creating StreamPlot3D..."];

p = StreamPlot3D[rhN,
  {i12, -0.05, 0.4},
  {i2, -0.05, 0.6},
  {s, -0.05, 1.3},
  StreamMarkers -> "Arrow3D",
  StreamStyle -> Directive[Opacity[0.7], Thickness[0.005]],
  StreamColorFunction -> "Rainbow",
  StreamPoints -> 12,
  BoxRatios -> {1, 1, 1},
  StreamScale -> 0.05
];

Print["Head[p]: ", Head[p]];
Print["Structure (first 500 chars): ", StringTake[ToString[InputForm[p]], UpTo[500]]];
Print["Done"];
