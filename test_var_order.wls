#!/usr/bin/env wolframscript

rhN = {0.1 - 0.1*s - 0.27*i12*s - 0.15*i2*s,
       -0.1*i2 - 0.1*i12*i2 + 0.15*i12*s + 0.15*i2*s,
       -0.1*i12 + 0.1*i12*i2 + 0.12*i12*s};

Print["rhN equations (in order):"];
Print["Eq 1 (should be for variable 1): ", rhN[[1]]];
Print["Eq 2 (should be for variable 2): ", rhN[[2]]];
Print["Eq 3 (should be for variable 3): ", rhN[[3]]];

Print["\nVariables[rhN]: ", Variables[rhN]];
Print["Order: alphabetical or based on RHS?"];

Print["\n The RHS structure suggests:"];
Print["Equation 1 has d(s)/dt on left"];
Print["Equation 2 has d(i2)/dt on left"];
Print["Equation 3 has d(i12)/dt on left"];

Print["\nSo StreamPlot3D needs ranges in order: {s, ...}, {i2, ...}, {i12, ...}"];
Print["But Variables gives: ", Variables[rhN]];
