#!/usr/bin/env wolframscript

Get["EpidCRN.wl"];

RHS = {
  -be1*i12*s - be2*i12*s - mu0*s - (i1*mu1*R1*s)/s0 - (i2*mu2*R2*s)/s0 - (i12*mu3*R3*s)/s0 + mu0*s0,
  -et1*i1*i12 - ga1*i1*i2 - i1*mu1 + be1*i12*s + (i1*mu1*R1*s)/s0,
  -ga2*i1*i2 - et2*i12*i2 - i2*mu2 + be2*i12*s + (i2*mu2*R2*s)/s0,
  et1*i1*i12 + ga1*i1*i2 + ga2*i1*i2 + et2*i12*i2 - i12*mu3 + (i12*mu3*R3*s)/s0
};

RHSc1 = RHS /. {i1 -> 0};
cN = {s0 -> 1, mu0 -> 0.1, mu1 -> 0.1, mu2 -> 0.1, mu3 -> 0.1,
      be1 -> 0, be2 -> 0.15, et1 -> 0, et2 -> 0.1,
      ga1 -> 0, ga2 -> 0, R1 -> 1.8, R2 -> 1.5, R3 -> 1.2};

rhN3D = RHSc1[[{1, 3, 4}]] //. cN;

Print["=== Testing phase3 with CORRECT variable order ==="];
Print["Equation 1 (d(s)/dt):   ", rhN3D[[1]]];
Print["Equation 2 (d(i2)/dt):  ", rhN3D[[2]]];
Print["Equation 3 (d(i12)/dt): ", rhN3D[[3]]];
Print["Variable order: {s, i2, i12}"];

{fp3, eig3, stab3, pFP3, pStream3} = phase3[rhN3D, {s, i2, i12}, 50, 12, 0.01, StreamScale -> 0.05];

Print["\n=== Results ==="];
Print["Fixed points: ", fp3];
Print["Stabilities: ", stab3];

Print["\n=== Checking pStream3 ==="];
Print["Head: ", Head[pStream3]];
pStreamStr = ToString[InputForm[pStream3]];
If[StringContainsQ[pStreamStr, "Arrow"],
  Print["SUCCESS: pStream3 contains Arrow primitives!"],
  Print["FAIL: pStream3 is empty or missing arrows"]
];
