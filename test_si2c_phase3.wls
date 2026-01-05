#!/usr/bin/env wolframscript

Get["EpidCRN.wl"];

(* SI2Co model *)
RHS = {
  -be1*i12*s - be2*i12*s - mu0*s - (i1*mu1*R1*s)/s0 - (i2*mu2*R2*s)/s0 - (i12*mu3*R3*s)/s0 + mu0*s0,
  -et1*i1*i12 - ga1*i1*i2 - i1*mu1 + be1*i12*s + (i1*mu1*R1*s)/s0,
  -ga2*i1*i2 - et2*i12*i2 - i2*mu2 + be2*i12*s + (i2*mu2*R2*s)/s0,
  et1*i1*i12 + ga1*i1*i2 + ga2*i1*i2 + et2*i12*i2 - i12*mu3 + (i12*mu3*R3*s)/s0
};

(* i1=0 face *)
RHSc1 = RHS /. {i1 -> 0};

(* Parameters *)
cN = {s0 -> 1, mu0 -> 0.1, mu1 -> 0.1, mu2 -> 0.1, mu3 -> 0.1,
      be1 -> 0, be2 -> 0.15, et1 -> 0, et2 -> 0.1,
      ga1 -> 0, ga2 -> 0, R1 -> 1.8, R2 -> 1.5, R3 -> 1.2};

Print["RHSc1 length: ", Length[RHSc1]];
Print["RHSc1[[{1,3,4}]] length: ", Length[RHSc1[[{1, 3, 4}]]]];

rhN3D = RHSc1[[{1, 3, 4}]] //. cN;
Print["rhN3D: ", rhN3D];
Print["Variables in rhN3D: ", Variables[rhN3D]];
Print["Number of variables: ", Length[Variables[rhN3D]]];

(* Try phase3 *)
Print["\nTesting phase3:"];
{fp3, eig3, stab3, pFP3, pStream3} = phase3[rhN3D, 50, 12, 0.01, StreamScale -> 0.05];
