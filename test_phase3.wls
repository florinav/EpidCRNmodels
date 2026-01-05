#!/usr/bin/env wolframscript

(* Test phase3 with 3D SI2Co system *)
Get["EpidCRN.wl"];

Print["=== Testing phase3 with 3D system ==="];

(* SI2Co model on i1=0 face *)
RHS = {
  -be1*i12*s - be2*i12*s - mu0*s - (i1*mu1*R1*s)/s0 - (i2*mu2*R2*s)/s0 - (i12*mu3*R3*s)/s0 + mu0*s0,
  -ga2*i1*i2 - et2*i12*i2 - i2*mu2 + be2*i12*s + (i2*mu2*R2*s)/s0,
  et2*i12*i2 - i12*mu3 + (i12*mu3*R3*s)/s0
};

(* Set i1=0 for the i1=0 face *)
RHSc1 = RHS /. {i1 -> 0};

(* Extract 3D subsystem (s, i2, i12) *)
RHS3D = RHSc1[[{1, 2, 3}]];

(* Parameters *)
cN = {s0 -> 1, mu0 -> 0.1, mu1 -> 0.1, mu2 -> 0.1, mu3 -> 0.1,
      be1 -> 0, be2 -> 0.15, et1 -> 0, et2 -> 0.1,
      ga1 -> 0, ga2 -> 0, R1 -> 1.8, R2 -> 1.5, R3 -> 1.2};

(* Apply parameter substitutions *)
rhN = RHS3D //. cN;

Print["3D RHS: ", rhN];

(* Test phase3 with custom options *)
{fp3, eig3, stab3, pFP3, pStream3} = phase3[rhN, 50, 8, 0.012,
  StreamScale -> 0.05];

Print["\n=== Results ==="];
Print["Fixed points: ", fp3];
Print["Eigenvalues: ", eig3];
Print["Stabilities: ", stab3];

Print["\n=== Visualization ==="];
Print["Use Show[pFP3, pStream3] to visualize"];
Print["Example: Show[pFP3, pStream3, ImageSize -> Large]"];
Print["pFP3 and pStream3 created successfully"];
