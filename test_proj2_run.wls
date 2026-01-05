#!/usr/bin/env wolframscript

(* Test proj2 *)
Get["test_proj2.wl"];
Get["EpidCRN.wl"];

Print["=== Test proj2 with 3D system ==="];

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

(* Test proj2: project 3D (s,i2,i12) onto 2D (i2,i12) *)
{fpProj, eig, stab, pFP, pStream} = proj2[RHSc1[[{1, 3, 4}]], {s, i2, i12}, {i2, i12}, cN];

Print["\n=== Results ==="];
Print["Projected FPs: ", fpProj];
Print["Stabilities: ", stab];

Print["\n=== Plotting ==="];
Print["pFP and pStream created successfully"];
