#!/usr/bin/env wolframscript

Get["EpidCRN.wl"];

(* SI2Co model *)
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

{fp3, eig3, stab3, pFP3, pStream3} = phase3[rhN3D, 50, 12, 0.01, StreamScale -> 0.05];

Print["\n=== Debugging pStream3 ==="];
Print["Head[pStream3]: ", Head[pStream3]];
Print["Length[pStream3]: ", Length[pStream3]];
Print["pStream3: ", Short[pStream3, 5]];

Print["\n=== Debugging pFP3 ==="];
Print["Head[pFP3]: ", Head[pFP3]];
Print["pFP3: ", Short[pFP3, 3]];
