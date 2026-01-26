(* ::Package:: *)

(* Test that ClaudeFuckups compiles *)
ClearAll["Global`*"];
Print["Testing ClaudeFuckups.wl compilation..."];

Check[
  <<EpidCRN`,
  Print["ERROR: Package failed to load"];
  Abort[],
  {ToExpression::notstrbox, Get::noopen}
];

Print["SUCCESS: Package loaded without errors"];

(* Test minSiph1 with RHS as third parameter *)
Print["\nTesting minSiph1[var, RN, RHS]:"];
Block[{RHS, var, RN, rts, result, mSi, cDFE, cE0, nonm},
  RHS = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
  var = {s, i};
  {RN, rts} = Take[EpidCRN`ODE2RN[RHS, var], 2];

  result = EpidCRN`minSiph1[var, RN, RHS];
  {mSi, cDFE, cE0, nonm} = result;

  Print["mSi=", mSi];
  Print["cDFE=", cDFE];
  Print["cE0=", cE0];
  Print["nonm=", nonm];
];

(* Test minSiph0 with RHS as third parameter *)
Print["\nTesting minSiph0[var, RN, RHS]:"];
Block[{RHS, var, RN, rts, result, mSi, cDFE, cE0, nonm},
  RHS = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
  var = {s, i};
  {RN, rts} = Take[EpidCRN`ODE2RN[RHS, var], 2];

  result = EpidCRN`minSiph0[var, RN, RHS];
  {mSi, cDFE, cE0, nonm} = result;

  Print["mSi=", mSi];
  Print["cDFE=", cDFE];
  Print["cE0=", cE0];
  Print["nonm=", nonm];
];

Print["\nAll tests passed."];
