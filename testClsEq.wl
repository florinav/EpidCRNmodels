(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

RN = {
  0 -> "x1",
  "x1" + "x2" -> 2*"x2",
  "x2" + "x3" -> 2*"x3",
  "x2" + "B1" -> "S1" + "x2" + "B1",
  "x2" + "B2" -> "S2" + "x2" + "B2",
  "S1" -> "B1",
  "S2" -> "B2",
  "B1" -> "R",
  "B2" -> "R",
  "R" -> "x3",
  "x1" -> 0,
  "x2" -> 0,
  "x3" -> 0
};

(* Extract species and variables *)
{spe, alp, bet, gam, Rv, RHS1, def} = extMat[RN];
var = ToExpression[spe];

(* Rates from OSN.nb *)
l4 = be1/(1 + a1*x2 + ep1*B1);
l5 = be2/(1 + a2*x2 + ep2*B2);

rts = {
  La,
  al*x1*x2,
  et*x2*x3,
  be1*B1*x2/(B1*ep1 + a1*x2 + 1),
  be2*B2*x2/(B2*ep2 + a2*x2 + 1),
  ga1*S1,
  ga2*S2,
  mu1*B1,
  mu2*B2,
  om*R,
  mu*x1,
  mun*x2,
  mu*x3
};
Print["RN has ", Length[RN], " reactions, var=", var];

(* Compute RHS *)
RHS = gam . rts;

(* Run bdAn *)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars} = bdAn[RN, rts, var];
Print["mSi=", mSi];
Print["E0=", E0];

(* Define formatting function *)
prF[ex_] := ex /.
  {x1 -> Subscript[x, 1], x2 -> Subscript[x, 2], x3 -> Subscript[x, 3],
   S1 -> Subscript[S, 1], S2 -> Subscript[S, 2],
   B1 -> Subscript[B, 1], B2 -> Subscript[B, 2],
   et -> \[Eta], de -> \[Delta],
   be -> \[Beta], la -> \[Lambda], La -> \[CapitalLambda], mu -> \[Mu],
   mu1 -> Subscript[\[Mu], 1], mu2 -> Subscript[\[Mu], 2],
   mup -> Subscript[\[Mu], P], mun -> Subscript[\[Mu], n],
   mue -> Subscript[\[Mu], e],
   ep1 -> Subscript[\[Epsilon], 1], ep2 -> Subscript[\[Epsilon], 2],
   ga1 -> Subscript[\[Gamma], 1], ga2 -> Subscript[\[Gamma], 2],
   a1 -> Subscript[\[Alpha], 1], a2 -> Subscript[\[Alpha], 2], al -> \[Alpha],
   be1 -> Subscript[\[Beta], 1], be2 -> Subscript[\[Beta], 2],
   om -> \[Omega]};

(* Test clsEq function for endemic equilibrium *)
Print[""];
Print["Testing clsEq for endemic equilibrium with tim=60:"];
endemic = clsEq[RHS, var, mSi, E0, prF, 60];
