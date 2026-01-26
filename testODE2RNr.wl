(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

prF[expr_] := expr /. {
  x1 -> Subscript[x, 1],
  x2 -> Subscript[x, 2],
  x3 -> Subscript[x, 3],
  x4 -> Subscript[x, 4],
  x5 -> Subscript[x, 5],
  x6 -> Subscript[x, 6],
  al3 -> Subscript[\[Alpha], 3],
  al4 -> Subscript[\[Alpha], 4],
  be3 -> Subscript[\[Beta], 3],
  be4 -> Subscript[\[Beta], 4],
  ep3 -> Subscript[\[Epsilon], 3],
  ep4 -> Subscript[\[Epsilon], 4],
  mu -> \[Mu],
  mu1 -> Subscript[\[Mu], 1],
  mu3 -> Subscript[\[Mu], 3],
  mu4 -> Subscript[\[Mu], 4],
  ga3 -> Subscript[\[Gamma], 3],
  ga4 -> Subscript[\[Gamma], 4],
  ga35 -> Subscript[\[Gamma], 31],
  ga46 -> Subscript[\[Gamma], 41],
  et -> \[Eta],
  La -> \[CapitalLambda]
};

RHS = {
 La - mu*x1 - et*x1*x2
           - (be3*x1*x3)/(1 + al3*x1 + ep3*x3)
           - (be4*x1*x4)/(1 + al4*x1 + ep4*x4)
           + ga3*x5*x5 + ga4*x6*x6 + mu1*x2,
 et*x1*x2 - (mu + mu1)*x2,
(be3*x1*x3)/(1 + al3*x1 + ep3*x3)
           - mu3*(x1*x3)/(1 + al3*x1 + ep3*x3)
           - ga35*x3*x5 - mu*x3,
be4*x1*x4/(1 + al4*x1 + ep4*x4)
           - mu4*(x1*x4)/(1 + al4*x1 + ep4*x4)
           - ga46*x4*x6 - mu*x4,
mu3*(x1*x3)/(1 + al3*x1 + ep3*x3)
           - ga3*x5*x5 - mu*x5 + ga35*x3*x5,
mu4*(x1*x4)/(1 + al4*x1 + ep4*x4)
           - ga4*x6*x6 - mu*x6 + ga46*x4*x6
};
var = {x1, x2, x3, x4, x5, x6};

Print["Testing ODE2RNr with rational rates"];
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNr[RHS, var, prF];
Print["verif RN=", (gam . rts - RHS) // Simplify];
