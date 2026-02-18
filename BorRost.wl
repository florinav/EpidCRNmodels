(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

(* Boros-Rost 3-species Hopf witness for SIRWS, with birth into S2 *)
var = {S2, i, W}; RHS = {La + ga*W - be*i*S2 - mu2*S2,
  i*(be*S2 - mu3) - nu*i*W,
  nu*i*W - ga*W - mu1*W};
prF[ex_] := ex /.
  {be -> \[Beta], La -> \[CapitalLambda],
   mu1 -> Subscript[\[Mu], 1], mu2 -> Subscript[\[Mu], 2],
   mu3 -> Subscript[\[Mu], 3],
   ga -> \[Gamma], nu -> \[Nu]};
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var, prF];
RHS - gam . rts // FullSimplify  (*Should be {0,0,0}*)

{RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, sipVar} = NGMRN[RN, rts, var];
jac = Grad[RHS, var];
Print["mSi=", mSi, " cE0=", Short[cE0]]

(* Cauchy-Binet child selection analysis *)
{n, m} = {Length[var], Length[RN]};
RM = symRM[alp, n, m];
G = gam . RM;
cpG = CharacteristicPolynomial[G, lam];
co = Rest@Reverse@CoefficientList[Expand[cpG], lam];

nCS = negCS[co, gam, n];
upf = Select[nCS, upfQ];
{gH, edges} = hasseCS[nCS];
Print["dim=", {n, m}, " #nCS=", Length[nCS], " #upf=", Length[upf],
      " #edges=", Length[edges]]
upfRes = upfCS[nCS, gH, spe, rts];

(* Hopf via Cor 5.21 *)
csN = Select[nCS, #[[1]] == n &];
stabN = Select[csN, hurQ[#[[5]]] &];
fs = fixSgn[nCS, n];
Print["fixed sign by order=", fs,
      " #", n, "x", n, " CS=", Length[csN], " dets=", #[[6]] & /@ csN]
Print["upfs exist=", Length[upf] > 0, " stable ", n, "x", n, " CS exist=",
      Length[stabN] > 0, ", Cor. 5.21 Hopf=", Length[stabN] > 0 && Length[upf] > 0]

(* Recipe I *)
riCores = oscRI[nCS, gH, spe, rts];
(* Recipe II *)
riiCores = oscRII[nCS, gH, spe, rts];
