(* ::Package:: *)

SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
<<EpidCRN`;

(* Example from Helton paper: A -> 2B, 2A + B -> 0 *)
RN = {a -> 2*b, 2*a + b -> 0};
{spe, alp, bet, gam, Rv, RHSorig, defInfo} = extMat[RN];
var=ToExpression[spe];
(* Jacobian and fixed points of original *)
jacOrig = Grad[RHSorig, var];
fpOrig = Solve[Thread[RHSorig == 0], var];

Print["Original RHS=", RHSorig];
Print["Original Jac=", jacOrig // MatrixForm];
Print["Original fixed pts=", fpOrig];

(* Apply sign pattern fix *)
Get["SignPattern.wl"];
rts = Table[k[i], {i, 1, Length[RN]}];
{gamFix, varFix, rtsFix, diag} = signPatFix[gam, var, rts];

(* Build flux: original fluxes + new flux k[3]*bp for reaction bp -> 2b *)
RvFix = Join[Rv, {rtsFix[[-1]] * varFix[[-1]]}];
RHSfix = gamFix . RvFix;

(* Jacobian and fixed points of fixed system *)
jacFix = Grad[RHSfix, varFix];
fpFix = Solve[Thread[RHSfix == 0], varFix];

Print["Fixed RHS=", RHSfix];
Print["Fixed Jac=", jacFix // MatrixForm];
Print["Fixed fixed pts=", fpFix];



