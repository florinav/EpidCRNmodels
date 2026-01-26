(* ::Package:: *)

SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
<<EpidCRN`;
Get["SignPattern.wl"];

(* Example from Helton paper: A -> 2B, 2A + B -> 0 *)
RN = {a -> 2*b, 2*a + b -> 0};
rts = {k[1], k[2]};

(* Get original system *)
{spe, alp, bet, gam, Rv, RHSorig, defInfo} = extMat[RN];
var = ToExpression[spe];
jacOrig = Grad[RHSorig, var];
fpOrig = Solve[Thread[RHSorig == 0], var];

Print["Original RHS=", RHSorig, " Jac=", jacOrig // MatrixForm, " fp=", fpOrig];

(* Apply Sfix *)
{varFix, RHSfix, gamFix, rtsFix, diag} = Sfix[RN, rts];
jacFix = Grad[RHSfix, varFix];
fpFix = Solve[Thread[RHSfix == 0], varFix];

Print["Fixed RHS=", RHSfix, " Jac=", jacFix // MatrixForm, " fp=", fpFix];



