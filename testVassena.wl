(* ::Package:: *)

SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
<<EpidCRN`;
Get["SignPattern.wl"];

(* Example 2 from Vassena paper (good child selection) *)
RN2 = {d + e -> f, e -> d, f -> 0};
rts2 = {k[1], k[2], k[3]};

{spe2, alp2, bet2, gam2, Rv2, RHS2, def2} = extMat[RN2];
var2 = ToExpression[spe2];
jac2 = Grad[RHS2, var2];
fp2 = Quiet[Solve[Thread[RHS2 == 0], var2]];

Print["Example 2 (good): RHS=", RHS2, " Jac=", jac2 // MatrixForm, " fp=", fp2];

{var2fix, RHS2fix, gam2fix, rts2fix, diag2} = Sfix[RN2, rts2];
jac2fix = Grad[RHS2fix, var2fix];
fp2fix = Quiet[Solve[Thread[RHS2fix == 0], var2fix]];

Print["Fixed: RHS=", RHS2fix, " Jac=", jac2fix // MatrixForm, " fp=", fp2fix];

(* Example 1 from Vassena paper (bad child selection) *)
RN1 = {a -> b + c, b -> c, c -> a};
rts1 = {k[1], k[2], k[3]};

{spe1, alp1, bet1, gam1, Rv1, RHS1, def1} = extMat[RN1];
var1 = ToExpression[spe1];
jac1 = Grad[RHS1, var1];
fp1 = Quiet[Solve[Thread[RHS1 == 0], var1]];

Print["Example 1 (bad): RHS=", RHS1, " Jac=", jac1 // MatrixForm, " fp=", fp1];

{var1fix, RHS1fix, gam1fix, rts1fix, diag1} = Sfix[RN1, rts1];
jac1fix = Grad[RHS1fix, var1fix];
fp1fix = Quiet[Solve[Thread[RHS1fix == 0], var1fix]];

Print["Fixed: RHS=", RHS1fix, " Jac=", jac1fix // MatrixForm, " fp=", fp1fix];




