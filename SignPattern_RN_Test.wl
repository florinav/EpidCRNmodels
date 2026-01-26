(* ::Package:: *)

(* Test SignPattern with complete RN from Helton paper *)

SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
<<EpidCRN`;

(* Example from Helton paper Section 2.1:
   Reaction 1: A -> 2B (rate k1)
   Reaction 2: 2A + B -> 0 (rate k2) *)

RN = {a[1] -> 2*a[2], 2*a[1] + a[2] -> 0};
Print["RN=", RN];

(* Extract stoichiometric matrix using extMat *)
{spe, alp, bet, gam, Rv, RHSorig, defInfo} = extMat[RN];
Print["spe=", spe, " alp=", alp, " bet=", bet];
Print["Original S (gam)=", gam // MatrixForm];
Print["Original RHS=", RHSorig];

(* Apply sign pattern fix *)
rts = Table[k[i], {i, 1, Length[RN]}];
Print["Original rates=", rts];

Get["SignPattern.wl"];
{gamFix, speFix, rtsFix, diag} = signPatFix[gam, spe, rts];
Print["Fixed S=", gamFix // MatrixForm];
Print["Fixed species=", speFix];
Print["Fixed rates=", rtsFix];

(* Compute extended ODE *)
RHSfix = gamFix . rtsFix;
Print["Extended ODE RHS=", RHSfix // MatrixForm];

(* Print as system *)
Print["Extended ODE system:"];
Do[Print["d/dt ", speFix[[i]], " = ", RHSfix[[i]]], {i, 1, Length[speFix]}];




