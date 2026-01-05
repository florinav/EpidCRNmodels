(* Test script for factoring det J1 *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<< EpidCRN`;

(* Model definition from RasKooi.wl *)
RHSfull =  {
    - s0/m * (be1 *v1 + be2 *v2) + mu * (n - s0),
    be1/m * s0 * v1 - (ga1 + mu1) * i1,
    be2/m * s0 * v2 - (ga2 + mu2) * i2,
    ga1 * i1 - (al1 + mu) * r1,
    ga2 * i2 - (al2 + mu) * r2,
    al1 * r1 - mu * s1 - be2 * si2/m * s1 * v2,
    al2 * r2 - mu * s2 - be1 * si1/m * s2 * v1,
    be2 * si2/m * s1 * v2 - (ga2 + mu2) * i12,
    be1 * si1/m * s2 * v1 - (ga1 + mu1) * i21,
    -bev/n * sv * (i1 + siv * i21 + i2 + siv * i12) + muv * (m - sv),
    bev/n * sv * (i1 + siv * i21) - muv * v1,
    bev/n * sv * (i2 + siv * i12) - muv * v2,
    ga1 * i21 + ga2 * i12 - mu * R
  };

RHS = Drop[RHSfull,-1];
var = {s0, i1, i2, r1, r2, s1, s2, i12, i21, sv, v1, v2};

{RN, rts, spe, alp, bet, gam,rnR} = ODE2RN[RHS, var];
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, Esys, EA} = bdCom[RN, rts, var];
E1=var/.EA[[2]][[1]]//FullSimplify;

(* Symmetry condition *)
csym={be2->be,be1->be,mu2->mu,mu1->mu,si2->si,si1->si,
      al2->al,al1->al,n->1,m->1};

(* Compute Jacobian at E1 with symmetry *)
jac=Grad[RHS,var];
j1=jac/.var -> E1;
j1s = j1/.csym;

Print["Computing det1..."];
det1raw = Det[j1s];
Print["det1 raw LeafCount=", LeafCount[det1raw]];

(* Try different factorizations *)
Print["\n--- Factor ---"];
det1a = det1raw//Factor;
Print["Length=", det1a//Length];
If[Head[det1a]===Times, Print["First 3 factors=", Take[List@@det1a,UpTo[3]]]];

Print["\n--- Collect by be ---"];
det1b = det1raw//Collect[#,be]&//Factor;
Print["Length=", det1b//Length];

Print["\n--- Collect by muv ---"];
det1c = det1raw//Collect[#,muv]&//Factor;
Print["Length=", det1c//Length];

Print["\n--- Simplify then Factor ---"];
det1d = det1raw//Simplify//Factor;
Print["Length=", det1d//Length];
If[Head[det1d]===Times, Print["First 3 factors=", Take[List@@det1d,UpTo[3]]]];

Print["\n--- Together then Factor ---"];
det1e = det1raw//Together//Factor;
Print["Length=", det1e//Length];

Print["\n--- Factor with Extension->Sqrt ---"];
det1f = det1raw//Factor[#,Extension->Automatic]&;
Print["Length=", det1f//Length];
