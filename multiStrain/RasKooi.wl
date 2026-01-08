(* ::Package:: *)

(* ::Section:: *)
(*Two-Strain SIR MODEL DEFINITION AND BOUNDARY ANALYSIS*)


ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<< EpidCRN`;
Get[FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels", "red.wl"}]];

prF[expr_] := expr /. {
    (* Variables *)
    s0 -> Subscript[s, 0],
    i1 -> Subscript[i, 1],
    i2 -> Subscript[i, 2],
    r1 -> Subscript[r, 1],
    r2 -> Subscript[r, 2],
    sv -> Subscript[s, v],
    s1 -> Subscript[s, 1],
    s2 -> Subscript[s, 2],
    i12 -> Subscript[i, 12],
    i21 -> Subscript[i, 21],
    v1 -> Subscript[v, 1],
    v2 -> Subscript[v, 2],
    (* Parameters - Greek letters *)
    be1 -> Subscript[\[Beta], 1],
    be2 -> Subscript[\[Beta], 2],
    bev -> Subscript[\[Beta], v],
    ga1 -> Subscript[\[Gamma], 1],
    ga2 -> Subscript[\[Gamma], 2],
    al1 -> Subscript[\[Alpha], 1],
    al2 -> Subscript[\[Alpha], 2],
    mu1 -> Subscript[\[Mu], 1],mu2 -> Subscript[\[Mu], 2],
    mu -> \[Mu],
    muv -> Subscript[\[Mu], v],
    si1 -> Subscript[\[Sigma], 1],
    si2 -> Subscript[\[Sigma], 2],
    siv -> Subscript[\[Sigma], v],
    (* Population sizes *)
    n -> N,
    m -> M
  };
(* Conservation check: With R class, Total[RHS/.cmu] should simplify to 0 *)
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
cmu={mu2->mu,mu1->mu};csym={be2->be,be1->be,mu2->mu,mu1->mu,si2->si,si1->si,
al2->al,al1->al,n->1,m->1};
Print["Conservation check: ",Total[RHSfull/.cmu]//FullSimplify//prF];

(* Working model: Drop R equation *)
RHS =  Drop[RHSfull,-1]//FullSimplify;
var = {s0, i1, i2, r1, r2, s1, s2, i12, i21, sv, v1, v2};

{RN, rts, spe, alp, bet, gam,rnR} = ODE2RN[RHS, var,prF];
(*Print[RN // Length," reactions and transitions: ",Transpose[{RN,rts}]//prF
//MatrixForm];*)

Print[var//Column//prF,"'= ",RHS//prF//MatrixForm];
Print["verif RHS-gam.rts=",gam . rts-RHS//Simplify];

mSi = minSiph[var, RN][[1]];
Print["minSiph=", mSi//prF];
(* MAXIMAL ELIMINABLE SETS *)
{maxId, maxSets,qss} = maxElim[RN, rts, var];
Print[maxSets//Length, "reduced models are",maxSets//prF];
Print["qss1=",(FullSimplify/@qss[[1]])//prF//MatrixForm];
Print["qss2=",(FullSimplify/@qss[[2]])//prF//MatrixForm];
Print["qss3=",(FullSimplify/@qss[[3]])//prF//MatrixForm];
(* Boundary analysis *)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0,F,V, K, R0A, infV} = bdAn[RN, rts, var];
R0A= R0A//remN;
Print[" cDFE=", E0,
" reproduction funcs=",R0A//prF ];
(* REPRODUCTION   NUMBERS *)
R01 = R0A[[1]] /. E0//FullSimplify ;
R02 = R0A[[2]] /. E0//FullSimplify ;
Print["R1(E0)=", R01//prF, " R2(E0)=", R02//prF];


per={1,4,2,3,5,6};K=K[[per,per]];
Print["NGM=",K//MatrixForm//prF];
Print["NGM^2=",K . K//MatrixForm//prF];
K //Eigenvalues//prF


(* Boundary analysis when EA rational*)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, Esys, EA} = bdCom[RN, rts, var];

Print[  EA // Length," non DFE bd pts:"];
E1=var/.EA[[2]][[1]]//FullSimplify;E2=var/.EA[[1]][[1]]//FullSimplify;
Print["E1=",E1//prF];


(* INVASION  NUMBERS *)
R12 = R0A[[1]] /. Thread[var -> E2]//FullSimplify ;
R21 = R0A[[2]] /. Thread[var -> E1]//FullSimplify ;
Print["R1(E0)=", R01//prF, " R2(E0)=", R02//prF, " R1(E2)=", R12//prF, " R2(E1)=", R21//prF," the inters. is at",cb2=
Solve[R01==R02,be2]//Flatten]


jac=Grad[RHS,var];j1=jac/.var -> E1;j1//MatrixForm
det1=Det[j1/.csym]//Factor
det1//Length


det1[[3]]


ch1=CharacteristicPolynomial[(j1/.csym),u];


chs=ch1//Factor;chs//Length


chs[[1]]


chs[[2]]


Print["conjectured coexistence when"];
Timing[BIC=Reduce[And@@cp&&R01>1&&R12>1&&1<R21&&R02>1];]
BICs=FullSimplify[BIC,cp];BICs//prF
gr1=D[R21,{{be1,be2}}]/.cb2//FullSimplify;
gr2=D[R12,{{be1,be2}}]/.cb2//FullSimplify;
Print["example of coexistence:"];
coP=FindInstance[(BIC),par][[1]];
Print["R01=",R01/.coP//N," R12=",R12/.coP//N," R02=",R02/.coP//N," R21=",R21/.coP//N]


(* ::Section:: *)
(*TEST REDUCTIONS*)


T = 2.5;
eps = 0.001;

(* MODEL 1: U1 - approximates near DFE *)
If[Length[maxId] >= 1,
  Print["\n========== MODEL 1: Eliminate ", var[[#]]& /@ maxId[[1]], " =========="];
  params1 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.3, be2 -> 0.25, bev -> 0.4,
            ga1 -> 0.1, ga2 -> 0.08,
            al1 -> 0.01, al2 -> 0.01,
            si1 -> 0.5, si2 -> 0.5, siv -> 0.5};
  ic01 = {95, eps, eps, eps, eps, eps, eps, eps, eps, 950, eps, eps};

  {qss1, Xvars1, err1, solF1, solR1} = testRed[RHS, var, maxId[[1]], params1, T, ic01];
  Print["QSS: ", qss1//Short];
  Print["Slow vars: ", Xvars1];
  Print["L2 errors: ", Thread[Xvars1 -> err1]];
  Print["Mean error: ", Mean[err1]];

  plots1 = Table[
    Plot[{Evaluate[Xvars1[[i]][t] /. solF1], Evaluate[Xvars1[[i]][t] /. solR1]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars1[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars1]]}
  ];
  Print[GraphicsRow[plots1, ImageSize -> 700]];
];

(* MODEL 2: U2 *)
If[Length[maxId] >= 2,
  Print["\n========== MODEL 2: Eliminate ", var[[#]]& /@ maxId[[2]], " =========="];
  params2 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.4, be2 -> 0.3, bev -> 0.5,
            ga1 -> 0.5, ga2 -> 0.4,
            al1 -> 0.02, al2 -> 0.02,
            si1 -> 0.5, si2 -> 0.5, siv -> 0.5};
  ic02 = {50, 5, 4, 3, 2, 10, 8, 2, 2, 500, 100, 80};

  {qss2, Xvars2, err2, solF2, solR2} = testRed[RHS, var, maxId[[2]], params2, T, ic02];
  Print["QSS: ", qss2//Short];
  Print["Slow vars: ", Xvars2];
  Print["L2 errors: ", Thread[Xvars2 -> err2]];
  Print["Mean error: ", Mean[err2]];

  plots2 = Table[
    Plot[{Evaluate[Xvars2[[i]][t] /. solF2], Evaluate[Xvars2[[i]][t] /. solR2]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars2[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars2]]}
  ];
  Print[GraphicsRow[plots2, ImageSize -> 700]];
];

(* MODEL 3: U3 - near E2 (strain 2 endemic) *)
If[Length[maxId] >= 3,
  Print["\n========== MODEL 3: Eliminate ", var[[#]]& /@ maxId[[3]], " =========="];
  params3 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.3, be2 -> 0.5, bev -> 0.6,
            ga1 -> 0.5, ga2 -> 0.4,
            al1 -> 0.02, al2 -> 0.02,
            si1 -> 0.5, si2 -> 0.5, siv -> 0.5};
  ic03 = {50, eps, 8, eps, 5, eps, 10, 3, eps, 500, eps, 150};

  {qss3, Xvars3, err3, solF3, solR3} = testRed[RHS, var, maxId[[3]], params3, T, ic03];
  Print["QSS: ", qss3//Short];
  Print["Slow vars: ", Xvars3];
  Print["L2 errors: ", Thread[Xvars3 -> err3]];
  Print["Mean error: ", Mean[err3]];

  plots3 = Table[
    Plot[{Evaluate[Xvars3[[i]][t] /. solF3], Evaluate[Xvars3[[i]][t] /. solR3]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars3[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars3]]}
  ];
  Print[GraphicsRow[plots3, ImageSize -> 700]];
];

(* MODEL 4: U4 - near E1 (strain 1 endemic) *)
If[Length[maxId] >= 4,
  Print["\n========== MODEL 4: Eliminate ", var[[#]]& /@ maxId[[4]], " =========="];
  params4 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.5, be2 -> 0.3, bev -> 0.6,
            ga1 -> 0.4, ga2 -> 0.5,
            al1 -> 0.02, al2 -> 0.02,
            si1 -> 0.5, si2 -> 0.5, siv -> 0.5};
  ic04 = {50, 8, eps, 5, eps, 10, eps, eps, 3, 500, 150, eps};

  {qss4, Xvars4, err4, solF4, solR4} = testRed[RHS, var, maxId[[4]], params4, T, ic04];
  Print["QSS: ", qss4//Short];
  Print["Slow vars: ", Xvars4];
  Print["L2 errors: ", Thread[Xvars4 -> err4]];
  Print["Mean error: ", Mean[err4]];

  plots4 = Table[
    Plot[{Evaluate[Xvars4[[i]][t] /. solF4], Evaluate[Xvars4[[i]][t] /. solR4]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars4[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars4]]}
  ];
  Print[GraphicsRow[plots4, ImageSize -> 700]];
];

(* MODEL 5: U5 *)
If[Length[maxId] >= 5,
  Print["\n========== MODEL 5: Eliminate ", var[[#]]& /@ maxId[[5]], " =========="];
  params5 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.35, be2 -> 0.35, bev -> 0.5,
            ga1 -> 0.3, ga2 -> 0.3,
            al1 -> 0.02, al2 -> 0.02,
            si1 -> 0.5, si2 -> 0.5, siv -> 0.5};
  ic05 = {60, 3, 3, 2, 2, 8, 8, 2, 2, 600, 80, 80};

  {qss5, Xvars5, err5, solF5, solR5} = testRed[RHS, var, maxId[[5]], params5, T, ic05];
  Print["QSS: ", qss5//Short];
  Print["Slow vars: ", Xvars5];
  Print["L2 errors: ", Thread[Xvars5 -> err5]];
  Print["Mean error: ", Mean[err5]];

  plots5 = Table[
    Plot[{Evaluate[Xvars5[[i]][t] /. solF5], Evaluate[Xvars5[[i]][t] /. solR5]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars5[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars5]]}
  ];
  Print[GraphicsRow[plots5, ImageSize -> 700]];
];

(* MODEL 6: U6 *)
If[Length[maxId] >= 6,
  Print["\n========== MODEL 6: Eliminate ", var[[#]]& /@ maxId[[6]], " =========="];
  params6 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.4, be2 -> 0.4, bev -> 0.55,
            ga1 -> 0.35, ga2 -> 0.35,
            al1 -> 0.015, al2 -> 0.015,
            si1 -> 0.6, si2 -> 0.6, siv -> 0.6};
  ic06 = {55, 4, 4, 3, 3, 9, 9, 2, 2, 550, 90, 90};

  {qss6, Xvars6, err6, solF6, solR6} = testRed[RHS, var, maxId[[6]], params6, T, ic06];
  Print["QSS: ", qss6//Short];
  Print["Slow vars: ", Xvars6];
  Print["L2 errors: ", Thread[Xvars6 -> err6]];
  Print["Mean error: ", Mean[err6]];

  plots6 = Table[
    Plot[{Evaluate[Xvars6[[i]][t] /. solF6], Evaluate[Xvars6[[i]][t] /. solR6]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars6[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars6]]}
  ];
  Print[GraphicsRow[plots6, ImageSize -> 700]];
];

(* MODEL 7: U7 *)
If[Length[maxId] >= 7,
  Print["\n========== MODEL 7: Eliminate ", var[[#]]& /@ maxId[[7]], " =========="];
  params7 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.45, be2 -> 0.35, bev -> 0.6,
            ga1 -> 0.4, ga2 -> 0.35,
            al1 -> 0.02, al2 -> 0.02,
            si1 -> 0.55, si2 -> 0.55, siv -> 0.55};
  ic07 = {52, 5, 3, 3, 2, 10, 7, 2, 2, 520, 110, 75};

  {qss7, Xvars7, err7, solF7, solR7} = testRed[RHS, var, maxId[[7]], params7, T, ic07];
  Print["QSS: ", qss7//Short];
  Print["Slow vars: ", Xvars7];
  Print["L2 errors: ", Thread[Xvars7 -> err7]];
  Print["Mean error: ", Mean[err7]];

  plots7 = Table[
    Plot[{Evaluate[Xvars7[[i]][t] /. solF7], Evaluate[Xvars7[[i]][t] /. solR7]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars7[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars7]]}
  ];
  Print[GraphicsRow[plots7, ImageSize -> 700]];
];

(* MODEL 8: U8 *)
If[Length[maxId] >= 8,
  Print["\n========== MODEL 8: Eliminate ", var[[#]]& /@ maxId[[8]], " =========="];
  params8 = {n -> 100, m -> 1000, mu -> 0.02, muv -> 0.5,
            be1 -> 0.35, be2 -> 0.45, bev -> 0.6,
            ga1 -> 0.35, ga2 -> 0.4,
            al1 -> 0.02, al2 -> 0.02,
            si1 -> 0.55, si2 -> 0.55, siv -> 0.55};
  ic08 = {52, 3, 5, 2, 3, 7, 10, 2, 2, 520, 75, 110};

  {qss8, Xvars8, err8, solF8, solR8} = testRed[RHS, var, maxId[[8]], params8, T, ic08];
  Print["QSS: ", qss8//Short];
  Print["Slow vars: ", Xvars8];
  Print["L2 errors: ", Thread[Xvars8 -> err8]];
  Print["Mean error: ", Mean[err8]];

  plots8 = Table[
    Plot[{Evaluate[Xvars8[[i]][t] /. solF8], Evaluate[Xvars8[[i]][t] /. solR8]},
      {t, 0, T}, PlotStyle -> {{Blue, Thick}, {Red, Dashed, Thick}},
      PlotLegends -> If[i == 1, {"Full", "Reduced"}, None],
      PlotLabel -> ToString[Xvars8[[i]]], ImageSize -> 200],
    {i, 1, Min[3, Length[Xvars8]]}
  ];
  Print[GraphicsRow[plots8, ImageSize -> 700]];
];


 ClearAll["Global`*"];
  <<EpidCRN`;
  ?bdCom
