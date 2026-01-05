(* ::Package:: *)

(* ::Section:: *)
(*Two-Strain SIR MODEL - Clean Testing Version*)


ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<< EpidCRN`;
Get[FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels", "red.wl"}]];
var = {s, i1, i2, r1, r2};
RHS = {
  la - (s)*(be1*i1 + be2*i2) + th1*r1 + th2*r2 - mu*s,
  (be1)*i1*(s + si1*r2) - ga1*i1 - mu1*i1,
  (be2)*i2*s + (be2)*si2*i2*r1 - ga2*i2 - mu2*i2,
  ga1*i1 - th1*r1 - (be2)*si2*r1*i2 - mu*r1,
  ga2*i2 - th2*r2 - (be1)*si1*r2*i1 - mu*r2
};
cmu={mu2->mu,mu1->mu};
Print["verif sum=",Total[RHS/.cmu]//FullSimplify];
{RN, rts, spe, alp, bet, gam,rnR} = ODE2RN[RHS, var];
Print["verif RHS-gam.rts=",gam . rts-RHS//Simplify];
Print[RN // Length," reactions and transitions: ",Transpose[{RN,rts}]
//MatrixForm];
mSi = minSiph[var, RN][[1]];
Print["minSiph=", mSi];
(* MAXIMAL ELIMINABLE SETS *)
{maxId, maxSets,qss} = maxElim[RN, rts, var];
Print["reduced models are",maxSets];
Print["qss1=",(FullSimplify/@qss[[1]])//MatrixForm];
(* Boundary analysis *)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F,V,K, R0A, infVars} = 
bdAn[RN, rts, var];
Print["RHS=",RHS//FullSimplify//MatrixForm, " cDFE=", E0, 
" reproduction funcs=", R0A];
{co,nc}=findCores[RN,3];co



{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, EA} = bdCo[RN, rts, var];

Print[  EA // Length," non DFE bd pts:"];
Print["E1=",E1=var/.EA[[2]][[1]]];Print[" E2=",E2=var/.EA[[1]][[1]]];


(* REPRODUCTION and INVASION  NUMBERS *)
params = {la -> 0.05, mu -> 0.002, mu1 -> 0.01, mu2 -> 0.01,
          be1 -> 2, be2 -> 1.5, ga1 -> 1, ga2 -> 0.8,
          th1 -> 0.05, th2 -> 0.04, si1 -> 0.6, si2 -> 0.5};
R10 = R0A[[1]] /. E0 ;
R20 = R0A[[2]] /. E0 ;
R12 = R0A[[1]] /. Thread[var -> E2] ;
R21 = R0A[[2]] /. Thread[var -> E1] ;
Print["R1(E0)=", R10, " R2(E0)=", R20, " R1(E2)=", R12, " R2(E1)=", R21]
cp=Thread[par>0];

Print["conjectured coexistence when"];
re=FullSimplify[Reduce[And@@cp&&R10>1&&R20>1&&R12>1&&R21>1],cp]


(* ::Section:: *)
(*TEST  REDUCTIONS*)


T = 5;
eps = 0.001;

(* MODEL 1: U1 - endemic coexistence parameters *)
If[Length[maxId] >= 1,
  Print["\n========== MODEL 1: Eliminate ", var[[#]]& /@ maxId[[1]], " =========="];
  (* Parameters chosen so R1(E0)>1, R2(E0)>1, ensuring endemic equilibrium *)
  params1 = {la -> 0.05, mu -> 0.002, mu1 -> 0.002, mu2 -> 0.002,
            be1 -> 0.12, be2 -> 0.1,
            ga1 -> 0.5, ga2 -> 0.4,
            th1 -> 0.02, th2 -> 0.02, si1 -> 0.5, si2 -> 0.5};

  (* Compute all invasion reproduction numbers *)
  R1E0 = (be1*(la/mu))/(ga1 + mu1) /. params1;
  R2E0 = (be2*(la/mu))/(ga2 + mu2) /. params1;
  R1E2 = (be1*(E2[[1]] + E2[[5]]*si1))/(ga1 + mu1) /. params1;
  R2E1 = (be2*(E1[[1]] + E1[[4]]*si2))/(ga2 + mu2) /. params1;

  Print["R1(E0)=", R1E0, " R2(E0)=", R2E0];
  Print["R1(E2)=", R1E2, " R2(E1)=", R2E1];
  Print["Coexistence: ", And[R1E0 > 1, R2E0 > 1, R1E2 > 1, R2E1 > 1]];

  (* Find endemic equilibrium and perturb *)
  endemicEqs = Join[Thread[RHS == 0], Thread[var > 0]];
  endemicSol = FindInstance[endemicEqs /. params1, var, Reals, 1];
  If[Length[endemicSol] > 0,
    ic01 = (var /. endemicSol[[1]]) + eps;
    Print["Starting near endemic: ", ic01];,
    (* Fallback if no endemic found *)
    ic01 = {0.7, 0.1, 0.05, 0.1, 0.05};
    Print["Using default IC: ", ic01];
  ];

  {qss1, Xvars1, err1, solF1, solR1} = testRed[RHS, var, maxId[[1]], params1, T, ic01];
  Print["QSS: ", qss1];
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

(* MODEL 2: U2 - approximates near DFE *)
If[Length[maxId] >= 2,
  Print["\n========== MODEL 2: Eliminate ", var[[#]]& /@ maxId[[2]], " =========="];
  params2 = {la -> 0.05, mu -> 0.002, mu1 -> 0.01, mu2 -> 0.01,
            be1 -> 2, be2 -> 1.5, ga1 -> 0.1, ga2 -> 0.08,
            th1 -> 0.05, th2 -> 0.04, si1 -> 0.6, si2 -> 0.5};
  ic02 = {.21, 0.001, 0.005, 0.001, 0.005};

  {qss2, Xvars2, err2, solF2, solR2} = testRed[RHS, var, maxId[[2]], params2, T, ic02];
  Print["QSS: ", qss2];
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

(* MODEL 3: U3 -  approximates near E2*)
If[Length[maxId] >= 3,
  Print["\n========== MODEL 3: Eliminate ", var[[#]]& /@ maxId[[3]], " =========="];
  params3 = {la -> 0.05, mu -> 0.002, mu1 -> 0.01, mu2 -> 0.01,
            be1 -> 2, be2 -> 1.5, ga1 -> 1, ga2 -> 0.8,
            th1 -> 0.05, th2 -> 0.04, si1 -> 0.6, si2 -> 0.5};
  ic03 = {0.8, 0, 0.15, 0, 0.05};

  {qss3, Xvars3, err3, solF3, solR3} = testRed[RHS, var, maxId[[3]], params3, T, ic03];
  Print["QSS: ", qss3];
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

(* MODEL 4: U4 approximates near E1*)
If[Length[maxId] >= 4,
  Print["\n========== MODEL 4: Eliminate ", var[[#]]& /@ maxId[[4]], " =========="];
  params4 = {la -> 0.05, mu -> 0.002, mu1 -> 0.01, mu2 -> 0.01,
            be1 -> 2, be2 -> 1.5, ga1 -> 1, ga2 -> 0.8,
            th1 -> 0.05, th2 -> 0.04, si1 -> 0.6, si2 -> 0.5};
  ic04 = {0.7, 0.1, 0.005, 0.1, 0.005};

  {qss4, Xvars4, err4, solF4, solR4} = testRed[RHS, var, maxId[[4]], params4, T, ic04];
  Print["QSS: ", qss4];
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
