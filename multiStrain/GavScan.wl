(* ::Package:: *)

(* GavScan - Two-strain model with reinfection analysis *)
ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession,NotebookAutoSave->True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;
Get[FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels", "red.wl"}]];
prF[ex_]:= ex/.{i1->Subscript[i,1],i2->Subscript[i,2],r1->Subscript[r,1],r2->Subscript[r,2],y1->Subscript[i,21],y2->Subscript[i,12],be1->Subscript[\[Beta],1],be2->Subscript[\[Beta],2],ga1->Subscript[\[Gamma],1],ga2->Subscript[\[Gamma],2],th1->Subscript[\[Theta],1],th2->Subscript[\[Theta],2],th3->Subscript[\[Theta],3],et1->Subscript[\[Eta],1],et2->Subscript[\[Eta],2],mu1->Subscript[\[Mu],1],mu2->Subscript[\[Mu],2],si1->Subscript[\[Sigma],1],
si2->Subscript[\[Sigma],2],La->\[CapitalLambda],mu->\[Mu]};
RHS = {
  La - be1*i1*s - be2*i2*s - mu*s + r1*th1 + r2*th2 + R*th3 - be1*et1*s*y1 - be2*et2*s*y2,
  -ga1*i1 - i1*mu + be1*i1*s + be1*et1*s*y1,
  -ga2*i2 - i2*mu + be2*i2*s + be2*et2*s*y2,
  ga1*i1 - mu*r1 - be2*i2*r1*si2 - r1*th1 - be2*et2*r1*si2*y2,
  ga2*i2 - mu*r2 - be1*i1*r2*si1 - r2*th2 - be1*et1*r2*si1*y1,
  be2*i2*r1*si2 - ga2*y2 - mu*y2 + be2*et2*r1*si2*y2,
  be1*i1*r2*si1 - ga1*y1 - mu*y1 + be1*et1*r2*si1*y1,
  -mu*R - R*th3 + ga1*y1 + ga2*y2
};
var = {s, i1, i2, r1, r2, y2, y1, R};

cmu={mu2->mu,mu1->mu};cmu3={mu2->mu,mu1->mu,th3->0};
cg={ga1->2 La/mu -mu,ga2->2 La/mu-mu};
Print["verif sum=",Total[RHS/.cmu]//FullSimplify];
{RN, rts, spe, alp, bet, gam,rnR} = ODE2RN[RHS, var,prF];
Print["verif RHS-gam.rts=",gam . rts-RHS//Simplify];

mSi = minSiph[var, RN][[1]];
Print["minSiph=", mSi];
(* MAXIMAL ELIMINABLE SETS *)
{maxId, maxSets,qss} = maxElim[RN, rts, var];
Print["reduced models are",maxSets];
Print["qss1=",(FullSimplify/@qss[[1]])//prF//MatrixForm];
Print["qss2=",(FullSimplify/@qss[[2]])//prF//MatrixForm];
Print["qss3=",(FullSimplify/@qss[[3]])//prF//MatrixForm];
(*{co,nc}=findCores[RN,3];co*)

(* Boundary analysis *)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0,F,V, K, R0A, infV} = bdAn[RN, rts, var];

infV={1,3,2,4};
Print["RHS=",RHS//FullSimplify//prF//MatrixForm, " cDFE=", E0//prF, 
" reproduction funcs=", R0A//prF," NGM=",K[[infV,infV]]//prF//MatrixForm];


(* Boundary analysis when EA rational*)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, EA} = bdCom[RN, rts, var];
Print[  EA // Length," non DFE bd pts:"];
Print["E1=",E1=var/.EA[[2]][[1]]];Print[" E2=",E2=var/.EA[[1]][[1]]];
(* REPRODUCTION and INVASION  NUMBERS *)

R01 = R0A[[1]] /. E0 ;
R02 = R0A[[2]] /. E0 ;
R12 = R0A[[1]] /. Thread[var -> E2] ;
R21 = R0A[[2]] /. Thread[var -> E1] ;
Print["R1(E0)=", R01, " R2(E0)=", R02, " R1(E2)=", R12, " R2(E1)=", R21," the inters. is at",cb2=
Solve[R01==R02,be2]//Flatten]
Print["conjectured coexistence when"];
BIC=Reduce[And@@cp&&R01>1&&R12>1&&1<R21&&R02>1];
BICs=FullSimplify[BIC,cp];BICs//prF
gr1=D[R21,{{be1,be2}}]/.cb2//FullSimplify;
gr2=D[R12,{{be1,be2}}]/.cb2//FullSimplify;
th=Thread[gr1+gr2==0]//FullSimplify;
Print["opposite grads when"];so=Reduce[Join[th,cp]]
fu1=FullSimplify[Reduce[And@@cp&&gr1[[1]]>0&&gr1[[2]]<0],cp];
fu2=FullSimplify[Reduce[And@@cp&&gr2[[1]]<0&&gr1[[2]]>0],cp];
Print["example of coexistence:"];
co=FindInstance[(BIC/.cg),Drop[par,{5,6}]][[1]];
coP=Join[co,cg/.co];
Print["R01=",R01/.coP//N," R12=",R12/.coP//N," R02=",R02/.coP//N," R21=",R21/.coP//N]





gridRes=80;plotInd={1,2};steTol=10^(-8);staTol=10^(-10);choTol=10^(-13);rangeExt=1.5;
(*tMax=300;nIc=8;*)
Timing[(*Step 3:Scanning-now with mSi from step 1*)
{plot,noSol,results}=scan[RHS,var,par,coP,plotInd,mSi,(*mSi auto-provided*)
gridRes,steTol,staTol,choTol,R01,R02,R12,R21,rangeExt];]
plot


Export["GavSca.pdf",plot]


(* ::Section:: *)
(*TEST  REDUCTIONS*)


T = 2.5;
eps = 0.001;

(* MODEL 1: U1 - approximates near DFE *)
Print["\n========== MODEL 1: Eliminate ", var[[#]]& /@ maxId[[1]], " =========="];
params1 = {La -> 0.1, mu -> 0.02,
          be1 -> 0.8, be2 -> 0.6,
          ga1 -> 0.1, ga2 -> 0.08,
          th1 -> 0.01, th2 -> 0.01, th3 -> 0.01,
          et1 -> 0.3, et2 -> 0.3,
          si1 -> 0.5, si2 -> 0.5};
ic01 = {4.5, eps, eps, eps, eps, eps, eps, eps};

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

(* MODEL 2: U2 - endemic coexistence parameters *)
Print["\n========== MODEL 2: Eliminate ", var[[#]]& /@ maxId[[2]], " =========="];
(* Parameters chosen so R1(E0)>1, R2(E0)>1, ensuring endemic equilibrium *)
params2 = {La -> 0.1, mu -> 0.002,
          be1 -> 0.15, be2 -> 0.12,
          ga1 -> 0.5, ga2 -> 0.4,
          th1 -> 0.02, th2 -> 0.02, th3 -> 0.01,
          et1 -> 0.3, et2 -> 0.3,
          si1 -> 0.5, si2 -> 0.5};

(* Compute all invasion reproduction numbers *)
R1E0 = R0A[[1]] /. E0 /. params2;
R2E0 = R0A[[2]] /. E0 /. params2;
R1E2 = R0A[[1]] /. Thread[var -> E2] /. params2;
R2E1 = R0A[[2]] /. Thread[var -> E1] /. params2;

Print["R1(E0)=", R1E0, " R2(E0)=", R2E0];
Print["R1(E2)=", R1E2, " R2(E1)=", R2E1];
Print["Coexistence: ", And[R1E0 > 1, R2E0 > 1, R1E2 > 1, R2E1 > 1]];

(* Find endemic equilibrium and perturb *)
endemicEqs = Join[Thread[RHS == 0], Thread[var > 0]];
endemicSol = FindInstance[endemicEqs /. params2, var, Reals, 1];
If[Length[endemicSol] > 0,
  ic02 = (var /. endemicSol[[1]]) + eps;
  Print["Starting near endemic: ", ic02];,
  (* Fallback if no endemic found *)
  ic02 = {0.6, 0.1, 0.08, 0.05, 0.04, 0.03, 0.03, 0.07};
  Print["Using default IC: ", ic02];
];

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

(* MODEL 3: U3 -  approximates near E2*)
Print["\n========== MODEL 3: Eliminate ", var[[#]]& /@ maxId[[3]], " =========="];
params3 = {La -> 0.1, mu -> 0.02,
          be1 -> 0.8, be2 -> 0.6,
          ga1 -> 0.5, ga2 -> 0.4,
          th1 -> 0.01, th2 -> 0.01, th3 -> 0.01,
          et1 -> 0.3, et2 -> 0.3,
          si1 -> 0.5, si2 -> 0.5};
ic03 = {4.0, eps, 0.15, eps, 0.05, 0.03, eps, 0.07};

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

(* MODEL 4: U4 approximates near E1*)
Print["\n========== MODEL 4: Eliminate ", var[[#]]& /@ maxId[[4]], " =========="];
params4 = {La -> 0.1, mu -> 0.02,
          be1 -> 0.8, be2 -> 0.6,
          ga1 -> 0.5, ga2 -> 0.4,
          th1 -> 0.01, th2 -> 0.01, th3 -> 0.01,
          et1 -> 0.3, et2 -> 0.3,
          si1 -> 0.5, si2 -> 0.5};
ic04 = {4.0, 0.15, eps, 0.05, eps, eps, 0.03, 0.07};

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


(*bdFp is used in bdCom*)bdfp=bdFp[RHS,var,mSi];
Print["rat sols on  siphon faces are"]
bd1=bdfp[[1,1]]//FullSimplify
bd2=bdfp[[2,1]]//FullSimplify


(*second siphon fp cannot be positive*)
E12=remZ[var/.EA[[2]][[2]]];E12//Length
c12=Thread[E12>0];
re=FullSimplify[Reduce[Join[cp,c12]],cp]



?EpidCRN`*




