(* ::Package:: *)

(* GavScan - Two-strain model with reinfection analysis *)
ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession,NotebookAutoSave->True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;
(*Get[FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels", "red.wl"}]];*)
prF[ex_]:= ex/.{i1->Subscript[i,1],i2->Subscript[i,2],r1->Subscript[r,1],r2->Subscript[r,2],y1->Subscript[i,21],y2->Subscript[i,12],be1->Subscript[\[Beta],1],be2->Subscript[\[Beta],2],ga1->Subscript[\[Gamma],1],ga2->Subscript[\[Gamma],2],th1->Subscript[\[Theta],1],th2->Subscript[\[Theta],2],th3->Subscript[\[Theta],3],et1->Subscript[\[Eta],1],et2->Subscript[\[Eta],2],mu1->Subscript[\[Mu],1],mu2->Subscript[\[Mu],2],si1->Subscript[\[Sigma],1],
si2->Subscript[\[Sigma],2],La->\[CapitalLambda],mu->\[Mu]};
RHS = {
  La - be1*i1*s - be2*i2*s - mu*s + r1*th1 + r2*th2 + R*th3 - be1*et1*s*y1 - be2*et2*s*y2,
  -ga1*i1 - i1*mu + be1*i1*s + be1*et1*s*y1,
  be1*i1*r2*si1 - ga1*y1 - mu*y1 + be1*et1*r2*si1*y1,
  -ga2*i2 - i2*mu + be2*i2*s + be2*et2*s*y2,
  be2*i2*r1*si2 - ga2*y2 - mu*y2 + be2*et2*r1*si2*y2,
  ga1*i1 - mu*r1 - be2*i2*r1*si2 - r1*th1 - be2*et2*r1*si2*y2,
  ga2*i2 - mu*r2 - be1*i1*r2*si1 - r2*th2 - be1*et1*r2*si1*y1,
  -mu*R - R*th3 + ga1*y1 + ga2*y2
};
var = {s, i1, y1,  i2, y2,r1,r2,   R};

cth3={th3->0};
cg={ga1->2 La/mu -mu,ga2->2 La/mu-mu};cChu={th1->0,th2->0,th3->0,et1->1,et2->1};
cIm={th1->0,th2->0,th3->0};
Print["verif sum=",Total[RHS]//FullSimplify];
{RN, rts, spe, alp, bet, gam,rnR} = ODE2RNs[RHS, var];(*def and print of RN*)
Print["verif RHS-gam.rts=",gam . rts-RHS//Simplify];

mSi = minSiph[var, RN][[1]];
Print["minSiph=", mSi];

(* Boundary analysis *)
{RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, E0, F, V, K, R0A, infV} = NGMRN[RN, rts, var];
jac=Grad[RHS,var];jac//MatrixForm//prF
infV={1,3,2,4};
Print["RHS=",RHS//FullSimplify//prF//MatrixForm, " cDFE=", E0//prF, 
" reproduction funcs fo not depend on th!!", R0A//prF," NGM=",K[[infV,infV]]//prF//MatrixForm];
(* Boundary equilibria via bdFp *)
EA = bdFp[RHS, var, mSi];E1=var/.EA[[2,1]];E2=var/.EA[[1,1]];
Print[EA // Length," non DFE bd pts:"];
Print["E1=",E1//prF]; Print[" E2=",E2//prF];
(* REPRODUCTION and INVASION  NUMBERS *)

R01 = R0A[[1]] /. E0 ;
R02 = R0A[[2]] /. E0 ;
R12 = R0A[[1]] /. Thread[var -> E2] ;
R21 = R0A[[2]] /. Thread[var -> E1] ;
Print["R1(E0)=", R01//prF, " R2(E0)=", R02//prF, " R1(E2)=", R12//prF, " R2(E1)=", R21//prF," the inters. is at",cb2=
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
va1=Drop[var,1];vaI=Drop[va1,4];eq=Thread[Drop[RHS,1]==0];
(*eq=Join[Thread[Drop[RHS,1]==0],Thread[vaI>0]];*)
so=selZF[Solve[eq,va1]];so//Length
so1=so[[1]]//FullSimplify;so1//prF
pol=-Numerator[Together[RHS[[1]]/.so1]];
polc=FullSimplify/@CoefficientList[pol,s];polc//Length
po=(pol/.s->0)//Factor;po//prF
po/.cIm//prF
p2 = pol /. Thread[var -> E2] //Factor;
p1 = pol /. Thread[var -> E1]//Factor;p1//prF
p1/.cIm//prF





(* ::Section:: *)
(*RELAY TABLE*)

(* v_k = gamma_k + mu  (GavC.tex notation) *)
v1=ga1+mu; v2=ga2+mu;

(* Transversal blocks  M_k(s, r_{3-k})  (rank-one update of -v_k I, GavC.tex model) *)
M1fun[sv_,r2v_]:=be1*Outer[Times,{sv,si1*r2v},{1,et1}]-v1*IdentityMatrix[2];
M2fun[sv_,r1v_]:=be2*Outer[Times,{sv,si2*r1v},{1,et2}]-v2*IdentityMatrix[2];
(* Perron (dominant) eigenvalue of M_k *)
al1[sv_,r2v_]:=be1*(sv+et1*si1*r2v)-v1;
al2[sv_,r1v_]:=be2*(sv+et2*si2*r1v)-v2;

(* Explicit boundary equilibria (GavC.tex eq. eGa) *)
s0e=La/mu;
s1e=v1/be1;  s2e=v2/be2;
tg1e=ga1/(ga1+th1+mu);  tg2e=ga2/(ga2+th2+mu);   (* tilde gamma_k *)
r1st=tg1e*(s0e-s1e);    r2st=tg2e*(s0e-s2e);     (* r_k^* at E_k *)
i1st=(1-tg1e)*(s0e-s1e); i2st=(1-tg2e)*(s0e-s2e);
Print["E1 verify (s,i1,r1 -- all 0=OK): ",FullSimplify/@{s1e-E1[[1]],i1st-E1[[2]],r1st-E1[[6]]}];
Print["E2 verify (s,i2,r2 -- all 0=OK): ",FullSimplify/@{s2e-E2[[1]],i2st-E2[[4]],r2st-E2[[7]]}];

(* Four reproduction/invasion numbers with symbolic formulas and cross-check vs cell 1 *)
R01e=(be1/v1)*s0e;
R02e=(be2/v2)*s0e;
R21e=(be2/v2)*(s1e+et2*si2*r1st);   (* R_2(E_1), eq:R21Gav *)
R12e=(be1/v1)*(s2e+et1*si1*r2st);   (* R_1(E_2), eq:R12Gav *)
Print["Four numbers vs cell-1 (all 0=OK): ",FullSimplify/@{R01e-R01,R02e-R02,R21e-R21,R12e-R12}];
Print["R1(E0)=",R01e//prF,"  R2(E0)=",R02e//prF,
      "\nR2(E1)=",R21e//prF,"  R1(E2)=",R12e//prF];

(* Relay table: six regimes from the four numbers (tab:Gav / Thm t:Gav-DFE) *)
relayClass[r1_,r2_,r21_,r12_]:=Which[
  r1<1&&r2<1,       "DFE stable",
  r1>1&&r2<1,       "CEP strain 1 (E1 LAS)",
  r2>1&&r1<1,       "CEP strain 2 (E2 LAS)",
  r1>r2>1&&r21<1,   "CEP strain 1 dominant (E1 LAS)",
  r2>r1>1&&r12<1,   "CEP strain 2 dominant (E2 LAS)",
  r21>1&&r12>1,     "Persistence + coexistence (E12 exists)",
  True,             "boundary/equality case"];
{r01n,r02n,r21n,r12n}=N@({R01,R02,R21,R12}/.coP);
Print["coP: R1(E0)=",r01n,"  R2(E0)=",r02n,"  R2(E1)=",r21n,"  R1(E2)=",r12n];
Print["Regime: ",relayClass[r01n,r02n,r21n,r12n]];

(* ── Relay coefficient b at E1 (Sigma2={i2,y2} invades, unfolding param be2) ── *)
(* var={s,i1,y1,i2,y2,r1,r2,R}; resI1 positions: 1=s 2=i1 3=y1 4=r1 5=r2 6=R *)
resI1={1,2,3,6,7,8}; invI1={4,5};
Jyy1=FullSimplify[jac[[resI1,resI1]]/.Thread[var->E1]];
Jyx1=FullSimplify[jac[[resI1,invI1]]/.Thread[var->E1]];
(* Perron pair of M2 at E1, normalised ell.r=1 (ho.tex Rem. r:relay-E1) *)
q2e=s1e+et2*si2*r1st;            (* = v2/be2 at criticality R21=1 *)
rr2={s1e,si2*r1st};  ll2={1,et2}/q2e;
Print["Perron ll2.rr2 (=1): ",FullSimplify[ll2.rr2]];
(* Crossing coefficient a = ll2.(d_{be2}M2).rr2 = q2e > 0 *)
aE1=FullSimplify[ll2.(Outer[Times,{s1e,si2*r1st},{1,et2}].rr2)];
Print["Crossing coeff a at E1 (should equal q2e>0): ",aE1//prF,
      "  match: ",FullSimplify[aE1-q2e]===0];
(* eta = -Jyy1^{-1} Jyx1.rr2;   D_x M2 = 0 (M2 depends only on s,r1) *)
eta2=FullSimplify[-LinearSolve[Jyy1,Jyx1.rr2]];
(* b = be2*(eta_s + et2*si2*eta_{r1}); only s(pos 1) and r1(pos 4) enter via DyM2 *)
bE1=FullSimplify[be2*(eta2[[1]]+et2*si2*eta2[[4]]),cp];
Print["Branch coeff b at E1: ",bE1//prF];
Print["b<0 (forward relay E1->E12): ",FullSimplify[bE1<0,cp]];

(* ── Relay coefficient b at E2 (Sigma1={i1,y1} invades, unfolding param be1) ── *)
(* resI2 positions: 1=s 2=i2 3=y2 4=r1 5=r2 6=R *)
resI2={1,4,5,6,7,8}; invI2={2,3};
Jyy2=FullSimplify[jac[[resI2,resI2]]/.Thread[var->E2]];
Jyx2=FullSimplify[jac[[resI2,invI2]]/.Thread[var->E2]];
q1e=s2e+et1*si1*r2st;
rr1={s2e,si1*r2st};  ll1={1,et1}/q1e;
Print["Perron ll1.rr1 (=1): ",FullSimplify[ll1.rr1]];
aE2=FullSimplify[ll1.(Outer[Times,{s2e,si1*r2st},{1,et1}].rr1)];
Print["Crossing coeff a at E2 (should equal q1e>0): ",aE2//prF,
      "  match: ",FullSimplify[aE2-q1e]===0];
eta1=FullSimplify[-LinearSolve[Jyy2,Jyx2.rr1]];
(* b = be1*(eta_s + et1*si1*eta_{r2}); only s(pos 1) and r2(pos 5) enter via DyM1 *)
bE2=FullSimplify[be1*(eta1[[1]]+et1*si1*eta1[[5]]),cp];
Print["Branch coeff b at E2: ",bE2//prF];
Print["b<0 (forward relay E2->E12): ",FullSimplify[bE2<0,cp]];

(* Numerical verification at coP *)
Print["b at coP: bE1=",N[bE1/.coP],"  bE2=",N[bE2/.coP],
      "  both<0: ",N[bE1/.coP]<0&&N[bE2/.coP]<0];

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



gridRes=80;plotInd={1,2};steTol=10^(-8);staTol=10^(-10);choTol=10^(-13);rangeExt=1.5;
(*tMax=300;nIc=8;*)
Timing[(*Step 3:Scanning-now with mSi from step 1*)
{plot,noSol,results}=scan[RHS,var,par,coP,plotInd,mSi,(*mSi auto-provided*)
gridRes,steTol,staTol,choTol,R01,R02,R12,R21,rangeExt];]
plot
Export["GavSca.pdf",plot]


(* MAXIMAL ELIMINABLE SETS *)
{maxId, maxSets,qss} = maxElim[RN, rts, var];
Print["reduced models are",maxSets];
Print["qss1=",(FullSimplify/@qss[[1]])//prF//MatrixForm];
Print["qss2=",(FullSimplify/@qss[[2]])//prF//MatrixForm];
Print["qss3=",(FullSimplify/@qss[[3]])//prF//MatrixForm];
(*{co,nc}=findCores[RN,3];co*)
