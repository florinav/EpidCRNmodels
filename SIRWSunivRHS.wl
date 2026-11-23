(* ::Package:: *)

ClearAll["Global`*"]
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]]; <<EpidCRN`;

(* ===== MODEL-SPECIFIC: change var, RHS, prF, yvar, xvar, pars for each model ===== *)
var  = {S, i, R, W};
RHS  = {La + gaw*W - be*S*i - mus*S,
        i*(be*S - gai - mui),
        gai*i + nu*i*W - (mur + gar)*R,
        gar*R - nu*i*W - gaw*W - muw*W};
prF[ex_] := ex /. {be -> \[Beta], La -> \[CapitalLambda], nu -> \[Nu],
   mus -> Subscript[\[Mu], s], mur -> Subscript[\[Mu], r],
   muw -> Subscript[\[Mu], w], mui -> Subscript[\[Mu], i],
   gai -> Subscript[\[Gamma], i], gar -> Subscript[\[Gamma], r],
   gaw -> Subscript[\[Gamma], w]};
yvar = {S, R, W};   (* susceptible-like classes *)
xvar = {i};         (* infected classes *)
pars = {La->1, be->1/2, mus->1/10, mui->1/10, gai->1/5,
        mur->1/10, gar->1/10, gaw->3/10, muw->1/10, nu->1/2};
(* ===== END MODEL-SPECIFIC ===== *)

(* CRN extraction *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var, prF];
Print["RHS check gam.rts-RHS=", FullSimplify[gam . rts - RHS]];

(* Siphons, NGM, R0 *)
{RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} = NGMRN[RN, rts, var];
Print["mSi=", mSi, " cDFE=", cDFE, " R0A=", R0A // prF];
Print["K=", K // MatrixForm // prF, " cE0=", Short[cE0 // prF]];

(* ABa decomposition *)
{LaE, ASE, vbeE, baE, AE} = extABa[RHS, var, yvar, xvar];

(* Four Laws *)
{y0, mR, R0, Dw, ye, k, xe} = FLaws[AE, vbeE, baE, ASE, LaE];

(* General P,B decomposition *)
{LaP, ASP, BP, AP, PP, yvarP} = extABP[RHS, var, xvar];

(* Oscillation analysis *)
{n, m} = {Length[var], Length[RN]};
RM       = symRM[alp, n, m];
co       = Rest@Reverse@CoefficientList[Expand@CharacteristicPolynomial[gam . RM, lam], lam];
nCS      = negCS[co, gam, n];
upf      = Select[nCS, upfQ];
{gH, edges} = hasseCS[nCS];
Print["#nCS=", Length[nCS], " #upf=", Length[upf], " #edges=", Length[edges]];
upfRes   = upfCS[nCS, gH, spe, rts];
riCores  = oscRI[nCS, gH, spe, rts];
riiCores = oscRII[nCS, gH, spe, rts];
Print["riCores=", riCores, " riiCores=", riiCores];

(* Hopf analysis *)
csN   = Select[nCS, #[[1]] == n &];
stabN = Select[csN, hurQ[#[[5]]] &];
fs    = fixSgn[nCS, n];
Print["fixSgn=", fs, " Cor.5.21 Hopf=", Length[stabN] > 0 && Length[upf] > 0];
{coH, resH} = resHopf[Grad[RHS, var]];
Print["Hopf coeff=", coH // prF, " Hopf resultant=", Short[resH // prF]];

(* Numerical checks for Four Laws *)
{y0n, mRn, R0n, Dwn, yen, kn, xen} = N[{y0, mR, R0, Dw, ye, k, xe} /. pars];
Print["R0=", R0n, " Law1 y0.mR-R0=", Chop[y0n . mRn - R0n], " Law4 ye.mR-1=", Chop[yen . mRn - 1]];
Print["Law2 xe-k*Dw=", Chop[xen - kn*Dwn], " RHS at endemic eq=",
      Chop[N[FullSimplify[RHS /. Thread[var -> Flatten[{ye, xe}]]] /. pars]]];
