(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;

(* 3-strain SIR: scalar i1 + rank-one {i2,i3} block.
   Strains 2,3 share foi = pi2*i2 + pi3*i3 (rank-one NGM block). *)
var = {s, i1, i2, i3};
foi = pi2*i2 + pi3*i3;
RHS = {
  lam - mu*s - be1*s*i1 - s*foi,
  be1*s*i1 - v1*i1,
  s*w2*foi - v2*i2,
  s*(1 - w2)*foi - v3*i3
};
toDisp = {
  lam   -> \[Lambda],
  mu    -> \[Mu],
  be1   -> Subscript[\[Beta],  1],
  v1    -> Subscript[\[Nu],    1],
  v2    -> Subscript[\[Nu],    2],
  v3    -> Subscript[\[Nu],    3],
  pi2   -> Subscript[\[Pi],    2],
  pi3   -> Subscript[\[Pi],    3],
  w2    -> Subscript[\[Omega], 2],
  foi   -> \[Phi],
  rho23 -> Subscript[\[Rho],  23],
  s0    -> Subscript[s,         0],
  s1    -> Subscript[s,         1],
  s23   -> Subscript[s,        23],
  i1    -> Subscript[i,         1],
  i2    -> Subscript[i,         2],
  i3    -> Subscript[i,         3],
  i1st  -> SuperStar[Subscript[i, 1]],
  i2st  -> SuperStar[Subscript[i, 2]],
  i3st  -> SuperStar[Subscript[i, 3]],
  R1    -> Subscript[R,  1],
  R23   -> Subscript[R, 23],
  K23   -> Subscript[K, 23]
};
prF[ex_] := ex /. toDisp;
Column[MapThread[TraditionalForm[Row[{#1, "' = ", #2}]] &,
  {var /. toDisp, RHS /. toDisp}]]



(* ::Section:: *)
(*CRN decomposition*)


{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var];
Print["Reactions (", Length[RN], ") verif gam.rts-RHS=", (gam . rts - RHS) // Simplify];



(* ::Section:: *)
(*NGM and reproduction numbers*)


{RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, sipVar} =
  NGMRN[RN, rts, var];
s0    = lam/mu;
R1    = K[[1, 1]] // FullSimplify;
K23   = K[[2 ;; 3, 2 ;; 3]] // FullSimplify;
R23   = Tr[K23]/.s->s0 // FullSimplify;
rho23 = w2*pi2/v2 + (1 - w2)*pi3/v3;
Print["K=", K // MatrixForm, " R1=", R1 // prF,
      " K23=", K23 // MatrixForm, " R23=", R23 // prF,
      " R23=s0*rho23?: ", FullSimplify[s0*rho23 - R23] ];



(* ::Section:: *)
(*Minimal siphons*)


mSI = {{"i1"}, {"i2", "i3"}};
Print["Minimal siphons=", mSI];



(* ::Section:: *)
(*Boundary equilibria*)


cE0 = {s -> lam/mu, i1 -> 0, i2 -> 0, i3 -> 0};
s1   = v1/be1;
i1st = FullSimplify[mu*(R1 - 1)/be1];
cE1  = {s -> s1, i1 -> i1st, i2 -> 0, i3 -> 0};
s23   = FullSimplify[1/rho23];
Amat  = {{s23*w2*pi2 - v2, s23*w2*pi3}, {s23*(1 - w2)*pi2, s23*(1 - w2)*pi3 - v3}};
{a2, a3} = NullSpace[Amat // FullSimplify][[1]];
foiVal  = FullSimplify[(lam - mu*s23)/s23];
normFoi = FullSimplify[a2*pi2 + a3*pi3];
i2st    = FullSimplify[foiVal*a2/normFoi];
i3st    = FullSimplify[foiVal*a3/normFoi];
cE23    = {s -> s23, i1 -> 0, i2 -> i2st, i3 -> i3st};
Print["E0 verif=", (RHS /. cE0) // FullSimplify,
      "\nE1=", cE1 // prF, " verif=", (RHS /. cE1) // FullSimplify, " exists iff R1>1",
      "\nE23: s23=", s23 // prF, " verif=", (RHS /. cE23) // FullSimplify, " exists iff R23>1"];



(* ::Section:: *)
(*Invasion numbers and CEP relay table*)


R1atE23 = FullSimplify[be1*s23/v1];
R23atE1 = FullSimplify[s1*rho23];
Print["R1(E0)=", R1 // prF, " R23(E0)=", R23 // prF,
      "\nR1(E23)=", R1atE23 // prF, " =R1/R23: ", FullSimplify[R1atE23 - R1/R23] === 0,
      "\nR23(E1)=", R23atE1 // prF, " =R23/R1: ", FullSimplify[R23atE1 - R23/R1] === 0,
      "\nCEP: ", Grid[{
        {"Region", "Conditions", "LAS"},
        {"\[CapitalOmega]0",  "R1<1, R23<1",  "E0"},
        {"\[CapitalOmega]1",  "R1>R23, R1>1", "E1"},
        {"\[CapitalOmega]23", "R23>R1, R23>1","E23"}
      }, Frame -> All]];
