(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

(* 2-patch SEIR graph-diffusion test.
   Task 4: verify R0^patch = R0^local for homogeneous coefficients.
   Both I and E have d=0 (no diffusion); result is trivially R0^patch = R0^local.
   Then repeat with I diffusing (d > 0); same equality must hold. *)

(* ── Local GraphDiffuseRN (self-contained copy) ─────────────────────────── *)
GraphDiffuseRN[RHS_, var_, diff_, L_] := Module[
  {nS, nP, diffA, pV, sub, varPatch, rhsPatch},
  nS = Length[var]; nP = Length[L];
  diffA = If[AssociationQ[diff], diff, Association[Rule @@@ diff]];
  pV[a_, p_] := ToExpression[ToString[var[[a]]] <> ToString[p]];
  sub[p_]    := Table[var[[a]] -> pV[a, p], {a, nS}];
  varPatch   = Flatten[Table[pV[a, p], {p, nP}, {a, nS}]];
  rhsPatch   = Flatten[Table[
    (RHS[[a]] /. sub[p]) +
    Lookup[diffA, var[[a]], 0] * Sum[L[[p, q]] * pV[a, q], {q, nP}],
    {p, nP}, {a, nS}]];
  {varPatch, rhsPatch}
];

(* ── prF ────────────────────────────────────────────────────────────────── *)
prF[ex_] := ex /. {
  lam -> \[Lambda], mu -> \[Mu], be -> \[Beta], ga -> \[Gamma], si -> \[Sigma],
  s1 -> Subscript[S,1], e1 -> Subscript[E,1], i1 -> Subscript[I,1], r1 -> Subscript[R,1],
  s2 -> Subscript[S,2], e2 -> Subscript[E,2], i2 -> Subscript[I,2], r2 -> Subscript[R,2]
};

(* ── Local SEIR skeleton ─────────────────────────────────────────────────── *)
var = {s, e, i, r};
RHS = {
  lam - mu*s - be*s*i,
  be*s*i - (mu + si)*e,
  si*e - (mu + ga)*i,
  ga*i - mu*r
};
(* Local R0 = be*si / ((mu+si)*(mu+ga)) at DFE s0 = lam/mu *)
R0local = be*(lam/mu)*si / ((mu + si)*(mu + ga));

(* ── Test A: no diffusion (trivial) ─────────────────────────────────────── *)
L2 = {{-1, 1}, {1, -1}};
diffNone = <||>;
{vP0, RP0} = GraphDiffuseRN[RHS, var, diffNone, L2];
{RN0, rts0, spe0, alp0, bet0, gam0, rnR0} = ODE2RNs[RP0, vP0, prF];
Print["Test A verif: ", RP0 - gam0 . rts0 // FullSimplify];
{RHS2a, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} =
    NGMRN[RN0, rts0, vP0];
R0pA = Max[R0A] // FullSimplify;
Print["Test A (no diff): R0^patch = ", R0pA // prF, "  R0^local = ", R0local // prF,
      "  diff = ", FullSimplify[R0pA - R0local]];

(* ── Test B: I diffuses ─────────────────────────────────────────────────── *)
diffI = <|i -> d|>;
{vP1, RP1} = GraphDiffuseRN[RHS, var, diffI, L2];
{RN1, rts1, spe1, alp1, bet1, gam1, rnR1} = ODE2RNs[RP1, vP1, prF];
Print["Test B verif: ", RP1 - gam1 . rts1 // FullSimplify];
{RHS2b, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} =
    NGMRN[RN1, rts1, vP1];
R0pB = Max[R0A] // FullSimplify;
Print["Test B (I diff): R0^patch = ", R0pB // prF, "  R0^local = ", R0local // prF,
      "  diff = ", FullSimplify[R0pB - R0local]];

(* ── Siphons of 2-patch SEIR with I-diffusion ────────────────────────────── *)
allSi = minSiph[vP1, RN1];
Print["minSiphs = ", allSi // prF];
