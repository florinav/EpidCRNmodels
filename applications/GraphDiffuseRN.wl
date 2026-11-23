(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;

(* GraphDiffuseRN -- extend a local ODE RHS to a multi-patch graph-diffusion ODE.
   Convention: variable for species a in patch p = ToExpression[ToString[var[[a]]] <> ToString[p]].
   Example: var={s,i}, patch 1 -> s1,i1; patch 2 -> s2,i2.
   Demonstration: 2-patch SIS; Task 4 constant-mode test: R0^patch = R0^local. *)

(* ────────────────────────────────────────────────────────────────────────
   {varPatch, RHSPatch} = GraphDiffuseRN[RHS, var, diff, L]
   RHS      : local RHS vector, length nS (one-patch ODE skeleton)
   var      : local variable list, length nS
   diff     : Association {species -> diffusion_coefficient}; omitted species get d=0
   L        : nP x nP graph Laplacian (row sums = 0, diagonal <= 0)
   varPatch : flattened patch variables ordered (patch1 all species, patch2, ...)
   RHSPatch : patch ODE RHS incorporating local kinetics + graph-diffusion
   ──────────────────────────────────────────────────────────────────────── *)
GraphDiffuseRN[RHS_, var_, diff_, L_] := Module[
  {nS, nP, diffA, pV, sub, varPatch, rhsPatch},
  nS = Length[var];
  nP = Length[L];
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

(* ── prF for 2-patch SIS ────────────────────────────────────────────────── *)
prF[ex_] := ex /. {
  be -> \[Beta], ga -> \[Gamma], mu -> \[Mu], d -> Subscript[d, i],
  s1 -> Subscript[S, 1], i1 -> Subscript[I, 1],
  s2 -> Subscript[S, 2], i2 -> Subscript[I, 2]
};

(* ── Local (one-patch) SIS skeleton ────────────────────────────────────── *)
var = {s, i};
RHS = {
  mu - mu*s - be*s*i,
  be*s*i - (mu + ga)*i
};
(* Local R0 = be/(mu+ga) *)
R0local = be/(mu + ga);

(* ── 2-patch setup ──────────────────────────────────────────────────────── *)
(* L for 2 nodes: each entry multiplied by d in diff Association *)
L2 = {{-1, 1}, {1, -1}};   (* graph Laplacian for K_2; actual diffusion rate in diff *)
diff2 = <|i -> d|>;         (* only I diffuses *)

{varPatch, RHSPatch} = GraphDiffuseRN[RHS, var, diff2, L2];
Print["varPatch = ", varPatch, "  RHSPatch = ", RHSPatch // prF];

(* ── CRN decomposition ──────────────────────────────────────────────────── *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNs[RHSPatch, varPatch, prF];
Print["verif RHS = ", RHSPatch - gam . rts // FullSimplify];

(* ── NGM of patch model ─────────────────────────────────────────────────── *)
{RHS2, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} =
    NGMRN[RN, rts, varPatch];
Print["E0 = ", cE0 // prF, "  infVars = ", infV, "  R0 = ", R0A // prF];
Print["K = ", K // prF // MatrixForm];

(* ── Task 4: constant-mode test ─────────────────────────────────────────── *)
(* For homogeneous coefficients the spectral radius of K_patch equals R0_local.
   NGMRN returns spectral radius as the first (largest) entry of R0A. *)
R0patch = R0A[[1]] // FullSimplify;
diff4   = FullSimplify[R0patch - R0local];
Print["Task 4 constant-mode test:  R0^patch - R0^local = ", diff4 // prF,
      "  (expected 0)"];
