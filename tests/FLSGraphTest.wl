(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

(* FLS 3-patch graph-diffusion test.
   Single-strain FLS (Frog-Chytrid) / cholera-type model with B diffusing.
   Task 4: verify R0^patch = R0^local for homogeneous coefficients.
   Uses NeumannLaplacian1D[3, h] for a 3-patch 1D grid. *)

(* ── Local helpers (self-contained) ────────────────────────────────────── *)
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

NeumannLaplacian1D[n_, hh_] := Module[{entries},
  entries = Flatten[{
    {{1,1} -> -1/hh^2, {1,2} -> 1/hh^2},
    Table[{{j,j-1} -> 1/hh^2, {j,j} -> -2/hh^2, {j,j+1} -> 1/hh^2}, {j, 2, n-1}],
    {{n,n-1} -> 1/hh^2, {n,n} -> -1/hh^2}
  }, 1];
  Normal[SparseArray[entries, {n, n}]]
];

(* ── prF ────────────────────────────────────────────────────────────────── *)
prF[ex_] := ex /. {
  lam -> \[Lambda], mu -> \[Mu], be -> \[Beta], ga -> \[Gamma],
  xi -> \[Xi], de -> \[Delta], dB -> Subscript[d, B]
};

(* ── Local FLS model: S, I, B ───────────────────────────────────────────── *)
(* Environmental (cholera-type): S susceptible hosts, I infected, B pathogen in water *)
var = {s, i, b};
RHS = {
  lam - mu*s - be*s*b,
  be*s*b - (mu + ga)*i,
  xi*i - de*b
};
(* DFE: s0 = lam/mu.  R0 = be*xi*s0 / ((mu+ga)*de) *)
R0local = be*xi*(lam/mu) / ((mu + ga)*de);

(* ── 3-patch Neumann Laplacian (h=1) ───────────────────────────────────── *)
L3 = NeumannLaplacian1D[3, 1];
Print["L3 = ", L3 // MatrixForm, "  row sums = ", Total[L3, {2}]];

(* B diffuses; S, I do not *)
diff3 = <|b -> dB|>;

{varPatch, RHSPatch} = GraphDiffuseRN[RHS, var, diff3, L3];
Print["varPatch = ", varPatch];
Print["verif dim: nS=", Length[var], " nP=3  len(varPatch)=", Length[varPatch]];

(* ── CRN decomposition ──────────────────────────────────────────────────── *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNs[RHSPatch, varPatch, prF];
Print["verif RHS = ", RHSPatch - gam . rts // FullSimplify];

(* ── NGM ────────────────────────────────────────────────────────────────── *)
{RHS2, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} =
    NGMRN[RN, rts, varPatch];
Print["E0 = ", cE0 // prF, "  infVars = ", infV];
Print["R0 from NGMRN = ", R0A // prF];

(* ── Task 4: constant-mode test ─────────────────────────────────────────── *)
R0patch = Max[R0A] // FullSimplify;
test4   = FullSimplify[R0patch - R0local];
Print["Task 4 FLS 3-patch:  R0^patch - R0^local = ", test4 // prF,
      "  (expected 0)"];

(* ── Minimal siphons ────────────────────────────────────────────────────── *)
allSi = minSiph[varPatch, RN];
Print["minSiphs = ", allSi // prF];
