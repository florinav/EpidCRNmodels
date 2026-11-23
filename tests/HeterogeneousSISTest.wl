(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;

(* Task 5: Heterogeneous 2-patch SIS counterexample.
   Demonstrates that R0^patch != R0^local,1 and R0^patch != R0^local,2
   when beta varies across patches.
   Uses PatchNGM directly (Kronecker-product construction).
   Symbolic result + numeric example. *)

(* ── Local PatchNGM (self-contained) ───────────────────────────────────── *)
PatchNGM[FF_, VV_, infVarsLocal_, diff_, L_] := Module[
  {nI, nP, diffA, dVec, Dinf, FP, VP, KP},
  nI    = Length[infVarsLocal];  nP = Length[L];
  diffA = If[AssociationQ[diff], diff, Association[Rule @@@ diff]];
  dVec  = Table[Lookup[diffA, infVarsLocal[[i]], 0], {i, nI}];
  Dinf  = DiagonalMatrix[dVec];
  FP    = KroneckerProduct[IdentityMatrix[nP], FF];
  VP    = KroneckerProduct[IdentityMatrix[nP], VV] - KroneckerProduct[L, Dinf];
  KP    = Simplify[FP . Inverse[VP]];
  {FP, VP, KP}
];

(* ── prF ────────────────────────────────────────────────────────────────── *)
prF[ex_] := ex /. {
  mu -> \[Mu], ga -> \[Gamma], d -> Subscript[d, i],
  be1 -> Subscript[\[Beta], 1], be2 -> Subscript[\[Beta], 2]
};

(* ── Symbolic 2-patch heterogeneous SIS ──────────────────────────────────── *)
(* At DFE (s_p* = 1 for each patch):
   F_patch = DiagonalMatrix[{be1, be2}]  (heterogeneous beta)
   V_patch local = {{mu+ga}}  (same in both patches)
   I diffuses with rate d between patches. *)
FFhet = DiagonalMatrix[{be1, be2}];   (* 2x2 because 2 patches, 1 infectious var each *)
VV1   = {{mu + ga}};                   (* 1x1 local V *)
infV1 = {i};
diff1 = <|i -> d|>;
L2    = {{-1, 1}, {1, -1}};

{FP2, VP2, KP2} = PatchNGM[{{be1}}, VV1, infV1, diff1, L2];
(* Note: PatchNGM builds F_patch = I_2 (x) {{be1}} which gives diag(be1,be1).
   For heterogeneous beta we must supply the full F_patch manually. *)
(* Correct heterogeneous F_patch: *)
FPhet = FFhet;   (* diag(be1,be2) *)
VPhet = KroneckerProduct[IdentityMatrix[2], VV1] - KroneckerProduct[L2, {{d}}];
(* VPhet = {{mu+ga+d, -d}, {-d, mu+ga+d}} *)
KPhet = Simplify[FPhet . Inverse[VPhet]];
Print["Heterogeneous KP = ", KPhet // prF // MatrixForm];

(* Spectral radius (characteristic polynomial of 2x2) *)
{ev1, ev2} = Eigenvalues[KPhet] // FullSimplify;
R0patch = Max[ev1, ev2] // FullSimplify;
Print["R0^patch (symbolic) = ", R0patch // prF];

(* Local R0 values *)
R0p1 = be1/(mu + ga);
R0p2 = be2/(mu + ga);
Print["R0^local,1 = ", R0p1 // prF, "  R0^local,2 = ", R0p2 // prF];

(* Symbolic inequality: R0^patch strictly between min and max local R0 *)
diff12 = FullSimplify[R0patch - be1/(mu + ga), {be1 > 0, be2 > 0, mu > 0, ga > 0, d > 0, be1 > be2}];
Print["R0^patch - R0^local,1 (should be < 0 when be1 > be2) = ", diff12 // prF];

(* ── Numerical example ──────────────────────────────────────────────────── *)
(* Parameters: mu=0.01, ga=0.1, be1=0.5, be2=0.3, d=0.05 *)
numRules = {mu -> 1/100, ga -> 1/10, be1 -> 1/2, be2 -> 3/10, d -> 1/20};
KPnum    = KPhet /. numRules // N;
R0num    = Max[Eigenvalues[KPnum]];
R0p1num  = (be1/(mu + ga)) /. numRules // N;
R0p2num  = (be2/(mu + ga)) /. numRules // N;
Print["Numerical example:  R0^patch = ", R0num,
      "  R0^local,1 = ", R0p1num, "  R0^local,2 = ", R0p2num];
Print["Check:  R0^local,2 < R0^patch < R0^local,1 ? ",
      R0p2num < R0num < R0p1num];

(* ── Limiting cases ─────────────────────────────────────────────────────── *)
(* d -> 0: R0^patch -> max(be1,be2)/(mu+ga) *)
R0dTo0 = Limit[R0patch, d -> 0] // FullSimplify;
Print["d->0:  R0^patch -> ", R0dTo0 // prF,
      "  (expected max of R0^local,1 and R0^local,2)"];

(* d -> Inf: R0^patch -> (be1+be2)/(2*(mu+ga)) (arithmetic mean) *)
R0dInf = Limit[R0patch, d -> Infinity] // FullSimplify;
Print["d->Inf:  R0^patch -> ", R0dInf // prF,
      "  (expected arithmetic mean of local R0s)"];
