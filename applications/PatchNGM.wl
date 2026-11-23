(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;

(* PatchNGM -- next-generation matrix for a multi-patch epidemic system
   using the Kronecker-product structure of the patch model.
   Assumes homogeneous coefficients (same F,V at every patch) and an arbitrary
   graph Laplacian L for the diffusing infectious variables.
   Reference: Diekmann-Heesterbeek-Metz (1990) formulation for patch models. *)

(* ────────────────────────────────────────────────────────────────────────
   {FP, VP, KP} = PatchNGM[FF, VV, infVarsLocal, diff, L]
   FF           : k x k infection (F) matrix at DFE for one patch
   VV           : k x k transition (V) matrix at DFE for one patch
   infVarsLocal : list of k infectious variable symbols (local, 1-patch)
   diff         : Association {species -> diffusion_coefficient}; d=0 if absent
   L            : nP x nP graph Laplacian (row sums=0, diag <= 0)
   FP           : (nP*k) x (nP*k) patch infection matrix  = I_nP (x) FF
   VP           : (nP*k) x (nP*k) patch transition matrix = I_nP (x) VV - L (x) Dinf
   KP           : FP . Inverse[VP]  (patch next-generation matrix)
   ──────────────────────────────────────────────────────────────────────── *)
PatchNGM[FF_, VV_, infVarsLocal_, diff_, L_] := Module[
  {nI, nP, diffA, dVec, Dinf, FP, VP, KP},
  nI    = Length[infVarsLocal];
  nP    = Length[L];
  diffA = If[AssociationQ[diff], diff, Association[Rule @@@ diff]];
  dVec  = Table[Lookup[diffA, infVarsLocal[[i]], 0], {i, nI}];
  Dinf  = DiagonalMatrix[dVec];
  FP    = KroneckerProduct[IdentityMatrix[nP], FF];
  VP    = KroneckerProduct[IdentityMatrix[nP], VV] - KroneckerProduct[L, Dinf];
  KP    = Simplify[FP . Inverse[VP]];
  {FP, VP, KP}
];

(* ── prF for SIS ────────────────────────────────────────────────────────── *)
prF[ex_] := ex /. {
  be -> \[Beta], ga -> \[Gamma], mu -> \[Mu], d -> Subscript[d, i]
};

(* ── Symbolic 2-patch SIS: Task 4 constant-mode test via PatchNGM ────────── *)
(* Local SIS at DFE (s* = mu/mu = 1):
   F  = {{be}},  V = {{mu+ga}}  (1x1 matrices, single infectious variable i) *)
FF1 = {{be}};
VV1 = {{mu + ga}};
infV1 = {i};
diff1 = <|i -> d|>;
L2    = {{-1, 1}, {1, -1}};   (* 2-patch K_2 Laplacian *)

{FP2, VP2, KP2} = PatchNGM[FF1, VV1, infV1, diff1, L2];
Print["2-patch SIS  FP = ", FP2 // MatrixForm];
Print["2-patch SIS  VP = ", VP2 // prF // MatrixForm];
Print["2-patch SIS  KP = ", KP2 // prF // MatrixForm];

(* Spectral radius of KP2 *)
R0patch2 = Max[Eigenvalues[KP2]] // FullSimplify;
R0local1 = be/(mu + ga);
test4    = FullSimplify[R0patch2 - R0local1];
Print["Task 4:  R0^patch - R0^local = ", test4 // prF, "  (expected 0)"];

(* ── 3-patch cholera/FLS (B diffusing): verify same result ─────────────── *)
(* Local model: S -> I (via B), I -> B shedding, B is infectious variable.
   Infectious variables: {i, b}
   F = {{0, be*s0}, {0, 0}}  (only S+B -> I term is new infection)
   V = {{mu+ga, 0}, {-xi, de}}
   s0 = lam/mu (DFE for S) *)
FF2 = {{0, be*lam/mu}, {0, 0}};
VV2 = {{mu + ga, 0}, {-xi, de}};
infV2 = {i, b};
diff2 = <|b -> dB|>;   (* only B diffuses *)
L3    = {{-1, 1, 0}, {1, -2, 1}, {0, 1, -1}};   (* NeumannLaplacian1D[3,1] *)

{FP3, VP3, KP3} = PatchNGM[FF2, VV2, infV2, diff2, L3];
Print["3-patch FLS  size KP = ", Dimensions[KP3]];

(* Local R0 (rank-1 model) *)
R0localFLS = be*lam*xi / (mu*(mu + ga)*de);
R0patchFLS = Max[Eigenvalues[KP3]] // FullSimplify;
testFLS    = FullSimplify[R0patchFLS - R0localFLS];
Print["Task 4 FLS 3-patch:  R0^patch - R0^local = ", testFLS // prF,
      "  (expected 0)"];
