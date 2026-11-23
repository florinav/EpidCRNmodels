(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;

(* NeumannLaplacian1D -- finite-difference Laplacian on a 1D uniform grid
   with homogeneous Neumann (zero-flux) boundary conditions.
   This is the standard metapopulation discretization of the 1D diffusion PDE
   u_t = d * u_xx on [0,L] with u_x(0)=u_x(L)=0. *)

(* ────────────────────────────────────────────────────────────────────────
   L = NeumannLaplacian1D[n, h]
   n : number of grid points (patches), n >= 2
   h : grid spacing (physical distance between adjacent patches)
   L : n x n dense matrix; row sums = 0; negative semi-definite
       Interior row j (1 < j < n): [0,...,1/h^2, -2/h^2, 1/h^2,...,0]
       Boundary row 1:              [-1/h^2,  1/h^2, 0, ..., 0]
       Boundary row n:              [0, ..., 0,  1/h^2, -1/h^2]
   ──────────────────────────────────────────────────────────────────────── *)
NeumannLaplacian1D[n_, h_] := Module[{entries},
  entries = Flatten[{
    (* left boundary *)
    {{1, 1} -> -1/h^2, {1, 2} -> 1/h^2},
    (* interior *)
    Table[{{j, j-1} -> 1/h^2, {j, j} -> -2/h^2, {j, j+1} -> 1/h^2}, {j, 2, n-1}],
    (* right boundary *)
    {{n, n-1} -> 1/h^2, {n, n} -> -1/h^2}
  }, 1];
  Normal[SparseArray[entries, {n, n}]]
];

(* ── Demonstration: 3-patch Neumann Laplacian ──────────────────────────── *)
L3 = NeumannLaplacian1D[3, h];
Print["NeumannLaplacian1D[3,h] = ", L3 // MatrixForm];
Print["Row sums (all 0): ", Total[L3, {2}]];

(* ── Eigenvalues (should be 0, -1/h^2, -3/h^2 for n=3) ────────────────── *)
evals3 = Eigenvalues[L3] // FullSimplify;
Print["Eigenvalues[L3] = ", evals3];

(* ── 4-patch example ────────────────────────────────────────────────────── *)
L4 = NeumannLaplacian1D[4, h];
Print["NeumannLaplacian1D[4,h] = ", L4 // MatrixForm];
Print["Eigenvalues[L4] = ", Eigenvalues[L4] // FullSimplify];

(* ── Relationship to Neumann eigenvalues ────────────────────────────────── *)
(* The eigenvalues of NeumannLaplacian1D[n,h] are
   lambda_k = -(2/h^2)(1 - Cos[k*Pi/(n-1)]) ... actually for this discretization
   they match the continuous spectrum -lambda_k of -d/dx^2 with Neumann BC. *)
Print["Continuous Neumann eigenvalues for n=3: ",
      Table[-(2/h^2)*(1 - Cos[k*Pi/2]), {k, 0, 2}] // FullSimplify];
