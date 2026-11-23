(* ::Package:: *)

(* CEP.wl — Competitive Exclusion Principle: siphon-lattice decision algorithm.
   See alg:CEP (alg.tex), Theorem t:CEP (Lsif.tex).

   covers format: {{parentNode, childNode, invaderBlock, faceR-expression}, ...}
   Threshold for cover k: T_k = covers[[k,4]] - 1  (i.e. R_sigma - 1).
   All routines return lists; never bare atoms.

   Exported: nStbQ, hotCh, cepDec, cepDecX, cepDecSym

   cepDec    — numerical, strict chain walk (Theorem t:CEP / alg:CEP numerical mode)
   cepDecX   — numerical, BFS extension (non-strict / coexistence case)
   cepDecSym — symbolic, sign-pattern enumeration (alg:CEP symbolic mode)
*)

ClearAll[nStbQ, hotCh, cepLat, cepCovers, cepDec, cepDecX, cepDecSym];

(* --- Lattice construction ---------------------------------------------- *)

(* cepLat[mSi]: cover skeleton of the siphon lattice.
   mSi: list of minimal siphons, each a sorted list of species strings.
   Returns: {{parentSiphon, childSiphon, invaderBlock}, ...}
   parentSiphon = union of some minimal siphons (sorted);
   invaderBlock = one minimal siphon si ⊆ parentSiphon;
   childSiphon  = parentSiphon \ invaderBlock.
   One cover per (node, minimal-siphon) pair with si ⊆ node.          *)
cepLat[mSi_] :=
  Module[{nodes, skeleton = {}},
    nodes = Sort /@ Union[Sort /@ (Union @@ # & /@ Subsets[mSi])];
    Do[
      Do[
        If[SubsetQ[nd, si],
          AppendTo[skeleton, {nd, Sort[Complement[nd, si]], si}]
        ],
        {si, mSi}
      ],
      {nd, nodes}
    ];
    skeleton
  ];

(* --- Cover assembly ---------------------------------------------------- *)

(* cepCovers[mSi, RHS, var]
   Builds the complete covers list from model data alone; no model-specific
   code needed beyond mSi, RHS, var (all available from NGMRN).
   For each cover (parentSip, childSip, invader) from cepLat[mSi]:
     1. Face equilibrium: solve RHS=0 with parentSip species = 0.
     2. Transversal Jacobian block: differentiate RHS[invader] w.r.t.
        invader at invader=0 (symbolic in susceptible variables).
     3. Regular splitting via posM: Fsig = posM[Jsig], Vsig = posM[-Jsig].
     4. R_sigma = rho(Fsig . Inverse[Vsig]) evaluated at face equilibrium.
   Node identifiers are sorted siphon lists (absent-species sets).
   Returns: {{parentSip, childSip, invaderStr, R_sigma}, ...}              *)
cepCovers[mSi_, RHS_, var_] :=
  Module[{skel = cepLat[mSi]},
    Table[
      Module[{parentSip, childSip, invader, absSym, invSyms, invIdx,
              JsigSym, FsigSym, VsigSym, faceEq, Fsig, Vsig, Rsig},
        {parentSip, childSip, invader} = cv;
        absSym  = ToExpression /@ parentSip;
        invSyms = ToExpression /@ invader;
        invIdx  = Flatten[Position[var, #] & /@ invSyms];
        (* Transversal Jacobian block (symbolic in susceptible vars) *)
        JsigSym = Outer[D, RHS[[invIdx]], invSyms] /. Thread[invSyms -> 0];
        (* Regular splitting: F = positive part, V = positive part of -J *)
        FsigSym = posM[ JsigSym];
        VsigSym = posM[-JsigSym];
        (* Face equilibrium: set parentSip species to 0, solve for others *)
        faceEq = First @ Solve[
          Thread[(RHS /. Thread[absSym -> 0]) == 0],
          Complement[var, absSym]];
        Fsig = FsigSym /. faceEq;
        Vsig = VsigSym /. faceEq;
        (* Face reproduction number *)
        Rsig = FullSimplify[Max @ Eigenvalues[Fsig . Inverse[Vsig]]];
        {parentSip, childSip, StringJoin[Sort[invader]], Rsig}
      ],
      {cv, skel}
    ]
  ];

Options[cepDec]  = {tol -> 10^-10, strict -> True};
Options[cepDecX] = {tol -> 10^-10};

(* --- Shared primitives ------------------------------------------------- *)

(* nStbQ[nd, cv, ps]  — True iff all outgoing T_k = R_k - 1 < -tol.
   Vacuously True at bottom nodes (no outgoing covers).               *)
nStbQ[nd_, cv_, ps_, tol_:10^-10] :=
  AllTrue[
    N[(#[[4]] /. ps) - 1] & /@ Select[cv, #[[1]] === nd &],
    # < -tol &];

(* hotCh[nd, cv, ps]  — outgoing covers with T_k = R_k - 1 > tol.   *)
hotCh[nd_, cv_, ps_, tol_:10^-10] :=
  Select[cv, #[[1]] === nd && N[(#[[4]] /. ps) - 1] > tol &];

(* --- Numerical: strict chain walk (alg:CEP numerical mode) ------------ *)
(*
   Returns one of (alg.tex output format):
     {"stable",       nd}            — nd is the unique LAS node
     {"threshold",    nd}            — no hot child and nd not stable
                                       (T within tol of 0)
     {"cepViolation", nd, blocks}    — multiple hot children at nd;
                                       blocks = list of invader block names
   strict->False skips the violation check (first hot child taken).   *)
cepDec[cv_, top_, ps_, OptionsPattern[]] :=
  Module[{nd = top, hot, tl = OptionValue[tol], st = OptionValue[strict]},
    While[True,
      If[nStbQ[nd, cv, ps, tl], Return[{"stable", nd}]];
      hot = hotCh[nd, cv, ps, tl];
      If[hot === {}, Return[{"threshold", nd}]];
      If[st && Length[hot] > 1,
        Return[{"cepViolation", nd, hot[[All, 3]]}]];
      nd = hot[[1, 2]]
    ]
  ];

(* --- Numerical: BFS extension for non-strict / coexistence case ------- *)
(*
   Returns one of:
     {"stable", nd}    — nd is the unique LAS node
     {"none",   seen}  — no stable node reached (diagnostic; should not
                         occur under the CEP hypotheses)               *)
cepDecX[cv_, top_, ps_, OptionsPattern[]] :=
  Module[{q = {top}, seen = {}, nd, hot, tl = OptionValue[tol]},
    While[q =!= {},
      nd = First[q]; q = Rest[q];
      If[MemberQ[seen, nd], Continue[]];
      seen = Join[seen, {nd}];
      If[nStbQ[nd, cv, ps, tl], Return[{"stable", nd}]];
      hot = hotCh[nd, cv, ps, tl][[All, 2]];
      q = Join[q, hot]
    ];
    {"none", seen}
  ];

(* --- Symbolic: sign-pattern relay walk (internal helper) --------------- *)
(*
   cepDecPat[cv, top, pat]
   pat: integer list, one entry per cover, +1 = hot (T>0), -1 = cold (T<0).
   Runs the strict relay walk on the sign pattern.
   Returns {"stable", nd}, {"threshold", nd}, or {"cepViolation", nd, blocks}.
*)
cepDecPat[cv_, top_, pat_] :=
  Module[{nd = top, outIdx, hotIdx},
    While[True,
      outIdx = Flatten[Position[cv[[All, 1]], nd]];
      If[AllTrue[pat[[outIdx]], # < 0 &], Return[{"stable", nd}]];
      hotIdx = Select[outIdx, pat[[#]] > 0 &];
      If[hotIdx === {}, Return[{"threshold", nd}]];
      If[Length[hotIdx] > 1,
        Return[{"cepViolation", nd, cv[[hotIdx, 3]]}]];
      nd = cv[[hotIdx[[1]], 2]]
    ]
  ];

(* --- Symbolic: sign-pattern enumeration (alg:CEP symbolic mode) -------- *)
(*
   cepDecSym[cv, top]
   Enumerates all 2^n sign patterns for the n cover thresholds T_k = R_k - 1.
   For each pattern yielding {"stable", nd}, records the region:
     {"region", conditions, nd}
   where conditions is a conjunction of symbolic inequalities T_k > 0 or < 0.
   Patterns yielding "cepViolation" or "threshold" are omitted; the caller
   must handle those cases separately (e.g. via cepDecX).
   Multiple patterns with the same LAS node appear as separate entries;
   simplification / merging of conditions is left to the caller.           *)
cepDecSym[cv_, top_] :=
  Module[{n = Length[cv], thrs, results = {}},
    thrs = cv[[All, 4]] - 1;
    Do[
      Module[{res = cepDecPat[cv, top, pat], conds},
        If[res[[1]] === "stable",
          conds = And @@ MapThread[If[#1 > 0, #2 > 0, #2 < 0] &, {pat, thrs}];
          AppendTo[results, {"region", conds, res[[2]]}]
        ]
      ],
      {pat, Tuples[{-1, 1}, n]}
    ];
    results
  ];
