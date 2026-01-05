(* ::Package:: *)

(* ::Section:: *)
(*Maximal linear-eliminable reduction (FSW-style) for mass-action CRNs*)


(* Requires: extSpe, comp2Asso from Core.wl *)

ClearAll[
  ComplexAssoc, AssocToExpr, RxnSideAssoc, StoichVector,
  RHSFromCRN, NonInteractingQ, ULinearMassActionQ, BuildGU,
  SpanningTreeConditionQ, MaxLinearEliminableSet,
  LinearEliminationRules, LinearReduceRHS, ManualLinearReduce, StoichMatrix,subsetQ,
properSubsetQ,maximalByInclusion,maxElim];

subsetQ[A_List, B_List] := Complement[Union[A], Union[B]] === {};
properSubsetQ[A_List, B_List] := subsetQ[A, B] && Length[Union[A]] < Length[Union[B]];

(* Keep only inclusion-maximal sets *)
maximalByInclusion[sets_List] := Module[
  {sorted, kept = {}, s},
  sorted = SortBy[sets, -Length[Union[#]] &];   (* largest first *)
  Do[
    If[!AnyTrue[kept, subsetQ[s, #] &],
      AppendTo[kept, s]
    ],
    {s, sorted}
  ];
  kept
];

(* Main routine: maximal linear eliminations (FSW conditions) *)
maxElim[RN_, rts_, var_List] := Module[
  {n, allSubsets, validIdx, Ustr, g, maximalIdx, maximalSets, replacementRules,
   spe, alp, bet, gam, RHS, U, X, qssEqs, qssRules},

  n = Length[var];
  allSubsets = Subsets[Range[n], {1, n - 1}];

  validIdx = Reap[
      Do[
        Ustr = var[[idx]];
        g = BuildGU[RN, Ustr];
        If[
          NonInteractingQ[RN, Ustr] &&
          ULinearMassActionQ[rts, Ustr] &&
          SpanningTreeConditionQ[g, Ustr],
          Sow[idx]
        ],
        {idx, allSubsets}
      ]
    ][[2]];
  validIdx = If[validIdx === {}, {}, First[validIdx]];

  maximalIdx = maximalByInclusion[validIdx];
  maximalSets = var[[#]] & /@ maximalIdx;

  (* Compute RHS from RN and rts *)
  {spe, alp, bet, gam} = Take[extMat[RN, var], 4];
  RHS = gam . rts;

  (* Compute replacement rules and filter trivial QSS *)
  Module[{nonTrivialIdx, nonTrivialSets, nonTrivialRules, qssEqs, qssRules, isTrivial},
    {nonTrivialIdx, nonTrivialSets, nonTrivialRules} = {{}, {}, {}};

    Do[
      U = maximalSets[[i]];
      X = Complement[var, U];
      (* QSS: set dU/dt = 0 and solve for U in terms of X *)
      qssEqs = Thread[RHS[[Flatten[Position[var, #] & /@ U]]] == 0];
      qssRules = Solve[qssEqs, U];

      If[Length[qssRules] > 0,
        (* Check if solution is trivial (all U variables -> 0) *)
        isTrivial = AllTrue[U, (# /. qssRules[[1]]) === 0 &];
        If[!isTrivial,
          AppendTo[nonTrivialIdx, maximalIdx[[i]]];
          AppendTo[nonTrivialSets, maximalSets[[i]]];
          AppendTo[nonTrivialRules, qssRules[[1]]];
        ];
      ];
    , {i, Length[maximalSets]}];

    {nonTrivialIdx, nonTrivialSets, nonTrivialRules}
  ]
];


(* parse complex like 2*"i1" + "s" into Association["i1"->2,"s"->1] *)
ComplexAssoc[expr_] := Module[{terms, addTerm, a = <||>, t, c, sp},
  If[expr === 0, Return[<||>]];
  terms = List @@ Expand[expr] /. Plus -> List;
  addTerm[term_] := Module[{coef = 1, sym},
    Which[
      term === 0, Null,
      MatchQ[term, _String], a[term] = Lookup[a, term, 0] + 1,
      MatchQ[term, Times[_Integer, _String]],
        coef = First@term; sym = Last@term;
        a[sym] = Lookup[a, sym, 0] + coef,
      MatchQ[term, Times[_Rational, _String]],
        coef = First@term; sym = Last@term;
        a[sym] = Lookup[a, sym, 0] + coef,
      True,
        (* allow 1*string, etc *)
        t = term /. Times -> List;
        If[Length[t] == 2 && IntegerQ[t[[1]]] && MatchQ[t[[2]], _String],
          a[t[[2]]] = Lookup[a, t[[2]], 0] + t[[1]],
          (* fallback: ignore unexpected *)
          Null
        ]
    ]
  ];
  Scan[addTerm, terms];
  KeySelect[a, (# =!= 0) &]
];

AssocToExpr[a_Association] := Total[KeyValueMap[#2*#1 &, a]] /. {0 -> 0};

RxnSideAssoc[rxn_] := Module[{lhs, rhs},
  {lhs, rhs} = List @@ rxn;
  {ComplexAssoc[lhs], ComplexAssoc[rhs]}
];

StoichVector[species_List, lhs_Association, rhs_Association] :=
  (Lookup[rhs, #, 0] - Lookup[lhs, #, 0]) & /@ species;

(* build RHS from RN,rts (mass-action or general rates in rts list) *)
RHSFromCRN[RN_List, rts_List] := Module[{sp, n, m, sto, rates, lhs, rhs, v},
  sp = extSpe[RN];
  n = Length[sp]; m = Length[RN];
  rates = rts;
  sto = ConstantArray[0, {n, m}];
  Do[
    {lhs, rhs} = RxnSideAssoc[RN[[k]]];
    v = StoichVector[sp, lhs, rhs];
    sto[[All, k]] = v;
  , {k, 1, m}];
  Thread[ToExpression /@ sp -> (sto . rates)]
];

(* FSW noninteracting: no complex contains two distinct U species; each U coeff is 0 or 1 *)
NonInteractingQ[RN_List, U_List] := Module[{u = AssociationThread[U -> True], ok = True, lhs, rhs, a, us, coeffs},
  Do[
    {lhs, rhs} = RxnSideAssoc[rxn];
    Do[
      a = side;
      us = Select[Keys[a], KeyExistsQ[u, #] &];
      If[Length[us] >= 2, ok = False; Break[]];
      coeffs = Lookup[a, us, {}];
      If[Or @@ ( (# > 1) & /@ coeffs), ok = False; Break[]];
    , {side, {lhs, rhs}}];
    If[!ok, Break[]];
  , {rxn, RN}];
  ok
];

(* U-linear for mass-action: each rate is at most linear in each Ui, and no rate contains product Ui*Uj *)
ULinearMassActionQ[rts_List, U_List] := Module[{uSyms = ToExpression /@ U, ok = True, rt, exps, totUexp},
  Do[
    rt = rts[[k]];
    exps = Exponent[rt, #] & /@ uSyms;
    If[Or @@ ( (# > 1) & /@ exps), ok = False; Break[]];
    totUexp = Total[exps];
    If[totUexp > 1, ok = False; Break[]];
  , {k, Length[rts]}];
  ok
];

(* Build G_U: nodes U plus "*" ; directed edges for reactions involving U in reactant/product *)
BuildGU[RN_List, U_List] := Module[
  {u = AssociationThread[U -> True], edges = {}, lhs, rhs, ul, ur, src, tgt},
  Do[
    {lhs, rhs} = RxnSideAssoc[rxn];
    ul = Select[Keys[lhs], KeyExistsQ[u, #] &];
    ur = Select[Keys[rhs], KeyExistsQ[u, #] &];
    If[ul == {} && ur == {}, Continue[]];
    src = If[ul == {}, "*", First@ul];
    tgt = If[ur == {}, "*", First@ur];
    AppendTo[edges, src -> tgt];
  , {rxn, RN}];
  Graph[Join[U, {"*"}], edges, VertexLabels -> "Name", DirectedEdges -> True]
];

(* spanning-tree condition: per undirected component H:
   if "*" in H then every node reaches "*";
   else exists a node r in H reachable from all nodes *)
SpanningTreeConditionQ[g_Graph, U_List] := Module[
  {ug, comps, ok = True, H, nodes, hasStar, subg, reachToStar, cand, allReach},
  ug = UndirectedGraph[g];
  comps = ConnectedComponents[ug];
  Do[
    nodes = comps[[k]];
    hasStar = MemberQ[nodes, "*"];
    subg = Subgraph[g, nodes];
    If[hasStar,
      reachToStar = And @@ (VertexInComponentQ[subg, #] && (Length@FindPath[subg, #, "*", Infinity] > 0 || # == "*") & /@ nodes);
      If[!reachToStar, ok = False; Break[]],
      cand = Select[nodes, # =!= "*" &];
      allReach = Or @@ Table[
        And @@ Table[Length@FindPath[subg, v, r, Infinity] > 0 || v == r, {v, cand}],
        {r, cand}
      ];
      If[!allReach, ok = False; Break[]]
    ];
  , {k, Length[comps]}];
  ok
];

(* greedy maximal eliminable subset starting from U0 *)
greedyElim[RN_List, rts_List, U0_List] := Module[
  {U = DeleteDuplicates@U0, changed = True, g},
  While[changed,
    changed = False;
    If[!NonInteractingQ[RN, U],
      U = Rest[U]; changed = True; Continue[]
    ];
    If[!ULinearMassActionQ[rts, U],
      U = Rest[U]; changed = True; Continue[]
    ];
    g = BuildGU[RN, U];
    If[!SpanningTreeConditionQ[g, U],
      U = Rest[U]; changed = True; Continue[]
    ];
  ];
  U
];

(* find maximal eliminable sets by trying different orderings *)
(* NOTE: This is a heuristic algorithm using greedy search from different
   starting permutations. It may not find ALL possible maximal sets, but
   finds several representative ones. For specific elimination sets, the
   user can verify FSW conditions and use directly. *)
(* OUTPUT: Returns indices in species list (var) instead of species strings *)
maxElim[RN_List, rts_List] := Module[
  {sp = extSpe[RN], maxSets = {}, U, perms, i, maxSetsIdx},

  (* Helper to check if set already in list (compare sorted) *)
  inList[set_, list_] := MemberQ[Sort /@ list, Sort[set]];

  (* Try greedy from all species *)
  U = greedyElim[RN, rts, sp];
  If[Length[U] > 0, AppendTo[maxSets, U]];

  (* Try greedy from different orderings *)
  perms = {Reverse[sp], Sort[sp], Sort[sp, Greater]};
  Do[
    U = greedyElim[RN, rts, perm];
    If[Length[U] > 0 && !inList[U, maxSets],
      AppendTo[maxSets, U]
    ],
  {perm, perms}];

  (* Try complement of each found set (if A is eliminable, maybe sp\A is too) *)
  Do[
    U = greedyElim[RN, rts, Complement[sp, maxSets[[i]]]];
    If[Length[U] > 0 && !inList[U, maxSets],
      AppendTo[maxSets, U]
    ],
  {i, Length[maxSets]}];

  (* Try starting from different random permutations *)
  Do[
    U = greedyElim[RN, rts, RandomSample[sp]];
    If[Length[U] > 0 && !inList[U, maxSets],
      AppendTo[maxSets, U]
    ],
  {i, 30}];  (* Increased from 10 to 30 for better coverage *)

  (* Filter to keep only truly maximal sets (not contained in any other) *)
  maxSets = DeleteDuplicates[maxSets];
  maxSets = Select[maxSets,
    Function[s, !AnyTrue[DeleteCases[maxSets, s], ContainsAll[#, s]&]]
  ];

  (* Convert species strings to indices in sp *)
  maxSetsIdx = Map[
    Function[elimSet, Map[Position[sp, #][[1,1]]&, elimSet]],
    maxSets
  ];

  maxSetsIdx
];

(* compute stoichiometric matrix gamma from RN *)
StoichMatrix[RN_List] := Module[{sp, m, sto, lhs, rhs, v, k},
  sp = extSpe[RN];
  m = Length[RN];
  sto = ConstantArray[0, {Length[sp], m}];
  Do[
    {lhs, rhs} = RxnSideAssoc[RN[[k]]];
    v = StoichVector[sp, lhs, rhs];
    sto[[All, k]] = v,
  {k, m}];
  sto
];

(* compute linear elimination rules u -> phi(x) using conservation laws *)
LinearEliminationRules[RN_List, rts_List, U_List, consLaws_:Automatic] := Module[
  {RHS, sp, uSyms, xSyms, eqs, A, b, sol, gamma, consVecs, uIdx, rank, indep},
  RHS = RHSFromCRN[RN, rts];
  sp = extSpe[RN];
  uSyms = ToExpression /@ U;
  xSyms = ToExpression /@ Complement[sp, U];
  eqs = (ToExpression[#] /. RHS) & /@ U;

  (* CoefficientArrays returns {b, A} where eqs == A.u + b *)
  {b, A} = CoefficientArrays[eqs, uSyms] // Normal;

  (* Check rank and handle singular systems *)
  rank = MatrixRank[A];
  If[rank < Length[uSyms],
    (* Singular system - use least-norm solution via pseudo-inverse *)
    sol = Quiet@PseudoInverse[A] . (-b),
    (* Non-singular - use direct solve *)
    sol = Quiet@LinearSolve[A, -b]
  ];
  Thread[uSyms -> sol]
];

(* reduced RHS after eliminating U via steady-state linear rules *)
LinearReduceRHS[RN_List, rts_List, Uspec_:Automatic] := Module[
  {allMaxElim, U, rules, RHS, sp, keep, redRHS},

  (* If U not specified, find all maximal eliminable sets and use first *)
  If[Uspec === Automatic,
    allMaxElim = maxElim[RN, rts];
    If[Length[allMaxElim] == 0,
      Print["No eliminable sets found"];
      Return[{{}, {}, {}}]
    ];
    If[Length[allMaxElim] > 1,
      Print["Found ", Length[allMaxElim], " maximal eliminable sets"];
      Print["Using first: ", allMaxElim[[1]]];
      Print["Alternatives: ", Rest[allMaxElim]];
    ];
    U = allMaxElim[[1]],
    (* else use specified U *)
    U = Uspec
  ];

  rules = LinearEliminationRules[RN, rts, U];
  RHS = RHSFromCRN[RN, rts];
  sp = extSpe[RN];
  keep = Complement[sp, U];
  redRHS = Thread[ToExpression /@ keep -> ((ToExpression /@ keep /. RHS) /. rules // Together)];
  {U, rules, redRHS}
];

(* Manual reduction with user-specified conservation law and steady-state expressions *)
ManualLinearReduce[RN_List, rts_List, U_List, uExpr_List] := Module[
  {RHS, sp, keep, uSyms, rules, redRHS},
  RHS = RHSFromCRN[RN, rts];
  sp = extSpe[RN];
  uSyms = ToExpression /@ U;
  keep = Complement[sp, U];
  rules = Thread[uSyms -> uExpr];
  redRHS = Thread[ToExpression /@ keep -> ((ToExpression /@ keep /. RHS) /. rules // Together // Simplify)];
  {rules, redRHS}
];



(* ::Section:: *)
(*NOTES ON LIMITATIONS*)


(* The automatic LinearReduceRHS works for networks WITHOUT conservation laws among U species.
   For networks WITH conservation laws (e.g., enzyme-substrate systems), the proper FSW algorithm requires:
   1. Spanning tree enumeration in G_U
   2. Cycle-based construction of reduced reactions
   3. Weight calculations using tree paths

   For such cases, use ManualLinearReduce with explicitly computed steady-state expressions,
   or compute conservation laws using cons[StoichMatrix[RN]] and solve manually. *)



(* ::Section:: *)
(*QSS Solution Utilities*)


(* solveRP: solve for rational positive solutions (solveRatPos abbreviated) *)
(* Excludes boundary equilibria (solutions with zeros), returns interior solutions only *)
(* Output: {status, solutions} where status in {"Success", "OnlyBoundary", "TimedOut", "NotRational", "NoSolution"} *)
solveRP[eqs_, vars_, timeout_: 10] := Module[
  {sol, filtered, isRational, hasZero},

  (* Helper: check if solution has any zero component *)
  hasZero[s_] := AnyTrue[vars /. s, # === 0 || PossibleZeroQ[#] &];

  (* Helper: check if solution is rational (can be tested with onlyP) *)
  isRational[s_] := AllTrue[vars /. s,
    Head[#] =!= Root && FreeQ[#, _Root] && FreeQ[#, _Complex] &];

  (* Solve with timeout *)
  sol = TimeConstrained[
    Quiet[Solve[eqs, vars]],
    timeout,
    $TimedOut
  ];

  (* Check for timeout *)
  If[sol === $TimedOut,
    Return[{"TimedOut", {}}]
  ];

  (* Check for no solutions *)
  If[Length[sol] == 0,
    Return[{"NoSolution", {}}]
  ];

  (* Filter out boundary equilibria (solutions with zeros) *)
  filtered = Select[sol, !hasZero[#] &];

  (* Check if filtered solutions are empty *)
  If[Length[filtered] == 0,
    Return[{"OnlyBoundary", sol}]  (* Return original for reference *)
  ];

  (* Check if solutions are rational *)
  If[!AllTrue[filtered, isRational],
    Return[{"NotRational", filtered}]
  ];

  (* Success: return interior rational solutions *)
  {"Success", filtered}
];

(* Helper: show reduced model with reaction network *)
showRedMod[RN_, rts_, U_, qss_] := Module[{rtsRed},
  rtsRed = rts /. qss;
  Print[RN // Length, " reactions, eliminating ", U, ":"];
  Print[Transpose[{RN, rtsRed}] // MatrixForm]
];

(* testEP: test endemic positivity for eliminable sets (FLWW Lemma 11 blanket condition i) *)
(* Shows each reduced model and tests positivity *)
testEP[RN_, rts_, spe_, RHS_, allMaxElim_, parReg1_, parReg2_] := Module[
  {sampleStates},

  Print["\n==== FLWW Lemma 11 Positivity Check ===="];
  Print["Test: reduced rates positive (blanket condition i)\n"];

  Do[
    Module[{U, Uvar, Xvar, qssEqs, qssSol, status, hasViolation, regime1Neg, regime2Neg},
      U = allMaxElim[[idx]];
      Print["[", idx, "] U = ", U];

      (* Convert to symbols *)
      Uvar = ToExpression /@ U;
      Xvar = Complement[ToExpression /@ spe, Uvar];

      (* Sample states for slow variables *)
      sampleStates = Which[
        Sort[U] == Sort[{"i1", "i2"}],
        {{s -> 0.7, r1 -> 0.1, r2 -> 0.05}, {s -> 0.5, r1 -> 0.2, r2 -> 0.15}, {s -> 0.3, r1 -> 0.25, r2 -> 0.2}},
        Sort[U] == Sort[{"r2", "r1", "s"}],
        {{i1 -> 0.1, i2 -> 0.05}, {i1 -> 0.15, i2 -> 0.1}, {i1 -> 0.2, i2 -> 0.15}},
        Sort[U] == Sort[{"i2", "r2"}],
        {{s -> 0.6, i1 -> 0.1, r1 -> 0.15}, {s -> 0.5, i1 -> 0.15, r1 -> 0.2}, {s -> 0.4, i1 -> 0.2, r1 -> 0.25}},
        Sort[U] == Sort[{"r1", "i1"}],
        {{s -> 0.6, i2 -> 0.05, r2 -> 0.1}, {s -> 0.5, i2 -> 0.1, r2 -> 0.15}, {s -> 0.4, i2 -> 0.15, r2 -> 0.2}},
        True,
        {{s -> 0.5, i1 -> 0.1, i2 -> 0.05, r1 -> 0.15, r2 -> 0.1}}
      ];

      (* QSS equations: set U' = 0 *)
      Module[{UindRHS},
        UindRHS = Position[ToExpression /@ spe, #][[1, 1]] & /@ Uvar;
        qssEqs = Thread[(RHS[[#]] & /@ UindRHS) == 0];
      ];

      (* Solve for rational positive solutions *)
      {status, qssSol} = solveRP[qssEqs, Uvar, 10];

      Which[
        status == "Success",
          hasViolation = False;
          Print["  Found ", Length[qssSol], " endemic QSS solution(s)"];
          Do[
            Module[{sol},
              sol = qssSol[[solIdx]];
              (* Show reduced model *)
              showRedMod[RN, rts, U, sol];
              (* Test positivity *)
              regime1Neg = Table[
                Count[rts /. sol /. parReg1 /. state, _?Negative, Infinity],
                {state, sampleStates}
              ];
              regime2Neg = Table[
                Count[rts /. sol /. parReg2 /. state, _?Negative, Infinity],
                {state, sampleStates}
              ];
              If[Max[regime1Neg] > 0 || Max[regime2Neg] > 0,
                Print["  \[Cross] Sol ", solIdx, " VIOLATES Lem11: reg1 neg=", Max[regime1Neg], " reg2 neg=", Max[regime2Neg]];
                hasViolation = True,
                Print["  \[Checkmark] Sol ", solIdx, " positive at test points"]
              ]
            ],
            {solIdx, Length[qssSol]}
          ];
          If[hasViolation, Print["  INADMISSIBLE"]],
        status == "OnlyBoundary",
          Print["  Only DFE (skip)"],
        status == "TimedOut",
          Print["  Timed out"],
        status == "NotRational",
          Print["  Not rational"],
        True,
          Print["  No QSS solution"]
      ];
    ],
    {idx, Length[allMaxElim]}
  ];

  Print["\n==== Symbolic Analysis (onlyP) ====\n"];

  Do[
    Module[{U, Uvar, qssEqs, qssSol, status},
      U = allMaxElim[[idx]];
      Print["[", idx, "] U = ", U];

      Uvar = ToExpression /@ U;
      Module[{UindRHS},
        UindRHS = Position[ToExpression /@ spe, #][[1, 1]] & /@ Uvar;
        qssEqs = Thread[(RHS[[#]] & /@ UindRHS) == 0];
      ];

      {status, qssSol} = solveRP[qssEqs, Uvar, 10];

      Which[
        status == "Success",
          Print["  Testing ", Length[qssSol], " symbolic solution(s)"];
          Do[
            Module[{rtsReduced, symbolicTests, allSymbolicPositive},
              rtsReduced = rts /. qssSol[[solIdx]];
              symbolicTests = Table[
                Module[{num, den},
                  {num, den} = NumeratorDenominator[Together[Expand[rtsReduced[[k]]]]];
                  onlyP[num] && onlyP[den]
                ],
                {k, Length[rtsReduced]}
              ];
              allSymbolicPositive = And @@ symbolicTests;
              If[allSymbolicPositive,
                Print["  Sol ", solIdx, ": \[Checkmark] Symbolically positive"],
                Print["  Sol ", solIdx, ": \[Cross] ", Count[symbolicTests, False], " rates with negative coefficients"]
              ]
            ],
            {solIdx, Length[qssSol]}
          ],
        status == "OnlyBoundary",
          Print["  Only DFE"],
        True,
          Print["  Cannot solve"]
      ];
    ],
    {idx, Length[allMaxElim]}
  ];

  Print["\nRemark: Lemma 11 requires all rates positive for bio-relevant parameters.\n"];
];



(* ::Section:: *)
(*END PACKAGE*)


(* Run testRed.nb for examples *)

(* Test reduced model approximation *)
testRed[RHS_, var_, Uc_, params_, T_, ic0_] := Module[
  {varSym, RHSsym, eqsFull, icFull, solFull, moc, cec, qssc,
   UvarSym, XvarSym, RHSred, eqsRed, icRed, solRed, lvVars, cekEnd, err},

  varSym = ToExpression /@ var;
  RHSsym = RHS /. Thread[var -> varSym];

  (* Full system *)
  eqsFull = Table[D[varSym[[i]][t], t] == (RHSsym[[i]] /.
    Thread[varSym -> (# [t] & /@ varSym)]), {i, Length[varSym]}];
  icFull = Thread[(# [0] & /@ varSym) == ic0];
  solFull = NDSolve[Join[eqsFull /. params, icFull], varSym, {t, 0, T}][[1]];

  (* Compute QSS *)
  moc = Complement[Range[Length[var]], Uc];
  cec = Quiet[Solve[Thread[RHS[[Uc]]==0], var[[Uc]]]];

  (* Select endemic QSS *)
  lvVars = Intersection[Uc, {2,3}];
  cekEnd = If[Length[lvVars] > 0,
    Select[cec, And @@ Table[(var[[i]] /. #) =!= 0, {i, lvVars}] &], cec];
  qssc = If[Length[cekEnd] > 0, cekEnd[[1]], If[Length[cec] > 0, cec[[1]], {}]];

  (* Reduced system *)
  UvarSym = ToExpression /@ (var[[#]]& /@ Uc);
  XvarSym = ToExpression /@ (var[[#]]& /@ moc);
  RHSred = (RHSsym[[#]]& /@ moc) /. qssc /. params;
  eqsRed = Table[D[XvarSym[[i]][t], t] == (RHSred[[i]] /.
    Thread[XvarSym -> (# [t] & /@ XvarSym)]), {i, Length[XvarSym]}];
  icRed = Table[XvarSym[[i]][0] == ic0[[Position[varSym, XvarSym[[i]]][[1,1]]]],
    {i, Length[XvarSym]}];
  solRed = NDSolve[Join[eqsRed, icRed], XvarSym, {t, 0, T}][[1]];

  (* L2 error *)
  err = Table[NIntegrate[((XvarSym[[i]][t] /. solFull) - (XvarSym[[i]][t] /. solRed))^2,
    {t, 0, T}], {i, Length[XvarSym]}];

  {qssc, XvarSym, err, solFull, solRed}
];
