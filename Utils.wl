(* ::Package:: *)

(* ========================================================================= *)
(* PHASE 3: ENHANCED USER EXPERIENCE *)
(* ========================================================================= *)

(* validateNetwork: Validates reaction network and rate consistency
   RN: reaction network
   rts: rates (optional)

   Returns: Association with validation results
*)
validateNetwork[RN_, rts_:Automatic] := Module[{
  issues, warnings, speciesCheck, rateCheck, stoichCheck,
  species, numReactions, stringCheck, coeffCheck},

  issues = {};
  warnings = {};

  (* Check 1: Verify reaction format *)
  If[!ListQ[RN],
    AppendTo[issues, "RN must be a list of reactions"];
    Return[<|"Valid" -> False, "Issues" -> issues|>];
  ];

  numReactions = Length[RN];

  (* Check 2: Verify each reaction is properly formatted *)
  Do[
    If[!MatchQ[RN[[i]], {_, _}] && !MatchQ[RN[[i]], _ -> _],
      AppendTo[issues, "Reaction " <> ToString[i] <> " is not in proper format {lhs, rhs} or lhs -> rhs"];
    ];
  , {i, numReactions}];

  If[issues =!= {},
    Return[<|"Valid" -> False, "Issues" -> issues|>];
  ];

  (* Check 3: Verify species are strings (critical note from Section 8) *)
  stringCheck = Check[
    species = extSpe[RN];
    True,
    False
  ];

  If[!stringCheck,
    AppendTo[issues, "Failed to extract species - check that all species are strings"];
  ];

  (* Check 4: Verify positive stoichiometric coefficients *)
  Do[
    Module[{lhs, rhs, lhsAsso, rhsAsso},
      {lhs, rhs} = If[MatchQ[RN[[i]], _ -> _], List @@ RN[[i]], RN[[i]]];

      lhsAsso = Check[comp2Asso[lhs], $Failed];
      rhsAsso = Check[comp2Asso[rhs], $Failed];

      If[lhsAsso === $Failed || rhsAsso === $Failed,
        AppendTo[issues, "Reaction " <> ToString[i] <> " has invalid stoichiometry"];
      , (* else *)
        (* Check for negative coefficients *)
        If[AnyTrue[Values[lhsAsso], # < 0 &] || AnyTrue[Values[rhsAsso], # < 0 &],
          AppendTo[warnings, "Reaction " <> ToString[i] <> " contains negative stoichiometric coefficients"];
        ];
      ];
    ];
  , {i, numReactions}];

  (* Check 5: Validate rates if provided *)
  If[rts =!= Automatic,
    If[Length[rts] != numReactions,
      AppendTo[issues, "Number of rates (" <> ToString[Length[rts]] <>
        ") does not match number of reactions (" <> ToString[numReactions] <> ")"];
    ];
  ];

  (* Return validation result *)
  <|
    "Valid" -> (issues === {}),
    "Issues" -> issues,
    "Warnings" -> warnings,
    "NumReactions" -> numReactions,
    "NumSpecies" -> If[stringCheck, Length[species], "unknown"],
    "Species" -> If[stringCheck, species, "extraction failed"]
  |>
];

(* estimateComplexity: Predicts computation time before analysis
   RN: reaction network

   Returns: Association with complexity estimates
*)
estimateComplexity[RN_] := Module[{
  species, numReactions, numSpecies, numSiphons, deficiency,
  expectedTime, complexity, recommendation, matResults},

  (* Extract basic metrics *)
  matResults = Check[extMat[RN], $Failed];

  If[matResults === $Failed,
    Return[<|"Complexity" -> "unknown", "Recommendation" -> "Network validation failed"|>];
  ];

  species = matResults[[1]];
  numSpecies = Length[species];
  numReactions = Length[RN];

  (* Estimate number of siphons (exponential in worst case) *)
  numSiphons = Min[2^numSpecies - 1, 1000]; (* upper bound for estimation *)

  (* Complexity heuristics *)
  complexity = Which[
    numSpecies <= 3 && numReactions <= 5, "Low",
    numSpecies <= 5 && numReactions <= 10, "Medium",
    numSpecies <= 7 && numReactions <= 15, "High",
    True, "Very High"
  ];

  expectedTime = Which[
    complexity === "Low", "< 10 seconds",
    complexity === "Medium", "10 seconds - 1 minute",
    complexity === "High", "1-5 minutes",
    complexity === "Very High", "> 5 minutes, timeouts likely"
  ];

  recommendation = Which[
    complexity === "Low", "Full symbolic analysis recommended",
    complexity === "Medium", "Symbolic analysis should complete; consider caching",
    complexity === "High", "Consider using bdAn instead of bd2; expect some timeouts",
    complexity === "Very High", "Recommend preliminary numerical analysis; symbolic timeouts expected"
  ];

  <|
    "Complexity" -> complexity,
    "NumSpecies" -> numSpecies,
    "NumReactions" -> numReactions,
    "EstimatedSiphons" -> "≤ " <> ToString[numSiphons],
    "ExpectedTime" -> expectedTime,
    "Recommendation" -> recommendation,
    "Tips" -> {
      "Use clearEpidCRNCache[] to free memory if needed",
      "bdAn is faster than bd2 for preliminary analysis",
      "Consider fixing some parameters to reduce symbolic complexity"
    }
  |>
];

(* explainTimeout: Suggests model simplifications when timeouts occur
   bdfpT: boundary fixed points result from bdFp
   RN: reaction network (optional)

   Returns: Association with diagnostic information
*)
explainTimeout[bdfpT_, RN_:Automatic] := Module[{
  frozenFacets, totalFacets, frozenRatio, suggestions, diagnosis},

  (* Count frozen facets *)
  frozenFacets = Count[bdfpT, {"froze", "froze"}];
  totalFacets = Length[bdfpT];
  frozenRatio = N[frozenFacets / totalFacets];

  diagnosis = Which[
    frozenRatio == 0, "No timeouts - analysis completed successfully",
    frozenRatio < 0.3, "Partial timeout - " <> ToString[frozenFacets] <> " of " <> ToString[totalFacets] <> " facets",
    frozenRatio < 0.7, "Significant timeout - over half of facets failed",
    True, "Complete timeout - most facets failed to solve"
  ];

  suggestions = {};

  If[frozenRatio > 0,
    AppendTo[suggestions, "Consider using FindInstance for numerical solutions"];
    AppendTo[suggestions, "Try fixing some parameters to reduce symbolic complexity"];
    AppendTo[suggestions, "Use bdAn instead of bd2 for preliminary analysis"];
  ];

  If[frozenRatio > 0.5,
    AppendTo[suggestions, "Model may be too complex for symbolic analysis"];
    AppendTo[suggestions, "Consider simplifying the reaction network"];
    AppendTo[suggestions, "Numerical continuation methods may be more appropriate"];
  ];

  If[RN =!= Automatic,
    Module[{complexity},
      complexity = estimateComplexity[RN];
      AppendTo[suggestions, "Model complexity: " <> complexity["Complexity"]];
      AppendTo[suggestions, complexity["Recommendation"]];
    ];
  ];

  <|
    "Diagnosis" -> diagnosis,
    "FrozenFacets" -> frozenFacets,
    "TotalFacets" -> totalFacets,
    "FrozenRatio" -> frozenRatio,
    "Suggestions" -> suggestions
  |>
];

(* ========================================================================= *)
(* END PHASE 3: ENHANCED USER EXPERIENCE *)
(* ========================================================================= *)

posM= 
Replace[#,{_?Negative->0,e_:>Replace[Expand[e],
{Times[_?Negative,_]->0,
t_Plus:>Replace[t,_?Negative|Times[_?Negative,_]->0,1]}]},{2}]&;
(*posL = Replace[#, {
  _?Negative -> 0,
  e_ :> Replace[Expand[e], {
    Times[_?Negative, _] -> 0,
    t_Plus :> Replace[t, _?Negative | Times[_?Negative, _] -> 0, 1]
  }]
}, {1}]&;*)


reCL[re_] :=DeleteCases[re, _Symbol > 0 | Subscript[_, __] > 0, Infinity];

 remZ[li_]:=Select[li, # =!= 0 &];

whenP=Function[exprs,
  Module[{params,reducedConditions},
    params=Variables[exprs];
    reducedConditions=Reduce[And @@ (#>0 & /@ exprs)&&And @@ (#>0 & /@ params),params,Reals];
    Simplify[reducedConditions,Assumptions->And @@ (#>0 & /@ params)]
  ]
];

(*whenP[{\[Mu]3/\[Alpha]3,r (K \[Alpha]3-\[Mu]3)/(K \[Alpha]3^2)}]*)

Hur4M[mat_] := Module[{coefs, det, checks},
  coefs = CharacteristicPolynomial[mat, \[Lambda]] // CoefficientList[#, \[Lambda]] &;
  det = Det[mat];
  checks = Positive @@@ {coefs, {det}};
  {coefs, det, checks}
];
(* Hur4M[{{1, 2}, {3, 4}}] *)

makeLPM[mat_] := Module[{minors},
  minors = Table[Det[Take[Take[mat, i], All, i]], {i, Length[mat]}];
  minors
];
(* makeLPM[{{2, 1}, {1, 3}}] *)

GetVec[mat_] := Module[{vals, vecs},
  {vals, vecs} = Eigensystem[mat];
  Transpose[{vals, vecs}]
];
(* GetVec[{{0, 1}, {-2, 0}}] *)

Deg[poly_, var_] := Exponent[poly, var];
(* Deg[x^3 + 2 x^2 + 1, x] *)

Stab[mat_] := Eigenvalues[mat];
(* Stab[{{0, 1}, {-2, 0}}] *)

(*
BadModule[] := Module[{mat = {{x + y, x - y}, {x, y}}},
  mat[[3]]
];
(* BadModule[] *)
*)

Grobpol[RHS_, var_, par_, ind_, cn_ : {}] := Module[{dyn, X, eq, elim, pol, ratsub},
  dyn = RHS; 
  X = var; 
  eq = Thread[dyn == 0]; 
  elim = Complement[Range[Length[X]], ind];
  pol = Collect[GroebnerBasis[Numerator[Together[dyn /. cn]], 
    Join[par, X[[ind]]], X[[elim]]], X[[ind]]];
  { pol}
]
