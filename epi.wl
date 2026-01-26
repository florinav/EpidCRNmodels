(* ::Package:: *)

(* ========================================================================= *)
(* EPI - MINIMAL EPIDEMIC MODELING SUBPACKAGE *)
(* Contains the essential functions for epidemic CRN modeling *)
(* IMPORTANT: This conains ODE2RN , extMat, minSiph, ... *)
(* Goal: Serve as foundation for future Python port *)
(* NOTE: Backup versions exist in Core1.wl and EpidCRNo.wl *)
(* ========================================================================= *)


posM=
Replace[#,{_?Negative->0,e_:>Replace[Expand[e],
{Times[_?Negative,_]->0,
t_Plus:>Replace[t,_?Negative|Times[_?Negative,_]->0,1]}]},{2}]&;


(* ========================================================================= *)
(* FROM CORE.WL *)
(* ========================================================================= *)

Par[RHS_,X_]:=Complement[Variables[RHS],X];

chk[RHS_] := Module[{vars, reserved, conflicts},
  vars = ToString /@ Variables[RHS];
  reserved = {"RHS","var", "rts", "RN","V", "F", "K", "Jx", "Jy", "Jxy", "Jyx", "Kd", "mSi",
              "alp", "bet", "gam", "spe", "E0", "cE0", "cDFE",
               "inf", "infc", "mod", "par", "cp","I","E"};
  conflicts = Intersection[vars, reserved];
  If[Length[conflicts] > 0,
    Print["WARNING: RHS variables conflict with package names: ", conflicts, " from ", reserved];
  ];
  ToExpression /@ conflicts
];

extSpe[reactions_] := Module[{allSpecies, reactants, products},
  allSpecies = {};
  Do[
   reactants = comp2Asso[reactions[[i, 1]]];
   products = comp2Asso[reactions[[i, 2]]];
   allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
   , {i, Length[reactions]}];
  DeleteDuplicates[allSpecies]
];



(* Main extMat function - now handles symbolic stoichiometry *)

allT[expr_] := Module[{expanded},
  expanded = Expand[expr];
  If[Head[expanded] === Plus, List @@ expanded, {expanded}]
];

(* isN: Detect if a term has negative sign (complement of posM logic) *)
isN[term_] := Module[{},
  Which[
    (* Numeric negative *)
    NumericQ[term] && Negative[term], True,
    (* Times expression starting with -1 *)
    Head[term] === Times && MemberQ[List @@ term, -1], True,
    (* Times expression with negative numeric coefficient *)
    Head[term] === Times && NumericQ[First[List @@ term]] && Negative[First[List @@ term]], True,
    (* Otherwise not negative *)
    True, False
  ]
];

ODE2RN[RHS_List, var_List, prF_:Identity] := Module[{
  n, allTermsByEq, allTermsFlat, sources, sinks, rtsRaw,
  gam, alp, bet, RN, spe, ii, jj, kk,
  alpRaw, gamRaw, uniqueAlpCols, mergeIdx, rts,
  coefToExpr, catCoef, newAlp, newBet, rnRed
  },

  n = Length[var];
  spe = ToString /@ var;  (* CRITICAL: spe must equal ToString[var], NOT lowercase *)

  (* STEP 1: Extract all terms from each equation using allT *)
  allTermsByEq = Table[allT[RHS[[ii]]], {ii, n}];
  allTermsFlat = Flatten[allTermsByEq];

  (* STEP 2: Separate into sources (positive-only) and sinks (negative) using isN *)
  sources = Select[allTermsFlat, !isN[#]&];
  sinks = Select[allTermsFlat, isN[#]&];
  sources = DeleteDuplicates[sources];
  sinks = DeleteDuplicates[sinks];

  (* STEP 3: Build rtsRaw from sources and negated sinks *)
  rtsRaw = Join[sources, -sinks];
  rtsRaw = DeleteDuplicates[rtsRaw];
  rtsRaw = Select[rtsRaw, # =!= 0 &];

  (* Build raw gamma matrix *)
  gamRaw = Table[
    Module[{eqTerms, coeff},
      eqTerms = allT[RHS[[ii]]];
      Table[
        coeff = 0;
        Do[
          Module[{term, termAbs},
            term = eqTerms[[kk]];
            termAbs = If[isN[term], -term, term];
            If[termAbs === rtsRaw[[jj]],
              coeff = If[isN[term], -1, 1];
              Break[];
            ];
          ],
          {kk, Length[eqTerms]}
        ];
        coeff
      , {jj, Length[rtsRaw]}]
    ],
    {ii, n}
  ];

  (* Build raw alpha matrix: exponents of variables in each rate *)
  alpRaw = Table[
    Module[{rate},
      rate = rtsRaw[[jj]];
      If[FreeQ[rate, var[[ii]]], 0, Exponent[rate, var[[ii]]]]
    ],
    {ii, n}, {jj, Length[rtsRaw]}
  ];

  (* STEP 4: Merge reactions with same source AND same gamma (same alp and gam columns) *)
  (* Create combined key: {alp column, gam column} *)
  Module[{alpGamPairs, uniquePairs},
    alpGamPairs = Table[{Transpose[alpRaw][[jj]], gamRaw[[All, jj]]}, {jj, Length[rtsRaw]}];
    uniquePairs = DeleteDuplicates[alpGamPairs];

    rts = {};
    gam = ConstantArray[0, {n, Length[uniquePairs]}];
    alp = ConstantArray[0, {n, Length[uniquePairs]}];

    Do[
      (* Find all indices that match this unique pair *)
      mergeIdx = Flatten[Position[alpGamPairs, uniquePairs[[jj]]]];
      (* Sum the rates *)
      AppendTo[rts, Total[rtsRaw[[mergeIdx]]]];
      (* Set alp and gam from the unique pair *)
      alp[[All, jj]] = uniquePairs[[jj, 1]];
      gam[[All, jj]] = uniquePairs[[jj, 2]],
      {jj, Length[uniquePairs]}
    ]
  ];

  (* Build beta matrix: bet = gam + alp *)
  bet = gam + alp;

  (* Fix negative stoichiometry: move negatives from bet to alp *)
  Do[
    If[bet[[ii, jj]] < 0,
      alp[[ii, jj]] = alp[[ii, jj]] - bet[[ii, jj]];
      bet[[ii, jj]] = 0
    ],
    {ii, n}, {jj, Length[rts]}
  ];

  (* Build reaction network RN using string format "s" + "i" -> 2*"i" *)
  RN = Table[
    Module[{left, right, leftTerms, rightTerms},
      (* Left side (reactants from alp) *)
      leftTerms = Table[
        If[alp[[kk, jj]] == 0, Nothing,
          If[alp[[kk, jj]] == 1, spe[[kk]], alp[[kk, jj]] * spe[[kk]]]
        ],
        {kk, n}
      ];
      left = If[Length[leftTerms] == 0, 0,
        If[Length[leftTerms] == 1, leftTerms[[1]], Plus @@ leftTerms]
      ];

      (* Right side (products from bet) *)
      rightTerms = Table[
        If[bet[[kk, jj]] == 0, Nothing,
          If[bet[[kk, jj]] == 1, spe[[kk]], bet[[kk, jj]] * spe[[kk]]]
        ],
        {kk, n}
      ];
      right = If[Length[rightTerms] == 0, 0,
        If[Length[rightTerms] == 1, rightTerms[[1]], Plus @@ rightTerms]
      ];

      left -> right
    ],
    {jj, Length[rts]}
  ];

  (* Compute reduced reactions (rnRed) with catalysts *)
  coefToExpr = Function[{coefs, vars},
    Module[{terms},
      terms = MapThread[If[#1 > 0, If[#1 == 1, #2, #1*#2], Nothing] &, {coefs, vars}];
      If[terms === {}, {}, terms]
    ]
  ];

  catCoef = MapThread[If[#1 > 0 && #2 > 0, Min[#1, #2], 0] &, {alp, bet}, 2];
  newAlp = alp - catCoef;
  newBet = bet - catCoef;

  rnRed = Table[
    Module[{catList, newL, newR, fixedL, fixedR, catExpr, lExpr, rExpr},
      catList = catCoef[[All, jj]];
      newL = newAlp[[All, jj]];
      newR = newBet[[All, jj]];

      (* Move negative stoichiometry from products to reactants *)
      fixedL = MapThread[If[#2 < 0, #1 - #2, #1] &, {newL, newR}];
      fixedR = MapThread[If[#2 < 0, 0, #2] &, {newL, newR}];

      catExpr = coefToExpr[catList, var];
      lExpr = coefToExpr[fixedL, var];
      rExpr = coefToExpr[fixedR, var];
      If[catExpr === {},
        lExpr -> rExpr,
        Row[{lExpr, Overscript[" \[Rule] ", catExpr], rExpr}]
      ]
    ],
    {jj, Length[rts]}
  ];

  Print[Length[RN], " reactions and rts=",
    Row[{
      Transpose[{prF[RN], prF[Factor /@ rts]}] // MatrixForm,
      Spacer[20],
      Transpose[{prF[rnRed], prF[Factor /@ rts]}] // MatrixForm
    }]
  ];

  {RN, rts, spe, alp, bet, gam, rnRed}
];

ODE2RNp[RHS_List, var_List, prF_:Identity] := Module[{
  n, RHSpol, RHSexp, allTerms, numList, denList, dens,
  RN, rts, spe, alp, bet, gam, rnRed, rtsTrue
  },

  n = Length[var];

  (* STEP 1: Expand RHS and extract all terms *)
  RHSexp = Expand /@ RHS;

  (* STEP 2: Process each equation to separate numerators and denominators *)
  RHSpol = Table[
    Module[{eq, terms, numPol},
      eq = RHSexp[[ii]];
      terms = If[Head[eq] === Plus, List @@ eq, {eq}];
      numPol = Total[Table[
        Module[{term, num, den},
          term = terms[[jj]];
          {num, den} = {Numerator[term], Denominator[term]};
          num
        ],
        {jj, Length[terms]}
      ]]
    ],
    {ii, n}
  ];

  (* STEP 3: Call ODE2RN on polynomial model *)
  {RN, rts, spe, alp, bet, gam, rnRed} =
    Block[{Print = (Null &)}, ODE2RN[RHSpol, var, Identity]];

  (* STEP 4: Build dens list by matching rts back to original terms *)
  dens = Table[
    Module[{rate, allDens, matchFound},
      rate = rts[[jj]];
      matchFound = False;
      allDens = 1;

      (* Search through all equations for terms whose numerator matches rate *)
      Do[
        Module[{eq, terms},
          eq = RHSexp[[ii]];
          terms = If[Head[eq] === Plus, List @@ eq, {eq}];
          Do[
            Module[{term, num, den},
              term = terms[[kk]];
              {num, den} = {Numerator[term], Denominator[term]};
              If[num === rate || num === -rate,
                allDens = 1/den;
                matchFound = True;
              ]
            ],
            {kk, Length[terms]}
          ];
          If[matchFound, Break[]];
        ],
        {ii, n}
      ];

      allDens
    ],
    {jj, Length[rts]}
  ];

  (* STEP 5: Compute true rates *)
  rtsTrue = MapThread[#1 * #2 &, {rts, dens}];

  (* STEP 6: Print always *)
  Print[Length[RN], " reactions and rts=",
    Row[{
      Transpose[{prF[RN], prF[Factor /@ rtsTrue]}] // MatrixForm,
      Spacer[20],
      Transpose[{prF[rnRed], prF[Factor /@ rtsTrue]}] // MatrixForm
    }]
  ];

  {RN, rtsTrue, spe, alp, bet, gam, rnRed}
];

defi[RN_] := Module[{complexes, adj, graph, l, s, gam, spe, al, be, Nc, ii, jj, numSpe, numRea, reactants, products},
  (* Get complexes *)
  complexes = DeleteDuplicates[Flatten[{First[#], Last[#]} & /@ RN]];
  Nc = Length[complexes];

  (* Build adjacency matrix for undirected graph *)
  adj = ConstantArray[0, {Nc, Nc}];
  Do[
    Module[{src, tgt, srcPos, tgtPos},
      src = First[RN[[ii]]];
      tgt = Last[RN[[ii]]];
      srcPos = Position[complexes, src][[1, 1]];
      tgtPos = Position[complexes, tgt][[1, 1]];
      adj[[srcPos, tgtPos]] = 1;
      adj[[tgtPos, srcPos]] = 1;
    ], {ii, Length[RN]}
  ];

  (* Linkage classes = connected components *)
  graph = AdjacencyGraph[adj];
  l = Length[ConnectedComponents[graph]];

  (* Compute gamma directly (avoid recursion with extMat) *)
  spe = extSpe[RN];
  numSpe = Length[spe];
  numRea = Length[RN];
  al = ConstantArray[0, {numSpe, numRea}];
  be = ConstantArray[0, {numSpe, numRea}];
  Do[
    reactants = comp2Asso[RN[[jj, 1]]];
    products = comp2Asso[RN[[jj, 2]]];
    Do[
      If[KeyExistsQ[reactants, spe[[ii]]], al[[ii, jj]] = reactants[spe[[ii]]]];
      If[KeyExistsQ[products, spe[[ii]]], be[[ii, jj]] = products[spe[[ii]]]];
    , {ii, numSpe}];
  , {jj, numRea}];
  gam = be - al;

  s = If[gam === {} || AllTrue[Flatten[gam], # == 0 &], 0,
    If[AllTrue[Flatten[gam], NumericQ], MatrixRank[gam], "symbolic"]
  ];

  (* Return {Nc, l, s, delta} *)
  If[s === "symbolic", {Nc, l, s, "symbolic"}, {Nc, l, s, Nc - l - s}]
];

extMat[reactions_, speOrder_List:{}] := Module[{
  spe, alp, bet, gamma, Rv, RHS,
  numReactions, numSpecies, reactants, products,
  var, rv, tk, defFormula, defTerms, defResult, Nc, l, s, mSi},

  (* Extract species preserving case, use provided order if given *)
  If[Length[speOrder] > 0,
    spe = ToString[#] & /@ speOrder;
  ,
    spe = extSpe[reactions];
  ];

  (* Check if species were found *)
  If[Length[spe] == 0,
    Print["Error: No species found in reactions"];
    Return[$Failed]
  ];

  numSpecies = Length[spe];
  numReactions = Length[reactions];

  (* Initialize stoichiometric matrices - can contain symbolic entries *)
  alp = ConstantArray[0, {numSpecies, numReactions}];
  bet = ConstantArray[0, {numSpecies, numReactions}];

  (* Build stoichiometric matrices *)
  Do[
    reactants = comp2Asso[reactions[[j, 1]]];
    products = comp2Asso[reactions[[j, 2]]];

    (* Fill matrices - match lowercase species *)
    Do[
      If[KeyExistsQ[reactants, spe[[i]]],
        alp[[i, j]] = reactants[spe[[i]]]];
      If[KeyExistsQ[products, spe[[i]]],
        bet[[i, j]] = products[spe[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];

  (* Calculate derived matrices *)
  gamma = bet - alp;

  (* Convert lowercase species to variables *)
  var = ToExpression /@ spe;

  (* Calculate reaction rates *)
  rv = Table[
    Product[
      If[alp[[i, j]] == 0,
        1,
        (* Handle symbolic exponents *)
        If[NumericQ[alp[[i, j]]],
          var[[i]]^alp[[i, j]],
          Power[var[[i]], al[[i, j]]]
        ]
      ],
      {i, numSpecies}
    ],
    {j, numReactions}
  ];

  (* Rate constants and RHS *)
  tk = Table[Symbol["k" <> ToString[j]], {j, numReactions}];
  Rv = tk*rv;
  RHS = gamma . Rv;

  (* Use defi for deficiency calculation *)
  {Nc, l, s, mSi} = defi[reactions];

  defFormula = "\[Delta] = Nc - \[ScriptL] - s";
  defTerms = "Nc = " <> ToString[Nc] <> " (complexes), " <>
             "\[ScriptL] = " <> ToString[l] <> " (linkage classes), " <>
             "s = " <> ToString[s] <> " (stoich dimension)";
  defResult = If[mSi === "symbolic",
    "\[Delta] = symbolic (requires numeric parameter values)",
    "\[Delta] = " <> ToString[Nc] <> " - " <> ToString[l] <>
      " - " <> ToString[s] <> " = " <> ToString[mSi]
  ];

  (* Return *)
  {spe, alp, bet, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
];



asoRea[RN_]:=Module[{parseSide,extractSpecies},
  extractSpecies[expr_]:=Module[{terms,species},
    If[expr===0,Return[{}]];
    terms=If[Head[expr]===Plus,List@@expr,{expr}];
    species={};
    Do[Which[
      StringQ[term]&&StringContainsQ[term," "],AppendTo[species,StringTrim[StringDrop[term,1]]],
      Head[term]===Times,AppendTo[species,Cases[term,_String][[1]]],
      StringQ[term],AppendTo[species,term]];,{term,terms}];
    DeleteDuplicates[species]];
  parseSide[expr_]:=extractSpecies[expr];
  Map[Function[r,Association["Substrates"->parseSide[r[[1]]],"Products"->parseSide[r[[2]]]]],RN]];



subsCon[str_String] := Module[{base, digits},
  If[StringMatchQ[str, RegularExpression["^([A-Za-z]+)([0-9]+)$"]],
    base = StringCases[str, RegularExpression["^([A-Za-z]+)([0-9]+)$"] -> "$1"][[1]];
    digits = StringCases[str, RegularExpression["^([A-Za-z]+)([0-9]+)$"] -> "$2"][[1]];
    Subscript[ToExpression[base], ToExpression[digits]],
    ToExpression[str]
  ]
];

subsCon[sym_Symbol] := Module[{str},
  str = ToString[sym];
  subsCon[str]
];

subsCon[expr_] := expr;


(* ========================================================================= *)
(* FROM SIPHONS.WL *)
(* ========================================================================= *)

(* Helper: Convert to string preserving case *)
toStr[s_String] := s;
toStr[s_Symbol] := SymbolName[s];
toStr[s_] := SymbolName[s];

pos[v_] := AllTrue[v, # >= 0 &] && Total[v] > 0;

isSiph[W_List, species_List, alpha_, beta_] := Module[
  {indices, alphaW, betaW, nReac, producingReactions},

  (* Convert species names to indices *)
  indices = Flatten[Position[species, #] & /@ W];
  If[Length[indices] != Length[W], Return[False]];

  nReac = Dimensions[alpha][[2]];
  alphaW = alpha[[indices]];
  betaW = beta[[indices]];

  (* Reactions where (S+)_W > 0 *)
  producingReactions = Select[Range[nReac], pos[betaW[[All, #]]] &];

  (* If no producing reactions, vacuously true *)
  If[Length[producingReactions] == 0, Return[True]];

  (* All producing reactions must also have (S-)_W > 0 *)
  AllTrue[producingReactions, pos[alphaW[[All, #]]] &]
];



minSiph[vars_, reactions_, RHS_:Null] := Module[{
  species, alpha, beta, n, m, reactionsAsso, siphons, minimal, nonm,
  cDFEspecies, cDFE, E0, nonSiphonVars, eqsAtDFE, sol},

  species = toStr /@ vars;
  reactionsAsso = asoRea[reactions];

  n = Length[species];
  m = Length[reactions];

  alpha = ConstantArray[0, {n, m}];
  beta = ConstantArray[0, {n, m}];

  Do[
    Module[{subs, prods},
      subs = Lookup[reactionsAsso[[j]], "Substrates", {}];
      prods = Lookup[reactionsAsso[[j]], "Products", {}];
      subs = toStr /@ subs;
      prods = toStr /@ prods;

      Do[
        If[MemberQ[subs, species[[i]]], alpha[[i, j]] = 1];
        If[MemberQ[prods, species[[i]]], beta[[i, j]] = 1];
      , {i, n}];
    ];
  , {j, m}];

  siphons = Select[Subsets[species, {1, n}], isSiph[#, species, alpha, beta] &];

  minimal = Select[siphons,
    Function[s, AllTrue[siphons, Function[t, t === s || !SubsetQ[s, t]]]]
  ];
  nonm = Complement[siphons, minimal];

  (* Compute cDFE and E0 if RHS is provided *)
  If[RHS =!= Null,
    (* cDFE: union of all minimal siphons set to zero *)
    cDFEspecies = DeleteDuplicates[Flatten[minimal]];
    cDFE = Thread[ToExpression[cDFEspecies] -> 0];

    (* cE0: solve RHS==0 under cDFE constraint *)
    nonSiphonVars = Complement[vars, ToExpression[cDFEspecies]];
    If[Length[nonSiphonVars] > 0,
      eqsAtDFE = Thread[(RHS /. cDFE) == 0];
      sol = Quiet[Solve[eqsAtDFE, nonSiphonVars], Solve::svars];
      If[sol === {} || sol === $Failed,
        cE0 = cDFE,
        cE0 = Join[cDFE, First[sol]]
      ],
      cE0 = cDFE
    ];

    {minimal, cDFE, cE0, nonm},

    (* Legacy output if RHS not provided *)
    {minimal, nonm}
  ]
];



comp2Asso[expr_] := Module[{terms, result, coeff, species},
  If[expr === 0 || expr === Null,
    Association[], (* Empty association for zero or null *)
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = Association[];
    Do[
      Which[
        (* Handle Times with coefficient * species *)
        Head[term] === Times && Length[term] == 2,
        (* Determine which part is species (String) and which is coefficient *)
        If[StringQ[First[term]],
          (* First is species (string), second is coefficient *)
          species = First[term];
          coeff = Last[term],
          (* Else: first is coefficient, second is species (string) *)
          coeff = First[term];
          species = Last[term]
        ];
        (* Store species preserving case *)
        If[StringQ[species],
          result[species] = coeff,
          (* Fallback if species is not a string (shouldn't happen) *)
          result[ToString[species]] = coeff
        ],

        (* Handle string species alone (coefficient 1) *)
        StringQ[term],
        result[term] = 1,

        (* Handle symbol species alone (coefficient 1) *)
        AtomQ[term],
        result[ToString[term]] = 1,

        (* Handle other cases - preserve case *)
        True,
        result[ToString[term]] = 1
      ];
    , {term, terms}];
    result
  ]
];


NGM[mod_, infVars_, Fuser_: {}] :=
  Module[{dyn, X, par, inf, infc, Jx, Jy, Jxy, Jyx, V1, F1, F, V, K, Kd, Vinv, allPos, 
  dfeRules},
   dyn = mod[[1]];
   X = mod[[2]];
   par = If[Length[mod] >= 3, mod[[3]], {}];
   inf = Flatten[Position[X, #] & /@ infVars];
   infc = Complement[Range[Length[X]], inf];
   Jx = Grad[dyn[[inf]], X[[inf]]];
   Jy = If[Length[infc] > 0, Grad[dyn[[infc]], X[[infc]]], {}];
   Jxy = If[Length[infc] > 0, Grad[dyn[[inf]], X[[infc]]], {}];
   Jyx = If[Length[infc] > 0, Grad[dyn[[infc]], X[[inf]]], {}];
   V1 = If[Length[infc] > 0, (-Jx) /. Thread[X[[infc]] -> 0], -Jx];
   F1 = (Jx + V1) /. Thread[X[[inf]] -> 0];
   If[Fuser === {},
     F = posM[F1];
     V = F - Jx;,
     V = Fuser - Jx;
     (* Check if V is nonsingular before inverting *)
     Module[{Vdfe = V /. Thread[X[[inf]] -> 0]},
       If[MatrixQ[Vdfe] && MatrixRank[Vdfe] == Length[Vdfe],
         Vinv = Inverse[V] /. Thread[X[[inf]] -> 0];
         allPos = AllTrue[Flatten[Vinv], onlyP];
         If[allPos,
           F = Fuser;
           Print["NGM: user F accepted (Inverse[V] all positive)"],
           Print["NGM: user F rejected (Inverse[V] not all positive), using posM"];
           F = posM[F1];
           V = F - Jx;
         ],
         (* V is singular, use posM *)
         Print["NGM: user F rejected (V singular at DFE), using posM"];
         F = posM[F1];
         V = F - Jx;
       ]
     ];
   ];

   (* Check if V is nonsingular at DFE before computing K *)
   Module[{Vdfe = V /. Thread[X[[inf]] -> 0]},
     If[MatrixQ[Vdfe] && MatrixRank[Vdfe] == Length[Vdfe],
       K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // Simplify;
       Kd = (Inverse[V] . F) /. Thread[X[[inf]] -> 0] // Simplify;
     ,
       (* V is singular - set dummy values *)
       K = 0;
       Kd = 0;
     ]
   ];
   {Jx, F, V, K, Jy, Jxy, Jyx, Kd}
  ]



NGMRN[RN_, rts_, var_] := Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, cv, ct, mS, mSi, inf, mod,
K, R0A, cDFE, RDFE, eq0, var0, cE0, Jx, Jy, F, VV, eigenvals, ngm, infVarsStr, infV},

  (* Extract stoichiometric matrices with var order *)
  {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
  RHS = gam . rts;
  par = Par[RHS, var];
  cp = Thread[par > 0];
  cv = Thread[var >= 0];
  ct = Join[cp, cv];
  mS = minSiph[var, RN, RHS];
  mSi = mS[[1]];
  cDFE = mS[[2]];
  cE0 = mS[[3]];

  (* infV: infection variable symbols (from minimal siphons) *)
  infVarsStr = Union[Flatten[mSi]];
  infV = ToExpression[infVarsStr];

  (* Call NGM with infection variables *)
  mod = {RHS, var, par};
  ngm = NGM[mod, infV];
  Jx = ngm[[1]];
  F = ngm[[2]];
  VV = ngm[[3]];  (* Renamed to avoid conflict with model variable v *)
  K = ngm[[4]];
  Jy = ngm[[5]];

  (* Check if VV is singular (K = 0 from NGM) *)
  If[K === 0,
    (* VV is singular - print informative warning and set dummy values *)
    Module[{Vdfe, rnk, dim},
      Vdfe = VV /. Thread[infV -> 0];
      rnk = MatrixRank[Vdfe];
      dim = Length[Vdfe];
      Print["VV with dim=", dim, " has rank=", rnk, ", so is singular"];
    ];
    R0A = {};
  ,
    (* VV is nonsingular - compute R0A *)
    eigenvals = Eigenvalues[K];
    R0A = Select[eigenvals, # =!= 0 &];
  ];

  (* Return with standard reserved names in signature *)
  V = VV;  (* Map internal VV to output V *)
  {RHS, var, par, cp, mSi, Jx, Jy, cDFE, cE0, F, V, K, R0A, infV}
];

Hur3[pol_, var_] := Module[{co, h3, inec, ineSys, polExpanded, deg},
  (* Expand polynomial *)
  polExpanded = Expand[pol];

  (* Check degree *)
  deg = Exponent[polExpanded, var];
  If[deg != 3,
    Print["Error: Expected degree-3 polynomial, got degree ", deg];
    Return[{}]
  ];

  (* Extract coefficients [a0, a1, a2, a3] using Coefficient *)
  co = {
    Coefficient[polExpanded, var, 0],
    Coefficient[polExpanded, var, 1],
    Coefficient[polExpanded, var, 2],
    Coefficient[polExpanded, var, 3]
  };

  (* Compute Hurwitz determinant h3 *)
  h3 = co[[2]]*co[[3]] - co[[1]]*co[[4]];

  (* Form stability inequalities *)
  inec = {co[[2]] > 0, co[[3]] > 0};
  ineSys = Append[inec, h3 > 0];

  ineSys
];

staPP[pol_, par_] := Module[{
   factors, processedFactors, linearConditions, quadraticConditions,
   higherFactors, linearFactors, quadraticFactors, deg, coeff, const, a, b, c, factor, var, lSta, qSta, cp
   },

  (* Define cp as positivity constraints for all parameters *)
  cp = Thread[par > 0];

  (* Find the variable - it's in Variables[pol] but not in par *)
  var = First[Complement[Variables[pol], par]];
  
  (* Factor the polynomial completely *)
  factors = Factor[pol];
  
  (* Convert to list of factors, handling both single factors and products *)
  processedFactors = If[Head[factors] === Times, List @@ factors, {factors}];
  
  (* Remove constant factors (those not involving var) *)
  processedFactors = Select[processedFactors, !FreeQ[#, var] &];
  
  (* Initialize condition and factor lists *)
  linearConditions = {};
  quadraticConditions = {};
  higherFactors = {};
  linearFactors = {};
  quadraticFactors = {};
  
  (* Process each factor by degree *)
  Do[
   deg = Exponent[factor, var];
   
   Which[
    deg == 1,
    (* Linear factors: a*var + b *)
    coeff = Coefficient[factor, var, 1]; 
    const = Coefficient[factor, var, 0]; 
    AppendTo[linearConditions, const > 0];
    AppendTo[linearFactors, factor],
    
    deg == 2,
    (* Quadratic factors: a*var^2 + b*var + c *)
    a = Coefficient[factor, var, 2];
    b = Coefficient[factor, var, 1]; 
    c = Coefficient[factor, var, 0];
    AppendTo[quadraticConditions, b > 0];
    AppendTo[quadraticConditions, c > 0];
    AppendTo[quadraticFactors, factor],
    
    deg >= 3,
    (* Higher degree factors *)
    AppendTo[higherFactors, factor]
    ],
   
   {factor, processedFactors}
   ];
  
  (* Combine conditions with And *)
  lSta = If[Length[linearConditions] > 0, Apply[And, linearConditions], True];
  qSta = If[Length[quadraticConditions] > 0, Apply[And, quadraticConditions], True];

  (* Simplify using positivity constraints *)
  If[cp =!= {},
    Module[{assumptions, red},
      assumptions = If[Length[cp] > 0, Apply[And, cp], True];
      (* red: remove cp conditions from Reduce output *)
      red[re_, cond_:{}] := re /. (# -> True & /@ cond);
      (* Use Reduce to eliminate factors, then remove cp from output *)
      lSta = red[Reduce[lSta && assumptions, Variables[lSta]], cp];
      qSta = red[Reduce[qSta && assumptions, Variables[qSta]], cp];
    ]
  ];

  (* Return as list {lSta, qSta, hDeg, linearFactors, quadraticFactors} *)
  {lSta, qSta, higherFactors, linearFactors, quadraticFactors}
 ]
cons[gamma_] := Module[{
  positiveConservationLaws, nullSpace, dims},

  dims = Dimensions[gamma];

  (* Find the left nullspace of gamma (vectors w such that w.gamma = 0) *)
  nullSpace = Quiet[NullSpace[Transpose[gamma]]];

  If[nullSpace === {} || nullSpace === $Failed,
    (* No conservation laws *)
    {},
    (* Filter for positive conservation laws *)
    positiveConservationLaws = Select[nullSpace,
      AllTrue[#, # >= 0 &] && AnyTrue[#, # > 0 &] &
    ];

    (* If no positive ones found, try to make them positive *)
    If[Length[positiveConservationLaws] == 0,
      positiveConservationLaws = Select[nullSpace,
        AllTrue[-#, # >= 0 &] && AnyTrue[-#, # > 0 &] &
      ];
      positiveConservationLaws = -# & /@ positiveConservationLaws;
    ];

    positiveConservationLaws
  ]
];

(* TEST posM - Commented out
Print["posM inp={{-a+2b,c-3d},{ef-g,2h}} out=", posM[{{-a + 2*b, c - 3*d}, {e*f - g, 2*h}}], " exp={{2b,c},{ef,2h}}"];
*)

(* ========================================================================= *)
(* EQUILIBRIA CLASSIFICATION *)
(* ========================================================================= *)

sipLat[mSi_, vars_] := Module[{lattice, unions, cDFE, filtered, rules},
  lattice = Subsets[mSi, {1, Length[mSi]}];
  unions = Map[Flatten[#, 1] &, lattice];
  unions = DeleteDuplicates[unions];
  cDFE = Flatten[mSi, 1];
  filtered = DeleteCases[unions, cDFE];
  rules = Map[Module[{sipVars},
    sipVars = Map[If[StringQ[#], ToExpression[#], #] &, #];
    Thread[sipVars -> 0]] &, filtered];
  rules = Reverse[SortBy[rules, Length]];
  rules
]

clsEqb[RHS_, vars_, mSi_, cE0_, prF_, tim_, keepVarUser_:Null] := Module[{
  eqs, sipRules, perSiphon, isDFEsol, isRational, printedRat, isTrulyRational, sameSol, alreadyPrinted,
  mSiSymbols, allSiphonVars, x1, x1Idx
  },
  eqs = Thread[RHS == 0];
  sipRules = sipLat[mSi, vars];

  (* Compute keep variable ONCE for all siphons *)
  mSiSymbols = Map[Map[If[StringQ[#], ToExpression[#], #]&, #]&, mSi];
  allSiphonVars = DeleteDuplicates[Flatten[mSiSymbols, 1]];
  x1 = If[keepVarUser =!= Null && MemberQ[vars, keepVarUser],
    keepVarUser,
    First[Select[vars, !MemberQ[allSiphonVars, #]&]]];
  x1Idx = Position[vars, x1][[1, 1]];

  isDFEsol[s_] := Module[{completed},
    completed = Join[s, Select[cE0, !MemberQ[s[[All, 1]], #[[1]]]&]];
    Sort[completed] === Sort[cE0]];

  isRational[s_] := FreeQ[s, Root | RootOf | AlgebraicNumber] &&
    FreeQ[s, _Power?((!IntegerQ[#[[2]]] && Head[#[[2]]] =!= Symbol)&)];

  isTrulyRational[s_] := isRational[s] && StringLength[ToString[s]] < 5000;
  sameSol[s1_, s2_] := Sort[s1] === Sort[s2];
  alreadyPrinted[s_] := AnyTrue[printedRat, sameSol[#, s]&];
  printedRat = {};

  perSiphon = Table[
    Module[{eqsWithRule, sol, rat, rur, siphonVars, allRat, allRur,
            ratS, rurEq, sipIdx, nBoundary, nIrrat, nDFE, sipName, solRed, nRur},

      siphonVars = rule[[All, 1]];
      sipIdx = Position[sipRules, rule][[1, 1]];
      sipName = "E" <> StringJoin[SymbolName /@ siphonVars];

      nonZat[sol_, varList_] := AllTrue[varList /. sol, # =!= 0 &];

      hasOtherSiphonZero[s_] := Module[{otherSiphons, antiSiph},
        otherSiphons = Select[mSiSymbols, !SubsetQ[siphonVars, #]&];
        antiSiph = DeleteDuplicates[Flatten[otherSiphons, 1]];
        !nonZat[s, antiSiph]];

      eqsWithRule = eqs /. rule;
      sol = MemoryConstrained[TimeConstrained[Quiet[Solve[eqsWithRule, vars], Solve::svars], tim, "limit"], 1024^3, "limit"];

      If[sol === "limit",
        Print["S", sipIdx, ":", sipName, ": timeout"];
        {siphonVars, 0, 0, "limit"},

        allRat = Select[sol, isRational];
        allRur = Map[Join[#, rule]&, Select[sol, !isRational[#]&]];

        ratS = {}; rurEq = {};
        If[Length[allRur] > 0,
          solRed = TimeConstrained[Quiet[Solve[Drop[eqsWithRule, x1Idx], Delete[vars, x1Idx]], Solve::svars], tim, "limit"];
          If[solRed =!= "limit" && solRed =!= {} && solRed =!= $Aborted,
            Module[{antiSiph, solComplete},
              antiSiph = Complement[Delete[vars, x1Idx], siphonVars];
              Do[solComplete = Join[solRed[[i]], rule];
                If[isTrulyRational[solRed[[i]]] && !alreadyPrinted[solRed[[i]]] && nonZat[solComplete, antiSiph],
                  ratS = solComplete;
                  rurEq = Collect[Numerator[Together[eqsWithRule[[x1Idx, 1]] /. solRed[[i]]]], x1, Factor];
                  Break[]],
                {i, Length[solRed]}]]]];

        rat = {};
        Do[Module[{solComplete},
            solComplete = Join[allRat[[i]], rule];
            If[!isDFEsol[allRat[[i]]] && !alreadyPrinted[allRat[[i]]] && !hasOtherSiphonZero[solComplete],
              AppendTo[rat, solComplete]]],
          {i, Length[allRat]}];
        rur = If[Length[ratS] > 0, {ratS, rurEq}, {}];

        nBoundary = Length[rat];
        nDFE = Count[allRat, _?(isDFEsol)];
        nIrrat = Length[allRur];

        If[nIrrat > 0 && Length[rur] == 0 && Length[allRur] > 0,
          Module[{irratSol, irratWithoutX1, poly},
            irratSol = allRur[[1]];
            irratWithoutX1 = Select[irratSol, #[[1]] =!= x1 &];
            poly = Collect[Numerator[Together[eqsWithRule[[x1Idx, 1]] /. irratWithoutX1]], x1, Factor];
            If[FreeQ[poly, Root | RootOf | AlgebraicNumber],
              ratS = irratWithoutX1; rurEq = poly; rur = {ratS, rurEq}]]];

        nRur = If[Length[rur] > 0, Ceiling[nIrrat/2], 0];
        If[nRur > 0, nIrrat = 0];

        (* Concise print: S#:sipName=solution *)
        If[Length[rat] > 0,
          Module[{short},
            short = Select[rat, StringLength[ToString[#]] < 3000 &];
            Do[Print["S", sipIdx, ":", sipName, "=", short[[j]] // prF], {j, Length[short]}];
            printedRat = Join[printedRat, rat]]];

        If[Length[rur] > 0,
          Print["S", sipIdx, ":", If[isRational[ratS], "RUR ", "nonRUR "], sipName, ": ratS=", ratS // prF, " poly=", rurEq // prF]];

        {nBoundary + nDFE, nRur, {rat, rur, allRur}}
      ]
    ],
    {rule, sipRules}];
  {perSiphon, Length[printedRat]}
]

clsEq[RHS_, vars_, mSi_, cE0_, prF_, tim_, keepVarUser_:Null] := Module[{
  eqs, sol, rat, irrat, perSiphon, isDFEsol, isTrulyRational, nonZeroSol, filtered, nSol, nS, emptySiphonResult,
  x1, x1Idx, solRed, rurEq, rurOK, nRat, nRUR, nNonRUR, eeNonRUR, clsEqbResult, nRatSiphon
  },
  clsEqbResult = clsEqb[RHS, vars, mSi, cE0, prF, tim, keepVarUser];
  perSiphon = clsEqbResult[[1]];
  nRatSiphon = clsEqbResult[[2]];
  nS = 2^Length[mSi] - 1;
  eqs = Thread[RHS == 0];
  eeNonRUR = 0;

  isDFEsol[s_] := Module[{completed},
    completed = Join[s, Select[cE0, !MemberQ[s[[All, 1]], #[[1]]]&]];
    Sort[completed] === Sort[cE0]];

  isTrulyRational[s_] := FreeQ[s, Root | RootOf | AlgebraicNumber] &&
    FreeQ[s, _Power?((!IntegerQ[#[[2]]] && Head[#[[2]]] =!= Symbol)&)] &&
    StringLength[ToString[s]] < 5000;

  hasSiphonZero[s_] := Module[{mSiSymbols},
    mSiSymbols = Map[Map[If[StringQ[#], ToExpression[#], #]&, #]&, mSi];
    AnyTrue[mSiSymbols, AllTrue[# /. s, # === 0 &]&]];

  sol = TimeConstrained[Quiet[Solve[eqs, vars], Solve::svars], tim, "timeout"];

  If[sol === "timeout",
    Print["EE: timeout"];
    emptySiphonResult = {0, 0, {{}, {}, {}}},

    rat = Select[sol, isTrulyRational];
    irrat = Select[sol, !isTrulyRational[#]&];
    filtered = Select[rat, !isDFEsol[#] && !hasSiphonZero[#]&];
    nonZeroSol = Select[filtered, !AllTrue[vars /. #, # === 0 &]&];
    nSol = Length[nonZeroSol];

    If[nSol > 0,
      If[nSol == 1,
        Print["EE=", prF[nonZeroSol[[1]]]],
        Do[Print["EE", i, "=", prF[nonZeroSol[[i]]]], {i, nSol}]];
      emptySiphonResult = {nSol, 0, {nonZeroSol, {}, irrat}},

      (* No rational endemic - check for irrational *)
      If[Length[irrat] > 0,
        Module[{nBefore = Length[irrat]},
        (* Filter irrational: no variable is zero *)
        irrat = Select[irrat, !MemberQ[#[[All, 2]], 0]&];
        If[Length[irrat] == 0,
          Print["EE: ", nBefore, " irrational but all have zeros"];
          emptySiphonResult = {0, 0, {{}, {}, {}}},

          (* Try RUR: keep x1 (non-siphon var or user choice), drop its eq, solve rest *)
          Module[{mSiSymbols, allSiphonVars, nonSiphonVars},
            mSiSymbols = Map[Map[If[StringQ[#], ToExpression[#], #]&, #]&, mSi];
            allSiphonVars = DeleteDuplicates[Flatten[mSiSymbols, 1]];
            nonSiphonVars = Select[vars, !MemberQ[allSiphonVars, #]&];
            x1 = If[keepVarUser =!= Null && MemberQ[vars, keepVarUser],
              keepVarUser,
              If[Length[nonSiphonVars] > 0, First[nonSiphonVars], vars[[1]]]];
            x1Idx = Position[vars, x1][[1, 1]]];
          rurOK = False;
          solRed = TimeConstrained[Quiet[Solve[Drop[eqs, {x1Idx}], Delete[vars, x1Idx]], Solve::svars], tim, "timeout"];
          If[solRed =!= "timeout" && solRed =!= {},
            Module[{ratRed, ratNoZero, poly},
              ratRed = Select[solRed, isTrulyRational];
              ratNoZero = Select[ratRed, !MemberQ[#[[All, 2]], 0]&];
              If[Length[ratNoZero] > 0,
                poly = Numerator[Together[eqs[[x1Idx, 1]] /. ratNoZero[[1]]]];
                poly = Collect[poly, x1, Factor];
                If[FreeQ[poly, Root | RootOf | AlgebraicNumber],
                  rurOK = True;
                  rurEq = poly;
                  Print["EE:RUR in ", x1, ": ratS=", prF[ratNoZero[[1]]], " poly=", prF[rurEq]]]]]];
          If[!rurOK,
            eeNonRUR = Length[irrat];
            Print["EE:nonRUR ", Length[irrat], " irrational"]];
          emptySiphonResult = {0, If[rurOK, Ceiling[Length[irrat]/2], 0], {{}, {}, irrat}}]],

        emptySiphonResult = {0, 0, {{}, {}, {}}}]]];

  perSiphon = Append[perSiphon, emptySiphonResult];

  (* Final summary: r=rational (unique), R=RUR, n=nonRUR *)
  nRat = nRatSiphon + nSol + 1;  (* unique siphon rational + EE rational + DFE *)
  nRUR = Total[perSiphon[[All, 2]]];
  (* nonRUR = irrational solutions where RUR failed (nRur=0 but has irrat) *)
  nNonRUR = Total[Map[If[#[[2]] == 0, Length[#[[3, 3]]], 0]&, Most[perSiphon]]] + eeNonRUR;
  Print["r=", nRat, " R=", nRUR, " n=", nNonRUR];

  {nS, perSiphon}
]

(* TEST - Commented out
RHS1 = {-be1 i12 s - be2 i12 s - mu0 s - (i1 mu1 R1 s)/s0 - (i2 mu2 R2 s)/s0 - (i12 mu3 R3 s)/s0 + mu0 s0,
        -et1 i1 i12 - ga1 i1 i2 - i1 mu1 + be1 i12 s + (i1 mu1 R1 s)/s0,
        -ga2 i1 i2 - et2 i12 i2 - i2 mu2 + be2 i12 s + (i2 mu2 R2 s)/s0,
        et1 i1 i12 + ga1 i1 i2 + ga2 i1 i2 + et2 i12 i2 - i12 mu3 + (i12 mu3 R3 s)/s0};
var = {s, i1, i2, i12};
{RN, rts, spe, alp, bet, gam} = ODE2RN[RHS1, var];
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars} = bdAn[RN, rts, var];
RHS - RHS1 // FullSimplify  (*Should be {0,0,0,0}*)
*)

(* TEST ODE2RNr with rational rates
RHSrat = {1 - x1 - (x1*x2)/(1 + x1), (x1*x2)/(1 + x1) - x2};
varRat = {x1, x2};
{RNrat, rtsRat, speRat, alpRat, betRat, gamRat, rnRedRat} =
  Block[{Print = (Null &)}, ODE2RNr[RHSrat, varRat, Identity]];
Print["ODE2RNr verif=", (gamRat . rtsRat - RHSrat) // Simplify, " exp={0,0}"];*)

(* TEST cons - find conservation laws
gamTest = {{1, -1, 0}, {0, 1, -1}, {-1, 0, 1}};
consLaws = cons[gamTest];
Print["cons gamTest=", gamTest, " consLaws=", consLaws, " verif=", consLaws.gamTest];*)

