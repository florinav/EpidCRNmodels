(* ::Package:: *)

(* ========================================================================= *)
(* CORE *)
(* ========================================================================= *)
(* MINIMAL EPIDEMIC MODELING SUBPACKAGE *)
(* Contains the essential functions for epidemic CRN modeling *)
(* IMPORTANT: This conains ODE2RN , extMat, minSiph, ... *)
(* Goal: Serve as foundation for future Python port *)
(* NOTE: Backup versions exist in Core1.wl and EpidCRNo.wl *)

cSo[RHS_, var_] := Module[
  {so, idx, irrPows, rootDegs, surdDegs, powDegs, algDeg, allNonzeroQ,
   ratIdx, sqrtIdx, deg3Idx, degGT3Idx, nzIdx, summary},

  so = Quiet[Solve[Thread[RHS == 0], var]];
  idx = Range[Length[so]];

  irrPows[sol_] := Cases[sol, Power[_, r_Rational] /; ! IntegerQ[r] :> r, Infinity];

  rootDegs[sol_] := Cases[
    sol,
    Root[f_, ___] :> Module[{t = Unique["t"]},
      Quiet @ Check[Exponent[f[t], t], Infinity]
    ],
    Infinity
  ];

  surdDegs[sol_] := Cases[sol, Surd[_, n_Integer] :> n, Infinity];

  powDegs[sol_] := Denominator /@ irrPows[sol];  (* q in (..)^(p/q); includes 2 for Sqrt *)

  algDeg[sol_] := Module[{d = Join[rootDegs[sol], surdDegs[sol], powDegs[sol]]},
    If[d === {}, 1, Max @ Replace[d, Infinity -> 10^9, {1}]]
  ];

  allNonzeroQ[sol_] := FreeQ[(var /. sol), 0];

  ratIdx    = Pick[idx, algDeg /@ so, 1];
  sqrtIdx   = Pick[idx, algDeg /@ so, 2];
  deg3Idx   = Pick[idx, algDeg /@ so, 3];
  degGT3Idx = Pick[idx, UnitStep[(algDeg /@ so) - 4], 1];
  nzIdx     = Pick[idx, allNonzeroQ /@ so, True];

  summary = {
    {"Total", Length[so]},
    {"Rational", Length@ratIdx, ratIdx},
    {"SquareRoots", Length@sqrtIdx, sqrtIdx},
    {"Degree3", Length@deg3Idx, deg3Idx},
    {"Degree>3", Length@degGT3Idx, degGT3Idx},
    {"AllNonzero", Length@nzIdx, nzIdx}
  };

  Join[summary, {{"Solutions", so}}]
];

cSol[RHS_, var_] := cSo[RHS, var]
cSol[RHS_, var_, pf_] := Module[{res, so},
  res = cSo[RHS, var];
  so = res[[-1, 2]];
  Append[Most[res], {"Solutions", pf /@ so}]
]

(*RHS={la-(mx+be1 v1) x,-v1 (mv+be1 x)+k1 y1,be1 v1 x-y1 (my+de w+ga1 z1),-mz z1+ga1 xiz y1 z1,-mw w+de w xiw y1};var={x,v1,y1,z1,w};
cSol[RHS,var] yields 4 quadratic irrational, etc*)
(*chkSol[sol_] := Module[
  {p, rootsurd, non2},
  p = Cases[sol, Power[_, r_Rational] /; ! IntegerQ[r] :> r, Infinity];
  rootsurd = Not@FreeQ[sol, Root | Surd];
  non2 = Select[p, Denominator[#] =!= 2 &];
  {
    "hasRootOrSurd" -> rootsurd,
    "allNonIntegerRationalPowers" -> p,
    "powersWithDenomNot2" -> non2,
    "containsSqrtExplicitly" -> Not@FreeQ[sol, Sqrt]
  }
];*)
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
               "inf", "infc", "mod", "par", "cp","I","E","rt"};
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
      Prepend[Transpose[{prF[RN], prF[Factor /@ rts]}], {"RN", "rts"}] // MatrixForm,
      Spacer[20],
      Prepend[Transpose[{prF[rnRed], prF[Factor /@ rts]}], {"RNr", "rts"}] // MatrixForm
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
      Prepend[Transpose[{prF[RN], prF[Factor /@ rtsTrue]}], {"RN", "rts"}] // MatrixForm,
      Spacer[20],
      Prepend[Transpose[{prF[rnRed], prF[Factor /@ rtsTrue]}], {"RNr", "rts"}] // MatrixForm
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

comp2Asso[expr_] := Module[{terms, result, coeff, species, stringParts, otherParts},
  If[expr === 0 || expr === Null,
    Association[],
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = Association[];
    Do[
      Which[
        (* Handle Times with any number of elements - extracts string as species *)
        Head[term] === Times,
        stringParts = Select[List @@ term, StringQ];
        otherParts = Select[List @@ term, !StringQ[#] &];
        If[Length[stringParts] == 1,
          species = First[stringParts];
          coeff = If[Length[otherParts] == 0, 1,
                    If[Length[otherParts] == 1, First[otherParts], Times @@ otherParts]];
          result[species] = coeff,
          result[ToString[term]] = 1
        ],

        (* Handle string species alone (coefficient 1) *)
        StringQ[term],
        result[term] = 1,

        (* Handle symbol species alone (coefficient 1) *)
        AtomQ[term],
        result[ToString[term]] = 1,

        (* Handle other cases *)
        True,
        result[ToString[term]] = 1
      ];
    , {term, terms}];
    result
  ]
];
NGM[mod_, infVars_, Fuser_: {}] :=
  Module[{dyn, X, par, inf, infc, Jx, Jy, Jxy, fx, V1, F1, F, V, K, Kd, Vinv, allPos,
  dfeRules, n, JxDFE, Fdiag},
   dyn = mod[[1]];
   X = mod[[2]];
   par = If[Length[mod] >= 3, mod[[3]], {}];
   inf = Flatten[Position[X, #] & /@ infVars];
   infc = Complement[Range[Length[X]], inf];
   Jx = Grad[dyn[[inf]], X[[inf]]];
   Jy = If[Length[infc] > 0, Grad[dyn[[infc]], X[[infc]]], {}];
   Jxy = If[Length[infc] > 0, Grad[dyn[[inf]], X[[infc]]], {}];
   fx = dyn[[inf]] /. Thread[X[[inf]] -> 0];  (* RHS for inf vars, evaluated on siphon facet *)
   V1 = If[Length[infc] > 0, (-Jx) /. Thread[X[[infc]] -> 0], -Jx];
   F1 = (Jx + V1) /. Thread[X[[inf]] -> 0];
   (* Extract positive diagonal terms from Jx at DFE *)
   n = Length[inf];
   JxDFE = Jx /. Thread[X[[inf]] -> 0];
   Fdiag = DiagonalMatrix[Table[posM[{{JxDFE[[i, i]]}}][[1, 1]], {i, n}]];
   If[Fuser === {},
     F = posM[F1] + Fdiag;
     V = F - Jx;,
     V = Fuser - Jx;
     (* Check if V is nonsingular - accept user F if Inverse[V] exists *)
     Module[{Vdfe = V /. Thread[X[[inf]] -> 0]},
       If[MatrixQ[Vdfe] && MatrixRank[Vdfe] == Length[Vdfe],
         F = Fuser;
         Print["NGM: user F accepted (V nonsingular)"],
         (* V is singular, use posM *)
         Print["NGM: user F rejected (V singular at DFE), using posM"];
         F = posM[F1] + Fdiag;
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
   {Jx, F, V, K, Jy, Jxy, fx, Kd}
  ]



NGMRN[RN_, rts_, var_] := Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, cv, ct, mS, mSi, inf, mod,
K, R0A, cDFE, RDFE, eq0, var0, cE0, Jx, Jy, Jxy, fx, F, VV, eigenvals, ngm, infVarsStr, infV},

  (* Extract stoichiometric matrices with var order *)
  {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
  RHS = gam . rts;  (* rts = full rates (monomials * constants) *)
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
  Jxy = ngm[[6]];
  fx = ngm[[7]];

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
  {RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, F, V, K, R0A, infV}
];

Hur3[pol_, var_] := Hur3[pol, var, 5];
Hur3[pol_, var_, tmo_] := Module[{co, h3, inec, ineSys, polExpanded, deg, par, cp, red},
  polExpanded = Expand[pol];
  deg = Exponent[polExpanded, var];
  If[deg != 3, Print["Hur3: deg ", deg, " != 3"]; Return[{}]];
  co = Table[Coefficient[polExpanded, var, j], {j, 0, 3}];
  h3 = co[[2]]*co[[3]] - co[[1]]*co[[4]];
  inec = {co[[2]] > 0, co[[3]] > 0};
  ineSys = Append[inec, h3 > 0];
  par = Complement[Variables[pol], {var}];
  cp = Thread[par > 0];
  red = TimeConstrained[
    FullSimplify[Reduce[FullSimplify[Join[cp, ineSys]]], cp], tmo, $Failed];
  {red, ineSys}
];

Hur4[pol_, var_] := Hur4[pol, var, 5];
Hur4[pol_, var_, tmo_] := Module[{co, H4, h4, inec, ineSys, polExpanded, deg, par, cp, red},
  polExpanded = Expand[pol];
  deg = Exponent[polExpanded, var];
  co = Table[Coefficient[polExpanded, var, j], {j, 0, deg}];
  H4 = {
    {co[[4]], co[[2]], 0,       0},
    {co[[5]], co[[3]], co[[1]], 0},
    {0,       co[[4]], co[[2]], 0},
    {0,       co[[5]], co[[3]], co[[1]]}
  };
  h4 = Det[H4];
  inec = {co[[1]] > 0, co[[2]] > 0, co[[3]] > 0};
  ineSys = Append[inec, h4 > 0];
  par = Complement[Variables[pol], {var}];
  cp = Thread[par > 0];
  red = TimeConstrained[
    FullSimplify[Reduce[FullSimplify[Join[cp, ineSys]]], cp], tmo, $Failed];
  {red, ineSys}
];

Hur4M[mat_] := Module[{lm, ch, cot, co, H4, h4, ine},
  lm = mat // Length;
  ch = ((-1)^lm*CharacteristicPolynomial[mat, \[Lambda]] // Factor);
  cot = CoefficientList[ch, \[Lambda]];
  co = Reverse[Drop[cot, -1]];
  H4 = {{co[[1]], 1, 0, 0},
        {co[[3]], co[[2]], co[[1]], 1},
        {0, co[[4]], co[[3]], co[[2]]},
        {0, 0, 0, co[[4]]}};
  h4 = Det[H4];
  ine = Thread[co > 0];
  {co, h4, ine}
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
    AppendTo[quadraticConditions,a b > 0];
    AppendTo[quadraticConditions,a c > 0];
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
(* cons moved to Core.wl *)

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
  mSiSymbols, allSiphonVars, x1, x1Idx, ratSiphons, rurSiphons, nonRurSiphons, printQueue
  },
  ratSiphons = {}; rurSiphons = {}; nonRurSiphons = {}; printQueue = {};
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

      hasOtherSiphonZero[s_] := Module[{otherSiphons},
        otherSiphons = Select[mSiSymbols, !SubsetQ[siphonVars, #]&];
        AnyTrue[otherSiphons, AllTrue[# /. s, # === 0 &]&]];

      eqsWithRule = eqs /. rule;
      sol = MemoryConstrained[TimeConstrained[Quiet[Solve[eqsWithRule, vars], Solve::svars], tim, "limit"], 1024^3, "limit"];

      If[sol === "limit",
        AppendTo[printQueue, {"timeout", sipIdx, sipName}];
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
        Do[Module[{solComplete, isDFE, alrPr, hasOther, solVars, ruleVarsNotInSol},
            solVars = allRat[[i]][[All, 1]];
            ruleVarsNotInSol = Select[rule, !MemberQ[solVars, #[[1]]]&];
            solComplete = Join[allRat[[i]], ruleVarsNotInSol];
            isDFE = isDFEsol[allRat[[i]]];
            alrPr = alreadyPrinted[solComplete];
            hasOther = hasOtherSiphonZero[solComplete];
            If[!isDFE && !alrPr && !hasOther,
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

        (* Collect print: S#:sipName=solution *)
        If[Length[rat] > 0,
          AppendTo[ratSiphons, sipIdx];
          Module[{short},
            short = Select[rat, StringLength[ToString[#]] < 3000 &];
            Do[AppendTo[printQueue, {"rat", sipIdx, sipName, short[[j]]}], {j, Length[short]}];
            printedRat = Join[printedRat, rat]]];

        If[Length[rur] > 0,
          If[isRational[ratS],
            AppendTo[rurSiphons, sipIdx];
            AppendTo[printQueue, {"rur", sipIdx, sipName, ratS, rurEq}],
            AppendTo[nonRurSiphons, {sipIdx, StringDrop[sipName, 1]}]]];

        {nBoundary + nDFE, nRur, {rat, rur, allRur}}
      ]
    ],
    {rule, sipRules}];
  {perSiphon, ratSiphons, rurSiphons, nonRurSiphons, printQueue}
]

(* clsEqR: Filter clsEqb output to exclude siphons containing S1 or S2 *)
clsEqR[RHS_, vars_, mSi_, cE0_, prF_, tim_, S1_, S2_:{}] := Module[{
  S1syms, S2syms, clsEqbResult, perSiphon, ratSiphons, rurSiphons, nonRurSiphons, printQueue,
  sipRules, containsExcluded, filteredQueue, filteredRat, filteredRur
  },
  Print["clsEqR: RUR poly res[[i,3,2,2]]; rat of Si res[[i,3,1]]"];
  (* Convert S1, S2 to symbols *)
  S1syms = Map[If[StringQ[#], ToExpression[#], #]&, S1];
  S2syms = Map[If[StringQ[#], ToExpression[#], #]&, S2];

  (* Call clsEqb with full mSi *)
  clsEqbResult = clsEqb[RHS, vars, mSi, cE0, prF, tim];
  perSiphon = clsEqbResult[[1]];
  ratSiphons = clsEqbResult[[2]];
  rurSiphons = clsEqbResult[[3]];
  nonRurSiphons = clsEqbResult[[4]];
  printQueue = clsEqbResult[[5]];

  (* Get sipRules to check which siphons to exclude *)
  sipRules = sipLat[mSi, vars];

  (* Check if a siphon (by index) contains ALL vars from S1 or ALL vars from S2 *)
  containsExcluded[sipIdx_] := Module[{siphonVars},
    siphonVars = sipRules[[sipIdx]][[All, 1]];
    SubsetQ[siphonVars, S1syms] || (S2syms =!= {} && SubsetQ[siphonVars, S2syms])];

  (* Filter printQueue *)
  filteredQueue = Select[printQueue, !containsExcluded[#[[2]]]&];
  filteredRat = Select[ratSiphons, !containsExcluded[#]&];
  filteredRur = Select[rurSiphons, !containsExcluded[#]&];

  (* Print filtered summary and details *)
  Print["r=", Length[filteredRat], " ", filteredRat, " R=", Length[filteredRur], " ", filteredRur];
  Do[Module[{msg = filteredQueue[[i]]},
    Switch[msg[[1]],
      "timeout", Print["S", msg[[2]], ":", msg[[3]], ": timeout"],
      "rat", Print["S", msg[[2]], ":", msg[[3]], "=", msg[[4]] // prF],
      "rur", Print["S", msg[[2]], ":RUR ", msg[[3]], ": ratS=", msg[[4]] // prF, " poly=", msg[[5]] // prF]
    ]], {i, Length[filteredQueue]}];

  {Length[sipRules], perSiphon}
]

clsEq[RHS_, vars_, mSi_, cE0_, prF_, tim_, keepVarUser_:Null] := Module[{
  eqs, sol, rat, irrat, perSiphon, isDFEsol, isTrulyRational, nonZeroSol, filtered, nSol, nS, emptySiphonResult,
  x1, x1Idx, solRed, rurEq, rurOK, nRat, nRUR, nNonRUR, eeNonRUR, clsEqbResult, ratSiphons, rurSiphons, nonRurSiphons, printQueue, eeQueue
  },
  Print["clsEq: use res[[i,3,2,2]] for RUR poly of siphon i, res[[i,3,1]] for rat solutions of siphon i"];
  clsEqbResult = clsEqb[RHS, vars, mSi, cE0, prF, tim, keepVarUser];
  perSiphon = clsEqbResult[[1]];
  ratSiphons = clsEqbResult[[2]];
  rurSiphons = clsEqbResult[[3]];
  nonRurSiphons = clsEqbResult[[4]];
  printQueue = clsEqbResult[[5]];
  nS = 2^Length[mSi] - 1;
  eqs = Thread[RHS == 0];
  eeNonRUR = 0;
  eeQueue = {};
  emptySiphonResult = Null;

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
    AppendTo[eeQueue, {"EE: timeout"}];
    emptySiphonResult = {0, 0, {{}, {}, {}}},

    rat = Select[sol, isTrulyRational];
    irrat = Select[sol, !isTrulyRational[#]&];
    filtered = Select[rat, !isDFEsol[#] && !hasSiphonZero[#]&];
    nonZeroSol = Select[filtered, !AllTrue[vars /. #, # === 0 &]&];
    nSol = Length[nonZeroSol];

    If[nSol > 0,
      If[nSol == 1,
        AppendTo[eeQueue, {"ee", nonZeroSol[[1]]}],
        Do[AppendTo[eeQueue, {"ee", i, nonZeroSol[[i]]}], {i, nSol}]];
      emptySiphonResult = {nSol, 0, {nonZeroSol, {}, irrat}},

      (* No rational endemic - check for irrational *)
      If[Length[irrat] > 0,
        Module[{nBefore = Length[irrat]},
        (* Filter irrational: no variable is zero *)
        irrat = Select[irrat, !MemberQ[#[[All, 2]], 0]&];
        If[Length[irrat] == 0,
          AppendTo[eeQueue, {"noEE", Length[sol]}];
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
            Module[{ratRed, ratNoZero, poly, deg, x1val, fullSol},
              ratRed = Select[solRed, isTrulyRational];
              ratNoZero = Select[ratRed, !MemberQ[#[[All, 2]], 0]&];
              If[Length[ratNoZero] > 0,
                poly = Numerator[Together[eqs[[x1Idx, 1]] /. ratNoZero[[1]]]];
                poly = Collect[poly, x1, Factor];
                deg = Exponent[poly, x1];
                If[FreeQ[poly, Root | RootOf | AlgebraicNumber],
                  If[deg == 1,
                    (* Linear poly - solve directly as rational, substitute x1 into other vars *)
                    x1val = Simplify[Solve[poly == 0, x1][[1, 1, 2]]];
                    fullSol = Join[{x1 -> x1val}, Simplify[ratNoZero[[1]] /. x1 -> x1val]];
                    nonZeroSol = {fullSol};
                    nSol = 1;
                    AppendTo[eeQueue, {"ee", fullSol}];
                    emptySiphonResult = {1, 0, {{fullSol}, {}, irrat}},
                    (* Higher degree - treat as RUR *)
                    rurOK = True;
                    rurEq = poly;
                    AppendTo[eeQueue, {"eerur", x1, ratNoZero[[1]], rurEq}]]]]]];
          If[!rurOK && emptySiphonResult === Null,
            eeNonRUR = Length[irrat];
            AppendTo[eeQueue, {"eeNonRUR", Length[irrat]}];
            emptySiphonResult = {0, 0, {{}, {}, irrat}}];
          If[rurOK, emptySiphonResult = {0, Ceiling[Length[irrat]/2], {{}, {}, irrat}}]]],

        emptySiphonResult = {0, 0, {{}, {}, {}}}]]];

  perSiphon = Append[perSiphon, emptySiphonResult];

  (* Final summary: rat=siphons with rational, RUR=siphons with RUR poly, nonRUR=siphons with irrational *)
  nRUR = Length[rurSiphons];
  nNonRUR = Length[nonRurSiphons] + If[eeNonRUR > 0, 1, 0];

  (* Build rat string with counts: 4(2) means siphon 4 has 2 sols, just 4 if 1 sol *)
  Module[{ratStr, nonRurStr, totRat, allRatSip, eeRatCnt},
    (* Include empty siphon if it has rational solutions *)
    eeRatCnt = Length[emptySiphonResult[[3, 1]]];
    allRatSip = If[eeRatCnt > 0, Append[ratSiphons, nS + 1], ratSiphons];
    nRat = Length[allRatSip];

    (* Build string with counts *)
    ratStr = StringJoin[Riffle[Map[Module[{idx = #, cnt},
      cnt = If[idx <= nS, Length[perSiphon[[idx, 3, 1]]], eeRatCnt];
      If[cnt > 1, ToString[idx] <> "(" <> ToString[cnt] <> ")", ToString[idx]]
    ]&, allRatSip], ","]];

    (* Total rational solutions *)
    totRat = Total[Map[Length[perSiphon[[#, 3, 1]]]&, ratSiphons]] + eeRatCnt;

    nonRurStr = StringJoin[Riffle[Map[ToString[#[[1]]] <> ":" <> #[[2]] &, nonRurSiphons], ","]];
    If[eeNonRUR > 0, nonRurStr = If[nonRurStr == "", "EE", nonRurStr <> ",EE"]];
    Print[totRat, " rat sol on siphons {", ratStr, "}, ", nRUR, " RUR on siphons ", rurSiphons, ", ", nNonRUR, " nonRUR on siphons {", nonRurStr, "}"]];

  (* Print queued siphon messages *)
  Do[Module[{msg = printQueue[[i]]},
    Switch[msg[[1]],
      "timeout", Print["S", msg[[2]], ":", msg[[3]], ": timeout"],
      "rat", Print["S", msg[[2]], ":", msg[[3]], "=", msg[[4]] // prF],
      "rur", Print["S", msg[[2]], ":RUR ", msg[[3]], ": ratS=", msg[[4]] // prF, " poly=", msg[[5]] // prF]
    ]], {i, Length[printQueue]}];

  (* Print queued EE messages *)
  Do[Module[{msg = eeQueue[[i]]},
    Switch[msg[[1]],
      "EE: timeout", Print["EE: timeout"],
      "ee", If[Length[msg] == 2, Print["EE=", prF[msg[[2]]]], Print["EE", msg[[2]], "=", prF[msg[[3]]]]],
      "noEE", Print["no EE (all ", msg[[2]], " solutions found by Solve have zeros)"],
      "eerur", Print["EE:RUR in ", msg[[2]], ": ratS=", prF[msg[[3]]], " poly=", prF[msg[[4]]]],
      "eeNonRUR", Print["EE:nonRUR ", msg[[2]], " irrational"]
    ]], {i, Length[eeQueue]}];

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



(* ========================================================================== *)
(* SYMBOL TO STRING CONVERSION *)
(* ========================================================================== *)

(* Helper: Convert to string preserving case - handles both symbols and strings *)
toStr[s_String] := s;
toStr[s_Symbol] := SymbolName[s];
toStr[s_] := SymbolName[s]; (* fallback for other types *)

symbToStr[complex_] := 
  Module[{}, 
   Which[complex === 0, "0", Head[complex] === Plus, 
    StringJoin[
     Riffle[Map[
       Which[Head[#] === Times, 
         Module[{parts, strings, nonStrings}, parts = List @@ #;
          strings = Select[parts, StringQ];
          nonStrings = Select[parts, ! StringQ[#] &];
          
          If[Length[nonStrings] == 1 && nonStrings[[1]] == 1, 
           strings[[1]], 
           ToString[nonStrings[[1]]] <> "*" <> strings[[1]]]], 
         StringQ[#], #, True, ToString[#]] &, List @@ complex], 
      " + "]], Head[complex] === Times, 
    Module[{parts, strings, nonStrings}, parts = List @@ complex;
     strings = Select[parts, StringQ];
     nonStrings = Select[parts, ! StringQ[#] &];
     If[Length[nonStrings] == 1 && nonStrings[[1]] == 1, strings[[1]],
       ToString[nonStrings[[1]]] <> "*" <> strings[[1]]]], 
    StringQ[complex], complex, True, ToString[complex]]];

(* Test symbToStr *)
(*
symbToStr[2*"X" + 3*"Y"]  (* Should return "2*X + 3*Y" *)
symbToStr["X"*"Y"]        (* Should return "X*Y" *)
symbToStr["X"]            (* Should return "X" *)
*)

(* Enhanced string to symbolic conversion *)
strToSymb[reactions_] := Module[{convertSide},
  convertSide[side_String] := Module[{result, terms, converted},
    If[side == "0" || side == "", Return[0]];
    terms = StringSplit[side, "+"];
    converted = Map[Module[{term, cleanTerm, coeff, species},
      term = StringTrim[#];
      Which[
        StringMatchQ[term, RegularExpression["^\\d+[A-Za-z0-9]+$"]],
        Module[{digits, letters},
          digits = StringTake[term, 1];
          letters = StringDrop[term, 1];
          ToExpression[digits]*letters
        ],
        StringMatchQ[term, RegularExpression["^[A-Za-z0-9]+$"]],
        term,
        True, term
      ]
    ] &, terms];
    If[Length[converted] == 1, converted[[1]], Apply[Plus, converted]]
  ];
  Map[Function[reaction,
    If[Head[reaction] === Rule,
      convertSide[First[reaction]] -> convertSide[Last[reaction]],
      {convertSide[First[reaction]], convertSide[Last[reaction]]}
    ]
  ], reactions]
];



cons[gamma_, cp_: {}] := Module[{
  nullSpace, positiveConservationLaws, constraints},
  
  (* Find the left nullspace of gamma (vectors w such that w.gamma = 0) *)
  nullSpace = Quiet[NullSpace[Transpose[gamma]]];
  
  If[nullSpace === {} || nullSpace === $Failed,
    (* No conservation laws *)
    {},
    
    (* If gamma is numeric, filter for positive conservation laws *)
    If[MatrixQ[gamma, NumericQ],
      (* Numeric case: select positive vectors *)
      positiveConservationLaws = Select[nullSpace,
        AllTrue[#, # >= 0 &] && AnyTrue[#, # > 0 &] &
      ];
      
      (* If no positive ones found, try negation *)
      If[Length[positiveConservationLaws] == 0,
        positiveConservationLaws = Select[nullSpace,
          AllTrue[-#, # >= 0 &] && AnyTrue[-#, # > 0 &] &
        ];
        positiveConservationLaws = -# & /@ positiveConservationLaws;
      ];
      
      positiveConservationLaws,
      
      (* Symbolic case: return nullspace basis, user applies constraints cp *)
      (* The constraints cp can be used with Reduce/Solve to find positive parametrization *)
      If[cp === {},
        nullSpace,
        (* Apply constraints if provided *)
        constraints = Join[Flatten[Thread[# >= 0] & /@ nullSpace], cp];
        {nullSpace, constraints}
      ]
    ]
  ]
];
(* Numeric matrix - automatic positive parametrization 
cons[{{-1, 1, 0}, {0, -1, 1}}]*)
(* Returns positive conservation vectors *)

(* Symbolic matrix - returns basis 
cons[{{-a, b, 0}, {0, -c, d}}]*)
(* Returns nullspace basis for user to parametrize *)

(* Symbolic matrix with constraints 
cons[{{-a, b, 0}, {0, -c, d}}, {a > 0, b > 0, c > 0, d > 0}]*)
(* Returns {nullspace, constraints} for Reduce/Solve *)


(* extSpe extracts species names preserving case *)
extSpe[reactions_] := Module[{allSpecies, reactants, products}, 
  allSpecies = {};
  Do[
   reactants = comp2Asso[reactions[[i, 1]]];
   products = comp2Asso[reactions[[i, 2]]];
   allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
   , {i, Length[reactions]}];
  DeleteDuplicates[allSpecies]
];



(* extMat moved to epi.wl to avoid duplication *)
(* NOTE: Backup versions exist in Core1.wl and EpidCRNo.wl *)


arrow2pairReac[reactions_] := Module[{converted},
  If[Length[reactions] == 0, Return[{}]];
  If[Head[reactions[[1]]] === Rule,
    converted = Table[
      {First[reactions[[i]]], Last[reactions[[i]]]}, 
      {i, Length[reactions]}
    ];
    Return[converted],
    Return[reactions]
  ]
];


asoRea[RN_]:=Module[{parseSide,extractSpecies},
  extractSpecies[expr_]:=Module[{terms,species},
    If[expr===0 || expr==="0",Return[{}]];
    terms=If[Head[expr]===Plus,List@@expr,{expr}];
    species={};
    Do[Which[
      StringQ[term]&&StringContainsQ[term," "],AppendTo[species,StringTrim[StringDrop[term,1]]],
      Head[term]===Times,AppendTo[species,Cases[term,_String][[1]]],
      StringQ[term]&&term=!="0",AppendTo[species,term]];,{term,terms}];
    DeleteDuplicates[species]];
  parseSide[expr_]:=extractSpecies[expr];
  Map[Function[r,Association["Substrates"->parseSide[r[[1]]],"Products"->parseSide[r[[2]]]]],RN]];



stoichiometricMatrices[reactions_] := Module[{
  species, numReactions, numSpecies, 
  alpha, beta, gamma, reactants, products, leftSide, rightSide},
  
  species = extSpe[reactions];
  numReactions = Length[reactions];
  numSpecies = Length[species];
  
  alpha = Table[0, {numSpecies}, {numReactions}];
  beta = Table[0, {numSpecies}, {numReactions}];
  
  Do[
    If[Head[reactions[[j]]] === Rule,
      leftSide = First[reactions[[j]]];
      rightSide = Last[reactions[[j]]],
      leftSide = reactions[[j, 1]];
      rightSide = reactions[[j, 2]]
    ];
    reactants = comp2Asso[leftSide];
    products = comp2Asso[rightSide];
    Do[
      If[KeyExistsQ[reactants, species[[i]]], 
        alpha[[i, j]] = reactants[species[[i]]]];
      If[KeyExistsQ[products, species[[i]]], 
        beta[[i, j]] = products[species[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  gamma = beta - alpha;
  {alpha, beta, gamma, species}
];

reaToRHS[reactions_] := Module[{alpha, beta, gamma, species, var, rv, tk, Rv, RHS, convertedReactions},
  convertedReactions = If[Length[reactions] > 0 && Head[reactions[[1]]] === Rule,
    Table[{First[reactions[[i]]], Last[reactions[[i]]]}, {i, Length[reactions]}],
    reactions
  ];
  {alpha, beta, gamma, species} = stoichiometricMatrices[convertedReactions];
  var = ToExpression[species];
  rv = expM[var, alpha // Transpose];
  tk = Array[Symbol["k" <> ToString[#]] &, alpha // Transpose // Length];
  Rv = tk*rv;
  RHS = gamma . Rv;
  {RHS, species, Rv}
];

expM=Inner[OperatorApplied[Power],#2,#1,Times]&;
(*Basic parameter extraction*)
Par[RHS_,X_]:=Complement[Variables[RHS],X];
(*Par[{a*x+b*y,c*x+d*y},{x,y}] returns {a,b,c,d}*)

(* ========================================================================= *)
(* TERM EXTRACTION AND SIGN DETECTION *)
(* ========================================================================= *)

(* allT: Extract all terms from expanded polynomial *)
allT[expr_] := If[Head[#] === Plus, List @@ #, {#}] &[Expand[expr]];

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

(* Tests for allT and isN *)
(*
Print["=== allT and isN Tests ==="];
Print["allT[a*x+b*y-c*z]=", allT[a*x + b*y - c*z]];
Print["allT[La-be*i*s-mu*s]=", allT[La - be*i*s - mu*s]];
Print["allT[-i*(ga+mu)+be*i*s]=", allT[-i*(ga+mu) + be*i*s]];
Print["isN[a*x]=", isN[a*x], " isN[-be*i*s]=", isN[-be*i*s], " isN[La]=", isN[La]];
Print["isN[-mu*s]=", isN[-mu*s], " isN[be*i*s]=", isN[be*i*s], " isN[-ga*i]=", isN[-ga*i]];
*)

(* ========================================================================= *)
(* ODE TO REACTION NETWORK CONVERSION *)
(* ========================================================================= *)

(* ODE2RN moved to epi.wl to avoid duplication *)

(* ========================================================================= *)
(* DEFICIENCY AND WEAK REVERSIBILITY *)
(* ========================================================================= *)

(* defi[RN] - computes deficiency: delta = Nc - l - s *)
(* Nc = number of complexes, l = linkage classes, s = rank of gamma *)
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

(* wrQ[RN] - checks if RN is weakly reversible *)
(* WR means each linkage class is strongly connected *)
wrQ[RN_] := Module[{complexes, edges, components, ii, jj, subEdges, subGraph, Nc, isWR, numSCC},
  (* Get complexes and build directed edges *)
  complexes = DeleteDuplicates[Flatten[{First[#], Last[#]} & /@ RN]];
  Nc = Length[complexes];
  If[Nc == 0, Return[False]];

  (* Build directed edge list *)
  edges = Table[
    Module[{srcPos, tgtPos},
      srcPos = Position[complexes, First[RN[[ii]]]][[1, 1]];
      tgtPos = Position[complexes, Last[RN[[ii]]]][[1, 1]];
      DirectedEdge[srcPos, tgtPos]
    ], {ii, Length[RN]}
  ];

  (* Get linkage classes from undirected version *)
  components = WeaklyConnectedComponents[Graph[Range[Nc], edges]];

  (* Check: WR iff each linkage class has exactly 1 SCC *)
  isWR = True;
  Do[
    subEdges = Select[edges, MemberQ[components[[jj]], #[[1]]] && MemberQ[components[[jj]], #[[2]]] &];
    If[Length[subEdges] == 0, isWR = False; Break[]];
    subGraph = Graph[components[[jj]], subEdges];
    numSCC = Length[ConnectedGraphComponents[subGraph, Method -> "StronglyConnected"]];
    If[numSCC > 1, isWR = False; Break[]],
    {jj, Length[components]}
  ];
  isWR
];

(* OLD CODE TO DELETE *)
(*
  coeffsByEq = Table[
    Module[{expanded, termList, monomialList},
      expanded = Expand[RHS[[ii]]];
      Print["Equation ", ii, " expanded: ", expanded];
      termList = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
      Print["termList: ", termList];

      monomialList = Map[Function[{term},
        Module[{coeff, mono, varPart, coeffPart},
          Which[
            (* Constant term *)
            FreeQ[term, Alternatives @@ var],
            mono = 1;
            coeff = term,

            (* Single variable *)
            MemberQ[var, term],
            mono = term;
            coeff = 1,

            (* Product with coefficient *)
            Head[term] === Times,
            varPart = Select[List @@ term, !FreeQ[#, Alternatives @@ var]&];
            coeffPart = Select[List @@ term, FreeQ[#, Alternatives @@ var]&];
            mono = If[Length[varPart] == 0, 1,
                     If[Length[varPart] == 1, varPart[[1]], Times @@ varPart]];
            coeff = If[Length[coeffPart] == 0, 1,
                      If[Length[coeffPart] == 1, coeffPart[[1]], Times @@ coeffPart]],

            (* Power or other *)
            True,
            mono = term;
            coeff = 1
          ];
          {mono, coeff}
        ]
      ], termList];

      monomialList
    ],
    {ii, n}
  ];

  (* Find all unique monomials across all equations *)
  uniqueTerms = DeleteDuplicates[Flatten[coeffsByEq[[All, All, 1]]]];

  Print["coeffsByEq=", coeffsByEq];
  Print["uniqueTerms=", uniqueTerms];

  (* STEP 1: Identify monomials that appear ONLY with positive sign *)
  (* These are production terms: gam=1, alp=0 *)
  monomialSigns = Association[];
  Do[
    mono = uniqueTerms[[mm]];
    (* Check sign of this monomial across all equations *)
    hasPositive = False;
    hasNegative = False;
    Do[
      (* Find this monomial in equation ii *)
      pos = Position[coeffsByEq[[ii, All, 1]], mono];
      If[Length[pos] > 0,
        coeff = coeffsByEq[[ii, pos[[1,1]], 2]];
        (* Check if coefficient has negative sign *)
        (* Check for -1 in Times, or Negative number, or leading minus *)
        If[Negative[coeff] ||
           (Head[coeff] === Times && MemberQ[List @@ coeff, -1]) ||
           (Head[coeff] === Times && Length[coeff] > 0 && Negative[First[List @@ coeff]]),
          hasNegative = True,
          hasPositive = True
        ];
      ];
    , {ii, n}];
    monomialSigns[mono] = {hasPositive, hasNegative};
    Print["mono=", mono, " signs=", monomialSigns[mono]];
  , {mm, Length[uniqueTerms]}];

  (* Separate positive-only terms from negative-containing terms *)
  posOnlyTerms = {};
  negTerms = {};

  Do[
    mono = uniqueTerms[[mm]];
    {hasPos, hasNeg} = monomialSigns[mono];

    If[hasPos && !hasNeg,
      (* Positive only: production term *)
      (* Find which equation it appears in and with what coefficient *)
      Do[
        pos = Position[coeffsByEq[[ii, All, 1]], mono];
        If[Length[pos] > 0,
          coeff = coeffsByEq[[ii, pos[[1,1]], 2]];
          AppendTo[posOnlyTerms, {mono, coeff, ii}];
        ];
      , {ii, n}];
    ];

    If[hasNeg,
      (* Has negative: consumption term *)
      Do[
        pos = Position[coeffsByEq[[ii, All, 1]], mono];
        If[Length[pos] > 0,
          coeff = coeffsByEq[[ii, pos[[1,1]], 2]];
          If[Negative[coeff],
            AppendTo[negTerms, {mono, -coeff, ii}];
          ];
        ];
      , {ii, n}];
    ];
  , {mm, Length[uniqueTerms]}];

  (* Combine all reaction terms: positive-only first, then negative *)
  allReactionTerms = Join[posOnlyTerms, negTerms];

  (* Build rates vector *)
  rts = allReactionTerms[[All, 2]] * allReactionTerms[[All, 1]];

  (* Build gamma matrix *)
  gam = Table[0, {n}, {Length[allReactionTerms]}];
  Do[
    (* Find where this reaction term appears in the RHS *)
    mono = allReactionTerms[[kk, 1]];
    Do[
      pos = Position[coeffsByEq[[ii, All, 1]], mono];
      If[Length[pos] > 0,
        coeff = coeffsByEq[[ii, pos[[1,1]], 2]];
        gam[[ii, kk]] = coeff / allReactionTerms[[kk, 2]];
      ];
    , {ii, n}];
  , {kk, Length[allReactionTerms]}];

  (* Build alpha matrix (reactants) *)
  alp = Table[0, {n}, {Length[allReactionTerms]}];
  Do[
    mono = allReactionTerms[[jj, 1]];

    (* Check if this is a positive-only term *)
    If[MemberQ[posOnlyTerms[[All, 1]], mono],
      (* Positive-only: alpha = 0 (source term) *)
      Continue[],
      (* Consumption term: extract exponents from monomial *)
      If[mono === 1,
        (* Constant term: alpha = 0 *)
        Continue[],
        (* Extract stoichiometric coefficients *)
        Do[
          If[FreeQ[mono, var[[ii]]],
            alp[[ii, jj]] = 0,
            alp[[ii, jj]] = Exponent[mono, var[[ii]]];
          ];
        , {ii, n}];
      ];
    ];
  , {jj, Length[allReactionTerms]}];

  (* Build beta matrix *)
  bet = gam + alp;

  (* Construct reaction network *)
  RN = Table[
    Module[{left, right, kk2},
      (* Left side (reactants) *)
      left = If[Total[alp[[All, jj]]] == 0,
        0,
        Plus @@ Table[
          If[alp[[kk2, jj]] == 0,
            Nothing,
            If[alp[[kk2, jj]] == 1,
              spe[[kk2]],
              alp[[kk2, jj]] * spe[[kk2]]
            ]
          ],
          {kk2, n}
        ]
      ];

      (* Right side (products) *)
      right = If[Total[bet[[All, jj]]] == 0,
        0,
        Plus @@ Table[
          If[bet[[kk2, jj]] == 0,
            Nothing,
            If[bet[[kk2, jj]] == 1,
              spe[[kk2]],
              bet[[kk2, jj]] * spe[[kk2]]
            ]
          ],
          {kk2, n}
        ]
      ];

      left -> right
    ],
    {jj, Length[allReactionTerms]}
  ];
*)
(* END OLD CODE TO DELETE *)

(* ========================================================================= *)
(* ALPHA/BETA TO REACTION NETWORK CONVERSION *)
(* ========================================================================= *)
(* RN = ab2RN[alp, bet, spe] - converts stoichiometric matrices to RN format *)
(* alp: reactant matrix (species x reactions), bet: product matrix *)
(* spe: list of species names (strings) *)

ab2RN[alp_, bet_, spe_] := Module[{n, m, jj, kk, left, right, leftTerms, rightTerms},
  n = Length[spe];
  m = Length[alp[[1]]];

  Table[
    (* Left side (reactants from alp): "s" + "i" format *)
    leftTerms = Table[
      If[alp[[kk, jj]] == 0, Nothing,
        If[alp[[kk, jj]] == 1, spe[[kk]], alp[[kk, jj]] * spe[[kk]]]
      ],
      {kk, n}
    ];
    left = If[Length[leftTerms] == 0, 0,
      If[Length[leftTerms] == 1, leftTerms[[1]], Plus @@ leftTerms]
    ];

    (* Right side (products from bet): 2*"i" format *)
    rightTerms = Table[
      If[bet[[kk, jj]] == 0, Nothing,
        If[bet[[kk, jj]] == 1, spe[[kk]], bet[[kk, jj]] * spe[[kk]]]
      ],
      {kk, n}
    ];
    right = If[Length[rightTerms] == 0, 0,
      If[Length[rightTerms] == 1, rightTerms[[1]], Plus @@ rightTerms]
    ];

    left -> right,
    {jj, m}
  ]
];

(* ========================================================================= *)
(* ODE TO W,Y MATRICES CONVERSION *)
(* ========================================================================= *)
(* {W, Y} = ODE2WY[RHS, var] - RHS == W . x^Y *)

ODE2WY[RHS_List, var_List] := Module[{
  n, m, allTermsByEq, monos, W, Y, ii, jj, kk, allT, getMono, getCoef
  },

  n = Length[var];

  (* Helper: Extract all terms from expression *)
  allT[expr_] := If[Head[expr] === Plus, List @@ expr, {expr}];

  (* Helper: Extract monomial (variable part) from term *)
  getMono[t_] := Module[{varFactors},
    If[FreeQ[t, Alternatives @@ var], 1,
      If[Head[t] === Times,
        varFactors = Select[List @@ t, !FreeQ[#, Alternatives @@ var]&];
        If[Length[varFactors] == 0, 1,
          If[Length[varFactors] == 1, varFactors[[1]], Times @@ varFactors]
        ],
        t
      ]
    ]
  ];

  (* Helper: Extract coefficient from term *)
  getCoef[t_] := Module[{coefFactors},
    If[FreeQ[t, Alternatives @@ var], t,
      If[Head[t] === Times,
        coefFactors = Select[List @@ t, FreeQ[#, Alternatives @@ var]&];
        If[Length[coefFactors] == 0, 1,
          If[Length[coefFactors] == 1, coefFactors[[1]], Times @@ coefFactors]
        ],
        1
      ]
    ]
  ];

  (* Extract all terms from all equations *)
  allTermsByEq = Table[allT[Expand[RHS[[ii]]]], {ii, n}];

  (* Get unique monomials across all equations *)
  monos = DeleteDuplicates[Flatten[Table[getMono /@ allTermsByEq[[ii]], {ii, n}]]];
  monos = Select[monos, # =!= 0 &];
  m = Length[monos];

  (* Build Y matrix: exponents of each variable in each monomial *)
  Y = Table[
    If[FreeQ[monos[[jj]], var[[ii]]], 0, Exponent[monos[[jj]], var[[ii]]]],
    {ii, n}, {jj, m}
  ];

  (* Build W matrix: coefficient of monos[[j]] in RHS[[i]] *)
  W = Table[
    Module[{eqTerms, c},
      eqTerms = allTermsByEq[[ii]];
      c = 0;
      Do[
        If[getMono[eqTerms[[kk]]] === monos[[jj]],
          c = c + getCoef[eqTerms[[kk]]];
        ],
        {kk, Length[eqTerms]}
      ];
      c
    ],
    {ii, n}, {jj, m}
  ];

  {W, Y}
];

(* TEST ODE2WY
RHS = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
var = {s, i};
{W, Y} = ODE2WY[RHS, var];
xY = Table[Product[var[[ii]]^Y[[ii, jj]], {ii, Length[var]}], {jj, Length[Y[[1]]]}];
Print["ODE2WY: ", Simplify[W . xY - RHS] == {0, 0}];
*)

(* ========================================================================= *)
(* ODE TO WR0 REALIZATION *)
(* ========================================================================= *)
(* ODE2WR0[RHS, var] - Finds WR0 realization if it exists *)
(* Input: RHS (ODE right-hand sides), var (variables) *)
(* Output: {success, RN, rts} or {False, reason, {}} *)

ODE2WR0[RHS_List, var_List] := Module[{
  n, m, RNtmp, rts, spe, alp, bet, gam,
  kerGam, kerGamPos, extremeRays, supports, partition, Vp, l,
  NN0, affineIndep, yi, coneGens, decomp, kij,
  failReason, srcIdx, yVecs, ii, jj, p, failed
  },

  n = Length[var];
  spe = toStr /@ var;

  (* Use ODE2RN to get proper reaction decomposition *)
  {RNtmp, rts, spe, alp, bet, gam} = ODE2RN[RHS, var];
  m = Length[rts];

  (* gam is n x m, alp is n x m *)
  (* For WR0 algorithm, we need columns to be source complexes *)
  (* alp columns are the source complexes (Y in paper notation) *)
  (* gam columns are the net reaction vectors (W in paper notation) *)

  failReason = "NotChecked";
  NN0 = {};

  (* Step 2: Find extreme rays of ker(gam) \:2229 R>=0^m *)
  kerGam = NullSpace[gam];

  If[Length[kerGam] == 0,
    failReason = "NoKernel";
    Return[{False, failReason, {}, {}}]
  ];

  kerGamPos = Select[kerGam, AllTrue[#, # >= 0 &] &];

  If[Length[kerGamPos] == 0,
    failReason = "NoNonNegKernel";
    Return[{False, failReason, {}, {}}]
  ];

  extremeRays = kerGamPos;
  (* Get indices where each extreme ray has positive entries *)
  supports = Table[
    Flatten[Position[extremeRays[[ii]], _?(# > 0 &)]],
    {ii, Length[extremeRays]}
  ];
  l = Length[extremeRays];

  (* Step 3: Check if supports partition {1,...,m} *)
  partition = Flatten[supports];
  If[Sort[DeleteDuplicates[partition]] != Range[m] || Length[partition] != m,
    failReason = "NoPartition";
    Return[{False, failReason, {}, {}}]
  ];

  (* Step 6-20: Build WR0 network *)
  failed = False;
  Do[
    If[failed, Break[]];
    Vp = supports[[p]];

    (* Step 7-8: Check affine independence *)
    yVecs = Table[alp[[All, ii]], {ii, Vp}];
    If[Length[Vp] > 1,
      affineIndep = MatrixRank[
        Transpose[Table[yVecs[[ii]] - yVecs[[1]], {ii, 2, Length[yVecs]}]]
      ] == Min[Length[Vp] - 1, n],
      affineIndep = True
    ];

    If[!affineIndep,
      failReason = "NotAffineIndep";
      failed = True;
      Break[]
    ];

    (* Step 9-20: For each complex in linkage class *)
    Do[
      If[failed, Break[]];
      yi = alp[[All, ii]];

      (* Build cone generators yj - yi for j != i in Vp *)
      coneGens = Table[
        If[jj != ii, alp[[All, jj]] - yi, Nothing],
        {jj, Vp}
      ];

      If[Length[coneGens] == 0, Continue[]];

      (* Check if gam[[All,i]] in cone{yj - yi} *)
      decomp = Quiet[LinearProgramming[
        ConstantArray[0, Length[coneGens]],
        Transpose[coneGens],
        gam[[All, ii]],
        Table[{0, Infinity}, {Length[coneGens]}]
      ]];

      If[!VectorQ[decomp],
        failReason = "NotInCone";
        failed = True;
        Break[]
      ];

      (* Add reactions yi -> yj with rate kij *)
      kij = decomp;
      srcIdx = 1;
      Do[
        If[kij[[srcIdx]] > 10^-10,
          AppendTo[NN0, {ii, jj, kij[[srcIdx]]}]
        ];
        srcIdx++,
        {jj, Select[Vp, # != ii &]}
      ],
      {ii, Vp}
    ],
    {p, l}
  ];

  If[failed,
    Return[{False, failReason, {}, {}}]
  ];

  (* Build RN from NN0 *)
  (* NN0 contains {source_idx, target_idx, rate_multiplier} *)
  Module[{RN, finalRts, srcAlp, tgtAlp},
    RN = Table[
      srcAlp = alp[[All, NN0[[ii, 1]]]];
      tgtAlp = alp[[All, NN0[[ii, 2]]]];
      Module[{left, right, leftTerms, rightTerms, kk},
        leftTerms = Table[
          If[srcAlp[[kk]] == 0, Nothing,
            If[srcAlp[[kk]] == 1, spe[[kk]], ToString[srcAlp[[kk]]] <> spe[[kk]]]
          ],
          {kk, n}
        ];
        left = If[Length[leftTerms] == 0, "0",
          If[Length[leftTerms] == 1, leftTerms[[1]], Plus @@ leftTerms]
        ];
        rightTerms = Table[
          If[tgtAlp[[kk]] == 0, Nothing,
            If[tgtAlp[[kk]] == 1, spe[[kk]], ToString[tgtAlp[[kk]]] <> spe[[kk]]]
          ],
          {kk, n}
        ];
        right = If[Length[rightTerms] == 0, "0",
          If[Length[rightTerms] == 1, rightTerms[[1]], Plus @@ rightTerms]
        ];
        left -> right
      ],
      {ii, Length[NN0]}
    ];
    (* Final rates = original rate * decomposition coefficient *)
    finalRts = Table[
      NN0[[ii, 3]] * rts[[NN0[[ii, 1]]]],
      {ii, Length[NN0]}
    ];
    {True, RN, finalRts}
  ]
];

(* ========================================================================= *)
(* SUBSCRIPT CONVERTER *)
(* ========================================================================= *)

(* subsCon - Convert variable names to subscripted form

   Purpose: Converts variable names like "i12", "i1", "x23" to subscripted forms
   where all trailing digits become subscripts: Subscript[i, 12], Subscript[i, 1], etc.

   Input:
   - str: String or Symbol to convert

   Output:
   - Subscripted expression if digits found, original otherwise

   Examples:
   - subsCon["i12"] -> Subscript[i, 12]
   - subsCon["i1"] -> Subscript[i, 1]
   - subsCon["s"] -> s
   - subsCon[i12] -> Subscript[i, 12]
*)

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

(* Test subsCon *)
(*Print["subsCon[\"i12\"]=", subsCon["i12"], " subsCon[\"i1\"]=", subsCon["i1"], " subsCon[\"s\"]=", subsCon["s"]];*)

(* ========================================================================= *)
(* CHILD SELECTION ANALYSIS (Cauchy-Binet)
   ALGORITHM (called from scripts after computing co from char poly of G=gam.RM):
   1. symRM builds symbolic reactivity matrix RM (m x n, entries rt[j,k] or 0).
      G = gam.RM is the symbolic Jacobian. Its char poly cpG has coefficients
      co[[k]] (k=1..n) for lam^(n-k), each a sum of k-minors of G (Cauchy-Binet).
   2. negCS extracts negative summands from each co[[k]] via allT+isN,
      then csInfo maps each to {k, term, sIdx, rIdx, Sk, det} where
      Sk = gam[[sIdx,rIdx]] is the INTEGER stoichiometric submatrix.
      Result: nCS list.
   3. upfQ keeps entries with Sign[det]==(-1)^(k-1), det!=0
      (unstable positive feedback).
   4. hasseCS builds directed graph on nCS by r-set inclusion:
      edge a->b when rSet(a) is strict subset of rSet(b) and k(a)<k(b).
      TransitiveReductionGraph gives Hasse diagram.
   5. oscCS: for each UPF node u, finds descendants via VertexOutComponent,
      tests hurQ (Hurwitz stability of the integer Sk). If any descendant
      is stable, {u, stableDescendants} is an oscillatory core.
   6. Hopf via CS inheritance (Cor 5.21):
      csN = n-CS (full-size child selections from coeff of lam^0 = det G).
      stabN = Hurwitz-stable n-CS. If stabN nonempty AND upf nonempty -> Hopf. *)
(* ========================================================================= *)

(* symRM: symbolic reactivity matrix RM (m x n), entry rt[j,k] if species k
   appears as reactant in reaction j (alp[k,j]>0), else 0. *)
symRM[alp_, n_, m_] := Table[If[alp[[k, j]] > 0, Global`rt[j, k], 0], {j, m}, {k, n}];

(* rSet: extract all rt[j,k] symbols from a term (reactivity set). *)
rSet[term_] := Cases[term, Global`rt[_, _], Infinity];

(* csInfo: given order k, a negative term from co[[k]], and gam,
   extract the CS bijection J: species -> reaction from the rt[j,k] symbols,
   build Sk = gam[[sIdx, rBij]] with columns in bijection order,
   return {k, term, sIdx, rBij, Sk, Det[Sk]}. *)
csInfo[k_, term_, gam_] := Module[{rs, pairs, sIdx, rIdx, rBij, Sk, SkB},
  rs = rSet[term];
  pairs = {#[[1]], #[[2]]} & /@ (rs /. Global`rt -> List);
  sIdx = Sort[Last /@ pairs];
  rIdx = Sort[First /@ pairs];
  If[Length[Union[sIdx]] != Length[Union[rIdx]], Return[Nothing]];
  (* rBij: reactions in bijection order (matching species order) for display *)
  rBij = Table[First[Select[pairs, #[[2]] == s &]][[1]], {s, sIdx}];
  (* Sk: sorted columns for det/sign computations; SkB: bijection order for display *)
  Sk = gam[[sIdx, rIdx]];
  SkB = gam[[sIdx, rBij]];
  {k, term, sIdx, rBij, Sk, Det[Sk], SkB}];

(* upfQ: True if child selection has unstable positive feedback,
   i.e. Sign[det] == (-1)^(k-1) and det != 0. *)
upfQ[{k_, _, sIdx_, _, Sk_, detSk_, _}] :=
  Sign[detSk] == (-1)^(k - 1) && detSk =!= 0;

(* unfQ: True if child selection has unstable negative feedback,
   i.e. Sign[det] == (-1)^k and det != 0. *)
unfQ[{k_, _, sIdx_, _, Sk_, detSk_, _}] :=
  Sign[detSk] == (-1)^k && detSk =!= 0;

(* hurQ: Hurwitz stability test for k x k integer matrix Sk.
   k=1,2: direct conditions. k=3,4: via Hur3/Hur4. k>=5: numeric eigenvalues.
   For symbolic Sk, Hur3/Hur4 return inequality lists; And@@{..} stays
   unevaluated, so Select treats it as non-True (entry skipped). *)
hurQ[Sk_] := Module[{k = Length[Sk], cp1},
  Which[
    k == 1, Sk[[1, 1]] < 0,
    k == 2, Tr[Sk] < 0 && Det[Sk] > 0,
    k == 3, And @@ Hur3[CharacteristicPolynomial[Sk, lam], lam][[2]],
    k == 4, And @@ Hur4[CharacteristicPolynomial[Sk, lam], lam][[2]],
    True, cp1 = Eigenvalues[N[Sk]]; And @@ (Re[#] < 0 & /@ cp1)
  ]];

(* negCS: from coefficients co[[1..n]] of char poly, extract all negative
   summands (child selections with negative sign) via allT+isN,
   then map csInfo to get {k, term, sIdx, rIdx, Sk, det} for each. *)
negCS[co_, gam_, n_] := Module[{summ, neg},
  summ = Table[allT[Expand[co[[k]]]], {k, n}];
  neg = Table[Select[summ[[k]], isN], {k, n}];
  Flatten[Table[csInfo[k, #, gam] & /@ neg[[k]], {k, n}], 1]];

(* hasseCS: Hasse diagram on nCS by r-set inclusion.
   Edge a->b when rSet(a) strict subset of rSet(b) and k(a)<k(b).
   Returns {TransitiveReductionGraph, edges}. *)
hasseCS[nCS_] := Module[{edges = {}, gH},
  Do[Do[
    If[nCS[[a, 1]] < nCS[[b, 1]] &&
       SubsetQ[rSet[nCS[[b, 2]]], rSet[nCS[[a, 2]]]] &&
       rSet[nCS[[a, 2]]] =!= rSet[nCS[[b, 2]]],
      AppendTo[edges, a -> b]],
    {b, Length[nCS]}], {a, Length[nCS]}];
  gH = If[edges === {},
    Graph[Range[Length[nCS]], {}],
    TransitiveReductionGraph[Graph[Range[Length[nCS]], edges]]];
  {gH, edges}];

(* upfCS: UPF feedbacks with stable supersets. For each UPF node u, find
   descendants via VertexOutComponent, test hurQ on each descendant's Sk.
   If any descendant is Hurwitz-stable, {u, stableDescendants} is a UPF feedback.
   Prints summary and each feedback with species names and rates. *)
upfCS[nCS_, gH_, spe_, rts_] := Module[{res = {}, desc, stab, ui, di},
  Do[
    desc = Complement[Normal@VertexOutComponent[gH, u], {u}];
    stab = Select[desc, hurQ[nCS[[#, 5]]] &];
    If[stab =!= {}, AppendTo[res, {u, stab}]],
    {u, Select[Range[Length[nCS]], upfQ[nCS[[#]]] &]}];
  Print["upfCS found ", Length[res], " UPF feedbacks with stable supersets"];
  Do[ui = res[[j, 1]]; di = res[[j, 2]];
    Print["UPF#", j,
      " species=", spe[[nCS[[ui, 3]]]],
      " rates=", rts[[nCS[[ui, 4]]]],
      " stable supersets=", spe[[nCS[[#, 3]]]] & /@ di],
    {j, Length[res]}];
  res];

(* oscRI: Recipe I — Oscillatory Cores of Class I.
   A core is a stable CS that is MINIMAL with the property of being
   Hurwitz-stable and containing a UPF as a principal submatrix.
   Algorithm:
   1. Find all stable CS (hurQ True) that have at least one UPF ancestor.
   2. Among those, keep minimal ones: d is minimal if no other stable CS
      with a UPF ancestor sits strictly between a UPF and d in the Hasse graph.
   3. For each core, report which UPFs it contains. *)
oscRI[nCS_, gH_, spe_, rts_] := Module[
  {upfIdx, allStab, stabWithUPF, cores, ui, coreUPFs},
  upfIdx = Select[Range[Length[nCS]], upfQ[nCS[[#]]] &];
  (* all stable CS that are descendants of at least one UPF *)
  allStab = Select[Range[Length[nCS]], hurQ[nCS[[#, 5]]] &];
  stabWithUPF = Select[allStab, Function[d,
    AnyTrue[upfIdx, MemberQ[Normal@VertexOutComponent[gH, #], d] &]]];
  (* minimal: no other stabWithUPF member is an ancestor of d *)
  cores = Select[stabWithUPF, Function[d,
    !AnyTrue[Complement[stabWithUPF, {d}],
      MemberQ[Normal@VertexOutComponent[gH, #], d] &]]];
  If[Length[cores] == 0,
    Print["Recipe I does not hold: no oscillatory cores of Class I found"],
    Print["Recipe I holds: ", Length[cores], " oscillatory cores of Class I"];
    Do[Module[{c = cores[[j]], cS, cR},
      cS = nCS[[c, 3]]; cR = nCS[[c, 4]];
      coreUPFs = Select[upfIdx,
        MemberQ[Normal@VertexOutComponent[gH, #], c] &];
      Print["RI#", j, " stable ", nCS[[c, 1]], "-CS: species=", spe[[cS]],
        " reactions=", rts[[cR]], " det=", nCS[[c, 6]],
        " Sk=", nCS[[c, 7]] // MatrixForm];
      Do[Print["  contains UPF ", nCS[[u, 1]], "-CS: species=", spe[[nCS[[u, 3]]]],
          " reactions=", rts[[nCS[[u, 4]]]], " det=", nCS[[u, 6]],
          " Sk=", nCS[[u, 7]] // MatrixForm],
        {u, coreUPFs}]],
      {j, Length[cores]}]];
  cores];

(* oscRII: Recipe II -- Oscillatory Cores of Class II [BSV25, Def 4].
   A core is a UNF (unstable negative feedback) of size k whose Sk has
   a (k-1)x(k-1) Hurwitz-stable principal submatrix.
   Keep only minimal such UNFs (no smaller UNF sub-CS also qualifies).
   Prints only the smallest core; returns all cores. *)
oscRII[nCS_, gH_, spe_, rts_] := Module[
  {unfIdx, cands, k, Sk, cores, sorted},
  unfIdx = Select[Range[Length[nCS]], unfQ[nCS[[#]]] &];
  cands = Select[unfIdx, Function[u,
    k = nCS[[u, 1]]; Sk = nCS[[u, 5]];
    k >= 2 && AnyTrue[
      Table[Sk[[Delete[Range[k], j], Delete[Range[k], j]]], {j, k}],
      hurQ]]];
  cores = Select[cands, Function[d,
    !AnyTrue[Complement[cands, {d}],
      MemberQ[Normal@VertexOutComponent[gH, #], d] &]]];
  If[Length[cores] == 0,
    Print["Recipe II: no cores of Class II found"],
    sorted = SortBy[cores, nCS[[#, 1]] &];
    Print["Recipe II: ", Length[cores], " cores (showing smallest)"];
    Module[{c = sorted[[1]], cK, cSk, cSpe, stSub},
      cK = nCS[[c, 1]]; cSk = nCS[[c, 5]]; cSpe = nCS[[c, 3]];
      stSub = Select[Range[cK],
        hurQ[cSk[[Delete[Range[cK], #], Delete[Range[cK], #]]]] &];
      Print["RII UNF ", cK, "-CS: spe=", spe[[cSpe]],
        " rts=", rts[[nCS[[c, 4]]]], " det=", nCS[[c, 6]]];
      Do[Print["  stable ", (cK - 1), "-sub: spe=",
          spe[[cSpe[[Delete[Range[cK], s]]]]]],
        {s, stSub}]]];
  cores];

(* fixSgn: fixed-sign determinant test [Vas23, Cor 5.7/5.8].
   For each order k, checks if all k-CS have sign-modified determinants
   of the same sign. Returns list of {k, isFixedSign, signs} for k=1..n. *)
fixSgn[nCS_, n_] := Table[Module[{csK, dets},
  csK = Select[nCS, #[[1]] == k &];
  dets = #[[6]] & /@ csK;
  {k, And @@ (# >= 0 & /@ dets) || And @@ (# <= 0 & /@ dets),
   Union[Sign /@ dets]}],
  {k, n}];

(* csAna: full Child Selection analysis pipeline.
   Input: gam (stoichiometric), alp (reactant), spe (species), rts (rates).
   Runs: symRM, negCS, hasseCS, upfCS, oscRI, oscRII, fixSgn.
   Prints summary with unstable orders and Hopf(5.21) conclusion.
   Returns {G, co, nCS, upf, gH, fs, riC, riiC}. *)
csAna[gam_, alp_, spe_, rts_] := Module[
  {n, m, RM, G, cpG, co, nCS, upf, gH, edges,
   csN, stabN, fs, nf, riC, riiC},
  {n, m} = Dimensions[alp];
  RM = symRM[alp, n, m];
  G = gam . RM;
  cpG = CharacteristicPolynomial[G, lam];
  co = Rest@Reverse@CoefficientList[Expand[cpG], lam];
  nCS = negCS[co, gam, n]; upf = Select[nCS, upfQ];
  {gH, edges} = hasseCS[nCS];
  Print["CS: dim=", {n, m}, " #nCS=", Length[nCS], " #upf=", Length[upf]];
  upfCS[nCS, gH, spe, rts];
  riC = oscRI[nCS, gH, spe, rts];
  riiC = oscRII[nCS, gH, spe, rts];
  csN = Select[nCS, #[[1]] == n &];
  stabN = Select[csN, hurQ[#[[5]]] &];
  fs = fixSgn[nCS, n];
  nf = Select[fs, !#[[2]] &];
  Print["fixSgn: unstable orders=", #[[1]] & /@ nf,
        " Hopf(5.21)=", Length[stabN] > 0 && Length[upf] > 0];
  {G, co, nCS, upf, gH, fs, riC, riiC}];

(* ========================================================================= *)
(* ODE TO RN WITH SYMBOLIC RATES *)
(* ========================================================================= *)
(* ODE2RNs: like ODE2RN but handles opaque function calls f[S,i], T[i] etc.
   Uses Coefficient to extract integer stoichiometry, Level to get reactant species.
   Polynomial terms go to ODE2RN; function reactions are appended. *)

ODE2RNs[RHS_, var_] := Module[
  {n, spe, fC, nf, gamF, alpF, betF, pRHS, pO, funRN, funRNs},
  n = Length[var];
  spe = ToString /@ var;
  (* Detect opaque function calls h[x..] where h is a user symbol
     (not arithmetic, not NumericFunction like Sin/Exp/Log)
     and at least one argument is a state variable *)
  fC = DeleteDuplicates[Cases[RHS,
    h_[x__] /; AtomQ[h] &&
      !MemberQ[{Times, Plus, Power, List, Rational, Complex}, h] &&
      !MemberQ[Attributes[h], NumericFunction] &&
      !FreeQ[{x}, Alternatives @@ var] :> h[x], {0, Infinity}]];
  nf = Length[fC];
  If[nf == 0, Return[ODE2RN[RHS, var]]];
  gamF = Table[Coefficient[RHS[[j]], fC[[k]]], {j, n}, {k, nf}];
  alpF = Table[Boole[MemberQ[Level[fC[[k]], {1}], var[[j]]]], {j, n}, {k, nf}];
  betF = gamF + alpF;
  (* Fix negative stoichiometry: move negatives from bet to alp *)
  Do[If[TrueQ[betF[[ii, jj]] < 0],
    alpF[[ii, jj]] = alpF[[ii, jj]] - betF[[ii, jj]];
    betF[[ii, jj]] = 0], {ii, n}, {jj, nf}];
  pRHS = Expand[RHS - gamF . fC];
  pO = Block[{Print = (Null &)}, ODE2RN[pRHS, var]];
  (* String-species RN format (consistent with ODE2RN element 1) *)
  funRN = Table[Module[{lT, rT, left, right},
    lT = Table[If[alpF[[kk, jj]] === 0, Nothing,
      If[alpF[[kk, jj]] === 1, spe[[kk]], alpF[[kk, jj]] spe[[kk]]]], {kk, n}];
    left = If[lT === {}, 0, If[Length[lT] == 1, lT[[1]], Plus @@ lT]];
    rT = Table[If[betF[[kk, jj]] === 0, Nothing,
      If[betF[[kk, jj]] === 1, spe[[kk]], betF[[kk, jj]] spe[[kk]]]], {kk, n}];
    right = If[rT === {}, 0, If[Length[rT] == 1, rT[[1]], Plus @@ rT]];
    left -> right], {jj, nf}];
  (* Symbol-species reduced reactions (element 7) *)
  funRNs = Table[Total[alpF[[All, k]] var] -> Total[betF[[All, k]] var], {k, nf}];
  Module[{aRN, aRts, aRnR},
    aRN = Join[pO[[1]], funRN]; aRts = Join[pO[[2]], fC]; aRnR = Join[pO[[7]], funRNs];
    Print[Length[aRN], " reactions and rts=",
      Row[{Prepend[Transpose[{aRN, Factor /@ aRts}], {"RN", "rts"}] // MatrixForm,
        Spacer[20],
        Prepend[Transpose[{aRnR, Factor /@ aRts}], {"RNr", "rts"}] // MatrixForm}]]];
  {Join[pO[[1]], funRN], Join[pO[[2]], fC], pO[[3]],
   Join[pO[[4]], alpF, 2], Join[pO[[5]], betF, 2], Join[pO[[6]], gamF, 2],
   Join[pO[[7]], funRNs]}];

ODE2RNs[RHS_, var_, prF_] := Module[
  {n, spe, fC, nf, gamF, alpF, betF, pRHS, pO, funRN, funRNs},
  n = Length[var];
  spe = ToString /@ var;
  fC = DeleteDuplicates[Cases[RHS,
    h_[x__] /; AtomQ[h] &&
      !MemberQ[{Times, Plus, Power, List, Rational, Complex}, h] &&
      !MemberQ[Attributes[h], NumericFunction] &&
      !FreeQ[{x}, Alternatives @@ var] :> h[x], {0, Infinity}]];
  nf = Length[fC];
  If[nf == 0, Return[ODE2RN[RHS, var, prF]]];
  gamF = Table[Coefficient[RHS[[j]], fC[[k]]], {j, n}, {k, nf}];
  alpF = Table[Boole[MemberQ[Level[fC[[k]], {1}], var[[j]]]], {j, n}, {k, nf}];
  betF = gamF + alpF;
  Do[If[TrueQ[betF[[ii, jj]] < 0],
    alpF[[ii, jj]] = alpF[[ii, jj]] - betF[[ii, jj]];
    betF[[ii, jj]] = 0], {ii, n}, {jj, nf}];
  pRHS = Expand[RHS - gamF . fC];
  pO = Block[{Print = (Null &)}, ODE2RN[pRHS, var, prF]];
  funRN = Table[Module[{lT, rT, left, right},
    lT = Table[If[alpF[[kk, jj]] === 0, Nothing,
      If[alpF[[kk, jj]] === 1, spe[[kk]], alpF[[kk, jj]] spe[[kk]]]], {kk, n}];
    left = If[lT === {}, 0, If[Length[lT] == 1, lT[[1]], Plus @@ lT]];
    rT = Table[If[betF[[kk, jj]] === 0, Nothing,
      If[betF[[kk, jj]] === 1, spe[[kk]], betF[[kk, jj]] spe[[kk]]]], {kk, n}];
    right = If[rT === {}, 0, If[Length[rT] == 1, rT[[1]], Plus @@ rT]];
    left -> right], {jj, nf}];
  funRNs = Table[Total[alpF[[All, k]] var] -> Total[betF[[All, k]] var], {k, nf}];
  Module[{aRN, aRts, aRnR},
    aRN = Join[pO[[1]], funRN]; aRts = Join[pO[[2]], fC]; aRnR = Join[pO[[7]], prF /@ funRNs];
    Print[Length[aRN], " reactions and rts=",
      Row[{Prepend[Transpose[{prF[aRN], prF[Factor /@ aRts]}], {"RN", "rts"}] // MatrixForm,
        Spacer[20],
        Prepend[Transpose[{prF[aRnR], prF[Factor /@ aRts]}], {"RNr", "rts"}] // MatrixForm}]]];
  {Join[pO[[1]], funRN], Join[pO[[2]], fC], pO[[3]],
   Join[pO[[4]], alpF, 2], Join[pO[[5]], betF, 2], Join[pO[[6]], gamF, 2],
   Join[pO[[7]], prF /@ funRNs]}];

resHopf[M_] := Module[{lam, ch, lc, co, pIw, reP, imP, res, w},
  ch = Expand[CharacteristicPolynomial[M, lam]];
  lc = Coefficient[ch, lam, Length[M]];
  ch = Expand[ch/lc];
  co = Reverse[CoefficientList[ch, lam]];
  pIw = Expand[ch /. lam -> I*w];
  reP = ComplexExpand[Re[pIw]] // Expand;
  imP = ComplexExpand[Im[pIw]] // Expand;
  reP = Expand[reP/w^Exponent[reP, w, Min]];
  imP = Expand[imP/w^Exponent[imP, w, Min]];
  res = Resultant[reP, imP, w];
  {co, res}];

(* ========================================================================= *)
(* TEST ODE2RN *)
(* ====================================================

RHS0 = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
var0 = {s, i};
{RN0, rts0, spe0, alp0, bet0, gam0} = ODE2RN[RHS0, var0];
Print["RN0=", RN0];
Print["ODE2RN: ", Simplify[gam0 . rts0 - RHS0] == {0, 0}]*)

