(* ::Package:: *)

(* ========================================================================= *)
(* EPI - MINIMAL EPIDEMIC MODELING SUBPACKAGE *)
(* Contains the essential functions for epidemic CRN modeling *)
(* IMPORTANT: This is the canonical version of ODE2RN and extMat *)
(* Goal: Serve as foundation for future Python port *)
(* NOTE: Backup versions exist in Core1.wl and EpidCRNo.wl *)
(* ========================================================================= *)


(* ========================================================================= *)
(* FROM UTILS.WL *)
(* ========================================================================= *)

posM=
Replace[#,{_?Negative->0,e_:>Replace[Expand[e],
{Times[_?Negative,_]->0,
t_Plus:>Replace[t,_?Negative|Times[_?Negative,_]->0,1]}]},{2}]&;


(* ========================================================================= *)
(* FROM CORE.WL *)
(* ========================================================================= *)

Par[RHS_,X_]:=Complement[Variables[RHS],X];

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
    Module[{catList, newL, newR, catExpr, lExpr, rExpr},
      catList = catCoef[[All, jj]];
      newL = newAlp[[All, jj]];
      newR = newBet[[All, jj]];
      catExpr = coefToExpr[catList, var];
      lExpr = coefToExpr[newL, var];
      rExpr = coefToExpr[newR, var];
      If[catExpr === {},
        lExpr -> rExpr,
        Row[{lExpr, Overscript[" \[Rule] ", catExpr], rExpr}]
      ]
    ],
    {jj, Length[rts]}
  ];

  If[prF =!= Identity,
    Print[Length[RN], " reactions and rts=",
      Row[{
        Transpose[{prF[RN], prF[Factor /@ rts]}] // MatrixForm,
        Spacer[20],
        Transpose[{prF[rnRed], prF[Factor /@ rts]}] // MatrixForm
      }]]
  ];

  {RN, rts, spe, alp, bet, gam, rnRed}
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



minSiph[vars_, reactions_] := Module[{
  species, alpha, beta, n, m, reactionsAsso, siphons, minimal, nonm},

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
  {minimal, nonm}
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


(* ========================================================================= *)
(* FROM BOUNDARY.WL *)
(* ========================================================================= *)


NGM[mod_, infVars_, Fuser_: {}] :=
  Module[{dyn, X, par, inf, infc, Jx, Jy, Jxy, Jyx, V1, F1, F, V, K, Kd, Vinv, allPos, dfeRules},
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



bdAn[RN_, rts_, var_] := Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, cv, ct, mS, mSi, inf, mod,
K, R0A, cDFE, RDFE, eq0, var0, E0, Jx, Jy, F, V, eigenvals, ngm, infVarsStr, infVars},

  (* Extract stoichiometric matrices with var order *)
  {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
  RHS = gam . rts;
  par = Par[RHS, var];
  cp = Thread[par > 0];
  cv = Thread[var >= 0];
  ct = Join[cp, cv];
  mS = minSiph[var, RN];
  mSi = mS[[1]];

  (* infVars as strings for internal use *)
  infVarsStr = Union[Flatten[mSi]];

  (* Convert string species to symbol variables for substitution *)
  cDFE = Thread[ToExpression[infVarsStr] -> 0];
  RDFE = RHS /. cDFE;
  eq0 = Thread[RDFE == 0];
  var0 = Complement[var, ToExpression[infVarsStr]];
  E0 = Join[Solve[eq0, var0] // Flatten, cDFE];

  (* NGM expects symbols, so convert infVarsStr to symbols *)
  mod = {RHS, var, par};
  ngm = NGM[mod, ToExpression[infVarsStr]];
  Jx = ngm[[1]];
  F = ngm[[2]];
  V = ngm[[3]];
  K = ngm[[4]];
  Jy = ngm[[5]];

  (* infVars are the variable symbols corresponding to rows/columns of K *)
  infVars = ToExpression[infVarsStr];

  (* Check if V is singular (K = 0 from NGM) *)
  If[K === 0,
    (* V is singular - print informative warning and set dummy values *)
    Module[{Vdfe, rnk, dim},
      Vdfe = V /. Thread[infVars -> 0];
      rnk = MatrixRank[Vdfe];
      dim = Length[Vdfe];
      Print["V with dim=", dim, " has rank=", rnk, ", so is singular"];
    ];
    R0A = {};
  ,
    (* V is nonsingular - compute R0A *)
    eigenvals = Eigenvalues[K];
    R0A = Select[eigenvals, # =!= 0 &];
  ];

  {RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars}
];

(* TEST posM - Commented out
Print["posM inp={{-a+2b,c-3d},{ef-g,2h}} out=", posM[{{-a + 2*b, c - 3*d}, {e*f - g, 2*h}}], " exp={{2b,c},{ef,2h}}"];
*)

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

