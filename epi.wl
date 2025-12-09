(* ::Package:: *)

(* ========================================================================= *)
(* EPI - MINIMAL EPIDEMIC MODELING SUBPACKAGE *)
(* Exact copy-paste from Core.wl, Boundary.wl, Siphons.wl, Utils.wl *)
(* ========================================================================= *)


(* ========================================================================= *)
(* FROM UTILS.WL *)
(* ========================================================================= *)

posM= 
Replace[#,{_?Negative->0,e_:>Replace[Expand[e],
{Times[_?Negative,_]->0,
t_Plus:>Replace[t,_?Negative|Times[_?Negative,_]->0,1]}]},{2}]&;
isNNe[m_] := m // Together // NumeratorDenominator //
  Map@CoefficientArrays //
  ReplaceAll[sa_SparseArray :> sa["NonzeroValues"]] // Flatten //
  AllTrue[#, NonNegative] &;


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
  ToLowerCase /@ DeleteDuplicates[allSpecies]
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

ODE2RN[RHS_List, var_List] := Module[{
  n, allTermsByEq, allTermsFlat, sources, sinks, rtsRaw,
  gam, alp, bet, RN, spe, ii, jj, kk,
  alpRaw, gamRaw, uniqueAlpCols, mergeIdx, rts
  },

  n = Length[var];
  spe = ToLowerCase[ToString[#]] & /@ var;

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

  {RN, rts, spe, alp, bet, gam}
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
  spe, al, be, gamma, Rv, RHS,
  numReactions, numSpecies, reactants, products,
  var, rv, tk, defFormula, defTerms, defResult, Nc, l, s, mSi},

  (* Extract species (lowercase), use provided order if given *)
  If[Length[speOrder] > 0,
    spe = ToLowerCase[ToString[#]] & /@ speOrder;
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
  al = ConstantArray[0, {numSpecies, numReactions}];
  be = ConstantArray[0, {numSpecies, numReactions}];
  
  (* Build stoichiometric matrices *)
  Do[
    reactants = comp2Asso[reactions[[j, 1]]];
    products = comp2Asso[reactions[[j, 2]]];
    
    (* Fill matrices - match lowercase species *)
    Do[
      If[KeyExistsQ[reactants, spe[[i]]], 
        al[[i, j]] = reactants[spe[[i]]]];
      If[KeyExistsQ[products, spe[[i]]], 
        be[[i, j]] = products[spe[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  (* Calculate derived matrices *)
  gamma = be - al;
  
  (* Convert lowercase species to variables *)
  var = ToExpression /@ spe;
  
  (* Calculate reaction rates *)
  rv = Table[
    Product[
      If[al[[i, j]] == 0, 
        1, 
        (* Handle symbolic exponents *)
        If[NumericQ[al[[i, j]]], 
          var[[i]]^al[[i, j]],
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
  {spe, al, be, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
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
  
  species = If[ListQ[vars] && AllTrue[vars, StringQ], vars, ToString /@ vars];
  reactionsAsso = asoRea[reactions];
  
  n = Length[species];
  m = Length[reactions];
  
  alpha = ConstantArray[0, {n, m}];
  beta = ConstantArray[0, {n, m}];
  
  Do[
    Module[{subs, prods},
      subs = Lookup[reactionsAsso[[j]], "Substrates", {}];
      prods = Lookup[reactionsAsso[[j]], "Products", {}];
      
      (* Convert to lowercase strings for comparison *)
      subs = ToLowerCase[ToString[#]] & /@ subs;
      prods = ToLowerCase[ToString[#]] & /@ prods;
      
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
        (* Convert species to lowercase and store *)
        If[StringQ[species],
          result[ToLowerCase[species]] = coeff,
          (* Fallback if species is not a string (shouldn't happen) *)
          result[ToLowerCase[ToString[species]]] = coeff
        ],
        
        (* Handle string species alone (coefficient 1) *)
        StringQ[term],
        result[ToLowerCase[term]] = 1,
        
        (* Handle symbol species alone (coefficient 1) *)
        AtomQ[term],
        result[ToLowerCase[ToString[term]]] = 1,
        
        (* Handle other cases - convert to lowercase *)
        True,
        result[ToLowerCase[ToString[term]]] = 1
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
     Vinv = Inverse[V] /. Thread[X[[inf]] -> 0];
     allPos = AllTrue[Flatten[Vinv], isNNe];
     If[allPos,
       F = Fuser;
       Print["NGM: user F accepted (Inverse[V] all positive)"],
       Print["NGM: user F rejected (Inverse[V] not all positive), using posM"];
       F = posM[F1];
       V = F - Jx;
     ];
   ];
   K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // Simplify;
   Kd = (Inverse[V] . F) /. Thread[X[[inf]] -> 0] // Simplify;
   {Jx, F, V, K, Jy, Jxy, Jyx, Kd}
  ]



bdAn[RN_, rts_, var_] := Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, cv, ct, mS, mSi, inf, mod,
K, R0A, cDFE, RDFE, eq0, var0, E0, Jx, Jy, eigenSystem, eigenvals, eigenvecs,
nonzeroIndices, relevantEigenvals, ngm, infVars},

  (* Extract stoichiometric matrices with var order *)
  {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
  RHS = gam . rts;
  par = Par[RHS, var];
  cp = Thread[par > 0];
  cv = Thread[var >= 0];
  ct = Join[cp, cv];
  mS = minSiph[spe, asoRea[RN]];
  mSi = mS[[1]];

  (* infVars is union of all siphon variables - keep as strings for mSi output *)
  infVars = Union[Flatten[mSi]];

  (* Convert string species to symbol variables for substitution *)
  cDFE = Thread[ToExpression[infVars] -> 0];
  RDFE = RHS /. cDFE;
  eq0 = Thread[RDFE == 0];
  var0 = Complement[var, ToExpression[infVars]];
  E0 = Join[Solve[eq0, var0] // Flatten, cDFE];

  (* NGM expects symbols, so convert infVars to symbols *)
  mod = {RHS, var, par};
  ngm = NGM[mod, ToExpression[infVars]];
  Jx = ngm[[1]];
  Jy = ngm[[5]];
  K = ngm[[4]];

  eigenvals = Eigenvalues[K];
  R0A = Select[eigenvals, # =!= 0 &];
  (*eigenSystem = Eigensystem[K];
  eigenvals = eigenSystem[[1]];
  nonzeroIndices = Flatten[Position[eigenvals, _?(# =!= 0 &)]];
  relevantEigenvals = eigenvals[[nonzeroIndices]];
  R0A = relevantEigenvals;*)

  {RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ngm}
];

(* TEST - Commented out
RHS1 = {-be1 i12 s - be2 i12 s - mu0 s - (i1 mu1 R1 s)/s0 - (i2 mu2 R2 s)/s0 - (i12 mu3 R3 s)/s0 + mu0 s0,
        -et1 i1 i12 - ga1 i1 i2 - i1 mu1 + be1 i12 s + (i1 mu1 R1 s)/s0,
        -ga2 i1 i2 - et2 i12 i2 - i2 mu2 + be2 i12 s + (i2 mu2 R2 s)/s0,
        et1 i1 i12 + ga1 i1 i2 + ga2 i1 i2 + et2 i12 i2 - i12 mu3 + (i12 mu3 R3 s)/s0};
var = {s, i1, i2, i12};
{RN, rts, spe, alp, bet, gam} = ODE2RN[RHS1, var];
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ng} = bdAn[RN, rts, var];
RHS - RHS1 // FullSimplify  (*Should be {0,0,0,0}*)
*)

