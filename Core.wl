(* ::Package:: *)

(* ========================================================================= *)
(* CORE IMPLEMENTATION *)
(* ========================================================================= *)
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

  (* Step 2: Find extreme rays of ker(gam) ∩ R>=0^m *)
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
(* TEST ODE2RN *)
(* ====================================================

RHS0 = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
var0 = {s, i};
{RN0, rts0, spe0, alp0, bet0, gam0} = ODE2RN[RHS0, var0];
Print["RN0=", RN0];
Print["ODE2RN: ", Simplify[gam0 . rts0 - RHS0] == {0, 0}];
bdAn[RN0, rts0]*)



