(* ::Package:: *)

(* ========================================================================= *)
(* CORE IMPLEMENTATION *)
(* ========================================================================= *)
(* ========================================================================== *)
(* SYMBOL TO STRING CONVERSION *)
(* ========================================================================== *)

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


(* Fixed extSpe to return lowercase species names *)
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
extMat[reactions_] := Module[{
  spe, al, be, gamma, Rv, RHS, 
  numReactions, numSpecies, reactants, products, 
  var, rv, tk, complexes,
  defFormula, defTerms, defResult, Nc, l, s, mSi},
  
  (* Extract species (lowercase) *)
  spe = extSpe[reactions];
  
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
  
  (* Complexes for deficiency *)
  complexes = DeleteDuplicates[
    Flatten[Table[{reactions[[i, 1]], reactions[[i, 2]]}, {i, Length[reactions]}]]
  ];
  
  (* Deficiency - for symbolic matrices, may need numeric substitution *)
  Nc = Length[complexes];
  l = 1;
  s = If[gamma === {} || AllTrue[Flatten[gamma], # == 0 &], 
    0, 
    If[AllTrue[Flatten[gamma], NumericQ],
      MatrixRank[gamma],
      "symbolic"
    ]
  ];
  
  defFormula = "\[Delta] = Nc - \[ScriptL] - s";
  defTerms = "Nc = " <> ToString[Nc] <> " (complexes), " <>
             "\[ScriptL] = " <> ToString[l] <> " (linkage classes), " <> 
             "s = " <> ToString[s] <> " (stoich dimension)";
  defResult = If[s === "symbolic",
    "\[Delta] = symbolic (requires numeric parameter values)",
    "\[Delta] = " <> ToString[Nc] <> " - " <> ToString[l] <>
      " - " <> ToString[s] <> " = " <> ToString[Nc - l - s]
  ];

  (* Return *)
  {spe, al, be, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
];


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

ODE2RN[RHS_List, var_List] := Module[{
  n, allTermsByEq, allTermsFlat, sources, sinks, rts,
  gam, alp, bet, RN, spe, ii, jj, kk, mono, coeff, varPart, coeffPart
  },

  n = Length[var];
  spe = ToLowerCase[ToString[#]] & /@ var;

  (* STEP 1: Extract all terms from each equation using allT *)
  allTermsByEq = Table[allT[RHS[[ii]]], {ii, n}];

  (* Flatten all terms across all equations *)
  allTermsFlat = Flatten[allTermsByEq];

  (* STEP 2: Separate into sources (positive-only) and sinks (negative) using isN *)
  sources = Select[allTermsFlat, !isN[#]&];  (* positive terms *)
  sinks = Select[allTermsFlat, isN[#]&];     (* negative terms *)

  (* Remove duplicates *)
  sources = DeleteDuplicates[sources];
  sinks = DeleteDuplicates[sinks];

  (* STEP 3: Build rts from sources and negated sinks *)
  (* For sinks, negate them to get positive rates *)
  rts = Join[sources, -sinks];

  (* Remove duplicates and any zero rates *)
  rts = DeleteDuplicates[rts];
  rts = Select[rts, # =!= 0 &];

  (* Build gamma matrix: gam[[eq, rate]] = coefficient of rts[[rate]] in RHS[[eq]] *)
  gam = Table[
    Module[{eqTerms, coeff},
      eqTerms = allT[RHS[[ii]]];
      Table[
        (* Check if this rate appears in this equation *)
        coeff = 0;
        Do[
          Module[{term, termAbs},
            term = eqTerms[[kk]];
            termAbs = If[isN[term], -term, term];
            (* If this term matches the rate *)
            If[termAbs === rts[[jj]],
              coeff = If[isN[term], -1, 1];
              Break[];
            ];
          ],
          {kk, Length[eqTerms]}
        ];
        coeff
      , {jj, Length[rts]}]
    ],
    {ii, n}
  ];

  (* Build alpha matrix: alp[[species, rate]] = exponent of var[[species]] in rts[[rate]] *)
  alp = Table[
    Module[{rate},
      rate = rts[[jj]];
      (* Extract exponent of var[[ii]] from rate *)
      If[FreeQ[rate, var[[ii]]], 0, Exponent[rate, var[[ii]]]]
    ],
    {ii, n}, {jj, Length[rts]}
  ];

  (* Build beta matrix: bet = gam + alp *)
  bet = gam + alp;

  (* Build reaction network RN in String-Plus format (Format 2: "s" + "i") *)
  RN = Table[
    Module[{left, right, kk, leftTerms, rightTerms},
      (* Left side (reactants from alp) *)
      leftTerms = Table[
        If[alp[[kk, jj]] == 0,
          Nothing,
          If[alp[[kk, jj]] == 1,
            spe[[kk]],
            ToString[alp[[kk, jj]]] <> spe[[kk]]
          ]
        ],
        {kk, n}
      ];
      left = If[Length[leftTerms] == 0,
        "0",
        If[Length[leftTerms] == 1,
          leftTerms[[1]],
          Plus @@ leftTerms  (* Creates "s" + "i" Plus format *)
        ]
      ];

      (* Right side (products from bet) *)
      rightTerms = Table[
        If[bet[[kk, jj]] == 0,
          Nothing,
          If[bet[[kk, jj]] == 1,
            spe[[kk]],
            ToString[bet[[kk, jj]]] <> spe[[kk]]
          ]
        ],
        {kk, n}
      ];
      right = If[Length[rightTerms] == 0,
        "0",
        If[Length[rightTerms] == 1,
          rightTerms[[1]],
          Plus @@ rightTerms  (* Creates "s" + "i" Plus format *)
        ]
      ];

      left -> right
    ],
    {jj, Length[rts]}
  ];

  (* Print["gam=", gam]; *)
  (* Print["alp=", alp//MatrixForm]; *)
  (* Print["bet=", bet//MatrixForm]; *)
  (* Print["RN=", RN]; *)

  (* Return *)
  {RN, rts, spe, alp, bet, gam}
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
    (* Left side (reactants from alp) *)
    leftTerms = Table[
      If[alp[[kk, jj]] == 0, Nothing,
        If[alp[[kk, jj]] == 1, spe[[kk]], ToString[alp[[kk, jj]]] <> spe[[kk]]]
      ],
      {kk, n}
    ];
    left = If[Length[leftTerms] == 0, "0",
      If[Length[leftTerms] == 1, leftTerms[[1]], Plus @@ leftTerms]
    ];

    (* Right side (products from bet) *)
    rightTerms = Table[
      If[bet[[kk, jj]] == 0, Nothing,
        If[bet[[kk, jj]] == 1, spe[[kk]], ToString[bet[[kk, jj]]] <> spe[[kk]]]
      ],
      {kk, n}
    ];
    right = If[Length[rightTerms] == 0, "0",
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
(* TEST ODE2RN *)
(* ====================================================

RHS0 = {La - be*i*s - mu*s, be*i*s - (ga + mu)*i};
var0 = {s, i};
{RN0, rts0, spe0, alp0, bet0, gam0} = ODE2RN[RHS0, var0];
Print["RN0=", RN0];
Print["ODE2RN: ", Simplify[gam0 . rts0 - RHS0] == {0, 0}];
bdAn[RN0, rts0]*)



