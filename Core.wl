(* ::Package:: *)

BeginPackage["EpidCRN`"];
Begin["`Private`"];

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
  
  (* Minimal siphons *)
  mSi = minSiph[spe, reactions][[1]];
  
  (* Return *)
  {spe, al, be, gamma, Rv, RHS, mSi, {defFormula, defTerms, defResult}}
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
End[];
EndPackage[];

