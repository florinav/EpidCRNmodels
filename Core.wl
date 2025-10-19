(* ::Package:: *)

BeginPackage["EpidCRN`"];
Begin["`Private`"];

(* ========================================================================= *)
(* CORE IMPLEMENTATION *)
(* ========================================================================= *)

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

(* Convert compound expression to association with lowercase species names *)
compToAsso[expr_] := Module[{terms, result},
  If[expr === 0 || expr === Null,
    Association[], (* Empty association for zero or null *)
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = Association[];
    Do[
      Which[
        Head[term] === Times && Length[term] >= 2 && NumericQ[First[term]],
        (* Handle cases like 2*A - convert to lowercase *)
        result[ToLowerCase[ToString[Last[term]]]] = First[term],
        Head[term] === Times,
        (* Handle cases like A*B (coefficient 1) *)
        result[ToLowerCase[ToString[term]]] = 1,
        StringQ[term],
        (* Handle string species - convert to lowercase *)
        result[ToLowerCase[term]] = 1,
        True,
        (* Handle other cases - convert to lowercase *)
        result[ToLowerCase[ToString[term]]] = 1
      ];
    , {term, terms}];
    result
  ]
];

(* Test 
Print["=== Test compToAsso ==="];
Print[compToAsso[0]];
Print[compToAsso["S"]];
Print[compToAsso["S" + "I1"]];
Print[compToAsso[2*"I2"]];
Print[compToAsso["S" + 2*"I1" + 3*"I2"]];*)


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
   reactants = compToAsso[reactions[[i, 1]]];
   products = compToAsso[reactions[[i, 2]]];
   allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
   , {i, Length[reactions]}];
  ToLowerCase /@ DeleteDuplicates[allSpecies]
];

(* Main extMat function - simplified using extSpe *)
extMat[reactions_] := Module[{
  spe, al, be, gamma, Rv, RHS, 
  numReactions, numSpecies, reactants, products, 
  var, rv, tk, complexes, deficiency,
  defFormula, defTerms, defResult, Nc, l, s, mSi, speOrig},
  
  (* Get original species names for matching *)
  speOrig = Module[{allSpecies = {}, reactants, products}, 
    Do[
     reactants = compToAsso[reactions[[i, 1]]];
     products = compToAsso[reactions[[i, 2]]];
     allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
     , {i, Length[reactions]}];
    DeleteDuplicates[allSpecies]
  ];
  
  (* Extract species (now lowercase) *)
  spe = extSpe[reactions];
  
  (* Check if species were found *)
  If[Length[spe] == 0,
    Print["Error: No species found in reactions"];
    Return[$Failed]
  ];
  
  numSpecies = Length[spe];
  numReactions = Length[reactions];
  
  (* Initialize stoichiometric matrices *)
  al = ConstantArray[0, {numSpecies, numReactions}];
  be = ConstantArray[0, {numSpecies, numReactions}];
  
  (* Build stoichiometric matrices using original species names *)
  Do[
    reactants = compToAsso[reactions[[j, 1]]];
    products = compToAsso[reactions[[j, 2]]];
    
    (* Fill matrices - match original species *)
    Do[
      If[KeyExistsQ[reactants, speOrig[[i]]], 
        al[[i, j]] = reactants[speOrig[[i]]]];
      If[KeyExistsQ[products, speOrig[[i]]], 
        be[[i, j]] = products[speOrig[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  (* Calculate derived matrices *)
  gamma = be - al;
  
  (* Convert lowercase species to variables *)
  var = ToExpression /@ spe;
  
  (* Calculate reaction rates *)
  rv = Table[
    Product[
      If[al[[i, j]] == 0, 1, var[[i]]^al[[i, j]]], 
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
  
  (* Deficiency *)
  Nc = Length[complexes];
  l = 1;
  s = If[gamma === {} || AllTrue[Flatten[gamma], # == 0 &], 0, MatrixRank[gamma]];
  
  defFormula = "\[Delta] = Nc - \[ScriptL] - s";
  defTerms = "Nc = " <> ToString[Nc] <> " (complexes), " <>
             "\[ScriptL] = " <> ToString[l] <> " (linkage classes), " <> 
             "s = " <> ToString[s] <> " (stoich dimension)";
  defResult = "\[Delta] = " <> ToString[Nc] <> " - " <> ToString[l] <> 
              " - " <> ToString[s] <> " = " <> ToString[Nc - l - s];
  
  (* Minimal siphons *)
  mSi = minSiph[spe, asoRea[reactions]];
  
  (* Return with lowercase species *)
  {spe, al, be, gamma, Rv, RHS, mSi, {defFormula, defTerms, defResult}}
];

(* Test 
Print["=== Test with Capital Letters ==="];
RN = {"S" + "I1" -> 2*"I1", "S" + "I2" -> 2*"I2", 
      0 -> "S", "S" -> 0, "I1" -> 0, "I2" -> 0};
result = extMat[RN];
Print["Lowercase species: ", result[[1]]];
Print["Gamma matrix:"];
Print[MatrixForm[result[[4]], TableHeadings -> {result[[1]], Range[Length[RN]]}]];*)


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
    reactants = compToAsso[leftSide];
    products = compToAsso[rightSide];
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
(*Testing example-uncomment to test core package*)
(*RN={0->"S1","S1"+"I1"->2*"I1","I1"->0};
result=extMat[RN];
Print["Species: ",result[[1]]];
Print["Alpha matrix: ",result[[2]]//MatrixForm];
Print["Gamma matrix: ",result[[4]]//MatrixForm];*)


(* ---------------------------------------------------------------------
   findCores[RN, opts] : find autocatalytic cores (minimal autocatalytic subnetworks)
   Uses EpidCRN functions: extMat, compToAsso, extSpe
   --------------------------------------------------------------------- *)

ClearAll[findInternalReactions, isAutocatalyticBlockLP, findCores];

(* Find internal reactions for candidate T using compToAsso *)
findInternalReactions[RN_List, T_List] := Module[
  {Tstr = ToLowerCase /@ (ToString /@ T), intReacIdx},
  
  intReacIdx = Table[
    Module[{reactants, reactantSpe, isInternal},
      reactants = compToAsso[RN[[j, 1]]];
      reactantSpe = Keys[reactants];
      
      (* Reaction is internal if all reactants are in T (or no reactants) *)
      isInternal = (reactantSpe === {} || SubsetQ[Tstr, reactantSpe]);
      
      If[$DebugCores, 
        Print["Rxn ", j, ": ", RN[[j]], " -> reactants=", reactantSpe, 
          ", T=", Tstr, ", internal=", isInternal]
      ];
      
      If[isInternal, j, Nothing]
    ],
    {j, Length[RN]}
  ];
  
  intReacIdx
];

(* LP-based test: check if T is autocatalytic *)
isAutocatalyticBlockLP[gamma_, speciesIdx_List, intReacIdx_List, vmin_: 10^-9] := Module[
  {gammaT, m, n, A, b, c, sol, tsol, vsol},
  
  If[Length[intReacIdx] == 0 || Length[speciesIdx] == 0, 
    Return[<|"isCore" -> False, "reason" -> "empty block"|>]
  ];
  
  (* Extract submatrix *)
  gammaT = If[Length[speciesIdx] == 1, 
    {gamma[[speciesIdx[[1]], intReacIdx]]}, 
    gamma[[speciesIdx, intReacIdx]]
  ];
  
  m = Length[speciesIdx]; 
  n = Length[intReacIdx];
  
  (* LP: maximize t subject to gammaT\[CenterDot]v >= t\[CenterDot]1, v >= vmin, t >= 0 *)
  A = Join[
    Map[Join[#, {-1}] &, gammaT],
    Table[Join[UnitVector[n, j], {0}], {j, n}],
    {Join[ConstantArray[0, n], {1}]}
  ];
  
  b = Join[ConstantArray[0, m], ConstantArray[vmin, n], {0}];
  c = Join[ConstantArray[0, n], {-1}];
  
  sol = Quiet@LinearProgramming[Rationalize[c, 0], Rationalize[A, 0], Rationalize[b, 0]];
  
  If[sol === $Failed || sol === {} || !ListQ[sol] || Length[sol] != n + 1, 
    Return[<|"isCore" -> False, "reason" -> "LP failed"|>]
  ];
  
  tsol = N[sol[[n + 1]]];
  vsol = N[Take[sol, n]];
  
  If[NumericQ[tsol] && tsol > 10^-12,
    <|"isCore" -> True, "t" -> tsol, "flux" -> AssociationThread[intReacIdx -> vsol]|>,
    <|"isCore" -> False, "reason" -> "t <= 0"|>
  ]
];"flux" -> AssociationThread[intReacIdx -> vsol], 
      "gammaT" -> gammaT,
      "netProd" -> gammaT . vsol|>,
    <|"isCore" -> False, 
      "reason" -> "max t <= 0", 
      "t" -> tsol|>
  ]
];

(* Main function *)
Options[findCores] = {"CandidateSets" -> Automatic, "MaxSize" -> 6};

findCores[RN_List, opts : OptionsPattern[]] := Module[
  {RND, spe, gamma, nSpec, nReac, candidateSets, maxSize, 
   speciesIdxMap, results, cores, minimal},
  
  (* Use EpidCRN extMat *)
  RND = extMat[RN];
  spe = RND[[1]];  (* lowercase species from extSpe *)
  gamma = RND[[4]];
  
  nSpec = Length[spe]; 
  nReac = Length[RN];
  
  (* Map species names to indices *)
  speciesIdxMap = AssociationThread[spe -> Range[nSpec]];
  
  (* Determine candidate sets *)
  candidateSets = OptionValue["CandidateSets"];
  maxSize = OptionValue["MaxSize"];
  
  If[candidateSets === Automatic,
    (* All nonempty subsets up to maxSize *)
    candidateSets = Rest@Subsets[spe, {1, Min[maxSize, nSpec]}],
    (* Convert candidates to lowercase strings *)
    candidateSets = Map[ToLowerCase /@ (ToString /@ #) &, candidateSets]
  ];
  
  (* Test each candidate *)
  results = Table[
    Module[{T = candidateSets[[k]], intReacIdx, speciesIdx, lpres},
      
      (* Find internal reactions *)
      intReacIdx = findInternalReactions[RN, T];
      
      If[Length[intReacIdx] == 0, 
        <|"T" -> T, "isCore" -> False, "reason" -> "no internal reactions"|>,
        
        (* Get species indices *)
        speciesIdx = Lookup[speciesIdxMap, T, Nothing];
        
        If[Length[speciesIdx] != Length[T],
          <|"T" -> T, "isCore" -> False, "reason" -> "species not found"|>,
          
          (* LP test *)
          lpres = isAutocatalyticBlockLP[gamma, speciesIdx, intReacIdx];
          <|"T" -> T, "isCore" -> lpres["isCore"], 
            "intReacIdx" -> intReacIdx, "details" -> lpres|>
        ]
      ]
    ],
    {k, Length[candidateSets]}
  ];
  
  (* Extract cores *)
  cores = Select[results, #["isCore"] === True &];
  
  (* Find minimal cores (no proper subset is also a core) *)
  minimal = Select[cores, Function[c,
    Not[MemberQ[DeleteCases[cores, c], 
      d_ /; SubsetQ[c["T"], d["T"]] && d["T"] =!= c["T"]]]
  ]];
  
  <|"species" -> spe, 
    "gamma" -> gamma,
    "tested" -> Length[candidateSets],
    "allResults" -> results, 
    "cores" -> minimal|>
];

(* ==================== Test ==================== *)
RN1 = {"a" -> "a" + "b", "b" -> "a" + "b", "a" -> 0, "b" -> 0};
res1 = findCores[RN1, "CandidateSets" -> {{"a", "b"}}];
res1["cores"]
