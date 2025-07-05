(* ::Package:: *)

(* ::Package:: *)

(* Endotactic Reaction Network Analysis Package *)
(* Based on "A Geometric Approach to the Global Attractor Conjecture" *)
(* by Gopalkrishnan, Miller, and Shiu *)

BeginPackage["EndotacticAnalysis`"];

(* Public function declarations *)
endo::usage = "endo[reactions] analyzes a reaction network for endotactic properties. Options: \"ShowPlot\" -> True/False/Automatic (default: Automatic), \"Verbose\" -> True/False (default: False).";
diagnostics::usage = "diagnostics[reactions] provides detailed network analysis including deficiency.";
EucFHJ::usage = "EucFHJ[reactions] creates Euclidean Feinberg-Horn-Jackson graph visualization for 2D networks.";
extSpe::usage = "extSpe[reactions] extracts all species from a list of reactions.";
getComplexes::usage = "getComplexes[reactions] extracts all complexes from a list of reactions.";
isEndotactic::usage = "isEndotactic[reactions, speciesList] checks if network is endotactic.";
isStronglyEndotactic::usage = "isStronglyEndotactic[reactions, speciesList] checks if network is strongly endotactic.";

Begin["`Private`"];

(* Helper function to parse complex expressions *)
parseComplex[complex_String] := Module[{terms, result},
  If[complex == "0", Return[<||>]];
  terms = StringSplit[complex, "+"];
  result = <||>;
  Do[
    Module[{term, coeff, species},
      term = StringTrim[terms[[i]]];
      If[StringMatchQ[term, DigitCharacter.. ~~ LetterCharacter..],
        coeff = ToExpression[StringTake[term, {1, StringPosition[term, LetterCharacter][[1, 1]] - 1}]];
        species = StringDrop[term, StringPosition[term, LetterCharacter][[1, 1]] - 1];,
        If[StringMatchQ[term, LetterCharacter..],
          coeff = 1;
          species = term;,
          coeff = 0; species = "";
        ]
      ];
      If[species != "", result[species] = coeff];
    ], {i, Length[terms]}
  ];
  result
];

parseComplex[complex_] := parseComplex[ToString[complex]];

(* Extract all species from a list of reactions *)
extSpe[reactions_] := Module[{allSpecies},
  allSpecies = {};
  Do[
    Module[{reactant, product, reactantSpecies, productSpecies},
      reactant = parseComplex[reactions[[i]] /. Rule -> List // First];
      product = parseComplex[reactions[[i]] /. Rule -> List // Last];
      reactantSpecies = Keys[reactant];
      productSpecies = Keys[product];
      allSpecies = Union[allSpecies, reactantSpecies, productSpecies];
    ], {i, Length[reactions]}
  ];
  allSpecies
];

(* Convert complex to vector in species space *)
complexToVector[complex_, speciesList_] := Module[{parsed, vector},
  parsed = parseComplex[complex];
  vector = Table[If[KeyExistsQ[parsed, speciesList[[i]]], parsed[speciesList[[i]]], 0], {i, Length[speciesList]}];
  vector
];

(* Compute reaction vector (product - reactant) *)
reactionVector[reaction_, speciesList_] := Module[{reactantVec, productVec},
  reactantVec = complexToVector[reaction /. Rule -> List // First, speciesList];
  productVec = complexToVector[reaction /. Rule -> List // Last, speciesList];
  productVec - reactantVec
];

(* Compute stoichiometric subspace (span of reaction vectors) *)
stoichiometricSubspace[reactions_, speciesList_] := Module[{reactionVectors},
  reactionVectors = Table[reactionVector[reactions[[i]], speciesList], {i, Length[reactions]}];
  reactionVectors = DeleteDuplicates[Select[reactionVectors, # != Table[0, {Length[speciesList]}] &]];
  If[Length[reactionVectors] == 0, {}, reactionVectors]
];

(* Check if vector is orthogonal to stoichiometric subspace *)
orthogonalToStoichiometric[w_, stoichSubspace_] := Module[{},
  If[Length[stoichSubspace] == 0, True,
    AllTrue[stoichSubspace, w . # == 0 &]
  ]
];

(* Get w-maximal elements from a set of vectors *)
getMaximalElements[vectors_, w_] := Module[{dotProducts, maxValue},
  If[Length[vectors] == 0, Return[{}]];
  dotProducts = Table[w . vectors[[i]], {i, Length[vectors]}];
  maxValue = Max[dotProducts];
  Pick[vectors, dotProducts, maxValue]
];

(* Get all complexes from reactions *)
getComplexes[reactions_] := Module[{complexes},
  complexes = {};
  Do[
    Module[{reactant, product},
      reactant = reactions[[i]] /. Rule -> List // First;
      product = reactions[[i]] /. Rule -> List // Last;
      If[!MemberQ[complexes, reactant], AppendTo[complexes, reactant]];
      If[!MemberQ[complexes, product], AppendTo[complexes, product]];
    ], {i, Length[reactions]}
  ];
  complexes
];

(* Calculate number of linkage classes (connected components) *)
getLinkageClasses[reactions_] := Module[{complexes, adjMatrix, components},
  complexes = getComplexes[reactions];
  If[Length[complexes] == 0, Return[0]];
  
  (* Build adjacency matrix *)
  adjMatrix = ConstantArray[0, {Length[complexes], Length[complexes]}];
  Do[
    Module[{reactant, product, reactantPos, productPos},
      reactant = reactions[[i]] /. Rule -> List // First;
      product = reactions[[i]] /. Rule -> List // Last;
      reactantPos = Position[complexes, reactant][[1, 1]];
      productPos = Position[complexes, product][[1, 1]];
      adjMatrix[[reactantPos, productPos]] = 1;
      adjMatrix[[productPos, reactantPos]] = 1; (* Make undirected for linkage classes *)
    ], {i, Length[reactions]}
  ];
  
  (* Find connected components *)
  components = ConnectedComponents[AdjacencyGraph[adjMatrix]];
  Length[components]
];

(* Calculate deficiency: \[Delta] = n - \[ScriptL] - s *)
calculateDeficiency[reactions_, speciesList_] := Module[{n, l, s},
  n = Length[getComplexes[reactions]]; (* number of complexes *)
  l = getLinkageClasses[reactions]; (* number of linkage classes *)
  s = If[Length[stoichiometricSubspace[reactions, speciesList]] == 0, 0, 
        MatrixRank[stoichiometricSubspace[reactions, speciesList]]]; (* dimension of stoichiometric subspace *)
  n - l - s
];
getReactantVectors[reactions_, speciesList_] := Module[{reactants},
  reactants = Table[complexToVector[reactions[[i]] /. Rule -> List // First, speciesList], {i, Length[reactions]}];
  DeleteDuplicates[reactants]
];

(* Check if reaction is w-essential *)
isWEssential[reaction_, w_, speciesList_] := Module[{reacVec},
  reacVec = reactionVector[reaction, speciesList];
  w . reacVec != 0
];

(* Get w-essential reactions *)
getWEssentialReactions[reactions_, w_, speciesList_] := Module[{},
  Select[reactions, isWEssential[#, w, speciesList] &]
];

(* Get w-support (maximal reactants among w-essential reactions) *)
getWSupport[reactions_, w_, speciesList_] := Module[{essentialReactions, essentialReactants},
  essentialReactions = getWEssentialReactions[reactions, w, speciesList];
  If[Length[essentialReactions] == 0, Return[{}]];
  essentialReactants = Table[complexToVector[essentialReactions[[i]] /. Rule -> List // First, speciesList], {i, Length[essentialReactions]}];
  getMaximalElements[essentialReactants, w]
];

(* Check if network is w-endotactic *)
isWEndotactic[reactions_, w_, speciesList_] := Module[{essentialReactions, wSupport},
  essentialReactions = getWEssentialReactions[reactions, w, speciesList];
  wSupport = getWSupport[reactions, w, speciesList];
  
  AllTrue[essentialReactions, 
    Module[{reactantVec, reacVec, inSupport},
      reactantVec = complexToVector[# /. Rule -> List // First, speciesList];
      reacVec = reactionVector[#, speciesList];
      inSupport = MemberQ[wSupport, reactantVec];
      If[inSupport, w . reacVec < 0, True]
    ] &
  ]
];

(* Generate test directions for endotactic check *)
generateTestDirections[speciesList_] := Module[{n, directions},
  n = Length[speciesList];
  directions = {};
  
  (* Standard basis vectors and their negatives *)
  Do[
    AppendTo[directions, UnitVector[n, i]];
    AppendTo[directions, -UnitVector[n, i]];
  , {i, n}];
  
  (* Some random directions *)
  Do[
    AppendTo[directions, Normalize[RandomReal[{-1, 1}, n]]];
  , {20}];
  
  (* Some systematic combinations *)
  If[n >= 2,
    Do[
      Do[
        AppendTo[directions, Normalize[UnitVector[n, i] + UnitVector[n, j]]];
        AppendTo[directions, Normalize[UnitVector[n, i] - UnitVector[n, j]]];
      , {j, i + 1, n}];
    , {i, 1, n - 1}];
  ];
  
  DeleteDuplicates[directions, Norm[#1 - #2] < 10^(-10) &]
];

(* Check if network is endotactic *)
isEndotactic[reactions_, speciesList_] := Module[{testDirections},
  testDirections = generateTestDirections[speciesList];
  AllTrue[testDirections, isWEndotactic[reactions, #, speciesList] &]
];

(* Check strong endotactic condition *)
isStronglyEndotactic[reactions_, speciesList_] := Module[{stoichSubspace, testDirections},
  If[!isEndotactic[reactions, speciesList], Return[False]];
  
  stoichSubspace = stoichiometricSubspace[reactions, speciesList];
  testDirections = Select[generateTestDirections[speciesList], !orthogonalToStoichiometric[#, stoichSubspace] &];
  
  AllTrue[testDirections,
    Module[{w, reactantVectors, maximalReactants, hasGoodReaction},
      w = #;
      reactantVectors = getReactantVectors[reactions, speciesList];
      maximalReactants = getMaximalElements[reactantVectors, w];
      
      hasGoodReaction = False;
      Do[
        Module[{maxReactant, candidateReactions},
          maxReactant = maximalReactants[[i]];
          candidateReactions = Select[reactions, complexToVector[# /. Rule -> List // First, speciesList] == maxReactant &];
          Do[
            Module[{reacVec},
              reacVec = reactionVector[candidateReactions[[j]], speciesList];
              If[w . reacVec < 0, hasGoodReaction = True; Break[]];
            ], {j, Length[candidateReactions]}];
          If[hasGoodReaction, Break[]];
        ], {i, Length[maximalReactants]}];
      hasGoodReaction
    ] &
  ]
];

(* Newton Polytope function drawer *)
NewtPol[complexCoords_Association, opts___] := Module[{points, convexHull, polytope},
  points = Values[complexCoords];
  If[Length[points] >= 3 && Length[Dimensions[points]] == 2,
    convexHull = ConvexHullMesh[points];
    polytope = {LightGray, EdgeForm[{Black, Thin}], convexHull},
    polytope = {}
  ];
  polytope
];

(* EucFHJ: Euclidean Feinberg-Horn-Jackson Graph Drawer *)
EucFHJ[reactions_, opts___] := Module[{species, complexes, reactionPairs, coords, vertices, edges, labels, polytope, result},
  species = extSpe[reactions];
  If[Length[species] > 2, species = Take[species, 2]];
  If[Length[species] < 2, species = Join[species, Table["dummy" <> ToString[i], {i, 2 - Length[species]}]]];
  
  complexes = {};
  reactionPairs = {};
  Do[
    Module[{reactant, product},
      reactant = reactions[[i]] /. Rule -> List // First;
      product = reactions[[i]] /. Rule -> List // Last;
      If[!MemberQ[complexes, reactant], AppendTo[complexes, reactant]];
      If[!MemberQ[complexes, product], AppendTo[complexes, product]];
      AppendTo[reactionPairs, {reactant, product}];
    ], {i, Length[reactions]}
  ];
  
  coords = <||>;
  Do[
    Module[{prodAssoc, coord},
      prodAssoc = parseComplex[complex];
      coord = Table[If[KeyExistsQ[prodAssoc, species[[j]]], prodAssoc[species[[j]]], 0], {j, 2}];
      coords[complex] = coord;
    ], {complex, complexes}
  ];
  
  polytope = NewtPol[coords];
  vertices = Table[{PointSize[Large], Red, Point[coords[complex]]}, {complex, complexes}];
  labels = Table[Text[Style[ToString[complex], 12, Black], coords[complex], {-1.5, -1.5}], {complex, complexes}];
  edges = Table[
    Module[{start, end},
      start = coords[pair[[1]]];
      end = coords[pair[[2]]];
      If[start == end,
        {Blue, Red, Arrowheads[Large], Arrow[{start, start + {0.2, 0.2}}]},
        {Blue, Red, Arrowheads[Large], Arrow[{start, end}]}
      ]
    ], {pair, reactionPairs}
  ];
  
  result = Graphics[{polytope, edges, vertices, labels}, 
    PlotRange -> All, AspectRatio -> 1, Axes -> True, 
    AxesLabel -> {ToString[species[[1]]], ToString[species[[2]]]}, 
    PlotLabel -> "Euclidean Feinberg-Horn-Jackson Graph", 
    ImageSize -> 300, GridLines -> Automatic, GridLinesStyle -> LightGray, opts];
  result
];

(* Main function: check endotactic properties *)
Options[endo] = {"ShowPlot" -> Automatic, "Verbose" -> False};

endo[reactions_, OptionsPattern[]] := Module[{speciesList, isEndo, isStrongEndo, result, showPlot, verbose},
  speciesList = extSpe[reactions];
  showPlot = OptionValue["ShowPlot"];
  verbose = OptionValue["Verbose"];
  
  If[verbose,
    Print["Analyzing network with species: ", speciesList];
    Print["Reactions: ", reactions];
  ];
  
  isEndo = isEndotactic[reactions, speciesList];
  isStrongEndo = If[isEndo, isStronglyEndotactic[reactions, speciesList], False];
  
  result = <|
    "Species" -> speciesList,
    "Reactions" -> reactions,
    "Endotactic" -> isEndo,
    "StronglyEndotactic" -> isStrongEndo,
    "StoichiometricSubspace" -> stoichiometricSubspace[reactions, speciesList]
  |>;
  
  If[verbose,
    Print["Results:"];
    Print["  Endotactic: ", isEndo];
    Print["  Strongly Endotactic: ", isStrongEndo];
  ];
  
  (* Show plot for 2D cases *)
  If[(showPlot === True || (showPlot === Automatic && Length[speciesList] == 2)) && Length[speciesList] <= 2,
    If[verbose, Print[""]];
    Print[EucFHJ[reactions]];
  ];
  
  result
];

(* Diagnostic tools *)
diagnostics[reactions_] := Module[{species, stoichSubspace, reactants, complexes, linkageClasses, deficiency},
  species = extSpe[reactions];
  stoichSubspace = stoichiometricSubspace[reactions, species];
  reactants = getReactantVectors[reactions, species];
  complexes = getComplexes[reactions];
  linkageClasses = getLinkageClasses[reactions];
  deficiency = calculateDeficiency[reactions, species];
  
  Print["=== NETWORK DIAGNOSTICS ==="];
  Print["Species: ", species];
  Print["Number of reactions: ", Length[reactions]];
  Print["Number of complexes: ", Length[complexes]];
  Print["Number of linkage classes: ", linkageClasses];
  Print["Stoichiometric subspace dimension: ", If[Length[stoichSubspace] == 0, 0, MatrixRank[stoichSubspace]]];
  Print["Deficiency: ", deficiency];
  Print["Reactant vectors: ", reactants];
  Print[""];
];

End[];
EndPackage[];

