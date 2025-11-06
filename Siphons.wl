(* ::Package:: *)

(* Correct implementations following Gagrani's definitions *)

ClearAll[isSiph, isAutoC, minSiph, findCores];

(* Check if a reaction is an inflow (0 \[RightArrow] ...) or outflow (... \[RightArrow] 0) *)
isInflowOrOutflow[reaction_] := 
  reaction[[1]] === 0 || reaction[[2]] === 0;

(* Filter out inflows and outflows, keeping only internal reactions *)
filterInternalReactions[RN_List] := 
  Select[RN, !isInflowOrOutflow[#] &];

(* 
  NAME CHANGES FOR EXISTING EpidCRN FUNCTIONS:
  
  OLD: isSiph[W, RN] 
  NEW: isSiph[W, alpha, beta]
  - Now takes stoichiometric matrices directly instead of RN
  - Must call on filtered reactions: RNinternal = filterInternalReactions[RN]; RND = extMat[RNinternal]; isSiph[W, RND[[2]], RND[[3]]]
  
  OLD: isAutoC[W, R, RN]
  NEW: isAutoC[W, R, alpha, beta]
  - Now takes stoichiometric matrices directly instead of RN
  - Must call on filtered reactions
  
  OLD: minSiph[RN]
  NEW: minSiph[RN]
  - Signature unchanged but NOW AUTOMATICALLY FILTERS inflows/outflows
  
  OLD: findCores[RN, maxSize]
  NEW: findCores[RN, maxSize]
  - Signature unchanged but NOW AUTOMATICALLY FILTERS inflows/outflows
  - Returns {minimal, nonMinimal} as before
*)



isSiph[W_List, species_List, alpha_, beta_] := Module[
  {indices, alphaW, betaW, nReac, producingReactions, result},
  
  indices = Flatten[Position[species, #] & /@ W];
  If[Length[indices] != Length[W], Return[False]];
  
  nReac = Dimensions[alpha][[2]];
  alphaW = alpha[[indices]];
  betaW = beta[[indices]];
  
  (* Reactions where (S+)_W > 0 *)
  producingReactions = Select[Range[nReac], pos[betaW[[All, #]]] &];
  
  (* Debug output *)
  If[W == {"a"}, 
    Print["Testing W = ", W];
    Print["  Producing reactions: ", producingReactions];
    Print["  For each producing reaction, check alphaW:"];
    Do[
      Print["    Reaction ", r, ": alphaW = ", alphaW[[All, r]], 
            ", pos = ", pos[alphaW[[All, r]]]],
      {r, producingReactions}
    ];
  ];
  
  If[Length[producingReactions] == 0, Return[True]];
  
  result = AllTrue[producingReactions, pos[alphaW[[All, #]]] &];
  
  If[W == {"a"}, Print["  Result: ", result]];
  
  result
];

(* Siphon test - W is list of species names *)
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

(* Find minimal siphons *)
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


(* Helper: check if vector is >= 0 and != 0 *)
pos[v_] := AllTrue[v, # >= 0 &] && Total[v] > 0;

(* Check if W with reactions R forms an autocatalytic core *)
isAutoC[W_List, R_List, species_List, reactions_List, alpha_, beta_] := Module[
  {wIndices, rIndices, alphaWR, betaWR, p1, p2},
  
  If[Length[W] == 0 || Length[R] == 0, Return[False]];
  
  wIndices = Flatten[Position[species, #] & /@ W];
  rIndices = Flatten[Position[reactions, #] & /@ R];
  
  If[Length[wIndices] != Length[W] || Length[rIndices] != Length[R], 
    Return[False]
  ];
  
  alphaWR = alpha[[wIndices, rIndices]];
  betaWR = beta[[wIndices, rIndices]];
  
  (* P1: Every reaction consumes AND produces in W *)
  p1 = AllTrue[Range[Length[rIndices]], Function[j,
    pos[alphaWR[[All, j]]] && pos[betaWR[[All, j]]]
  ]];
  
  (* P2: Every species in W is consumed AND produced *)
  p2 = AllTrue[Range[Length[wIndices]], Function[i,
    pos[alphaWR[[i, All]]] && pos[betaWR[[i, All]]]
  ]];
  
  p1 && p2
];

(* Find autocatalytic cores *)
findCores[RN_, maxSize_] := Module[{
  species, alpha, beta, n, m, allSubsets, autocatalyticSets, minimal, nonMinimal},
  
  species = extSpe[RN];
  n = Length[species];
  m = Length[RN];
  
  alpha = ConstantArray[0, {n, m}];
  beta = ConstantArray[0, {n, m}];
  
  Do[
    Module[{subsAsso, prodsAsso},
      subsAsso = comp2Asso[RN[[j, 1]]];
      prodsAsso = comp2Asso[RN[[j, 2]]];
      
      Do[
        If[KeyExistsQ[subsAsso, species[[i]]], 
          alpha[[i, j]] = subsAsso[species[[i]]]
        ];
        If[KeyExistsQ[prodsAsso, species[[i]]], 
          beta[[i, j]] = prodsAsso[species[[i]]]
        ];
      , {i, n}];
    ];
  , {j, m}];
  
  allSubsets = Subsets[species, {1, Min[maxSize, n]}];
  
  autocatalyticSets = Select[allSubsets, Function[W,
    Module[{allReactionSubsets},
      allReactionSubsets = Subsets[RN, {1, m}];
      AnyTrue[allReactionSubsets, 
        isAutoC[W, #, species, RN, alpha, beta] &
      ]
    ]
  ]];
  
  minimal = Select[autocatalyticSets, 
    Function[s, AllTrue[autocatalyticSets, Function[t, t === s || !SubsetQ[s, t]]]]
  ];
  
  nonMinimal = Complement[autocatalyticSets, minimal];
  
  {minimal, nonMinimal}
];

(* Test with two-strain epidemic model 
Print["=== Two-strain epidemic with coinfection ==="];
RN = {0 -> "S", "S" + "I1" -> 2*"I1", "S" + "I2" -> 2*"I2", 
      "I1" + "I2" -> "I1" + "I12", "I1" + "I2" -> "I12" + "I2",
      "I1" + "I12" -> 2*"I12", "I12" + "I2" -> 2*"I12",
      "S" -> 0, "I1" -> 0, "I2" -> 0, "I12" -> 0};
result = findCores[RN, 3];
Print["Minimal cores: ", result[[1]]];
Print["Non-minimal cores: ", result[[2]]];*)


(* =================================================================== *)
(* IGMS with visible edge labels ABOVE edges (noninvasive to comp2Asso) *)
(* =================================================================== *)

ClearAll[safeLowerKeys, canActRea, edgIGMS, IGMS];

(* Map association keys to lowercase strings, preserving values *)
safeLowerKeys[asso_Association] := Association @ KeyValueMap[ToLowerCase[ToString[#1]] -> #2 &, asso];

(* reaction: a Rule lhs->rhs (your RN format)
   Si, Sj: lists of species names (strings) \[LongDash] can be mixed case.
   Uses your original comp2Asso (do NOT modify it elsewhere). *)
canActRea[reaction_, Si_List, Sj_List] := Module[
  {lhs, rhs, lhsComp, rhsComp, lhsL, rhsL, SiL, SjL, newSj},

  (* skip pure birth/death *)
  If[MatchQ[reaction, 0 -> _] || MatchQ[reaction, _ -> 0], Return[False]];

  {lhs, rhs} = List @@ reaction;

  (* build stoichiometric maps via your comp2Asso, then lowercase keys locally *)
  lhsComp = comp2Asso[lhs];  rhsComp = comp2Asso[rhs];
  lhsL = safeLowerKeys @ lhsComp;  rhsL = safeLowerKeys @ rhsComp;

  (* lowercase versions of siphon species *)
  SiL = ToLowerCase /@ ToString /@ Si;
  SjL = ToLowerCase /@ ToString /@ Sj;

  (* must consume something in Si *)
  If[Intersection[Keys[lhsL], SiL] === {}, Return[False]];

  (* net positive production of at least one new Sj \ Si species *)
  newSj = Complement[SjL, SiL];
  If[newSj === {}, Return[False]];

  AnyTrue[newSj, (Lookup[rhsL, #, 0] - Lookup[lhsL, #, 0]) > 0 &]
];

(* Build directed IGMS edges + aligned label rules for EdgeLabels *)
edgIGMS[RN_List, mSi_List] := Module[
  {n = Length[mSi], edges = {}, firstLabels = {}, allLabels = <||>, SiL, SjL, hits},
  Do[
    If[i =!= j,
      SiL = ToLowerCase /@ ToString /@ mSi[[i]];
      SjL = ToLowerCase /@ ToString /@ mSi[[j]];
      hits = Select[RN, canActRea[#, SiL, SjL] &];
      If[hits =!= {},
        AppendTo[edges, i \[DirectedEdge] j];
        AppendTo[firstLabels, ToString[First[hits], InputForm]];
        allLabels[i \[DirectedEdge] j] = ToString[#, InputForm] & /@ hits;
      ];
    ],
    {i, n}, {j, n}
  ];
  <|
    "Edges" -> edges,
    "EdgeLabelRules" -> MapThread[#1 -> Placed[Style[#2, Small], Above] &, {edges, firstLabels}],
    "EdgeAllLabels" -> allLabels
  |>
];

(* Draw graph; show nothing when empty *)
IGMS[RN_List, mSi_List] := Module[
  {n = Length[mSi], edata, edges, elabs, g, cycles, pairs, mSiStr, labText,gr,
  embedding ,
scaledEmbedding},
  Print["Minimal siphons: ", Table[Subscript["T", i] -> mSi[[i]], {i, n}]];
  edata = edgIGMS[RN, mSi];
  edges = edata["Edges"];  elabs = edata["EdgeLabelRules"];

  If[edges === {},
    Print["IGMS edges: {}"];
    Print["empty graph"];
    Return[{{}, {}, {}, None}]
  ];

  Print["IGMS edges: ", edges];

  (* build label text explicitly as in IGMSL *)
  pairs = (List @@ #) & /@ edges;
  mSiStr = Map[ToLowerCase@*ToString, mSi, {2}];
  labText = Table[
    Module[{i = pairs[[k, 1]], j = pairs[[k, 2]], hits},
      hits = Select[RN, canActRea[#, mSiStr[[i]], mSiStr[[j]]] &];
      If[hits === {}, "", ToString[First[hits], InputForm]]
    ],
    {k, Length[pairs]}
  ];
  elabs = MapThread[#1 -> Placed[Style[#2, Blue, Bold, FontSize -> 12], 1/2] &, {edges, labText}];

  
  (* Build your initial graph to get the embedding *)
g = Graph[
    Range[n], edges,
    DirectedEdges -> True,
    EdgeLabels -> elabs,
    VertexLabels -> Table[
      i -> Placed[
        Row[{Subscript["T", i], " = ", Row[mSi[[i]], ","]}],
        Center
      ],
      {i, n}
    ],
    VertexSize -> .4,
    VertexStyle -> LightBlue, EdgeStyle -> Black,
    VertexLabelStyle -> Directive[FontSize -> 12],
    GraphLayout -> "SpringElectricalEmbedding",
    ImageSize -> Medium
];

(* Obtain and scale the embedding *)
embedding = GraphEmbedding[g];
scaledEmbedding = 4 embedding; (* scale by factor, e.g., 3 *)

(* Now define the graph using only VertexCoordinates *)
gr=Graph[
    Range[n], edges,
    DirectedEdges -> True,
    EdgeLabels -> elabs,
    VertexLabels -> Table[
      i -> Placed[
        Row[{Subscript["T", i], " = ", Row[mSi[[i]], ","]}],
        Center
      ],
      {i, n}
    ],
    VertexSize -> .4,
    VertexStyle -> LightBlue, EdgeStyle -> Black,
    VertexLabelStyle -> Directive[FontSize -> 12],
    VertexCoordinates -> scaledEmbedding,
    ImageSize -> Medium
];


  cycles = FindCycle[g, Infinity, All];
  Print[If[cycles === {}, "IGMS is acyclic", Row[{"Cycles found: ", Length[cycles]}]]];
  If[cycles =!= {}, Print["Cycles: ", cycles]];

  Print[gr];
  {edges, elabs, cycles, gr}
];



persistenceAnalysis[reactions_, OptionsPattern[{"Verbose" -> False}]] := 
  Module[{siphonClass, drainableSiphons, criticalDrainableSiphons, isPersistent, threats, analysis, verbose}, 
   verbose = OptionValue["Verbose"];
   If[verbose, Print["=== PERSISTENCE ANALYSIS ==="]];
   
   siphonClass = siphonAnalysis[reactions, "Verbose" -> False];
   drainableSiphons = Select[Keys[siphonClass], siphonClass[#]["IsDrainable"] &];
   criticalDrainableSiphons = Select[drainableSiphons, siphonClass[#]["IsCritical"] &];
   isPersistent = Length[drainableSiphons] == 0;
   threats = Table[siphon -> siphonClass[siphon]["Significance"], {siphon, criticalDrainableSiphons}];
   
   analysis = <|"Persistent" -> isPersistent, 
     "AllSiphons" -> Keys[siphonClass], 
     "DrainableSiphons" -> drainableSiphons, 
     "CriticalDrainableSiphons" -> criticalDrainableSiphons, 
     "NumberOfDrainableSiphons" -> Length[drainableSiphons], 
     "ExtinctionThreats" -> threats, 
     "TheoreticalBasis" -> "Networks without drainable siphons are persistent (Deshpande & Gopalkrishnan, 2014)", 
     "SiphonClassification" -> siphonClass|>;
   
   If[verbose, 
    Print["Persistent: ", isPersistent];
    Print["Total siphons: ", Length[Keys[siphonClass]]];
    Print["Drainable siphons: ", Length[drainableSiphons]];
    Print["Critical drainable siphons: ", Length[criticalDrainableSiphons]];
    If[Length[criticalDrainableSiphons] > 0, 
     Print["Main extinction threats: ", criticalDrainableSiphons];];];
   analysis];

catalysisAnalysis[reactions_, OptionsPattern[{"Verbose" -> False}]] := 
  Module[{siphonClass, selfReplicableSiphons, catalyticSets, strictlyCatalyticSets, autocatalyticPotential, analysis, verbose}, 
   verbose = OptionValue["Verbose"];
   If[verbose, Print["=== CATALYSIS ANALYSIS ==="]];
   
   siphonClass = siphonAnalysis[reactions, "Verbose" -> False];
   selfReplicableSiphons = Select[Keys[siphonClass], siphonClass[#]["IsSelfReplicable"] &];
   catalyticSets = Select[selfReplicableSiphons, siphonClass[#]["IsCritical"] &];
   strictlyCatalyticSets = Select[catalyticSets, !siphonClass[#]["IsDrainable"] &];
   
   autocatalyticPotential = Which[
     Length[catalyticSets] == 0, "No autocatalytic behavior detected", 
     Length[strictlyCatalyticSets] > 0, "Strong autocatalytic potential (pure growth)", 
     True, "Autocatalytic potential with competition (growth vs extinction)"];
   
   analysis = <|"HasCatalyticSets" -> Length[catalyticSets] > 0, 
     "CatalyticSets" -> catalyticSets, 
     "StrictlyCatalyticSets" -> strictlyCatalyticSets, 
     "SelfReplicableSiphons" -> selfReplicableSiphons, 
     "AutocatalyticPotential" -> autocatalyticPotential, 
     "TheoreticalBasis" -> "Self-replicable critical siphons correspond to catalytic sets (Deshpande & Gopalkrishnan, 2014)", 
     "SiphonClassification" -> siphonClass|>;
   
   If[verbose, 
    Print["Catalytic sets found: ", Length[catalyticSets]];
    Print["Strictly catalytic sets: ", Length[strictlyCatalyticSets]];
    Print["Autocatalytic potential: ", autocatalyticPotential];
    If[Length[catalyticSets] > 0, Print["Catalytic sets: ", catalyticSets];];];
   analysis];

autocatalysisReport[reactions_, OptionsPattern[{"Verbose" -> True}]] := 
  Module[{matResults, persistenceResults, catalysisResults, minimalCriticalSiphons, summary, verbose}, 
   verbose = OptionValue["Verbose"];
   If[verbose, Print["=== COMPREHENSIVE AUTOCATALYSIS ANALYSIS ==="]];
   
   matResults = extMat[reactions];
   persistenceResults = persistenceAnalysis[reactions, "Verbose" -> False];
   catalysisResults = catalysisAnalysis[reactions, "Verbose" -> False];
   minimalCriticalSiphons = findMinimalCriticalSiphons[reactions];
   
   summary = <|"NetworkInfo" -> <|"Species" -> matResults[[1]], 
       "NumberOfSpecies" -> Length[matResults[[1]]], 
       "NumberOfReactions" -> Length[reactions], 
       "Deficiency" -> matResults[[7]]|>, 
     "Persistence" -> persistenceResults, 
     "Catalysis" -> catalysisResults, 
     "MinimalCriticalSiphons" -> minimalCriticalSiphons|>;
   
   If[verbose, 
    Print["Network: ", Length[matResults[[1]]], " species, ", Length[reactions], " reactions"];
    Print["Persistent: ", persistenceResults["Persistent"]];
    Print["Has catalytic sets: ", catalysisResults["HasCatalyticSets"]];
    Print["Minimal critical siphons: ", Length[minimalCriticalSiphons]];];
   summary];

checkPersistence[RN_] := 
  Module[{RND, spe, allSiphons, drainableSiphons, analysis, persistenceStatus, competing}, 
   RND = extMat[RN];
   spe = RND[[1]];
   allSiphons = minSiph[spe, asoRea[RN]];
   drainableSiphons = Select[allSiphons, isDrainable[RN, #] &];
   competing = Select[drainableSiphons, isSelfReplicable[RN, #] &];
   
   persistenceStatus = Which[
     Length[drainableSiphons] == 0, "Persistent", 
     Length[competing] == Length[drainableSiphons], "Unknown", 
     True, "Non-persistent"];
   
   analysis = <|"AllSiphons" -> allSiphons, 
     "DrainableSiphons" -> drainableSiphons, 
     "CompetingSiphons" -> competing, 
     "PersistenceStatus" -> persistenceStatus|>;
   {persistenceStatus, analysis}];

persistenceReport[RN_] := 
  Module[{RND, spe, persistenceCheck, siphonClassification, minimalCriticalSiphons, result}, 
   RND = extMat[RN];
   spe = RND[[1]];
   persistenceCheck = checkPersistence[RN];
   siphonClassification = siphonAnalysis[RN];
   minimalCriticalSiphons = findMinimalCriticalSiphons[RN];
   
   result = <|"Species" -> spe, 
     "NumberOfSpecies" -> Length[spe], 
     "NumberOfReactions" -> Length[RN], 
     "Deficiency" -> RND[[7]], 
     "PersistenceStatus" -> persistenceCheck[[1]], 
     "PersistenceAnalysis" -> persistenceCheck[[2]], 
     "SiphonClassification" -> siphonClassification, 
     "MinimalCriticalSiphons" -> minimalCriticalSiphons|>;
   
   Print["=== PERSISTENCE ANALYSIS ==="];
   Print["Persistence: ", persistenceCheck[[1]]];
   If[persistenceCheck[[1]] == "Unknown", 
    Print["WARNING: Competing forces - drainable siphons are also self-replicable"];
    Print["This represents the central challenge of the persistence conjecture"];];
   Print[""];
   Print["Siphon Analysis:"];
   Print["  Drainable siphons: ", persistenceCheck[[2]]["DrainableSiphons"]];
   Print["  Competing siphons: ", persistenceCheck[[2]]["CompetingSiphons"]];
   Print[""];
   Print["Siphon Classification:"];
   Do[Print["  ", siphon, " -> ", siphonClassification[siphon]["Type"]], {siphon, Keys[siphonClassification]}];
   result];



cons[gamma_, {}] := 
  Module[{leftKernel, positiveConservationLaws, nullSpace, dims}, 
   dims = Dimensions[gamma];
   nullSpace = Quiet[NullSpace[Transpose[gamma]]];
   
   If[nullSpace === {} || nullSpace === $Failed, {}, 
    positiveConservationLaws = Select[nullSpace, AllTrue[#, # >= 0 &] && AnyTrue[#, # > 0 &] &];
    If[Length[positiveConservationLaws] == 0, 
     positiveConservationLaws = Select[nullSpace, AllTrue[-#, # >= 0 &] && AnyTrue[-#, # > 0 &] &];
     positiveConservationLaws = -# & /@ positiveConservationLaws;];
    positiveConservationLaws]];

isInvariantFacet[facetSet_, reactions_] := 
  Module[{vf, vars, facetIndices, facetRules, isInvariant, derivative, varSymbols}, 
   {vf, vars} = Take[reaToRHS[reactions], 2];
   varSymbols = ToExpression[vars];
   facetIndices = Flatten[Position[vars, #] & /@ facetSet];
   facetRules = Table[varSymbols[[facetIndices[[j]]]] -> 0, {j, Length[facetIndices]}];
   
   isInvariant = True;
   Do[derivative = Simplify[vf[[facetIndices[[j]]]] /. facetRules];
    derivative = derivative /. Power[0, _] :> 0;
    derivative = Simplify[derivative];
    If[!(PossibleZeroQ[derivative] || TrueQ[Simplify[derivative <= 0]] || 
        TrueQ[Simplify[derivative <= 0, Assumptions -> And @@ (# >= 0 & /@ varSymbols)]] || 
        MatchQ[derivative, _?NonPositive] || MatchQ[derivative, _?Negative] || derivative === 0), 
     isInvariant = False; Break[]];, {j, Length[facetIndices]}];
   isInvariant];

invFace[reactions_, maxCodim_:10] := 
  Module[{species, n, invariantFacets, subsets}, 
   species = extSpe[reactions];
   n = Length[species];
   invariantFacets = {};
   
   Do[subsets = Subsets[species, {k}];
    Do[candidateSet = subsets[[i]];
     If[isInvariantFacet[candidateSet, reactions], AppendTo[invariantFacets, candidateSet]];, {i, Length[subsets]}];, {k, 1, Min[maxCodim, n]}];
   invariantFacets];


(*(* =================================================================== *)
(* IGMS with visible edge labels (exact RN reaction) shown ABOVE edges *)
(* =================================================================== *)

ClearAll[reacToString, edgIGMS, IGMS];

(* Convert a reaction (lhs -> rhs) to a compact RN-style string *)
reacToString[r_Rule] := ToString[r, InputForm];

(* Build directed IGMS edges and a label association for EdgeLabels *)
edgIGMS[RN_List, mSi_List] := Module[
  {n = Length[mSi], mSiStr, edgeLabelAssoc, edgeAllLabels, edges},

  (* Your canActRea expects lower-case species strings *)
  mSiStr = Map[ToLowerCase @* ToString, mSi, {2}];

  edgeLabelAssoc = <||>;              (* edge -> first label (string) *)
  edgeAllLabels  = <||>;              (* edge -> list of all labels  *)

  Do[
    If[i =!= j,
      Module[{hitRxns, e, lbls},
        hitRxns = Select[RN, canActRea[#, mSiStr[[i]], mSiStr[[j]]] &];
        If[hitRxns =!= {},
          e = i \[DirectedEdge] j;
          lbls = reacToString /@ hitRxns;
          If[! KeyExistsQ[edgeLabelAssoc, e], edgeLabelAssoc[e] = First[lbls]];
          edgeAllLabels[e] = lbls;  (* keeps all labels; useful for tooltips *)
        ];
      ]
    ],
    {i, n}, {j, n}
  ];

  edges = Keys[edgeLabelAssoc];

  <|"Edges" -> edges,
    "EdgeLabelAssoc" -> AssociationMap[Placed[Style[edgeLabelAssoc[#], Small], Above] &, edges],
    "EdgeAllLabels" -> edgeAllLabels|>
];

(* Main IGMS: draws graph with edge labels ABOVE; prints "empty graph" if none *)
IGMS[RN_List, mSi_List] := Module[
  {n = Length[mSi], edata, edges, elabs, g, cycles},

  Print["Minimal siphons: ",
    Table[Subscript["T", i] -> mSi[[i]], {i, n}]
  ];

  edata = edgIGMS[RN, mSi];
  edges = edata["Edges"];
  elabs = edata["EdgeLabelAssoc"];

  If[edges === {},
    Print["IGMS edges: {}"];
    Print["empty graph"];
    Return[{{}, {}, None}, Module]
  ];

  Print["IGMS edges: ", edges];

  g = Graph[
    Range[n],
    edges,
    DirectedEdges -> True,
    EdgeLabels -> elabs,                     (* <- label shown ABOVE each edge *)
    VertexLabels -> Table[
      i -> Placed[Row[{Subscript["T", i], " = ", Row[mSi[[i]], ","]}], Center],
      {i, n}
    ],
    VertexSize -> 0.3,
    VertexStyle -> LightBlue,
    EdgeStyle -> Arrowheads[0.03],
    VertexLabelStyle -> Directive[FontSize -> 10],
    GraphLayout -> "SpringElectricalEmbedding",
    ImageSize -> Medium
  ];

  cycles = FindCycle[g, Infinity, All];
  Print[If[cycles === {}, "IGMS is acyclic", Row[{"Cycles found: ", Length[cycles]}]]];
  If[cycles =!= {}, Print["Cycles: ", cycles]];

  Print[g];

  {edges, cycles, g}
];

(* -------------------------------------------------------------------
Notes:
- Edge labels now appear because:
  (1) Edges are DirectedEdge[i,j], and
  (2) EdgeLabels is given an association mapping each edge to
      Placed["exact RN reaction", Above].
- If you want tooltips listing all reactions per edge, you can add an
  EdgeShapeFunction which wraps the default shape in a Tooltip using
  edata["EdgeAllLabels"][edge].
------------------------------------------------------------------- *)
*)


(* Convert compound expression to association with lowercase species names *)
(* Enhanced to handle symbolic coefficients using * separator *)
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

(* Tests 
Print["=== Test comp2Asso ==="];
Print[comp2Asso[0]];                              (* <||> *)                       (* <|"i2" -> 2|> *)
Print[comp2Asso["S" + 2*"I1" + 3*"I2"]];         (* <|"s" -> 1, "i1" -> 2, "i2" -> 3|> *)
Print[comp2Asso[(1 + a)*"x1" + (1 + b)*"x2"]];  (* <|"x1" -> 1+a, "x2" -> 1+b|> *)*)


(* ========================================================================= *)
(* PHASE 2: INVASION GRAPH FRAMEWORK *)
(* ========================================================================= *)

(* findAdmissibleCommunities: Identifies admissible communities from siphon decomposition
   RHS: right-hand side of ODE system
   var: list of variables (symbols)
   mSi: list of minimal siphons (as variable names/strings)

   Returns: List of admissible communities, each represented as a list of siphon indices
*)
findAdmissibleCommunities[RHS_, var_, mSi_] := Module[{
  n, communities, singletons, pairs, triples, quadruples,
  allCommunities, admissible, testAdmissibility},

  n = Length[mSi];

  (* Helper function to test if a community is admissible *)
  (* A community is admissible if it corresponds to a valid boundary equilibrium *)
  testAdmissibility[communityIndices_] := Module[{
    activeSiphons, inactiveSiphons, activeVars, inactiveVars,
    constraints, feasible},

    (* Siphons in the community are active (non-zero) *)
    activeSiphons = mSi[[communityIndices]];
    (* Siphons not in the community are inactive (zero) *)
    inactiveSiphons = Delete[mSi, List /@ communityIndices];

    (* Variables corresponding to active and inactive siphons *)
    activeVars = Union[Flatten[ToExpression /@ activeSiphons]];
    inactiveVars = Union[Flatten[ToExpression /@ inactiveSiphons]];

    (* Check if there exists a positive steady state with inactive vars = 0 *)
    constraints = Join[
      Thread[inactiveVars -> 0],
      Thread[RHS /. Thread[inactiveVars -> 0] == 0]
    ];

    (* Try to find if such an equilibrium exists *)
    feasible = TimeConstrained[
      FindInstance[constraints, activeVars, Reals, 1],
      2,
      {}
    ];

    feasible =!= {}
  ];

  (* Generate all possible communities *)
  singletons = Table[{i}, {i, n}];

  (* Pairs *)
  pairs = If[n >= 2, Subsets[Range[n], {2}], {}];

  (* Triples *)
  triples = If[n >= 3, Subsets[Range[n], {3}], {}];

  (* Quadruples (for larger systems) *)
  quadruples = If[n >= 4, Subsets[Range[n], {4}], {}];

  allCommunities = Join[singletons, pairs, triples, quadruples];

  (* Filter to admissible communities *)
  admissible = Select[allCommunities, testAdmissibility];

  Print["Found ", Length[admissible], " admissible communities out of ", Length[allCommunities], " candidates"];

  admissible
];

(* computeInvasionRates: Computes invasion reproduction numbers for communities
   RHS: right-hand side of ODE system
   var: list of variables
   community: list of siphon indices representing resident community
   bdfpT: boundary fixed points from bdFp analysis

   Returns: Association mapping invading community -> invasion reproduction number
*)
computeInvasionRates[RHS_, var_, community_, bdfpT_] := Module[{
  residentSiphons, residentEquilibria, invasionRates, invader,
  residentEq, modAtResident, ngmResult, invadingVar, invasionR0},

  (* Get equilibria for resident community *)
  (* community is a list of siphon indices *)
  residentEquilibria = If[Length[community] > 0 && Length[bdfpT] >= Max[community],
    Flatten[bdfpT[[community]], 1],
    {}
  ];

  If[residentEquilibria === {} || residentEquilibria === {"froze"},
    Print["No valid resident equilibria for community ", community];
    Return[<||>];
  ];

  invasionRates = <||>;

  (* For each potential invader *)
  Do[
    If[!MemberQ[community, invader],
      (* This is a potential invading community *)
      Do[
        residentEq = eq;

        (* Evaluate system at resident equilibrium *)
        modAtResident = RHS /. residentEq;

        (* Compute NGM for invading strain *)
        ngmResult = Quiet[
          Check[NGM[{modAtResident, var}, var], $Failed],
          {Power::infy, Infinity::indet}
        ];

        If[ngmResult =!= $Failed,
          invadingVar = ngmResult[[7]]; (* K matrix *)
          If[MatrixQ[invadingVar],
            invasionR0 = Max[Eigenvalues[invadingVar]];
            invasionRates[invader] = invasionR0;
          ];
        ];
      , {eq, residentEquilibria}];
    ];
  , {invader, Length[bdfpT]}];

  invasionRates
];

(* rahmanInvasionGraph: Complete invasion graph construction following Rahman et al framework
   RHS: right-hand side of ODE system
   var: list of variables
   mSi: minimal siphons
   bdfpT: boundary fixed points

   Returns: Association with "Communities", "Edges", "InvasionNumbers", "Graph"
*)
rahmanInvasionGraph[RHS_, var_, mSi_, bdfpT_] := Module[{
  communities, edges, invasionNumbers, graph, vertices, edgeList,
  residentIdx, invaderIdx, invasionRate},

  Print["=== Rahman Invasion Graph Construction ==="];

  (* Step 1: Find admissible communities *)
  communities = findAdmissibleCommunities[RHS, var, mSi];

  If[communities === {},
    Print["No admissible communities found"];
    Return[<|"Communities" -> {}, "Edges" -> {}, "Graph" -> None|>];
  ];

  (* Step 2: Compute invasion rates between communities *)
  invasionNumbers = <||>;
  edges = {};

  Do[
    invasionRate = computeInvasionRates[RHS, var, communities[[residentIdx]], bdfpT];

    (* Create edges for successful invasions (R > 1) *)
    Do[
      If[KeyExistsQ[invasionRate, invaderIdx] && invasionRate[invaderIdx] > 1,
        AppendTo[edges, residentIdx -> invaderIdx];
        invasionNumbers[{residentIdx, invaderIdx}] = invasionRate[invaderIdx];
      ];
    , {invaderIdx, Length[communities]}];

  , {residentIdx, Length[communities]}];

  Print["Constructed invasion graph with ", Length[edges], " edges"];

  (* Step 3: Create graph visualization *)
  vertices = Range[Length[communities]];
  edgeList = DirectedEdge @@@ edges;

  graph = If[edgeList =!= {},
    Graph[vertices, edgeList,
      VertexLabels -> Table[
        i -> Placed[
          Row[{"C", i, ": ", Row[communities[[i]], ","]}],
          Center
        ],
        {i, Length[communities]}
      ],
      EdgeLabels -> Table[
        edge -> Placed[
          If[KeyExistsQ[invasionNumbers, List @@ edge],
            Style["R=" <> ToString[N[invasionNumbers[List @@ edge], 2]], Small],
            ""
          ],
          Above
        ],
        {edge, edgeList}
      ],
      VertexSize -> 0.3,
      VertexStyle -> LightGreen,
      EdgeStyle -> Directive[Black, Arrowheads[0.02]],
      GraphLayout -> "LayeredDigraphEmbedding",
      ImageSize -> Large
    ],
    Print["Empty invasion graph"];
    None
  ];

  <|
    "Communities" -> communities,
    "Edges" -> edges,
    "InvasionNumbers" -> invasionNumbers,
    "Graph" -> graph
  |>
];

(* ========================================================================= *)
(* END PHASE 2: INVASION GRAPH FRAMEWORK *)
(* ========================================================================= *)

(*=============================================================================
  MAS Module for EpidCRN: Minimal Autocatalytic Subnetworks

  Main functions:
  - isMAS[RN, T]  : Test if siphon T is self-replicable
  - MAS[RN]       : Find all MAS in network RN
  - testMAS[RN]   : Pretty-print MAS analysis
=============================================================================*)

(*-----------------------------------------------------------------------------
  Helper: Find reactions whose ALL reactants lie in siphon T
-----------------------------------------------------------------------------*)

(*-----------------------------------------------------------------------------
  isMAS[RN, T, eps]: Test if siphon T is self-replicable (an MAS)
-----------------------------------------------------------------------------*)
isMAS[RN_List, T_List, eps_ : 10^-6] := 
  Module[{RND, spe, gamma, speIdx, intReacIdx, gammaT, nReac, fluxVars, constr, sol},
   
   (* Get stoichiometry from EpidCRN *)
   RND = extMat[RN];
   spe = RND[[1]];
   gamma = RND[[4]];  (* rows=species, cols=reactions *)
   
   (* Species indices in T *)
   speIdx = Flatten[Position[spe, #] & /@ T];
   
   (* Internal reactions *)
   intReacIdx = findInternalReactions[RN, T];
   
   (* Early return if no internal reactions *)
   If[Length[intReacIdx] == 0,
    Return[<|"isMAS" -> False, 
      "reason" -> "no internal reactions",
      "siphon" -> T|>]
   ];
   
   (* Extract submatrix: species in T \[Times] internal reactions *)
   (* For single species, extract as list then wrap in list to make matrix *)
   If[Length[speIdx] == 1,
    gammaT = {gamma[[speIdx[[1]], intReacIdx]]};
    ,
    gammaT = gamma[[speIdx, intReacIdx]];
   ];
   
   nReac = Length[intReacIdx];
   
   (* Verify gammaT is non-empty *)
   If[gammaT === {} || gammaT === {{}},
    Return[<|"isMAS" -> False, 
      "reason" -> "empty submatrix",
      "siphon" -> T, "speIdx" -> speIdx, "intReacIdx" -> intReacIdx|>]
   ];
   
   (* Find v > eps with gammaT \[CenterDot] v > eps *)
   fluxVars = Table[c[i], {i, nReac}];
   constr = Flatten[{Thread[fluxVars >= eps], Thread[Flatten[gammaT . fluxVars] >= eps]}];
   
   sol = Quiet@FindInstance[constr, fluxVars, Reals, 1, WorkingPrecision -> 16];
   
   If[sol === {} || sol === {{}}, 
    <|"isMAS" -> False, 
      "reason" -> "no positive flux exists",
      "siphon" -> T,
      "intReac" -> intReacIdx,
      "gammaT" -> gammaT|>,
    <|"isMAS" -> True, 
      "siphon" -> T,
      "flux" -> (fluxVars /. First[sol]), 
      "intReac" -> intReacIdx,
      "netProd" -> Flatten[gammaT . (fluxVars /. First[sol])]|>]
  ];

(*-----------------------------------------------------------------------------
  MAS[RN]: Find all MAS in reaction network RN
-----------------------------------------------------------------------------*)
MAS[RN_List] := 
  Module[{RND, spe, gamma, minSiphQ, analysis, masQ},
   
   RND = extMat[RN];
   spe = RND[[1]];
   gamma = RND[[4]];
   
   (* Get minimal siphons *)
   minSiphQ = Quiet@minSiph[spe, asoRea[RN]];
   If[minSiphQ === $Failed || !ListQ[minSiphQ], minSiphQ = {}];
   
   (* Analyze each minimal siphon *)
   analysis = Table[
     Module[{T = minSiphQ[[k]], test},
      test = isMAS[RN, T];
      <|"siphon" -> T, 
        "isMAS" -> test["isMAS"],
        "drainable" -> isDrainable[RN, T],
        "critical" -> isCritical[RN, T],
        "details" -> test|>
     ],
     {k, Length[minSiphQ]}
   ];
   
   (* Extract only the MAS *)
   masQ = Select[minSiphQ, isMAS[RN, #]["isMAS"] &];
   
   <|"species" -> spe,
     "gamma" -> gamma,
     "minSiph" -> minSiphQ, 
     "MAS" -> masQ,
     "analysis" -> analysis|>
  ];

(*-----------------------------------------------------------------------------
  testMAS[RN]: Essential MAS analysis output
-----------------------------------------------------------------------------*)
testMAS[RN_List] := Module[{res},
  res = MAS[RN];
  
  Print["Species: ", res["species"]];
  Print["Minimal siphons: ", res["minSiph"]];
  Print["MAS: ", res["MAS"]];
  
  (* Analysis of each minimal siphon *)
  Do[
    Module[{a = res["analysis"][[i]], T},
      T = a["siphon"];
      Print[T, ": MAS=", a["isMAS"], ", Drain=", a["drainable"], 
        ", Crit=", a["critical"]];
      
      If[a["isMAS"],
        Print["  flux=", a["details"]["flux"], ", net=", a["details"]["netProd"]];
        ,
        Print["  ", a["details"]["reason"]];
      ];
    ],
    {i, Length[res["analysis"]]}
  ];
  
  res
];
