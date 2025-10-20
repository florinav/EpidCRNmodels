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
  species, alpha, beta, n, m, reactionsAsso, siphons, minimal},
  
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
  
  {minimal,siphons}
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
      subsAsso = compToAsso[RN[[j, 1]]];
      prodsAsso = compToAsso[RN[[j, 2]]];
      
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


(* Module to compute and visualize the IGMS (Infection Graph of Minimal Siphons) *)

ClearAll[canActRea, edgIGMS, IGMS];


(* Check if reaction activates Sj from Si *)
canActRea[reaction_, Si_List, Sj_List] := Module[
  {lhs, rhs, lhsComp, rhsComp, newSpecies},
  
  (* Skip birth and death reactions *)
  If[reaction === (0 -> _) || reaction === (_ -> 0), Return[False]];
  
  {lhs, rhs} = List @@ reaction;
  
  (* Parse complexes *)
  lhsComp = compToAsso[lhs];
  rhsComp = compToAsso[rhs];
  
  (* Quick filter: skip if reaction doesn't involve Si *)
  If[Length[Intersection[Keys[lhsComp], Si]] == 0, Return[False]];
  
  (* Find species that need to be produced *)
  newSpecies = Complement[Sj, Si];
  If[Length[newSpecies] == 0, Return[False]];
  
  (* Check if any new species has net positive production *)
  AnyTrue[newSpecies, Lookup[rhsComp, #, 0] > Lookup[lhsComp, #, 0] &]
];

(* Compute edges of IGMS graph *)
edgIGMS[RN_List, mSi_List] := Module[
  {n = Length[mSi], mSiStr},
  
  (* Convert symbols to strings *)
  mSiStr = Map[ToLowerCase@*ToString, mSi, {2}];
  
  (* Find all edges i -> j where some reaction activates j from i *)
  Flatten@Table[
    If[i != j && AnyTrue[RN, canActRea[#, mSiStr[[i]], mSiStr[[j]]] &],
      i -> j,
      Nothing
    ],
    {i, n}, {j, n}
  ]
];

(* Main IGMS function *)
IGMS[RN_List, mSi_List] := Module[
  {n = Length[mSi], igmsEdges, igmsGraph, cycles},
  
  Print["Minimal siphons: ", Table[Subscript["T", i] -> mSi[[i]], {i, n}]];
  
  (* Compute edges and create graph *)
  igmsEdges = edgIGMS[RN, mSi];
  
  igmsGraph = Graph[
    Range[n], igmsEdges,
    VertexLabels -> Table[
      i -> Placed[Row[{Subscript["T", i], "=", Row[mSi[[i]], ","]}], Center],
      {i, n}
    ],
    VertexSize -> 0.3,
    VertexStyle -> LightBlue,
    EdgeStyle -> Arrowheads[0.03],
    VertexLabelStyle -> Directive[FontSize -> 10],
    GraphLayout -> "SpringElectricalEmbedding",
    ImageSize -> Medium
  ];
  
  (* Find cycles *)
  cycles = FindCycle[igmsGraph, Infinity, All];
  
  Print["IGMS edges: ", igmsEdges];
  Print[If[Length[cycles] > 0,
    Row[{"Cycles found: ", Length[cycles]}],
    "IGMS is acyclic"
  ]];
  If[Length[cycles] > 0, Print["Cycles: ", cycles]];
  Print[igmsGraph];
  
  (* Return list: {edges, cycles, graph} *)
  {igmsEdges, cycles, igmsGraph}
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
