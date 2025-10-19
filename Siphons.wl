(* ::Package:: *)

minSiph[species_List, reactions_List, opt___] := 
  Module[{ns, sm, specs, constraints, solutions, siphons, minimal, reacAso, spe, returnAll},
   
   (* Check if "All" option is given *)
   returnAll = MemberQ[{opt}, "All"];
   
   (* Convert species to strings if needed *)
   spe = If[Length[species] > 0 && StringQ[species[[1]]], 
     species, 
     Map[If[StringQ[#], #, SymbolName[#]] &, species]
   ];
   ns = Length[spe];
   
   (* Convert reactions to association format if needed *)
   reacAso = If[Length[reactions] > 0 && Head[reactions[[1]]] === Rule,
     asoRea[reactions],
     reactions
   ];
   
   (* Create index mapping - using string keys *)
   sm = AssociationThread[spe -> Range[ns]];
   
   (* Create boolean variables - use non-conflicting prefix *)
   specs = Array[Symbol["sbo" <> ToString[#]] &, ns];
   constraints = {Or @@ specs};  (* At least one species in siphon *)
   
   (* Add constraints from reactions *)
   Do[
     Module[{subIdx, prodIdx, substrates, products},
       substrates = reacAso[[i]]["Substrates"];
       products = reacAso[[i]]["Products"];
       
       (* Get substrate indices - match lowercase *)
       subIdx = {};
       If[substrates =!= {} && substrates =!= {""}, 
         Do[
           Module[{subLower = ToLowerCase[sub]},
             If[KeyExistsQ[sm, subLower], AppendTo[subIdx, sm[subLower]]]
           ],
           {sub, substrates}
         ]
       ];
       
       (* Get product indices - match lowercase *)
       prodIdx = {};
       If[products =!= {} && products =!= {""}, 
         Do[
           Module[{prodLower = ToLowerCase[prod]},
             If[KeyExistsQ[sm, prodLower], AppendTo[prodIdx, sm[prodLower]]]
           ],
           {prod, products}
         ]
       ];
       
       (* Add siphon constraint for each product *)
       Do[
         Module[{constraint},
           If[Length[subIdx] == 0,
             (* Product from nothing - can't be in siphon *)
             constraint = Not[specs[[p]]],
             (* Product requires at least one substrate in siphon *)
             constraint = Implies[specs[[p]], 
               If[Length[subIdx] == 1, 
                 specs[[subIdx[[1]]]], 
                 Or @@ specs[[subIdx]]
               ]
             ]
           ];
           AppendTo[constraints, constraint]
         ],
         {p, prodIdx}
       ]
     ],
     {i, Length[reacAso]}
   ];
   
   (* Solve for siphons *)
   solutions = FindInstance[constraints, specs, Integers, 50];
   If[solutions === {}, Return[{}]];
   
   (* Extract siphon indices *)
   siphons = Map[Flatten@Position[specs /. #, True] &, solutions];
   siphons = DeleteDuplicates[siphons];
   siphons = Select[siphons, Length[#] > 0 &];
   
   (* Return all siphons or just minimal ones *)
   If[returnAll,
     (* Return all siphons as symbols *)
     Map[Symbol /@ spe[[#]] &, siphons],
     (* Find minimal siphons *)
     minimal = {};
     Do[
       If[Not[AnyTrue[siphons, 
           Function[other, other =!= siphon && SubsetQ[siphon, other]]]], 
         AppendTo[minimal, siphon]
       ],
       {siphon, siphons}
     ];
     (* Return minimal siphons as symbols - use Symbol[] to avoid subscript interpretation *)
     Map[Symbol /@ spe[[#]] &, minimal]
   ]
];

(* Test 
Print["=== Test with species 's' ==="];
RN = {"s" + "i1" -> 2*"i1", "s" + "i2" -> 2*"i2", 
      0 -> "s", "s" -> 0, "i1" -> 0, "i2" -> 0};
spe = extSpe[RN];
Print["Species: ", spe];
reac = asoRea[RN];
siphons = minSiph[spe, reac];
Print["Minimal siphons: ", siphons];*)


(* findCores: Find autocatalytic cores using EpidCRN functions *)

ClearAll[findInternalReactions, isAutocatalyticBlockLP, findCores];

findCores::usage = "findCores[RN, k] or findCores[RN, \"MaxSize\" -> k]
  findCores[RN, \"CandidateSets\" -> {list of candidate sets}]
  findCores[RN, \"MinimalOnly\" -> False]  (* return all cores, not just minimal ones *)
RETURNS: Association with keys:
  \"all\" - all cores found (isCore=True only)
  \"cores\" - list of minimal cores (or all if MinimalOnly->False)
  Each core contains:
    \"T\" - the species set
    \"intReacIdx\" - internal reaction indices
    \"details\" - LP results including:
      \"flux\" - reaction flux vector v
      \"growth\" - growth vector A\[CenterDot]v (should be > 0 for true cores)
      \"t\" - uniform lower bound on growth rates";

(* Find internal reactions for candidate T *)
(* A reaction is internal if: reactants \[SubsetEqual] T AND products \:2229 T \[NotEqual] \[EmptySet] *)
(* This excludes pure outflow reactions like a\[RightArrow]0 *)
findInternalReactions[RN_List, T_List] := Module[
  {Tstr = ToLowerCase /@ (ToString /@ T)},
  
  Select[Range[Length[RN]], 
    Module[{reactants = compToAsso[RN[[#, 1]]], products = compToAsso[RN[[#, 2]]]},
      (* Reactants from T or empty (inflow) *)
      (Keys[reactants] === {} || SubsetQ[Tstr, Keys[reactants]]) &&
      (* At least one product in T (excludes outflows to 0) *)
      Length[Intersection[Keys[products], Tstr]] > 0
    ] &
  ]
];

(* LP test for autocatalytic core *)
isAutocatalyticBlockLP[gamma_, speciesIdx_List, intReacIdx_List] := 
  Module[{gammaT, gammaExt, m, n, c, bounds, sol, tsol, vsol, growth, eps = 10^-8, vmin = 0, vmax = 1},
  
  If[Length[intReacIdx] == 0 || Length[speciesIdx] == 0, 
    Return[<|"isCore" -> False, "reason" -> "empty"|>]];
  
  (* Extract submatrix for species in T and internal reactions *)
  gammaT = If[Length[speciesIdx] == 1, 
    {gamma[[speciesIdx[[1]], intReacIdx]]}, 
    gamma[[speciesIdx, intReacIdx]]];
  
  m = Length[speciesIdx]; 
  n = Length[intReacIdx];
  
  (* Extend gamma by appending column of -1's: [\[CapitalGamma]_T | -1] *)
  gammaExt = Map[Append[#, -1] &, gammaT];
  
  (* Objective: minimize -t (i.e., maximize t) *)
  c = Append[ConstantArray[0, n], -1];
  
  (* Bounds: v_i >= 0, at least one v_i > 0 enforced by requiring sum(v) >= eps *)
  bounds = Join[ConstantArray[{vmin, vmax}, n], {{0, Infinity}}];
  
  (* Add constraint: sum of fluxes >= eps to ensure v \[NotEqual] 0 *)
  gammaExtWithSum = Join[gammaExt, {Append[ConstantArray[1, n], 0]}];
  rhsWithSum = Join[ConstantArray[0, m], {eps}];
  
  (* Solve LP: maximize t subject to [\[CapitalGamma]_T | -1]\[CenterDot][v; t] >= 0 and sum(v) >= eps *)
  (* Suppress warnings about infeasible LPs *)
  sol = Quiet[
    LinearProgramming[c, gammaExtWithSum, rhsWithSum, bounds],
    LinearProgramming::lpsnf
  ];
  
  If[sol === $Failed || !ListQ[sol] || MemberQ[sol, Indeterminate],
    Return[<|"isCore" -> False, "reason" -> "LP failed"|>]];
  
  tsol = sol[[n + 1]];
  vsol = Take[sol, n];
  
  (* Compute actual growth: \[CapitalGamma]_T \[CenterDot] v *)
  growth = gammaT . vsol;
  
  (* Check: t > eps means \[CapitalGamma]_T \[CenterDot] v >= t\[CenterDot]1 > 0 (strictly positive growth) *)
  If[tsol > eps && AllTrue[growth, # > eps &],
    <|"t" -> tsol, "flux" -> AssociationThread[intReacIdx -> vsol]|>,
    <|"isCore" -> False, "reason" -> "growth not strictly positive", 
      "t" -> tsol, "v" -> vsol, "growth" -> growth|>
  ]
];

(* Main function *)
Options[findCores] = {"CandidateSets" -> Automatic, "MaxSize" -> 6, "MinimalOnly" -> True};

findCores[RN_List, maxSize_Integer] := findCores[RN, "MaxSize" -> maxSize];

findCores[RN_List, opts : OptionsPattern[]] := Module[
  {RND, spe, gamma, candidateSets, maxSize, minimalOnly, speciesIdxMap, results, cores, minimal},
  
  RND = extMat[RN];
  spe = RND[[1]];
  gamma = RND[[4]];
  
  speciesIdxMap = AssociationThread[spe -> Range[Length[spe]]];
  
  candidateSets = OptionValue["CandidateSets"];
  maxSize = OptionValue["MaxSize"];
  minimalOnly = OptionValue["MinimalOnly"];
  
  If[candidateSets === Automatic,
    candidateSets = Rest@Subsets[spe, {1, Min[maxSize, Length[spe]]}],
    candidateSets = Map[ToLowerCase /@ (ToString /@ #) &, candidateSets]
  ];
  
  results = Table[
    Module[{T = candidateSets[[k]], intReacIdx, speciesIdx, lpres},
      intReacIdx = findInternalReactions[RN, T];
      
      If[Length[intReacIdx] == 0, 
        <|"T" -> T, "isCore" -> False, "reason" -> "no internal reactions"|>,
        
        speciesIdx = Lookup[speciesIdxMap, T, Nothing];
        
        If[Length[speciesIdx] != Length[T],
          <|"T" -> T, "isCore" -> False, "reason" -> "species not found"|>,
          
          lpres = isAutocatalyticBlockLP[gamma, speciesIdx, intReacIdx];
          If[KeyExistsQ[lpres, "isCore"] && lpres["isCore"] === False,
            <|"T" -> T, "isCore" -> False, "reason" -> lpres["reason"]|>,
            <|"T" -> T, "intReacIdx" -> intReacIdx, "details" -> lpres|>
          ]
        ]
      ]
    ],
    {k, Length[candidateSets]}
  ];
  
  cores = Select[results, !KeyExistsQ[#, "isCore"] || #["isCore"] === True &];
  
  minimal = If[minimalOnly,
    Select[cores, Function[c,
      Not[MemberQ[DeleteCases[cores, c], 
        d_ /; SubsetQ[c["T"], d["T"]] && d["T"] =!= c["T"]]]
    ]],
    cores
  ];
  
  <|"species" -> spe, "tested" -> Length[candidateSets], 
    "all" -> cores, "cores" -> minimal|>
];

(* TEST 
RN1 = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> 0};
res1 = findCores[RN1, "CandidateSets" -> {{"s", "i1"}}];
Print["Simple model - Cores: ", res1["cores"]];*)

(* Example with multiple independent cores *)
(* Two independent infection chains: S+I1->2I1 and S+I2->2I2 
RN3 = {0 -> "S", "S" -> 0, 
       "S" + "I1" -> 2*"I1", "I1" -> 0,
       "S" + "I2" -> 2*"I2", "I2" -> 0};
res3 = findCores[RN3, 3];
Print["\nMultiple cores example:"];
Print["  All cores found: ", res3["all"]];
Print["  Minimal cores: ", res3["cores"]];*)

(* SDAS 
RN2 = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> 0, 
   "I1" + "I2" -> 2*"I2", "I2" -> 0, "I2" + "I3" -> 2*"I3", "I3" -> 0};
res2 = findCores[RN2];
Print["\nFull epidemic model - minimal cores: ", res2["cores"]];*)


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
