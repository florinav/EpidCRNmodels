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


(* Module to compute and visualize the IGMS (Infection Graph of Minimal Siphons) *)

(* Helper: Check if one reaction activates Sj from Si *)
(* Si and Sj must be lists of strings *)
canActRea[reaction_, Si_, Sj_] := Module[
  {lhs, rhs, lhsComp, rhsComp, lhsSpecies, newSpecies, netProduction},
  
  (* Skip birth and death reactions *)
  If[reaction === (0 -> _) || reaction === (_ -> 0), Return[False]];
  
  lhs = reaction[[1]];
  rhs = reaction[[2]];
  
  (* Parse complexes - now returns lowercase from compToAsso *)
  lhsComp = compToAsso[lhs];
  rhsComp = compToAsso[rhs];
  
  (* Extract unique species as strings *)
  lhsSpecies = Keys[lhsComp];
  
  (* Quick filter: skip if reaction doesn't involve infection species *)
  If[Length[Intersection[lhsSpecies, Si]] == 0, Return[False]];
  
  (* Find species that need to be produced *)
  newSpecies = Complement[Sj, Si];
  If[Length[newSpecies] == 0, Return[False]];
  
  (* Check net production with stoichiometry *)
  netProduction = Table[
    Lookup[rhsComp, sp, 0] - Lookup[lhsComp, sp, 0],
    {sp, newSpecies}
  ];
  
  Max[netProduction] > 0
];

(* Compute edges of IGMS graph *)
edgIGMS[RN_, mSi_] := Module[{n, igmsEdges, mSiStr},
  n = Length[mSi];
  igmsEdges = {};
  
  (* Convert mSi from symbols to strings, using SymbolName *)
  mSiStr = Map[SymbolName /@ # &, mSi];
  
  Do[
    Do[
      If[i != j,
        Do[
          If[canActRea[reaction, mSiStr[[i]], mSiStr[[j]]],
            AppendTo[igmsEdges, i -> j];
            Break[]
          ],
          {reaction, RN}
        ]
      ],
      {j, n}
    ],
    {i, n}
  ];
  
  igmsEdges
];

(* Main IGMS function *)
IGMS[RN_, mSi_] := Module[
  {species, infSpecies, igmsEdges, igmsGraph, cycles, n},
  
  n = Length[mSi];
  
  (* Extract all species from reaction network *)
  species = DeleteDuplicates[Flatten[Cases[RN, _String, Infinity]]];
  
  (* Flatten minimal siphons to get all infection species *)
  infSpecies = DeleteDuplicates[Flatten[mSi]];
  
  Print["Minimal siphons: ", 
    Table[Row[{Subscript["T", i], "=", mSi[[i]]}], {i, n}]];
  
  (* Compute IGMS edges *)
  igmsEdges = edgIGMS[RN, mSi];
  
  (* Create graph *)
  igmsGraph = Graph[
    Range[n],
    igmsEdges,
    VertexLabels -> Table[
      i -> Placed[
        Row[{Subscript["T", i], "=", Row[mSi[[i]], ","]}],
        Center
      ],
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
  
  (* Print edges and cycles *)
  Print["IGMS edges: ", igmsEdges];
  If[Length[cycles] > 0,
    Print["Cycles in IGMS (", Length[cycles], " found): "];
    Do[Print["  Cycle ", i, ": ", cycles[[i]]], {i, Length[cycles]}],
    Print["IGMS is acyclic (no cycles found)"]
  ];
  Print[igmsGraph];
  
  (* Return list: edges, cycles, graph *)
  {igmsEdges, cycles, igmsGraph}
];

$DebugIGMS = False;


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
(* Robust findInternalReactions: accepts Symbol or String species; T may be symbols or strings *)
findInternalReactions[RN_List, T_List] := Module[
  {Tstr, reactants, reactantTerms, reactantSpe},
  (* canonicalize T to strings *)
  Tstr = ToString /@ T;
  Table[
    reactants = RN[[j, 1]];
    (* parse reactant complex into a list of terms, then convert each term to a string species name *)
    reactantSpe = If[reactants === 0,
      {}, (* null complex *)
      Module[{terms = If[Head[reactants] === Plus, List @@ reactants, {reactants}]},
        DeleteDuplicates[
          ToString /@ (terms /. {
            (m_Integer)*s_ :> s, 
            (m_?NumericQ)*s_ :> s,
            s_ :> s
          })
        ]
      ]
    ];
    (* choose reaction j if all reactants (as strings) are inside T (as strings) *)
    If[SubsetQ[reactantSpe, Tstr], j, Nothing],
    {j, Length[RN]}
  ]
];

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


isDrainable[RN_, speciesSet_List] := 
  Module[{RND, spe, gamma, speciesIndices, coeffs, constraints, solution}, 
   RND = extMat[RN];
   spe = RND[[1]];
   gamma = RND[[4]];
   speciesIndices = Flatten[Position[spe, #] & /@ speciesSet];
   
   coeffs = Array[c, Dimensions[gamma][[2]]];
   linearComb = coeffs . Transpose[gamma];
   constraints = Join[Thread[coeffs >= 0], 
     Thread[linearComb[[speciesIndices]] < 0], {Total[coeffs] > 0}];
   solution = FindInstance[constraints, coeffs, Reals];
   Length[solution] > 0];

isSelfReplicable[RN_, speciesSet_List] := 
  Module[{RND, spe, gamma, speciesIndices, coeffs, constraints, solution}, 
   RND = extMat[RN];
   spe = RND[[1]];
   gamma = RND[[4]];
   speciesIndices = Flatten[Position[spe, #] & /@ speciesSet];
 
   
   coeffs = Array[c, Dimensions[gamma][[2]]];
   linearComb = coeffs . Transpose[gamma];
   constraints = Join[Thread[coeffs >= 0], 
     Thread[linearComb[[speciesIndices]] > 0], {Total[coeffs] > 0}];
   solution = FindInstance[constraints, coeffs, Reals];
   Length[solution] > 0];

isCritical[RN_, speciesSet_List] := 
  Module[{RND, spe, conservationLaws, speciesIndices}, 
   RND = extMat[RN];
   spe = RND[[1]];
   speciesIndices = Flatten[Position[spe, #] & /@ speciesSet];
   conservationLaws = cons[RND[[4]], {}];
   
   If[Length[conservationLaws] == 0, True, 
    !AnyTrue[conservationLaws, 
     Function[law, 
      AllTrue[speciesIndices, law[[#]] >= 0 &] && 
       AnyTrue[speciesIndices, law[[#]] > 0 &] && 
       AllTrue[Complement[Range[Length[spe]], speciesIndices], law[[#]] == 0 &]]]]];

isSiph[species_List, reactions_List, siphon_List] := 
  Module[{ns, sm, siphonSet, isSiphonQ, subIdx, prodIdx, substrates, products}, 
   ns = Length[species];
   sm = AssociationThread[species -> Range[ns]];
   siphonSet = siphon;
   isSiphonQ = True;
   
   Do[substrates = reaction["Substrates"];
    products = reaction["Products"];
    subIdx = If[substrates === {} || substrates === {""}, {}, 
      Select[Lookup[sm, substrates, Nothing], IntegerQ]];
    prodIdx = If[products === {} || products === {""}, {}, 
      Select[Lookup[sm, products, Nothing], IntegerQ]];
    
    Do[If[MemberQ[siphonSet, species[[p]]], 
      If[Length[subIdx] == 0, isSiphonQ = False; Break[], 
       If[!AnyTrue[subIdx, MemberQ[siphonSet, species[[#]]] &], 
        isSiphonQ = False; Break[]]]];, {p, prodIdx}];
    If[!isSiphonQ, Break[]];, {reaction, reactions}];
   isSiphonQ];

siphonAnalysis[reactions_, OptionsPattern[{"Verbose" -> False}]] := 
  Module[{species, allSiphons, classification, verbose}, 
   verbose = OptionValue["Verbose"];
   species = extSpe[reactions];
   
   If[verbose, Print["Analyzing network with ", Length[species], " species: ", species]];
   allSiphons = minSiph[species, asoRea[reactions]];
   If[verbose, Print["Found ", Length[allSiphons], " siphons: ", allSiphons]];
   
   classification = Association[];
   Do[Module[{siphon, drainableResult, selfReplicableResult, criticalResult, siphonType, significance}, 
     siphon = allSiphons[[i]];
     drainableResult = isDrainable[reactions, siphon];
     selfReplicableResult = isSelfReplicable[reactions, siphon];
     criticalResult = isCritical[reactions, siphon];
     
     siphonType = Which[
       drainableResult && selfReplicableResult, "Competing (both drainable and self-replicable)", 
       drainableResult && !selfReplicableResult, "Extinction risk (drainable only)", 
       !drainableResult && selfReplicableResult, "Autocatalytic growth (self-replicable only)", 
       True, "Neutral (neither drainable nor self-replicable)"];
     
     significance = Which[
       !criticalResult, "Non-critical (protected by conservation)", 
       drainableResult && !selfReplicableResult, "High extinction risk", 
       drainableResult && selfReplicableResult, "Competition between growth and extinction", 
       !drainableResult && selfReplicableResult, "Pure autocatalytic potential", 
       True, "No direct persistence threat"];
     
     classification[siphon] = <|"IsDrainable" -> drainableResult, 
       "IsSelfReplicable" -> selfReplicableResult, 
       "IsCritical" -> criticalResult, 
       "Type" -> siphonType, 
       "Significance" -> significance, 
       "IsMinimal" -> True|>;
     
     If[verbose, 
      Print["Siphon ", siphon, ": ", siphonType];
      Print["  Drainable: ", drainableResult, ", Self-replicable: ", selfReplicableResult, ", Critical: ", criticalResult];
      Print["  Significance: ", significance];];], {i, Length[allSiphons]}];
   classification];
   
   isCatalytic[RN_] := 
  Module[{catalyticSets}, 
   catalyticSets = findCatalyticSets[RN];
   Length[catalyticSets] > 0];

findCatalyticSets[RN_] := 
  Module[{RND, spe, allSiphons, selfReplicableCriticalSiphons, reactionCatalysts, allCatalyticSets}, 
   RND = extMat[RN];
   spe = RND[[1]];
   allSiphons = Quiet[minSiph[spe, asoRea[RN]]];
   
   If[allSiphons === $Failed || !ListQ[allSiphons], allSiphons = {}];
   selfReplicableCriticalSiphons = Select[allSiphons, isSelfReplicable[RN, #] && isCritical[RN, #] &];
   
   reactionCatalysts = {};
   Do[Module[{reactants, products, catalysts}, 
     reactants = compToAsso[RN[[i, 1]]];
     products = compToAsso[RN[[i, 2]]];
     catalysts = Select[Keys[reactants], KeyExistsQ[products, #] && reactants[#] == products[#] &];
     If[Length[catalysts] > 0, reactionCatalysts = Union[reactionCatalysts, {catalysts}]];], {i, Length[RN]}];
   
   Union[selfReplicableCriticalSiphons, reactionCatalysts]];

findMinimalCriticalSiphons[RN_] := 
  Module[{RND, spe, allSiphons, criticalSiphons}, 
   RND = extMat[RN];
   spe = RND[[1]];
   allSiphons = minSiph[spe, asoRea[RN]];
   criticalSiphons = Select[allSiphons, isCritical[RN, #] &];
   criticalSiphons];




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
