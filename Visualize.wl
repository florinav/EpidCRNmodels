(* ::Package:: *)

phase2[RHS_, var_, cN_:{}, tMax_:50, nTraj_:15] := Module[{
  dyn, cFP, fp, jac, eigVals, stabilities, Xp, Xs, saddles, 
  xM, r1, r2, plotFP, plotStream, Gp, cP1, cP2},
  
  (* Apply parameter substitutions *)
  dyn = If[Length[cN] > 0, RHS /. cN, RHS];
  Print["dyn=", dyn];
  
  (* Find fixed points as conditions (rules) *)
  cFP = Quiet[NSolve[And @@ Thread[dyn == 0], var, Reals]];
  Print["cFP (conditions)=", cFP];
  
  (* Extract coordinate values from conditions *)
  fp = var /. cFP;
  Print["fp (all fixed points)=", fp];
  
  (* Check if any fixed points found *)
  If[Length[fp] == 0,
    Print["Warning: No fixed points found"];
    Return[{{}, {}, {}, {}, Graphics[], Graphics[]}]
  ];
  
  (* Jacobian *)
  jac = Outer[D, dyn, var];
  
  (* Analyze each fixed point *)
  eigVals = Table[
    Eigenvalues[jac /. Thread[var -> fp[[i]]]],
    {i, Length[fp]}
  ];
  
  stabilities = Table[
    Which[
      AllTrue[Re[eigVals[[i]]], # < 0 &], "Stable",
      AllTrue[Re[eigVals[[i]]], # > 0 &], "Unstable",
      True, "Saddle"
    ],
    {i, Length[fp]}
  ];
  
  Print["stabilities=", stabilities];
  
  (* Filter positive fixed points *)
  Xp = Select[fp, AllTrue[#, NonNegative] &];
  Print["Xp (positive fixed points)=", Xp];
  
  (* Identify saddle points among positive fps *)
  saddles = If[Length[Xp] > 0,
    Module[{indices},
      indices = Flatten[Position[fp, #] & /@ Xp];
      Select[
        MapThread[{#1, #2, #3} &, 
          {Xp, eigVals[[indices]], stabilities[[indices]]}],
        #[[3]] == "Saddle" &
      ]
    ],
    {}
  ];
  
  Print["saddles=", saddles];
  
  (* Sort fixed points *)
  Xs = If[Length[Xp] > 0, SortBy[Xp, Identity], {}];
  xM = If[Length[Xp] > 0, Max /@ Transpose[Xp], {1, 1}];
  
  (* Plotting ranges *)
  r1 = {var[[1]], -0.05, xM[[1]] + 0.3};
  r2 = {var[[2]], -0.05, xM[[2]] + 0.3};
  
  Print["Plot ranges: r1=", r1, ", r2=", r2];
  
  (* Fixed points with color coding *)
  Gp = If[Length[Xp] > 0,
    Module[{indices},
      indices = Flatten[Position[fp, #] & /@ Xp];
      Print["Plotting fixed points at: ", Xp];
      Graphics[{
        PointSize[0.05],
        EdgeForm[Directive[Black, Thick]],
        MapThread[
          {Switch[#2, "Stable", Green, "Unstable", Red, "Saddle", Orange, _, Black],
           Disk[#1, 0.03]} &,
          {Xp, stabilities[[indices]]}
        ]
      }]
    ],
    Graphics[]
  ];
  
  (* First nullcline - dx1/dt = 0 in Blue *)
  cP1 = ContourPlot[
    dyn[[1]],
    r1, r2,
    Contours -> {0},
    ContourStyle -> Directive[Blue, Thick],
    ContourShading -> None,
    Frame -> True,
    FrameLabel -> {ToString[var[[1]]], ToString[var[[2]]]}
  ];
  
  (* Second nullcline - dx2/dt = 0 in Red *)
  cP2 = ContourPlot[
    dyn[[2]],
    r1, r2,
    Contours -> {0},
    ContourStyle -> Directive[Red, Thick],
    ContourShading -> None,
    Frame -> True,
    FrameLabel -> {ToString[var[[1]]], ToString[var[[2]]]}
  ];
  
  (* Plot 1: Fixed points + both nullclines *)
  plotFP = Show[cP1, cP2, Gp, 
    PlotLabel -> Style["Nullclines and Fixed Points", Medium]];
  
  (* Plot 2: StreamPlot *)
  plotStream = StreamPlot[
    {dyn[[1]], dyn[[2]]},
    r1, r2,
    StreamStyle -> Arrowheads[0.02],
    ColorFunction -> "Rainbow",
    StreamPoints -> Fine,
    Frame -> True,
    FrameLabel -> {ToString[var[[1]]], ToString[var[[2]]]},
    PlotLabel -> Style["Stream Plot", Medium]
  ];
  
  (* Return results *)
  {Xs, eigVals, stabilities, saddles, plotFP, plotStream}
];


(* ========================================================================== *)
(* NEWTON POLYTOPE CONSTRUCTION *)
(* ========================================================================== *)

NewtPol[complexCoords_Association, opts___] := 
  Module[{points, convexHull, polytope}, 
   points = Values[complexCoords];
   If[Length[points] >= 3 && Length[Dimensions[points]] == 2, 
    convexHull = ConvexHullMesh[points];
    polytope = {LightGray, EdgeForm[{Black, Thin}], convexHull}, 
    polytope = {}];
   polytope];

(* Test NewtPol *)
(*
coords = <|"X" -> {0, 0}, "Y" -> {1, 0}, "Z" -> {0, 1}|>;
NewtPol[coords]  (* Should return a convex hull mesh graphics *)
*)


(* ========================================================================== *)
(* FEINBERG-HORN-JACKSON GRAPH VISUALIZATION *)
(* ========================================================================== *)

FHJ[comp_List, edges_List, rates_List, ver_ : {}, groups_List : {}] :=
   Module[{colorList, shapeList, vertexColors, options, vertexShapes, 
    defaultColor = Yellow}, 
   colorList = {Green, Red, Yellow, Purple, Orange};
   shapeList = {"Square", "Circle", "ConcaveHexagon", "Triangle", 
     "Hexagon", "Pentagon", "Star"};
   vertexColors = 
    Join[
     Flatten[MapIndexed[Thread[#1 -> colorList[[#2[[1]]]]] &, 
       groups]], # -> defaultColor & /@ 
      Complement[comp, Flatten[groups]]];
   vertexShapes = 
    Flatten[MapIndexed[Thread[#1 -> shapeList[[#2[[1]]]]] &, 
      groups]];
   options = {VertexShapeFunction -> vertexShapes, 
     VertexStyle -> vertexColors, VertexSize -> ver, 
     VertexLabels -> {_ -> Placed[Automatic, Center]}, 
     EdgeStyle -> {{Black, Thick}}, PerformanceGoal -> "Quality", 
     EdgeLabels -> Thread[edges -> rates], 
     EdgeLabelStyle -> Directive[Black, Bold, Background -> White]};
   LayeredGraphPlot[edges, Right, options]];

(* Test FHJ *)
(*
FHJ[{"X", "Y", "Z"}, {{"X", "Y"}, {"Y", "Z"}}, {1, 2}]  
(* Should return a layered graph plot *)
*)


(* ========================================================================== *)
(* EUCLIDEAN FEINBERG-HORN-JACKSON GRAPH *)
(* ========================================================================== *)

EucFHJ[reactions_] := EucFHJ[reactions, {}]

EucFHJ[reactions_, opts___] := 
  Module[{matResults, species, complexes, reactantComplexes, 
    productOnlyComplexes, reactionPairs, coords, vertices, edges, 
    labels, reactionPolygon, result},
   matResults = extMat[reactions];
   species = matResults[[1]];
   If[Length[species] > 2, species = Take[species, 2]];
   If[Length[species] < 2, 
    species = 
     Join[species, 
      Table["dummy" <> ToString[i], {i, 2 - Length[species]}]]];
   
   complexes = {};
   reactantComplexes = {};
   reactionPairs = {};
   Do[Module[{reactant, product},
     If[Head[reactions[[i]]] === Rule, 
      reactant = reactions[[i]] /. Rule -> List // First;
      product = reactions[[i]] /. Rule -> List // Last,
      reactant = reactions[[i, 1]];
      product = reactions[[i, 2]]];
     If[! MemberQ[complexes, reactant], AppendTo[complexes, reactant]];
     If[! MemberQ[complexes, product], AppendTo[complexes, product]];
     If[! MemberQ[reactantComplexes, reactant], 
      AppendTo[reactantComplexes, reactant]];
     AppendTo[reactionPairs, {reactant, product}];], {i, 
     Length[reactions]}];
   productOnlyComplexes = Complement[complexes, reactantComplexes];
   
   coords = <||>;
   Do[Module[{prodAssoc, coord}, prodAssoc = comp2Asso[complex];
     coord = 
      Table[If[KeyExistsQ[prodAssoc, species[[j]]], 
        prodAssoc[species[[j]]], 0], {j, 2}];
     coords[complex] = coord;], {complex, complexes}];
   
   reactionPolygon = 
    Module[{reactantPoints, convexHull}, 
     reactantPoints = Values[KeyTake[coords, reactantComplexes]];
     If[Length[reactantPoints] >= 3, 
      convexHull = ConvexHullMesh[reactantPoints];
      {Lighter[Blue, 0.8], EdgeForm[{Blue, Thick}], 
       convexHull}, {}]];
   
   vertices = 
    Module[{reactantVertices, productVertices}, 
     reactantVertices = 
      Table[{Red, EdgeForm[{Red, Thick}], 
        Disk[coords[complex], 0.03]}, {complex, reactantComplexes}];
     productVertices = 
      Table[{PointSize[0.015], Red, Point[coords[complex]]}, {complex,
         productOnlyComplexes}];
     Join[reactantVertices, productVertices]];
   
   labels = 
    Table[Text[Style[symbToStr[complex], 12, Black], 
      coords[complex], {-1.5, -1.5}], {complex, complexes}];
   edges = Table[Module[{start, end}, start = coords[pair[[1]]];
      end = coords[pair[[2]]];
      If[
       start == end, {Blue, Red, Arrowheads[Large], 
        Arrow[{start, start + {0.2, 0.2}}]}, {Blue, Red, 
        Arrowheads[Large], Arrow[{start, end}]}]], {pair, 
      reactionPairs}];
   
   result = 
    Graphics[{reactionPolygon, edges, vertices, labels}, 
     PlotRange -> All, AspectRatio -> 1, Axes -> True, 
     AxesLabel -> {ToString[species[[1]]], ToString[species[[2]]]}, 
     PlotLabel -> "Euclidean Feinberg-Horn-Jackson Graph", 
     ImageSize -> 400, GridLines -> Automatic, 
     GridLinesStyle -> LightGray, opts];
   result];

(* Test EucFHJ *)
(*
EucFHJ[{"X" -> "Y", "Y" -> "X"}]  (* Should return a 2D reaction graph *)
*)


(* ========================================================================== *)
(* HELPER FUNCTIONS FOR ENDOTACTIC ANALYSIS *)
(* ========================================================================== *)

complexToVector[complex_, speciesList_] := 
  Module[{parsed, vector}, parsed = comp2Asso[complex];
   vector = 
    Table[If[KeyExistsQ[parsed, speciesList[[i]]], 
      parsed[speciesList[[i]]], 0], {i, Length[speciesList]}];
   vector];

(* Test complexToVector *)
(*
complexToVector[2*"X" + "Y", {"X", "Y"}]  (* Should return {2, 1} *)
*)

reactionVector[reaction_, speciesList_] := 
  Module[{reactantVec, productVec}, 
   reactantVec = complexToVector[reaction[[1]], speciesList];
   productVec = complexToVector[reaction[[2]], speciesList];
   productVec - reactantVec];

(* Test reactionVector *)
(*
reactionVector[{"X", "Y"}, {"X", "Y"}]  (* Should return {-1, 1} *)
*)

stoichiometricSubspace[reactions_, speciesList_] := 
  Module[{reactionVectors}, 
   reactionVectors = 
    Table[reactionVector[reactions[[i]], speciesList], {i, 
      Length[reactions]}];
   reactionVectors = 
    DeleteDuplicates[
     Select[reactionVectors, # != 
        Table[0, {Length[speciesList]}] &]];
   If[Length[reactionVectors] == 0, {}, reactionVectors]];

(* Test stoichiometricSubspace *)
(*
stoichiometricSubspace[{{"X", "Y"}}, {"X", "Y"}]  (* Should return {{-1, 1}} *)
*)

getMaximalElements[vectors_, w_] := 
  Module[{dotProducts, maxValue, validVectors, result}, 
   If[Length[vectors] == 0, Return[{}]];
   validVectors = 
    Select[vectors, VectorQ[#] && Length[#] == Length[w] &];
   If[Length[validVectors] == 0, Return[{}]];
   dotProducts = Map[w . # &, validVectors];
   If[Length[dotProducts] == 0 || ! AllTrue[dotProducts, NumericQ], 
    Return[{}]];
   maxValue = Max[dotProducts];
   result = validVectors[[Flatten[Position[dotProducts, maxValue]]]];
   result];

(* Test getMaximalElements *)
(*
getMaximalElements[{{1, 0}, {0, 1}}, {1, 1}]  (* Should return {{1, 0}, {0, 1}} *)
*)

getReactantVectors[reactions_, speciesList_] := 
  Module[{reactants}, 
   reactants = 
    Table[complexToVector[reactions[[i, 1]], speciesList], {i, 
      Length[reactions]}];
   DeleteDuplicates[reactants]];

(* Test getReactantVectors *)
(*
getReactantVectors[{{"X", "Y"}}, {"X", "Y"}]  (* Should return {{1, 0}} *)
*)

isWEssential[reaction_, w_, speciesList_] := 
  Module[{reacVec}, reacVec = reactionVector[reaction, speciesList];
   w . reacVec != 0];

(* Test isWEssential *)
(*
isWEssential[{"X", "Y"}, {1, 0}, {"X", "Y"}]  (* Should return True *)
*)

getWEssentialReactions[reactions_, w_, speciesList_] := 
  Select[reactions, isWEssential[#, w, speciesList] &];

(* Test getWEssentialReactions *)
(*
getWEssentialReactions[{{"X", "Y"}}, {1, 0}, {"X", "Y"}]  (* Should return {{"X", "Y"}} *)
*)

getWSupport[reactions_, w_, speciesList_] := 
  Module[{essentialReactions, essentialReactants}, 
   essentialReactions = 
    getWEssentialReactions[reactions, w, speciesList];
   If[Length[essentialReactions] == 0, Return[{}]];
   essentialReactants = 
    Table[complexToVector[essentialReactions[[i, 1]], 
      speciesList], {i, Length[essentialReactions]}];
   getMaximalElements[essentialReactants, w]];

(* Test getWSupport *)
(*
getWSupport[{{"X", "Y"}}, {1, 0}, {"X", "Y"}]  (* Should return {{1, 0}} *)
*)

isWEndotactic[reactions_, w_, speciesList_] := 
  Module[{essentialReactions, wSupport}, 
   essentialReactions = 
    getWEssentialReactions[reactions, w, speciesList];
   wSupport = getWSupport[reactions, w, speciesList];
   AllTrue[essentialReactions, 
    Module[{reactantVec, reacVec, inSupport}, 
      reactantVec = complexToVector[#[[1]], speciesList];
      reacVec = reactionVector[#, speciesList];
      inSupport = MemberQ[wSupport, reactantVec];
      If[inSupport, w . reacVec < 0, True]] &]];

(* Test isWEndotactic *)
(*
isWEndotactic[{{"X", "Y"}}, {1, 0}, {"X", "Y"}]  (* Should return True *)
*)

generateTestDirections[speciesList_] := 
  Module[{n, directions}, n = Length[speciesList];
   directions = {};
   Do[AppendTo[directions, UnitVector[n, i]];
    AppendTo[directions, -UnitVector[n, i]];, {i, n}];
   Do[AppendTo[directions, 
      Normalize[RandomReal[{-1, 1}, n]]];, {20}];
   If[n >= 2, 
    Do[Do[AppendTo[directions, 
         Normalize[UnitVector[n, i] + UnitVector[n, j]]];
        AppendTo[directions, 
         Normalize[UnitVector[n, i] - UnitVector[n, j]]];, {j, i + 1, 
         n}];, {i, 1, n - 1}];];
   DeleteDuplicates[directions, Norm[#1 - #2] < 10^(-10) &]];

(* Test generateTestDirections *)
(*
generateTestDirections[{"X", "Y"}]  (* Should return list of 2D directions *)
*)

isEndotactic[reactions_, speciesList_] := 
  Module[{testDirections}, 
   testDirections = generateTestDirections[speciesList];
   AllTrue[testDirections, 
    isWEndotactic[reactions, #, speciesList] &]];

(* Test isEndotactic *)
(*
isEndotactic[{{"X", "Y"}}, {"X", "Y"}]  (* Should return True *)
*)

isStronglyEndotactic[reactions_, speciesList_] := 
  Module[{stoichSubspace, testDirections}, 
   If[! isEndotactic[reactions, speciesList], Return[False]];
   stoichSubspace = stoichiometricSubspace[reactions, speciesList];
   testDirections = 
    Select[generateTestDirections[
      speciesList], ! AllTrue[stoichSubspace, # . # == 0 &] &];
   AllTrue[testDirections, 
    Module[{w, reactantVectors, maximalReactants, hasGoodReaction}, 
      w = #;
      reactantVectors = getReactantVectors[reactions, speciesList];
      maximalReactants = getMaximalElements[reactantVectors, w];
      hasGoodReaction = False;
      Do[
       Module[{maxReactant, candidateReactions}, 
        maxReactant = maximalReactants[[i]];
        candidateReactions = 
         Select[reactions, 
          complexToVector[#[[1]], speciesList] == maxReactant &];
        Do[
         Module[{reacVec}, 
          reacVec = 
           reactionVector[candidateReactions[[j]], speciesList];
          If[w . reacVec < 0, hasGoodReaction = True; Break[]];], {j, 
          Length[candidateReactions]}];
        If[hasGoodReaction, Break[]];], {i, Length[maximalReactants]}];
      hasGoodReaction] &]];

(* Test isStronglyEndotactic *)
(*
isStronglyEndotactic[{{"X", "Y"}}, {"X", "Y"}]  (* Should return True *)
*)


(* ========================================================================== *)
(* ENDOTACTIC ANALYSIS WITH VISUALIZATION *)
(* ========================================================================== *)

Options[endo] = {"ShowPlot" -> Automatic, "Verbose" -> False};

endo[reactions_, OptionsPattern[]] := 
  Module[{matResults, speciesList, isEndo, isStrongEndo, result, 
    showPlot, verbose, processedReactions, axisDirections, 
    axisEndotacticResults, failingSpecies, oneEndotacticReport, 
    testDirections, endotacticDirections, 
    endotacticDirectionResults},
   matResults = extMat[reactions];
   speciesList = matResults[[1]];
   processedReactions = 
    Table[{reactions[[i, 1]], reactions[[i, 2]]}, {i, 
      Length[reactions]}];
   showPlot = OptionValue["ShowPlot"];
   verbose = OptionValue["Verbose"];
   
   axisDirections = 
    Table[UnitVector[Length[speciesList], i], {i, 
      Length[speciesList]}];
   axisEndotacticResults = 
    Table[
     isWEndotactic[processedReactions, axisDirections[[i]], 
      speciesList], {i, Length[speciesList]}];
   
   failingSpecies = {};
   Do[If[! axisEndotacticResults[[i]], 
     AppendTo[failingSpecies, speciesList[[i]]]], {i, 
     Length[speciesList]}];
   
   oneEndotacticReport = 
    Table[speciesList[[i]] -> axisEndotacticResults[[i]], {i, 
      Length[speciesList]}];
   
   testDirections = generateTestDirections[speciesList];
   endotacticDirectionResults = 
    Table[{testDirections[[i]], 
      isWEndotactic[processedReactions, testDirections[[i]], 
       speciesList]}, {i, Length[testDirections]}];
   endotacticDirections = 
    Select[endotacticDirectionResults, #[[2]] == True &][[All, 1]];
   
   isEndo = isEndotactic[processedReactions, speciesList];
   isStrongEndo = 
    If[isEndo, isStronglyEndotactic[processedReactions, speciesList], 
     False];
   result = <|"Species" -> speciesList, 
     "Reactions" -> processedReactions, "Endotactic" -> isEndo, 
     "StronglyEndotactic" -> isStrongEndo, 
     "OneEndotacticReport" -> oneEndotacticReport, 
     "FailingSpeciesDirections" -> failingSpecies, 
     "EndotacticDirections" -> endotacticDirections|>;
   
   If[(showPlot === 
        True || (showPlot === Automatic && 
         Length[speciesList] == 2)) && Length[speciesList] == 2, 
    Print[EucFHJ[reactions]];, 
    If[showPlot === True && Length[speciesList] != 2, 
     Print["Note: Visualization only available for 2-species \
networks. This network has ", Length[speciesList], " species."];]];
   result];

(* Test endo *)
(*
endo[{"X" -> "Y", "Y" -> "X"}]  (* Should return analysis results *)
*)


(* ========================================================================== *)
(* INVASION GRAPH CONSTRUCTION *)
(* ========================================================================== *)

(* invGr - Invasion Graph Construction and Visualization

   Purpose: Constructs invasion graphs for multi-strain epidemic models showing
   which communities can invade which boundary equilibria.

   Input:
   - mSi: minimal siphons (list of lists of species names as strings)
   - bdFps: boundary fixed points (list of solution lists, one per siphon)
   - R0A: basic reproduction numbers (list)
   - RHS: right-hand side of ODE
   - var: variable list
   - E0: disease-free equilibrium

   Output: {gra, ver, edg, col}
   - gra: Graph object with invasion edges
   - ver: vertex list (communities as sorted integer lists)
   - edg: edge list (rules showing invasion pathways)
   - col: vertex colors (invasion rates: positive=Green, negative=Yellow)

   Algorithm:
   1. Build vertices from siphon complements (communities)
   2. Compute invasion numbers Rj(ySc) for each strain at each community
   3. Create edge S->S\Tj when Rj(ySc)>1
   4. Color nodes by invasion success
*)

invGr[mSi_List, bdFps_List, R0A_List, RHS_List, var_List, E0_List] :=
  Module[{
    n, nSi, ver, edg, invNum, col, gra,
    ii, jj, allSiphonIndices, siphonIndicesList,
    fromCom, toCom, invasionRate, canInvade, vertexColors
  },

  (* Number of strains *)
  n = Length[R0A];
  nSi = Length[mSi];

  (* Helper: convert siphon species names to indices *)
  (* mSi[[i]] is a list of species names like {"x1", "x2"} *)
  siphonIndicesList = Table[
    Flatten[Position[var, ToExpression[#]] & /@ mSi[[i]]],
    {i, nSi}
  ];

  (* All infection indices (union of all siphons) *)
  allSiphonIndices = Sort[DeleteDuplicates[Flatten[siphonIndicesList]]];

  (* Build vertices: all possible subcommunities *)
  ver = Subsets[allSiphonIndices];
  ver = Sort[ver];

  (* Build edges and compute invasion numbers *)
  edg = {};
  invNum = Table[{}, {Length[ver]}];

  Do[
    fromCom = ver[[ii]];

    (* For each strain jj *)
    Do[
      (* Strain jj can invade if its siphon is contained in fromCom *)
      If[SubsetQ[fromCom, siphonIndicesList[[jj]]],
        (* Compute invasion number R0j at DFE for now *)
        invasionRate = R0A[[jj]] /. E0;

        (* Numerical check *)
        canInvade = TrueQ[N[invasionRate] > 1];

        (* Target community: remove strain jj's siphon *)
        toCom = Sort[Complement[fromCom, siphonIndicesList[[jj]]]];

        (* Add edge if invasion successful and target is a vertex *)
        If[canInvade && MemberQ[ver, toCom],
          AppendTo[edg, fromCom -> toCom];
        ];

        (* Store invasion number *)
        AppendTo[invNum[[ii]], {jj, invasionRate, canInvade}];
      ];
    , {jj, nSi}];
  , {ii, Length[ver]}];

  (* Color vertices by invasion success *)
  col = Table[
    If[Length[invNum[[ii]]] > 0,
      If[AnyTrue[invNum[[ii]], #[[3]] &],
        RGBColor[0, 0.7, 0],  (* Green: successful invasion *)
        RGBColor[1, 1, 0]     (* Yellow: no successful invasions *)
      ],
      RGBColor[0.8, 0.8, 0.8]  (* LightGray: no invasion attempts *)
    ],
    {ii, Length[ver]}
  ];

  (* Build graph with subscripted labels *)
  vertexColors = Thread[ver -> col];
  vertexLabels = Table[
    ver[[i]] -> If[Length[ver[[i]]] == 0,
      "",
      Module[{varNames, subsVars, labelPos},
        varNames = var[[ver[[i]]]];
        (* Convert each variable to subscripted form and use Row for display *)
        subsVars = Row[subsCon[ToString[#]] & /@ varNames];
        (* Position logic: use Center for intermediate levels, Above/Below for ends *)
        (* Also offset horizontally for vertically aligned nodes *)
        labelPos = Which[
          Length[ver[[i]]] == 0, Center,
          Length[ver[[i]]] == Length[allSiphonIndices], Below,
          Length[ver[[i]]] <= Length[allSiphonIndices]/3, Above,
          Length[ver[[i]]] >= 2*Length[allSiphonIndices]/3, Below,
          True, If[Mod[i, 2] == 0, {Right, Above}, {Left, Above}]
        ];
        Placed[subsVars, labelPos]
      ]
    ],
    {i, Length[ver]}
  ];
  gra = If[Length[edg] > 0,
    Graph[ver, edg,
      VertexLabels -> vertexLabels,
      VertexStyle -> vertexColors,
      VertexSize -> 0.4,
      EdgeStyle -> Directive[RGBColor[0, 0, 0], Thickness[0.005]],
      GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Top, "RootVertex" -> {}},
      ImageSize -> {300, 300},
      ImagePadding -> {{40, 40}, {30, 30}}
    ],
    Graph[ver, {},
      VertexLabels -> vertexLabels,
      VertexStyle -> vertexColors,
      VertexSize -> 0.4,
      GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Top, "RootVertex" -> {}},
      ImageSize -> {300, 300},
      ImagePadding -> {{40, 40}, {30, 30}}
    ]
  ];

  (* Return *)
  {gra, ver, edg, col}
];

(* Test: 4 invasion graphs like Figure 2 in AF.pdf *)
(* Example call:
Module[{mSi, RHS, var, E0, res1, res2, res3, res4},
  mSi = {{"x1"}, {"x2"}, {"x3"}};
  RHS = {r*x1*(1 - x1 - alp*x2 - alp*x3),
         r*x2*(1 - x2 - alp*x1 - alp*x3),
         r*x3*(1 - x3 - alp*x1 - alp*x2)};
  var = {x1, x2, x3};
  E0 = {x1 -> 0, x2 -> 0, x3 -> 0};

  res1 = invGr[mSi, {}, {1.2, 0.8, 0.8}, RHS, var, E0];
  res2 = invGr[mSi, {}, {1.5, 1.3, 0.9}, RHS, var, E0];
  res3 = invGr[mSi, {}, {1.8, 1.6, 1.4}, RHS, var, E0];
  res4 = invGr[mSi, {}, {2.0, 1.9, 1.7}, RHS, var, E0];

  GraphicsRow[{res1[[1]], res2[[1]], res3[[1]], res4[[1]]}, ImageSize -> Full, Spacings -> 10]
]
*)

