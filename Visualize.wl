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
   Do[Module[{prodAssoc, coord}, prodAssoc = compToAsso[complex];
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
  Module[{parsed, vector}, parsed = compToAsso[complex];
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
