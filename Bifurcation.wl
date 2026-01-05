(* ::Package:: *)

fpHopf[RHS_, var_, par_, p0val_] := 
  Module[{eqns, allSols, posSols, eq, jac, eigs, complexEigs, upperEig, angle, coP}, 
   (* Create substitution rules *)
   coP = Thread[par -> p0val];
   
   (* Create equations and find positive solutions *)
   eqns = Thread[(RHS /. coP) == 0];
   allSols = Quiet[NSolveValues[eqns, var, Reals]];
   allSols = Union[allSols, SameTest -> (Norm[#1 - #2] < 10^-10 &)];
   allSols = Select[allSols, AllTrue[#, # >= 0 &] &];
   posSols = Select[allSols, AllTrue[#, # > 0 &] &];
   
   If[Length[posSols] == 0, Return[{{}, {}, -90, {}}];];
   
   (* Use first positive solution for eigenvalue analysis *)
   eq = posSols[[1]];
   jac = D[RHS, {var}] /. coP /. Thread[var -> eq];
   eigs = Chop[Eigenvalues[jac] // N];
   complexEigs = Select[eigs, Im[#] != 0 &];
   
   (* Calculate angle *)
   angle = If[Length[complexEigs] > 0, 
     upperEig = First[Select[complexEigs, Im[#] > 0 &]];
     ArcTan[Re[upperEig]/Im[upperEig]]*180/Pi, -90];
   
   {posSols, complexEigs, angle, eigs}];
   


(*Equilibrium scanning with automatic varInd detection and performance optimizations*)
scan[RHS_, var_, par_, persRule_, plotInd_, mSi_,
  gridRes_: Automatic, steadyTol_: 10^(-5),
  stabTol_: 10^(-8), chopTol_: 10^(-10), R01_: Automatic,
  R02_: Automatic, R12_: Automatic, R21_: Automatic, rangeExtension_: 2.5] := 
 Block[{bifP1min, bifP1max, bifP2min, bifP2max, bifP1vals, bifP2vals, 
   totalPoints, res, outcomes, outcomeCounts, finalPlot, plotData, 
   useGridMode, bifParIdx1, bifParIdx2, bifP1Center, bifP2Center, 
   progressVar, currentProgress, numpar, conPar, delta, wRan, hRan, 
   rCurves, fixedParams, fInd, activeEquations, activeColors, 
   activeLabels, par1, par2, range1, range2, 
   varS, varI1, varI2, varInd, infPositions, susceptiblePositions, persVal, 
   intersectionPoint, intersectionSol, numSol, rootSol, startPt,
   zeroTol, posTol, eq, EE, jac, eigs, classification, R01val, 
   R02val, R12val, R21val, coexistencePoints, rCurveLabels, 
   combinedLegend, eqLegend, rLegend, intersectionPlot, legendItems},
  
  (* Extract numeric parameter values from rules *)
  persVal = par /. persRule;
  
  (*AUTO-DETERMINE varInd from mSi*)
  (* Convert variable names in mSi to their positions in var *)
  (* mSi contains strings, so convert to symbols before searching *)
  infPositions = Union[Flatten[Table[
    FirstPosition[var, ToExpression[mSi[[i]][[j]]]][[1]],
    {i, 1, Length[mSi]}, {j, 1, Length[mSi[[i]]]}
    ]]]; (* All infected compartment positions *)
  susceptiblePositions = Complement[Range[Length[var]], infPositions];

  (* Construct varInd automatically using positions *)
  varInd = {
    susceptiblePositions[[1]], (* First non-infected compartment *)
    FirstPosition[var, ToExpression[mSi[[1]][[1]]]][[1]], (* First infected compartment of strain 1 *)
    FirstPosition[var, ToExpression[mSi[[2]][[1]]]][[1]]  (* First infected compartment of strain 2 *)
    };
  
  (* Extract variables using auto-determined indices *)
  varS = var[[varInd[[1]]]]; (* Susceptible variable *)
  varI1 = var[[varInd[[2]]]]; (* Strain 1 variable *)
  varI2 = var[[varInd[[3]]]]; (* Strain 2 variable *)

  (*Compute endemic equilibrium at persVal*)
  Module[{endemicEqs, endemicSols, endemicEE},
   endemicEqs = Join[Thread[RHS == 0], Thread[var > 0]];
   endemicSols = Quiet[FindInstance[endemicEqs /. Thread[par -> persVal], var, Reals, 1],
     {Power::infy, Infinity::indet, Power::indet}];
   If[Length[endemicSols] > 0,
    endemicEE = N[var /. endemicSols[[1]]];
    Print["Endemic equilibrium: ", Thread[var -> endemicEE]];,
    Print["No endemic equilibrium found"];
   ];
  ];

  (*Set default scanning parameters*)
  delta = 1/10;
  wRan = 1;
  hRan = 1;
  
  (*Fixed color mapping*)
  colorMap = <|"DFE" -> RGBColor[0, 0, 1], 
    "E1" -> RGBColor[0, 1, 0], "E2" -> RGBColor[0.6, 0.2, 0.8], 
    "EE-Stable" -> RGBColor[1, 1, 0], 
    "EE-Unstable" -> RGBColor[0.9, 0.4, 0.4], 
    "NoSol" -> RGBColor[1, 0, 0]|>;
  
  (*Extract bifurcation parameter indices for plotting*)
  bifParIdx1 = plotInd[[1]];
  bifParIdx2 = plotInd[[2]];
  bifP1Center = persVal[[bifParIdx1]];
  bifP2Center = persVal[[bifParIdx2]];

  (*Parameters for R-curve plotting*)
  par1 = par[[bifParIdx1]];
  par2 = par[[bifParIdx2]];
  
  (*Set up fixed parameters for R-curves*)
  fInd = Complement[Range[Length[par]], plotInd];
  fixedParams = Thread[par[[fInd]] -> persVal[[fInd]]];
  
  (*Compute intersection of all 4 R curves to set ranges*)
  intersectionPoint = {};
  If[R01 =!= Automatic && R02 =!= Automatic && R12 =!= Automatic && R21 =!= Automatic,
   (* Try to find intersection of all 4 curves *)
   numSol = Quiet[NSolve[{(R01 /. fixedParams) == 1, (R02 /. fixedParams) == 1,
     (R12 /. fixedParams) == 1, (R21 /. fixedParams) == 1}, {par1, par2}, Reals],
     {Power::infy, Infinity::indet, Power::indet, Solve::ratnz}];

   If[Length[numSol] > 0 && AllTrue[numSol[[1]], NumericQ[#[[2]]] &],
    intersectionPoint = {par1, par2} /. numSol[[1]];,

    (* Fallback: find R01-R02 intersection *)
    numSol = Quiet[NSolve[{(R01 /. fixedParams) == 1, (R02 /. fixedParams) == 1}, {par1, par2}, Reals],
      {Power::infy, Infinity::indet, Power::indet, Solve::ratnz}];
    If[Length[numSol] > 0 && AllTrue[numSol[[1]], NumericQ[#[[2]]] &],
     intersectionPoint = {par1, par2} /. numSol[[1]];
     Print["4-curve intersection not found, using R01-R02"];
    ];
   ];

   If[Length[intersectionPoint] == 2 && AllTrue[intersectionPoint, NumericQ],
    (*Set ranges to exceed intersection point - use rangeExtension parameter*)
    bifP1min = Max[intersectionPoint[[1]] * 0.1, 0.001];
    bifP1max = intersectionPoint[[1]] * rangeExtension;
    bifP2min = Max[intersectionPoint[[2]] * 0.1, 0.001];
    bifP2max = intersectionPoint[[2]] * rangeExtension;

    Print["Intersection: ", par1, "=", N[intersectionPoint[[1]]], " ", par2, "=", N[intersectionPoint[[2]]],
      " Ranges: ", par1, "\[Element][", N[bifP1min], ",", N[bifP1max], "] ", par2, "\[Element][", N[bifP2min], ",", N[bifP2max], "]"];,

    (* Fallback to default ranges *)
    Print["Warning: Could not find intersection, using default ranges"];
    bifP1min = Max[bifP1Center*(1 - wRan), 0.001];
    bifP1max = bifP1Center*(1 + wRan);
    bifP2min = Max[bifP2Center*(1 - hRan), 0.001];
    bifP2max = bifP2Center*(1 + hRan);
    ];,

   (* If R curves not provided, use default ranges *)
   bifP1min = Max[bifP1Center*(1 - wRan), 0.001];
   bifP1max = bifP1Center*(1 + wRan);
   bifP2min = Max[bifP2Center*(1 - hRan), 0.001];
   bifP2max = bifP2Center*(1 + hRan);
   ];
  
  (*Update ranges after potential intersection adjustment*)
  range1 = {bifP1min, bifP1max};
  range2 = {bifP2min, bifP2max};
  
  (*Detect which mode we're using*)
  useGridMode = (gridRes =!= Automatic);
  
  If[useGridMode, 
   (*GRID MODE*)
   bifP1vals = Table[bifP1min + (bifP1max - bifP1min)*k/(gridRes - 1), {k, 0, gridRes - 1}];
   bifP2vals = Table[bifP2min + (bifP2max - bifP2min)*k/(gridRes - 1), {k, 0, gridRes - 1}];,
   (*RANGE MODE*)
   step1 = (bifP1max - bifP1min)*delta;
   step2 = (bifP2max - bifP2min)*delta;
   bifP1vals = Table[bifP1, {bifP1, bifP1min, bifP1max, step1}];
   bifP2vals = Table[bifP2, {bifP2, bifP2min, bifP2max, step2}];
   ];
  
  totalPoints = Length[bifP1vals]*Length[bifP2vals];
  Print["Scanning ", totalPoints, " parameter combinations..."];
  
  (*Initialize for scanning*)
  res = {};
  currentProgress = 0;
  progressVar = 0;
  coexistencePoints = 0;
  Print[ProgressIndicator[Dynamic[progressVar]]];
  
  (*Define tolerance for numerical zeros*)
  zeroTol = 10^(-10);
  posTol = 10^(-8);
  
  (*OPTIMIZED MAIN SCANNING LOOP*)
  Do[
   Do[currentProgress++;
    progressVar = N[currentProgress/totalPoints];

    (*Create parameter values for this grid point*)
    numpar = persVal;
    numpar[[bifParIdx1]] = bifP1;
    numpar[[bifParIdx2]] = bifP2;
    conPar = Thread[par -> numpar];

    (*Pre-evaluate reproduction numbers with error handling*)
    R01val = Quiet[N[R01 /. conPar], {Power::infy, Overflow::unfl, Power::indet}];
    R02val = Quiet[N[R02 /. conPar], {Power::infy, Overflow::unfl, Power::indet}];

    (*Early validation*)
    If[!NumericQ[R01val] || !NumericQ[R02val] || !FiniteQ[R01val] || !FiniteQ[R02val] ||
       R01val < 0 || R02val < 0,
     classification = "NoSol";
     AppendTo[res, {N[bifP1], N[bifP2], classification}];
     Continue[];
     ];
    
    (*Fast classification for simple cases*)
    If[R01val < 1 && R02val < 1,
     classification = "DFE";
     AppendTo[res, {N[bifP1], N[bifP2], classification}];
     Continue[];
     ];
    
    If[R01val < 1 && R02val > 1,
     classification = "E2";
     AppendTo[res, {N[bifP1], N[bifP2], classification}];
     Continue[];
     ];
    
    If[R01val > 1 && R02val < 1,
     classification = "E1";
     AppendTo[res, {N[bifP1], N[bifP2], classification}];
     Continue[];
     ];
    
    (*Complex case: both can invade - need invasion numbers*)
    R12val = Quiet[N[R12 /. conPar], {Power::infy, Overflow::unfl, Power::indet}];
    R21val = Quiet[N[R21 /. conPar], {Power::infy, Overflow::unfl, Power::indet}];
    
    If[!NumericQ[R12val] || !NumericQ[R21val] || !FiniteQ[R12val] || !FiniteQ[R21val] ||
       R12val < 0 || R21val < 0,
     classification = "NoSol";
     AppendTo[res, {N[bifP1], N[bifP2], classification}];
     Continue[];
     ];
    
    Which[
     (*Competitive exclusion cases*)
     R12val < 1 && R21val < 1,
     classification = If[R01val > R02val, "E1", "E2"],

     R12val > 1 && R21val < 1,
     classification = "E1",

     R12val < 1 && R21val > 1,
     classification = "E2",
     
     (*Potential coexistence*)
     R12val > 1 && R21val > 1,
     coexistencePoints++;
     
     (*Robust equilibrium solving with timeout*)
     eq = TimeConstrained[
       Quiet[NSolve[Join[Thread[(RHS /. conPar) == 0], Thread[var >= 0]], var, Reals], 
         {NSolve::ratnz, NSolve::precw, NSolve::incs}], 
       3 (* 3 second timeout *)
       ];
     
     If[eq === $Aborted || Head[eq] =!= List,
      classification = "NoSol",
      
      EE = Select[eq, (varS /. #) > 5*posTol && (varI1 /. #) > 5*posTol && (varI2 /. #) > 5*posTol &];
      
      If[Length[EE] >= 1,
       (*Check stability with timeout*)
       jac = TimeConstrained[
         Quiet[N[D[RHS, {var}] /. conPar /. EE[[1]]], {Power::infy}], 
         1 (* 1 second timeout *)
         ];
       
       If[jac === $Aborted || !MatrixQ[jac] || !AllTrue[Flatten[jac], NumericQ],
        classification = "EE-Unstable",
        
        eigs = Quiet[Eigenvalues[jac], {Eigenvalues::eivn0}];
        If[AllTrue[eigs, NumericQ] && AllTrue[Re[eigs], # < -stabTol &],
         classification = "EE-Stable",
         classification = "EE-Unstable"
         ];
        ];,
       
       classification = "NoSol"
       ];
      ];,
     
     True,
     classification = "NoSol"
     ];
    
    (*Store result*)
    AppendTo[res, {N[bifP1], N[bifP2], classification}];
    
    , {bifP2, bifP2vals}], {bifP1, bifP1vals}];
  
  (*Debug output*)
  Print["Checked ", coexistencePoints, " potential coexistence points"];

  (*Process results*)
  outcomes = DeleteDuplicates[Table[res[[i, 3]], {i, 1, Length[res]}]];
  outcomeCounts = Table[Count[res, {_, _, outcomes[[i]]}], {i, 1, Length[outcomes]}];

  (*Summary with percentages - print before plot generation*)
  Do[Print[outcomes[[i]], ": ", outcomeCounts[[i]], " (", Round[100.*outcomeCounts[[i]]/Length[res]], "%)"], {i, 1, Length[outcomes]}];
  
  (*Create plot with fixed colors*)
  plotData = Table[Select[res, #[[3]] == outcomes[[i]] &][[All, 1 ;; 2]], {i, 1, Length[outcomes]}];
  plotMarkers = Table[{Style["\[FilledSquare]", colorMap[outcomes[[i]]]], 12}, {i, 1, Length[outcomes]}];
  
  (*Create equilibrium classification plot*)
  finalPlot = ListPlot[plotData, PlotMarkers -> plotMarkers, AspectRatio -> 1, 
    PlotRange -> {range1, range2}, GridLines -> Automatic, PlotLegends -> None];
  
  (*Add intersection point if computed*)
  intersectionPlot = {};
  If[Length[intersectionPoint] == 2,
   intersectionPlot = ListPlot[{intersectionPoint}, 
     PlotStyle -> Red, PlotMarkers -> {Style["\[CircleDot]", Red, 20]}];
   finalPlot = Show[finalPlot, intersectionPlot];
   ];
  
  (*Create R-curve plots*)
  rCurves = {};
  rCurveLabels = {};
  If[R01 =!= Automatic || R02 =!= Automatic || R12 =!= Automatic || R21 =!= Automatic,
   activeEquations = {};
   activeColors = {};
   If[R01 =!= Automatic,
    AppendTo[activeEquations, Evaluate[Quiet[FullSimplify[(R01 /. fixedParams) == 1],
      {Power::infy, Infinity::indet, Power::indet}]]];
    AppendTo[activeColors, Directive[Red, Thick]];
    AppendTo[rCurveLabels, "R01=1"];];
   If[R02 =!= Automatic,
    AppendTo[activeEquations, Evaluate[Quiet[FullSimplify[(R02 /. fixedParams) == 1],
      {Power::infy, Infinity::indet, Power::indet}]]];
    AppendTo[activeColors, Directive[Blue, Thick]];
    AppendTo[rCurveLabels, "R02=1"];];
   If[R12 =!= Automatic,
    AppendTo[activeEquations, Evaluate[Quiet[FullSimplify[(R12 /. fixedParams) == 1],
      {Power::infy, Infinity::indet, Power::indet}]]];
    AppendTo[activeColors, Directive[Purple, Thick]];
    AppendTo[rCurveLabels, "R12=1"];];
   If[R21 =!= Automatic,
    AppendTo[activeEquations, Evaluate[Quiet[FullSimplify[(R21 /. fixedParams) == 1],
      {Power::infy, Infinity::indet, Power::indet}]]];
    AppendTo[activeColors, Directive[Green, Thick]];
    AppendTo[rCurveLabels, "R21=1"];];
   Print["Rij equations: ", activeEquations];
   If[Length[activeEquations] > 0, 
    rCurves = {ContourPlot[Evaluate[activeEquations], 
       Evaluate[{par1, range1[[1]], range1[[2]]}], 
       Evaluate[{par2, range2[[1]], range2[[2]]}], 
       ContourStyle -> activeColors, PlotPoints -> 50]};];
   ];
  
  (*Combine plots with comprehensive legend*)
  If[Length[rCurves] > 0, 
   (*With R-curves*)
   eqLegend = PointLegend[Table[colorMap[outcomes[[i]]], {i, Length[outcomes]}],
      outcomes, LegendFunction -> "Panel", LegendLabel -> "Equilibria",
      LegendMarkers -> Table["\[FilledSquare]", {i, Length[outcomes]}],
      LegendMarkerSize -> 15];
   rLegend = LineLegend[activeColors, rCurveLabels, 
     LegendFunction -> "Panel", LegendLabel -> "R-curves"];
   
   legendItems = {eqLegend, rLegend};
   If[Length[intersectionPoint] == 2,
    intersectionLegend = PointLegend[{Red}, {"R01=R02=1"}, 
      LegendFunction -> "Panel", LegendLabel -> "Intersection",
      LegendMarkers -> {"\[CircleDot]"}, LegendMarkerSize -> 20];
    AppendTo[legendItems, intersectionLegend];
    ];
   
   finalPlot = Show[finalPlot, rCurves[[1]], Frame -> True, 
     FrameLabel -> {ToString[par1], ToString[par2]}, 
     PlotLabel -> "Equilibrium Classification with R-curves", ImageSize -> 450];
   finalPlot = Legended[finalPlot, Placed[Column[legendItems, Spacings -> 0.5], Right]];, 
   
   (*No R-curves*)
   legendItems = {PointLegend[Table[colorMap[outcomes[[i]]], {i, Length[outcomes]}], outcomes,
      LegendFunction -> "Panel", LegendMarkers -> Table["\[FilledSquare]", {i, Length[outcomes]}],
      LegendMarkerSize -> 15]};
   
   If[Length[intersectionPoint] == 2,
    intersectionLegend = PointLegend[{Red}, {"R01=R02=1"}, 
      LegendFunction -> "Panel", LegendLabel -> "Intersection",
      LegendMarkers -> {"\[CircleDot]"}, LegendMarkerSize -> 20];
    AppendTo[legendItems, intersectionLegend];
    ];
   
   finalPlot = Show[finalPlot, Frame -> True,
     FrameLabel -> {ToString[par1], ToString[par2]},
     PlotLabel -> "Equilibrium Classification", ImageSize -> 450];
   finalPlot = Legended[finalPlot, Placed[Column[legendItems, Spacings -> 0.5], Right]];
   ];

  (*Return NoSol points as errors*)
  noSolPoints = Select[res, #[[3]] == "NoSol" &];
  {finalPlot, noSolPoints, res}
  ];


simpleOptHopf[RHS_, var_, par_, coP_, optInd_, numTries_ : 10] := 
  Module[{p0val, bestAngle, bestValues, i, testValues, newP0Val, 
    posSols, complexEigs, eigs, currentAngle}, 
   Print["Starting simple optimization - angle-based objective"];
   p0val = par /. coP;
   newP0Val = p0val;
   {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, newP0Val];
   
   bestAngle = If[Length[complexEigs] > 0, 
     Module[{upperEig}, 
      upperEig = First[Select[complexEigs, Im[#] > 0 &]];
      ArcTan[Re[upperEig]/Im[upperEig]]*180/Pi], 90];
   bestValues = p0val[[optInd]];
   Print["Initial angle = ", bestAngle, " degrees"];
   
   Do[testValues = bestValues*(1 + 0.5*RandomReal[{-1, 1}, Length[optInd]]);
    testValues = Max[#, 0.001] & /@ testValues;
    newP0Val = p0val;
    newP0Val[[optInd]] = testValues;
    {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, newP0Val];
    
    currentAngle = If[Length[complexEigs] > 0, 
      Module[{upperEig}, 
       upperEig = First[Select[complexEigs, Im[#] > 0 &]];
       ArcTan[Re[upperEig]/Im[upperEig]]*180/Pi], 90];
    
    (* Maximize angle: move from negative towards positive *)
    If[currentAngle > bestAngle, 
     bestAngle = currentAngle;
     bestValues = testValues;
     Print["Improvement: angle = ", currentAngle, " degrees at ", testValues];
     
     (* Stop if we reach positive (Hopf detected) *)
     If[currentAngle > 0, 
      Print["*** HOPF BIFURCATION DETECTED: angle = ", currentAngle, " degrees ***"];
      Break[];];];, {i, 1, numTries}];
   
   Print["Best angle = ", bestAngle, " degrees"];
   Print["Best values = ", bestValues];
   Module[{finalP0Val}, 
    finalP0Val = p0val;
    finalP0Val[[optInd]] = bestValues;
    {bestAngle, bestValues, finalP0Val}]];

optHopf[RHS_, var_, par_, coP_, optInd_, timeLimit_ : 60, method_ : "NelderMead", 
   accGoal_ : 4, precGoal_ : 4, maxIter_ : 200] := 
  Module[{(*objFunc,*) p0val, varList, constraints, lowerBounds, upperBounds, 
    result, startTime, initialAngle, elapsedTime}, 
   Print["Starting Mathematica optimization - angle-based objective"];
   (*objFunc = objHopf[RHS, var, par, coP, optInd];*) (* COMMENTED OUT - objHopf not defined *)
   p0val = par /. coP;
   (*initialAngle = objFunc[p0val[[optInd]]];*) (* COMMENTED OUT *)
   initialAngle = -90; (* Default fallback value *)
   Print["Initial angle = ", initialAngle, " degrees"];
   
   Print["ERROR: optHopf requires objHopf function to be defined"];
   Print["Returning fallback values"];
   {initialAngle, p0val[[optInd]], p0val}
   
   (* COMMENTED OUT ENTIRE OPTIMIZATION SECTION
   varList = par[[optInd]];
   lowerBounds = 0.5*p0val[[optInd]];
   upperBounds = 1.5*p0val[[optInd]];
   constraints = Table[lowerBounds[[i]] <= varList[[i]] <= upperBounds[[i]], {i, Length[optInd]}];
   
   startTime = AbsoluteTime[];
   
   Module[{var1, var2}, 
    {var1, var2} = varList;
    result = TimeConstrained[
      NMaximize[
       {objFunc[{var1, var2}], And @@ constraints}, {var1, var2}, 
       Method -> method, MaxIterations -> maxIter, 
       AccuracyGoal -> accGoal, PrecisionGoal -> precGoal], 
      timeLimit, $Failed];];
   
   elapsedTime = Round[AbsoluteTime[] - startTime, 0.1];
   Print["Elapsed time: ", elapsedTime, " seconds"];
   
   If[result =!= $Failed && Length[result] >= 2, 
    Module[{maxAngle, optParams, bestValues, finalP0Val, bestAngle}, 
     {maxAngle, optParams} = result;
     bestValues = Values[optParams];
     bestAngle = objFunc[bestValues];
     
     If[bestAngle > 0, 
      Print["*** HOPF BIFURCATION DETECTED: angle = ", bestAngle, " degrees ***"];, 
      Print["Best angle = ", bestAngle, " degrees"];];
     Print["Best values = ", bestValues];
     
     finalP0Val = p0val;
     finalP0Val[[optInd]] = bestValues;
     {bestAngle, bestValues, finalP0Val}], 
    
    Print["Optimization failed - returning initial values"];
    Print["Best angle = ", initialAngle, " degrees"];
    Print["Best values = ", p0val[[optInd]]];
    {initialAngle, p0val[[optInd]], p0val}]
    *)
   ];

cont[RHS_, var_, par_, p0val_, stepSize_, plotInd_ : {3, 4}, bifInd_ : 4, 
   analyticalPlot_ : None] := 
  Module[{posSols, bifParamIndex, currentVal, curve, newP0Val, paramVal, 
    complexEigs, eigs, eq, maxRealPart, complexPairs, bifPlot, plotC, 
    forwardCurve, backwardCurve, i, maxIter = 10, successCount = 0, 
    hopfPoints, testedPoints = {}, startingPoint, adaptiveStep, currentStep, 
    overlayPossible, beta1Index, beta2Index, beta1Val, beta2Val, startPosInCurve}, 
   
   (* Use parameter indices *)
   bifParamIndex = bifInd;
   beta1Index = plotInd[[1]];
   beta2Index = plotInd[[2]];
   currentVal = p0val[[bifParamIndex]];
   
   Print["Starting bifParam index ", bifParamIndex, " at value ", currentVal];
   Print["Plot parameters at indices: ", plotInd, " with values: ", p0val[[plotInd]]];
   
   (* Check overlay compatibility *)
   overlayPossible = analyticalPlot =!= None && MemberQ[plotInd, bifParamIndex];
   If[overlayPossible, 
    Print["Will overlay on analytical plot (bifParam in plotInd)"];, 
    If[analyticalPlot =!= None, 
     Print["Standalone plot only (bifParam not in plotInd)"];];];
   
   (* Get starting equilibrium *)
   {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, p0val];
   If[Length[posSols] == 0, 
    Print["ERROR: No equilibrium found at starting p0val!"];
    Return[{}];];
   
   (* Include starting equilibrium with beta coordinates *)
   {beta1Val, beta2Val} = p0val[[plotInd]];
   startingPoint = <|"Parameter" -> currentVal, "Values" -> posSols[[1]], 
     "MaxRealPart" -> If[Length[eigs] > 0, Max[Re[eigs]], 0], 
     "ComplexPairs" -> Length[complexEigs]/2, "Eigenvalues" -> eigs, 
     "ComplexEigs" -> complexEigs, "Beta1" -> beta1Val, "Beta2" -> beta2Val|>;
   
   curve = {};
   forwardCurve = {startingPoint};
   backwardCurve = {};
   
   (* Adaptive step size function *)
   adaptiveStep[stepNum_, baseStep_] := Which[
     stepNum <= 2, baseStep/5, 
     stepNum <= 4, baseStep/3, 
     stepNum <= 6, baseStep/2, 
     True, baseStep];
   
   (* Forward continuation *)
   i = 1;
   While[True, 
    currentStep = adaptiveStep[i, stepSize];
    paramVal = currentVal + i*currentStep;
    AppendTo[testedPoints, {"Forward", paramVal}];
    newP0Val = p0val;
    newP0Val[[bifParamIndex]] = paramVal;
    {beta1Val, beta2Val} = newP0Val[[plotInd]];
    {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, newP0Val];
    
    If[Length[posSols] > 0, 
     eq = posSols[[1]];
     maxRealPart = If[Length[eigs] > 0, Max[Re[eigs]], 0];
     complexPairs = Length[complexEigs]/2;
     successCount++;
     
     If[Length[complexEigs] > 0 && Abs[Max[Re[complexEigs]]] < 0.01, 
      Print["*** HOPF BIFURCATION FOUND at parameter = ", paramVal, " ***"];];
     
     AppendTo[forwardCurve, 
      <|"Parameter" -> paramVal, "Values" -> eq, "MaxRealPart" -> maxRealPart, 
        "ComplexPairs" -> complexPairs, "Eigenvalues" -> eigs, 
        "ComplexEigs" -> complexEigs, "Beta1" -> beta1Val, "Beta2" -> beta2Val|>];
     i++;
     If[i >= maxIter, Break[]];, 
     Break[];];];
   
   (* Backward continuation *)
   i = 1;
   While[True, 
    currentStep = adaptiveStep[i, stepSize];
    paramVal = currentVal - i*currentStep;
    If[paramVal <= 0, Break[]];
    AppendTo[testedPoints, {"Backward", paramVal}];
    newP0Val = p0val;
    newP0Val[[bifParamIndex]] = paramVal;
    {beta1Val, beta2Val} = newP0Val[[plotInd]];
    {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, newP0Val];
    
    If[Length[posSols] > 0, 
     eq = posSols[[1]];
     maxRealPart = If[Length[eigs] > 0, Max[Re[eigs]], 0];
     complexPairs = Length[complexEigs]/2;
     successCount++;
     
     If[Length[complexEigs] > 0 && Abs[Max[Re[complexEigs]]] < 0.01, 
      Print["*** HOPF BIFURCATION FOUND at parameter = ", paramVal, " ***"];];
     
     PrependTo[backwardCurve, 
      <|"Parameter" -> paramVal, "Values" -> eq, "MaxRealPart" -> maxRealPart, 
        "ComplexPairs" -> complexPairs, "Eigenvalues" -> eigs, 
        "ComplexEigs" -> complexEigs, "Beta1" -> beta1Val, "Beta2" -> beta2Val|>];
     i++;
     If[i > maxIter, Break[]];, 
     Break[];];];
   
   (* Combine curves *)
   curve = Join[backwardCurve, forwardCurve];
   Print["Found ", Length[curve], " points, tested ", Length[testedPoints], " parameters"];
   If[Length[curve] > 0, 
    Print["Parameter range: ", curve[[1]]["Parameter"], " to ", curve[[-1]]["Parameter"]];];
   curve];

hopfD[curve_] := 
  Module[{hopfPoints, stabilityPlot, hopfPlot}, 
   If[Length[curve] < 3, 
    Print["hopfD needs >=3 points, got ", Length[curve]];
    Return[{}]];
   
   (* Check if curve has required fields *)
   If[!KeyExistsQ[curve[[1]], "MaxRealPart"], 
    Print["Error: curve missing stability data. Run cont first."];
    Return[{}]];
   
   stabilityPlot = ListLinePlot[{#["Parameter"], #["MaxRealPart"]} & /@ curve, 
     AxesLabel -> {"Parameter", "Max Re(\[Lambda])"}, 
     PlotLabel -> "Stability Analysis", PlotStyle -> Blue, 
     GridLines -> {None, {0}}, GridLinesStyle -> Directive[Red, Dashed], 
     ImageSize -> 400];
   
   hopfPoints = Select[Table[
      Module[{prev, curr, next, prevReal, currReal, hasComplex}, 
       {prev, curr, next} = {curve[[i - 1]], curve[[i]], 
         If[i < Length[curve], curve[[i + 1]], curve[[i]]]};
       {prevReal, currReal} = {prev["MaxRealPart"], curr["MaxRealPart"]};
       hasComplex = curr["ComplexPairs"] > 0;
       If[hasComplex && prevReal*currReal < 0, 
        <|"parBif" -> curr["Parameter"], "Type" -> "Hopf"|>, Nothing]], 
      {i, 2, Length[curve] - 1}], # =!= Nothing &];
   hopfPoints];

scanPar[RHS_, var_, par_, p0val_, plotInd_, gridRes_ : Automatic, 
   plot_ : Automatic, hTol_ : 0.01, delta_ : 1/20, wRan_ : 1, hRan_ : 1, 
   R01_ : Automatic, R02_ : Automatic, R21_ : Automatic, R12_ : Automatic] := 
  Module[{betterColors, useGridMode, beta1Index, beta2Index, beta1Center, 
    beta2Center, beta1min, beta1max, beta2min, beta2max, beta1vals, beta2vals, 
    totalPoints, res, currentProgress, progressVar, tol, inP, numpar, conPar, 
    sol, I1val, I2val, currentType, jac, eigs, maxRealPart, outcomes, 
    outcomeCounts, finalPlot}, 
   
   betterColors = {RGBColor[0, 0, 1], RGBColor[0, 1, 0], RGBColor[0.8, 0.5, 0.2], 
     RGBColor[0.6, 0.3, 0.8], RGBColor[1, 0, 0], RGBColor[1.0, 0.0, 1.0]};
   
   useGridMode = (gridRes =!= Automatic);
   beta1Index = If[Length[plotInd] >= 1, plotInd[[1]], 1];
   beta2Index = If[Length[plotInd] >= 2, plotInd[[2]], 2];
   beta1Center = p0val[[beta1Index]];
   beta2Center = p0val[[beta2Index]];
   
   Print["Varying parameters at indices ", {beta1Index, beta2Index}, 
    " with center values: ", {beta1Center, beta2Center}];
   
   If[useGridMode, 
    Print["Using grid mode with ", gridRes, "x", gridRes, " points..."];
    beta1min = Max[beta1Center*(1 - wRan), 0.001];
    beta1max = beta1Center*(1 + wRan);
    beta2min = Max[beta2Center*(1 - hRan), 0.001];
    beta2max = beta2Center*(1 + hRan);
    beta1vals = Table[beta1min + (beta1max - beta1min)*k/(gridRes - 1), {k, 0, gridRes - 1}];
    beta2vals = Table[beta2min + (beta2max - beta2min)*k/(gridRes - 1), {k, 0, gridRes - 1}];, 
    
    Print["Using range mode with step size ", delta, "..."];
    beta1min = Max[beta1Center*(1 - wRan), 0.001];
    beta1max = beta1Center*(1 + wRan);
    beta2min = Max[beta2Center*(1 - hRan), 0.001];
    beta2max = beta2Center*(1 + hRan);
    beta1vals = Table[beta1, {beta1, beta1min, beta1max, delta}];
    beta2vals = Table[beta2, {beta2, beta2min, beta2max, delta}];];
   
   totalPoints = Length[beta1vals]*Length[beta2vals];
   Print["Total points to scan: ", totalPoints];
   res = {};
   currentProgress = 0;
   progressVar = 0;
   tol = 10^(-6);
   inP = Table[{var[[j]], 1/Length[var]}, {j, Length[var]}];
   Print[ProgressIndicator[Dynamic[progressVar]]];
   
   (* Main scanning loop *)
   Do[Do[currentProgress++;
     progressVar = N[currentProgress/totalPoints];
     numpar = p0val;
     numpar[[beta1Index]] = beta1;
     numpar[[beta2Index]] = beta2;
     conPar = Thread[par -> numpar];
     sol = Quiet[FindRoot[RHS /. conPar, inP]];
     
     If[Head[sol] === List, 
      {I1val, I2val} = {var[[2]], var[[3]]} /. sol; (* Assuming I1,I2 are at positions 2,3 *)
      currentType = Which[
        I1val < tol && I2val < tol, "DFE", 
        I1val >= tol && I2val < tol, "E1", 
        I1val < tol && I2val >= tol, "E2", 
        True, 
        jac = D[RHS, {var}] /. conPar /. sol;
        eigs = Chop[Eigenvalues[jac] // N];
        maxRealPart = Max[Re[eigs]];
        If[maxRealPart < 0, "EEstable", "EEunstable"]];
      AppendTo[res, {N[beta1], N[beta2], currentType}];, 
      AppendTo[res, {N[beta1], N[beta2], "NoSol"}];];, {beta2, beta2vals}];, {beta1, beta1vals}];
   
   Print["Scanning complete!"];
   outcomes = DeleteDuplicates[Table[res[[i, 3]], {i, 1, Length[res]}]];
   outcomeCounts = Table[Count[res, {_, _, outcomes[[i]]}], {i, 1, Length[outcomes]}];
   
   finalPlot = ListPlot[
     Table[Select[res, #[[3]] == outcomes[[i]] &][[All, 1 ;; 2]], {i, Length[outcomes]}], 
     PlotMarkers -> Table[{Style["\[FilledCircle]", betterColors[[i]]], 14}, {i, Min[Length[betterColors], Length[outcomes]]}], 
     PlotLegends -> Table[outcomes[[i]] <> " (" <> ToString[outcomeCounts[[i]]] <> ")", {i, Length[outcomes]}], 
     AspectRatio -> 1, 
     PlotLabel -> "Parameter Space Analysis (" <> ToString[Length[res]] <> " points)", 
     FrameLabel -> {"Parameter " <> ToString[beta1Index], "Parameter " <> ToString[beta2Index]}, 
     ImageSize -> 450];
   
   {finalPlot, {}, res}];

pertIC[equilibrium_, var_, factor_ : 0.1, minq_ : 0.001, n_ : 1] := 
  Module[{eqVals, pertVals}, 
   eqVals = var /. equilibrium;
   Table[
    pertVals = Max[#, minq] & /@ (eqVals + factor*RandomReal[{-1, 1}, Length[var]]*Abs[eqVals]);
    pertVals = pertVals (Total[eqVals]/Total[pertVals]);
    Thread[var -> pertVals], {n}]];

TS[RHS_, var_, par_, p0val_, tmax_ : 100] := 
  Module[{varT, initialConditions, rhsT, odeSystem, sol, solFun, coP, t}, 
   (* Create substitution rules *)
   coP = Thread[par -> p0val];
   
   (* Time-dependent variables *)
   varT = Through[var[t]];
   
   (* Initial conditions using parameter values *)
   initialConditions = Thread[(varT /. t -> 0) == (var /. coP)];
   
   (* RHS with substituted parameters and time variables *)
   rhsT = RHS /. Thread[var -> varT] /. coP;
   
   (* Full ODE system *)
   odeSystem = Join[Thread[D[varT, t] == rhsT], initialConditions];
   
   (* Solve ODE system *)
   sol = Quiet@NDSolve[odeSystem, varT, {t, 0, tmax}, Method -> {"BDF"}];
   If[!MatchQ[sol, {(_Rule | _List) ...}], 
    Print["NDSolve failed"];
    Return[$Failed];];
   
   (* Extract solution functions *)
   solFun = #[[0]] & /@ (varT /. sol[[1]]);
   
   (* Return full solution functions *)
   solFun];

intEq[RHS_, var_, par_, p0val_, att_] := 
  Module[{eqs, sols, cleanSols, positiveSols, threshold = 10^(-10), 
    perturbedICs, attRules, coP}, 
   (* Create substitution rules *)
   coP = Thread[par -> p0val];
   eqs = Thread[RHS == 0] /. coP;
   
   (* Convert att list to replacement rules *)
   attRules = Thread[var -> att];
   
   (* Use attRules as seed instead of calling TS *)
   perturbedICs = Join[{var /. attRules}, 
     var /. # & /@ pertIC[attRules, var, 0.1, 15], 
     var /. # & /@ pertIC[attRules, var, 0.3, 8]];
   
   sols = Table[Quiet[FindRoot[eqs, Evaluate[Thread[{var, perturbedICs[[i]]}]], Method -> "Newton"]], {i, 1, Length[perturbedICs]}];
   sols = Select[sols, Head[#] === List && Length[#] == Length[var] &];
   
   cleanSols = sols /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
   cleanSols = DeleteDuplicatesBy[cleanSols, Round[Sort[var /. #], 10^(-6)] &];
   positiveSols = Select[cleanSols, And @@ (# > threshold & /@ (var /. #)) &];
   
   Print["Found ", Length[positiveSols], " equilibria"];
   positiveSols];
