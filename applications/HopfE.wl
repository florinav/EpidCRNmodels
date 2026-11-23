(* ::Package:: *)

BeginPackage["Hopf`"]

fpHopf::usage = "pos=fpHopf[RHS, var, par, p0val] yields all positive solutions for p0val, and examines
the eigenvalues of the first for possible closeness to Hopf";
pertIC::usage = "pertIC[equilibrium, var, factor, n] creates perturbed initial conditions from an equilibrium."
TS::usage = "TS[RHS, var, par, p0val, tmax:100] computes time series and its limit from initial conditions 
1/var//Length."
intEqV::usage = "equilibria = intEqV[RHS, var, par, p0val] finds interior equilibrium points with verbose output and symbolic equilibria checking."

intEq::usage = "equilibria = intEq[RHS, var, par, p0val, att] finds interior equilibrium points using FindRoot."

cont::usage = "curve = cont2[RHS, var, att, par, p0val, range, inp] performs continuation in user specified range."

hopfD::usage = "hopfPoints = hopfD[curve] detects Hopf bifurcation points from continuation curve data."
phase::usage = "phasePlot = phase[RHS, var, par, p0val, tmax] creates phase portrait with multiple trajectories around equilibrium."
wF::usage = "wF[bifParam, p0val, RHS, var] performs workflow analysis at bifurcation parameter."
scanPar::usage = "{plot,errors,res} = scanPar[RHS, var, par, p0val, plotInd, plot] scans parameter space to classify stability regions."
pltRegions::usage = "pltRegions[scanResults, paramNames] plots parameter space regions from scanPar results."
ploR::usage = "ploR[p0val, R0A, E0, E1, E2, R12, R21] plots analytical boundaries from reproduction numbers."
objHopf;optHopf;pertIC;findH;simpleOptHopf;cont2;contS;
scan4::usage = "scanPR[RHS_, var_, par_, p0val_: Automatic, plotInd_: Automatic, plot_: Automatic,
  \[Beta]1range_: {1/10, 3}, \[Beta]2range_: {1/10, 3}, delta_: 1/20, 
  wRan_: 0.5, hRan_: 0.5, gridRes_: 8, hTol_: 0.01]";
  
Begin["`Private`"]
contS[RHS_, var_, par_, p0val_, step1_, step2_, initialBifInd_, 
tol_: 0.01, plot_: None] := 
 Module[{currentP0Val, currentAngle, globalBestAngle, globalBestP0Val, 
   usedParams, remainingParams, nextPair, cycleAngle, cycleP0Val, 
   cycle = 0, totalParams, currentBifInd, allResults = {}},
  
  Print["contS: Starting with ", Length[par], " parameters"];
  
  currentP0Val = p0val;
  globalBestP0Val = p0val;
  currentAngle = fpHopf[RHS, var, par, currentP0Val][[3]];
  globalBestAngle = currentAngle;
  usedParams = {};
  totalParams = Length[par];
  currentBifInd = initialBifInd;
  
  Print["Initial angle: ", currentAngle];
  Print["Total parameters: ", totalParams];
  
  While[Length[usedParams] < totalParams - 1,
   cycle++;
   Print["\n=== CYCLE ", cycle, " ==="];
   Print["Optimizing parameters: ", currentBifInd];
   
   Print["contS: Calling cont2..."];
   {cycleAngle, cycleP0Val} = cont2[RHS, var, par, currentP0Val, step1, step2, currentBifInd, tol, None];
   Print["contS: cont2 completed, angle = ", cycleAngle];
   
   currentP0Val = cycleP0Val;
   currentAngle = cycleAngle;
   
   AppendTo[allResults, <|"Cycle" -> cycle, "Parameters" -> currentBifInd, "EndAngle" -> cycleAngle|>];
   
   If[cycleAngle > globalBestAngle,
    globalBestAngle = cycleAngle;
    globalBestP0Val = cycleP0Val;
    Print["NEW GLOBAL BEST: angle = ", globalBestAngle];
    ];
   
   usedParams = Union[usedParams, currentBifInd];
   Print["Used parameters so far: ", usedParams];
   
   remainingParams = Complement[Range[totalParams], usedParams];
   
   If[Length[remainingParams] >= 2,
    nextPair = Take[remainingParams, 2];
    currentBifInd = nextPair;
    Print["Next pair: ", nextPair];,
    Print["Not enough unused parameters remaining"];
    Break[];
    ];
   ];
  
  Print["\n=== FINAL RESULTS ==="];
  Print["Total cycles completed: ", cycle];
  Print["Global best angle: ", globalBestAngle];
  
  Print["contS: Returning results..."];
  {globalBestAngle, globalBestP0Val, allResults}
  ];
  
  cont2[RHS_, var_, par_, p0val_, step1_, step2_, bifInd_, tol_: 0.01, plot_: None] := 
 Module[{currentP0Val, currentAngle, bestAngle, bestP0Val, 
   testP0Val, testResult, testAngle, maxIter = 100, iter = 0, 
   currentStep1, currentStep2, searchPath = {}, overlayPlot, 
   step, bestDirection, bestParam, bestImprovement, 
   param1Plus, param1Minus, param2Plus, param2Minus},
  
  (* Initialize *)
  currentP0Val = p0val;
  bestP0Val = p0val;
  
  (* Safe evaluation of fpHopf *)
  testResult = Quiet[fpHopf[RHS, var, par, currentP0Val]];
  If[!ListQ[testResult] || Length[testResult] < 3,
   Print["Error: fpHopf failed at starting point"];
   Return[{-1000, p0val}];
   ];
  
  currentAngle = testResult[[3]];
  bestAngle = currentAngle;
  currentStep1 = step1;
  currentStep2 = step2;
  
  (* Record starting point for search path *)
  AppendTo[searchPath, currentP0Val[[bifInd]]];
  
  Print["Starting optimization at angle = ", currentAngle];
  Print["Optimizing parameters at indices: ", bifInd, " with steps: {", currentStep1, ", ", currentStep2, "}"];
  
  (* Main optimization loop *)
  While[iter < maxIter && Max[currentStep1, currentStep2] > tol/100,
   iter++;
   
   (* Test all 4 directions *)
   bestImprovement = 0;
   bestDirection = 0;
   bestParam = 0;
   param1Plus = -1000; param1Minus = -1000; param2Plus = -1000; param2Minus = -1000;
   
   (* Test param 1 + direction *)
   testP0Val = currentP0Val;
   testP0Val[[bifInd[[1]]]] = currentP0Val[[bifInd[[1]]]] + currentStep1;
   
   If[testP0Val[[bifInd[[1]]]] > 0,
    testResult = Quiet[fpHopf[RHS, var, par, testP0Val]];
    If[ListQ[testResult] && Length[testResult] >= 3,
     testAngle = testResult[[3]];
     If[NumericQ[testAngle],
      param1Plus = testAngle - currentAngle;
      If[param1Plus > bestImprovement,
       bestImprovement = param1Plus;
       bestDirection = 1;
       bestParam = 1;
       ];
      ];
     ];
    ];
   
   (* Test param 1 - direction *)
   testP0Val = currentP0Val;
   testP0Val[[bifInd[[1]]]] = currentP0Val[[bifInd[[1]]]] - currentStep1;
   
   If[testP0Val[[bifInd[[1]]]] > 0,
    testResult = Quiet[fpHopf[RHS, var, par, testP0Val]];
    If[ListQ[testResult] && Length[testResult] >= 3,
     testAngle = testResult[[3]];
     If[NumericQ[testAngle],
      param1Minus = testAngle - currentAngle;
      If[param1Minus > bestImprovement,
       bestImprovement = param1Minus;
       bestDirection = -1;
       bestParam = 1;
       ];
      ];
     ];
    ];
   
   (* Test param 2 + direction *)
   testP0Val = currentP0Val;
   testP0Val[[bifInd[[2]]]] = currentP0Val[[bifInd[[2]]]] + currentStep2;
   
   If[testP0Val[[bifInd[[2]]]] > 0,
    testResult = Quiet[fpHopf[RHS, var, par, testP0Val]];
    If[ListQ[testResult] && Length[testResult] >= 3,
     testAngle = testResult[[3]];
     If[NumericQ[testAngle],
      param2Plus = testAngle - currentAngle;
      If[param2Plus > bestImprovement,
       bestImprovement = param2Plus;
       bestDirection = 1;
       bestParam = 2;
       ];
      ];
     ];
    ];
   
   (* Test param 2 - direction *)
   testP0Val = currentP0Val;
   testP0Val[[bifInd[[2]]]] = currentP0Val[[bifInd[[2]]]] - currentStep2;
   
   If[testP0Val[[bifInd[[2]]]] > 0,
    testResult = Quiet[fpHopf[RHS, var, par, testP0Val]];
    If[ListQ[testResult] && Length[testResult] >= 3,
     testAngle = testResult[[3]];
     If[NumericQ[testAngle],
      param2Minus = testAngle - currentAngle;
      If[param2Minus > bestImprovement,
       bestImprovement = param2Minus;
       bestDirection = -1;
       bestParam = 2;
       ];
      ];
     ];
    ];
   
   Print["4-way test: P1+:", param1Plus, " P1-:", param1Minus, " P2+:", param2Plus, " P2-:", param2Minus];
   
   (* If best improvement > tolerance, go in that direction *)
   If[bestImprovement > tol,
    step = bestDirection * If[bestParam == 1, currentStep1, currentStep2];
    
    Print["Best direction: Param ", bifInd[[bestParam]], " ", If[bestDirection > 0, "+", "-"], 
          " step=", If[bestParam == 1, currentStep1, currentStep2], " improvement=", bestImprovement];
    
    (* Keep going in this direction while improving *)
    While[True,
     testP0Val = currentP0Val;
     testP0Val[[bifInd[[bestParam]]]] = currentP0Val[[bifInd[[bestParam]]]] + step;
     
     If[testP0Val[[bifInd[[bestParam]]]] <= 0, Break[]];
     
     testResult = Quiet[fpHopf[RHS, var, par, testP0Val]];
     
     If[!ListQ[testResult] || Length[testResult] < 3, Break[]];
     
     testAngle = testResult[[3]];
     
     If[!NumericQ[testAngle] || testAngle <= currentAngle + tol, Break[]];
     
     (* Accept the step *)
     currentP0Val = testP0Val;
     currentAngle = testAngle;
     AppendTo[searchPath, currentP0Val[[bifInd]]];
     
     If[currentAngle > bestAngle,
      bestAngle = currentAngle;
      bestP0Val = currentP0Val;
      ];
     
     Print["Continue: angle = ", currentAngle];
     ];,
    
    (* No improvement > tolerance in any direction, reduce step sizes *)
    If[bestImprovement > 2*tol,
     Print["Best improvement ", bestImprovement, " < tolerance ", tol, " but > 2*tol, continuing"];,
     currentStep1 = currentStep1/2;
     currentStep2 = currentStep2/2;
     Print["All 4 directions failed (best=", bestImprovement, "], reducing step sizes to: {", currentStep1, ", ", currentStep2, "}"];
     ];
    ];
   ];
  
  Print["Final result: angle = ", bestAngle];
  Print["Best parameters: ", bestP0Val[[bifInd]]];
  Print["Final step sizes: {", currentStep1, ", ", currentStep2, "}"];
  Print["Iterations: ", iter];
  
  (* Create overlay plot if requested *)
  If[plot =!= None && Length[searchPath] > 1,
   Print["Creating overlay plot with search path..."];
   
   overlayPlot = Show[
    plot,
    Graphics[{
      Blue, PointSize[0.008],
      Point[searchPath],
      Blue, Thick,
      Line[searchPath],
      Red, PointSize[0.012],
      Point[searchPath[[1]]],
      Green, PointSize[0.012],
      Point[bestP0Val[[bifInd]]],
      Text[Style["Start", FontSize -> 10, FontWeight -> Bold, FontColor -> Red], 
           searchPath[[1]] + {0.02, 0.02}],
      Text[Style["Best", FontSize -> 10, FontWeight -> Bold, FontColor -> Green], 
           bestP0Val[[bifInd]] + {0.02, 0.02}]
      }],
    PlotLabel -> "4-Way Search with Independent Step Sizes"
    ];
   Print[overlayPlot];
   ];
  
  {bestAngle, bestP0Val}
  ];


(*
(* Full test example *)
var = {s, i, r1, r2, r3};
par = {\[Lambda], \[Mu], \[Beta], \[Beta]1, \[Beta]2, \[Beta]3, \[Gamma], \[Gamma]1, \[Gamma]2, \[Gamma]3};
RHS = {
  \[Lambda] - \[Mu] s - \[Beta] s i + \[Gamma]3 r3,
  i (\[Beta] s + \[Beta]1 r1 + \[Beta]2 r2 + \[Beta]3 r3 - \[Gamma] - \[Mu]),
  \[Gamma] i - r1 (\[Beta]1 i + \[Gamma]1 + \[Mu]),
  \[Gamma]1 r1 - r2 (\[Beta]2 i + \[Gamma]2 + \[Mu]),
  \[Gamma]2 r2 - r3 (\[Beta]3 i + \[Gamma]3 + \[Mu])
};
p0val = {0.02, 0.012, 0.3, 0.1, 0.05, 0.02, 0.1, 0.05, 0.03, 0.01};

(* Without plot *)
{bestAngle, bestP0Val} = cont2[RHS, var, par, p0val, 0.001, {2, 3}, 0.01];

(* With plot overlay *)
basePlot = ContourPlot[x + y, {x, 0.005, 0.02}, {y, 0.2, 0.4}];
{bestAngle, bestP0Val} = cont2[RHS, var, par, p0val, 0.001, {2, 3}, 0.01, basePlot];
*)
(* Simplified fpHopf that outputs angle as 3rd output, eigs as 4th *)
fpHopf[RHS_, var_, par_, p0val_] := 
 Module[{eqns, allSols, posSols, eq, jac, eigs, complexEigs, 
 upperEig, angle, coP},
  
  (* Create substitution rules *)
  coP = Thread[par -> p0val];
  
  (* Create equations and find positive solutions *)
  eqns = Thread[(RHS /. coP) == 0];
  allSols = Quiet[NSolveValues[eqns, var, Reals]];
  allSols = Union[allSols, SameTest -> (Norm[#1 - #2] < 10^-10 &)];
  allSols = Select[allSols, AllTrue[#, # >= 0 &] &]; 
  posSols = Select[allSols, AllTrue[#, # > 0 &] &];
  
  If[Length[posSols] == 0, 
   Return[{{}, {}, -90, {}}];
   ];
  
  (* Use first positive solution for eigenvalue analysis *)
  eq = posSols[[1]];
  jac = D[RHS, {var}] /. coP /. Thread[var -> eq];
  eigs = Chop[Eigenvalues[jac] // N];
  complexEigs = Select[eigs, Im[#] != 0 &];
  
  (* Calculate angle *)
  angle = If[Length[complexEigs] > 0,
    upperEig = First[Select[complexEigs, Im[#] > 0 &]];
    ArcTan[Re[upperEig]/Im[upperEig]] * 180/Pi,
    -90
    ];
   
  {posSols, complexEigs, angle, eigs}
  ];

(* Fixed simpleOptHopf function *)
simpleOptHopf[RHS_, var_, par_, coP_, optInd_, numTries_: 10] := 
 Module[{p0val, bestAngle, bestValues, i, testValues, newP0Val, posSols, complexEigs, eigs, currentAngle},
  Print["Starting simple optimization - angle-based objective"];
  
  p0val = par /. coP;
  newP0Val = p0val;
  {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, newP0Val];
  
  bestAngle = If[Length[complexEigs] > 0,
    Module[{upperEig},
     upperEig = First[Select[complexEigs, Im[#] > 0 &]];
     ArcTan[Re[upperEig]/Im[upperEig]] * 180/Pi
     ],
    90
    ];
  bestValues = p0val[[optInd]];
  
  Print["Initial angle = ", bestAngle, " degrees"];
  
  Do[
   testValues = bestValues * (1 + 0.5*RandomReal[{-1, 1}, Length[optInd]]);
   testValues = Max[#, 0.001] & /@ testValues;
   
   newP0Val = p0val;
   newP0Val[[optInd]] = testValues;
   {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, newP0Val];
   
   currentAngle = If[Length[complexEigs] > 0,
     Module[{upperEig},
      upperEig = First[Select[complexEigs, Im[#] > 0 &]];
      ArcTan[Re[upperEig]/Im[upperEig]] * 180/Pi
      ],
     90
     ];
   
   (* Maximize angle: move from negative towards positive *)
   If[currentAngle > bestAngle,
    bestAngle = currentAngle;
    bestValues = testValues;
    Print["Improvement: angle = ", currentAngle, " degrees at ", testValues];
    
    (* Stop if we reach positive (Hopf detected) *)
    If[currentAngle > 0,
     Print["*** HOPF BIFURCATION DETECTED: angle = ", currentAngle, " degrees ***"];
     Break[];
     ];
    ];
   
   , {i, 1, numTries}];
  
  Print["Best angle = ", bestAngle, " degrees"];
  Print["Best values = ", bestValues];
  
  Module[{finalP0Val},
   finalP0Val = p0val;
   finalP0Val[[optInd]] = bestValues;
   {bestAngle, bestValues, finalP0Val}
   ]
  ];

(* Fixed optHopf function *)
optHopf[RHS_, var_, par_, coP_, optInd_, timeLimit_: 60, method_: "NelderMead", 
        accGoal_: 4, precGoal_: 4, maxIter_: 200] := 
 Module[{objFunc, p0val, varList, constraints, lowerBounds, upperBounds, result, startTime, initialAngle, elapsedTime},
  Print["Starting Mathematica optimization - angle-based objective"];
  
  objFunc = objHopf[RHS, var, par, coP, optInd];
  p0val = par /. coP;
  
  initialAngle = objFunc[p0val[[optInd]]];
  Print["Initial angle = ", initialAngle, " degrees"];
  
  varList = par[[optInd]];
  lowerBounds = 0.5 * p0val[[optInd]];
  upperBounds = 1.5 * p0val[[optInd]];
  
  constraints = Table[
    lowerBounds[[i]] <= varList[[i]] <= upperBounds[[i]], 
    {i, Length[optInd]}
    ];
  
  startTime = AbsoluteTime[];
  
  (* FIX: Extract individual variables for NMaximize *)
  Module[{var1, var2},
   {var1, var2} = varList;
   
   result = TimeConstrained[
    NMaximize[  (* Maximize angle: from negative towards positive *)
     {objFunc[{var1, var2}], And @@ constraints},
     {var1, var2},
     Method -> method,
     MaxIterations -> maxIter,
     AccuracyGoal -> accGoal,
     PrecisionGoal -> precGoal
     ],
    timeLimit,
    $Failed
    ];
   ];
  
  (* Print why optimization stopped *)
  Print["NMaximize result structure: ", result];
  
  elapsedTime = Round[AbsoluteTime[] - startTime, 0.1];
  Print["Elapsed time: ", elapsedTime, " seconds"];
  
  If[result =!= $Failed && Length[result] >= 2,
   Module[{maxAngle, optParams, bestValues, finalP0Val, bestAngle},
    {maxAngle, optParams} = result;
    bestValues = Values[optParams];
    bestAngle = objFunc[bestValues];  (* Get the actual angle *)
    
    If[bestAngle > 0,
     Print["*** HOPF BIFURCATION DETECTED: angle = ", bestAngle, " degrees ***"];,
     Print["Best angle = ", bestAngle, " degrees"];
     ];
    Print["Best values = ", bestValues];
    
    finalP0Val = p0val;
    finalP0Val[[optInd]] = bestValues;
    {bestAngle, bestValues, finalP0Val}
    ],
   Print["Optimization failed - returning initial values"];
   Print["Best angle = ", initialAngle, " degrees"];
   Print["Best values = ", p0val[[optInd]]];
   {initialAngle, p0val[[optInd]], p0val}
   ]
  ];
  
  (*
myObjFunc = objHopf[RHS, var, par, coP, {3, 4}];
 myObjFunc[{0.6, 1.2}]

maxIt = 500;
accG = 2;
precG = 2;
tiLim = 120;
meth = "NelderMead";
{bestAngle, bestValues, finalP0Val} = 
optHopf[RHS, var, par, coP, {3, 4}, 120, "NelderMead", accG, precG, maxIt];

result1 = optHopf[RHS, var, par, coP, {3, 4}, 60, "DifferentialEvolution", 4, 4, 200];
result2 = optHopf[RHS, var, par, coP, {3, 4}, 90, "SimulatedAnnealing", 4, 4, 300];
*)
(*===PERTURBED INITIAL CONDITIONS===*)
pertIC[equilibrium_, var_, factor_: 0.1, minq_: 0.001, n_: 1] :=
   Module[{eqVals, pertVals},
    eqVals = var /. equilibrium;
    Table[
      pertVals = Max[#, minq] & /@ (eqVals + factor*RandomReal[{-1, 1},
Length[var]]*Abs[eqVals]);
      pertVals = pertVals(Total[eqVals]/Total[pertVals]);
      Thread[var -> pertVals],{n}]];

TS[RHS_, var_, par_, p0val_, tmax_: 100] := Module[
   {varT, initialConditions, rhsT, odeSystem, sol, solFun, coP, t},

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
   sol = Quiet@NDSolve[odeSystem, varT, {t, 0, tmax}, Method->{"BDF"}];

   If[!MatchQ[sol, {(_Rule | _List)...}],
     Print["NDSolve failed"];
     Return[$Failed];
   ];

   (* Extract solution functions *)
   solFun = #[[0]] & /@ (varT /. sol[[1]]);

   (* Return full solution functions *)
   solFun
];
  
scanPR[RHS_, var_, par_, coP_, plotInd_, 
  R01_: Automatic, R02_: Automatic, R21_: Automatic, R12_: Automatic,
  rangeB1_: Automatic, rangeB2_: Automatic, stepSize_: 0.1] := 
 Module[{par1, par2, pl1, pl2},
  
  (* Debug: show what we're working with *)
  Print["DEBUG: par = ", par];
  Print["DEBUG: coP = ", coP];
  Print["DEBUG: plotInd = ", plotInd];
  
  (* Extract the parameter symbols *)
  par1 = par[[plotInd[[1]]]];
  par2 = par[[plotInd[[2]]]];
  Print["DEBUG: par1 = ", par1, ", par2 = ", par2];
  
  (* Try to get numeric values *)
  pl1 = par1 /. coP;
  pl2 = par2 /. coP;
  Print["DEBUG: pl1 = ", pl1, ", pl2 = ", pl2];
  
  (* Check if extraction worked *)
  If[!NumericQ[pl1] || !NumericQ[pl2],
    Print["ERROR: Parameter extraction failed"];
    Print["pl1 = ", pl1, " (should be numeric)"];
    Print["pl2 = ", pl2, " (should be numeric)"];
    Return[$Failed];
  ];
  
  Print["SUCCESS: Extracted pl1 = ", pl1, ", pl2 = ", pl2];
  
  (* Simple test plot *)
  Graphics[{
    Red, PointSize[0.02], Point[{{pl1, pl2}}],
    Text["Test Point", {pl1, pl2}, {0, 1.5}]
  }, 
  Frame -> True,
  FrameLabel -> {ToString[par1], ToString[par2]},
  PlotRange -> {{0.5*pl1, 1.5*pl1}, {0.5*pl2, 1.5*pl2}}]
]
(*Fast parameter scanning using claP logic*)scanPar[RHS_,var_,par_,p0val_,plotInd_,gridRes_:Automatic,plot_:Automatic,hTol_:0.01,delta_:1/20,wRan_:1,hRan_:1,R01_:Automatic,R02_:Automatic,R21_:Automatic,R12_:Automatic]:=Module[{(*Parameter ranges and scanning*)\[Beta]1,\[Beta]2,\[Beta]1min,\[Beta]1max,\[Beta]2min,\[Beta]2max,stepSize,\[Beta]1vals,\[Beta]2vals,totalPoints,(*Results collection*)res,outcomes,outcomeCounts,(*Plotting and visualization*)finalPlot,plotData,betterColors,simpleTooltips,(*Mode detection*)useGridMode,\[Beta]1Index,\[Beta]2Index,\[Beta]1Center,\[Beta]2Center,(*Progress tracking*)progressVar,currentProgress,(*R curve variables*)R01F,R02F,R21F,R12F,rCurves,coPForR,(*Fast classification variables*)numpar,conPar,inP,sol,I1val,I2val,tol,currentType,errors,gridResults,nB1,nB2,jac,eigs,maxRealPart},(*Define better colors with old EE colors*)betterColors={RGBColor[0,0,1],(*Pure Blue for DFE*)RGBColor[0,1,0],(*Pure Green for E1*)RGBColor[0.8,0.5,0.2],(*Orange for E2*)RGBColor[0.6,0.3,0.8],(*Purple for EEstable*)RGBColor[1,0,0],(*Pure Red for EEunstable*)RGBColor[1.0,0.0,1.0](*Magenta for errors/NoSol*)};
(*Detect which mode we're using*)useGridMode=(gridRes=!=Automatic);
(*Extract parameter indices*)\[Beta]1Index=If[Length[plotInd]>=1,plotInd[[1]],1];
\[Beta]2Index=If[Length[plotInd]>=2,plotInd[[2]],2];
\[Beta]1Center=p0val[[\[Beta]1Index]];
\[Beta]2Center=p0val[[\[Beta]2Index]];
Print["Varying parameters at indices ",{\[Beta]1Index,\[Beta]2Index}," with center values: ",{\[Beta]1Center,\[Beta]2Center}];
If[useGridMode,(*GRID MODE*)Print["Using grid mode with ",gridRes,"\[Times]",gridRes," points..."];
\[Beta]1min=Max[\[Beta]1Center*(1-wRan),0.001];
\[Beta]1max=\[Beta]1Center*(1+wRan);
\[Beta]2min=Max[\[Beta]2Center*(1-hRan),0.001];
\[Beta]2max=\[Beta]2Center*(1+hRan);
\[Beta]1vals=Table[\[Beta]1min+(\[Beta]1max-\[Beta]1min)*k/(gridRes-1),{k,0,gridRes-1}];
\[Beta]2vals=Table[\[Beta]2min+(\[Beta]2max-\[Beta]2min)*k/(gridRes-1),{k,0,gridRes-1}];
stepSize="grid";,(*RANGE MODE*)Print["Using range mode with step size ",delta,"..."];
\[Beta]1min=Max[\[Beta]1Center*(1-wRan),0.001];
\[Beta]1max=\[Beta]1Center*(1+wRan);
\[Beta]2min=Max[\[Beta]2Center*(1-hRan),0.001];
\[Beta]2max=\[Beta]2Center*(1+hRan);
\[Beta]1vals=Table[\[Beta]1,{\[Beta]1,\[Beta]1min,\[Beta]1max,delta}];
\[Beta]2vals=Table[\[Beta]2,{\[Beta]2,\[Beta]2min,\[Beta]2max,delta}];
stepSize=delta;];
(*R curve preparation*)rCurves={};
If[R01=!=Automatic||R02=!=Automatic||R21=!=Automatic||R12=!=Automatic,coPForR=Thread[par->p0val];
coPForR=DeleteCases[coPForR,HoldPattern[par[[\[Beta]1Index]]->_]|HoldPattern[par[[\[Beta]2Index]]->_]];
If[R01=!=Automatic,R01F=R01/. coPForR;
AppendTo[rCurves,ContourPlot[R01F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Red,Thick,Opacity[0.8]],PlotPoints->30]]];
If[R02=!=Automatic,R02F=R02/. coPForR;
AppendTo[rCurves,ContourPlot[R02F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Blue,Thick,Opacity[0.8]],PlotPoints->30]]];
If[R21=!=Automatic,R21F=R21/. coPForR;
AppendTo[rCurves,ContourPlot[R21F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Green,Thick,Opacity[0.8]],PlotPoints->30]]];
If[R12=!=Automatic,R12F=R12/. coPForR;
AppendTo[rCurves,ContourPlot[R12F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Orange,Thick,Opacity[0.8]],PlotPoints->30]]];];
totalPoints=Length[\[Beta]1vals]*Length[\[Beta]2vals];
Print["Total points to scan: ",totalPoints];
(*Initialize for fast scanning*)res={};
currentProgress=0;
progressVar=0;
tol=10^(-6);
inP=Table[{var[[j]],1/Length[var]},{j,Length[var]}];
Print[ProgressIndicator[Dynamic[progressVar]]];
(*Fast scanning loop using claP logic*)Do[Do[currentProgress++;
progressVar=N[currentProgress/totalPoints];
(*Create parameter values for this grid point*)numpar=p0val;
numpar[[\[Beta]1Index]]=\[Beta]1;
numpar[[\[Beta]2Index]]=\[Beta]2;
conPar=Thread[par->numpar];
(*Single FindRoot call with equal initial conditions*)sol=Quiet[FindRoot[RHS/. conPar,inP]];
If[Head[sol]===List,(*Extract I1,I2 and classify*){I1val,I2val}={I1,I2}/. sol;
currentType=Which[I1val<tol&&I2val<tol,"DFE",I1val>=tol&&I2val<tol,"E1",I1val<tol&&I2val>=tol,"E2",True,(*Both strains present-check EE stability*)Module[{jac,eigs,maxRealPart},jac=D[RHS,{var}]/. conPar/. sol;
eigs=Chop[Eigenvalues[jac]//N];
maxRealPart=Max[Re[eigs]];
If[maxRealPart<0,"EEstable","EEunstable"]]];
res=Append[res,{N[\[Beta]1],N[\[Beta]2],currentType}];,(*FindRoot failed*)res=Append[res,{N[\[Beta]1],N[\[Beta]2],"NoSol"}];];,{\[Beta]2,\[Beta]2vals}];,{\[Beta]1,\[Beta]1vals}];
Print["Scanning complete!"];
(*Simple error detection for grid mode*)errors={};
If[useGridMode&&gridRes>2,gridResults=Table[Table[Null,{j,1,gridRes}],{i,1,gridRes}];
Do[Do[pos=Position[res,{\[Beta]1vals[[i]],\[Beta]2vals[[j]],_}];
If[Length[pos]>0,gridResults[[i,j]]=res[[pos[[1,1]],3]]];,{j,1,gridRes}];,{i,1,gridRes}];
Do[Do[If[gridResults[[i,j]]=!=Null,neighbors={};
If[i>1,AppendTo[neighbors,gridResults[[i-1,j]]]];
If[i<gridRes,AppendTo[neighbors,gridResults[[i+1,j]]]];
If[j>1,AppendTo[neighbors,gridResults[[i,j-1]]]];
If[j<gridRes,AppendTo[neighbors,gridResults[[i,j+1]]]];
neighbors=DeleteCases[neighbors,Null];
If[Length[neighbors]>=3&&Length[DeleteDuplicates[Join[{gridResults[[i,j]]},neighbors]]]>=4,AppendTo[errors,{\[Beta]1vals[[i]],\[Beta]2vals[[j]],gridResults[[i,j]]}];
pos=Position[res,{\[Beta]1vals[[i]],\[Beta]2vals[[j]],gridResults[[i,j]]}];
If[Length[pos]>0,res[[pos[[1,1]],3]]="error"];];];,{j,2,gridRes-1}];,{i,2,gridRes-1}];];
(*Process results*)outcomes=DeleteDuplicates[Table[res[[i,3]],{i,1,Length[res]}]];
outcomeCounts=Table[Count[res,{_,_,outcomes[[i]]}],{i,1,Length[outcomes]}];
Print["Found equilibrium types: ",Table[outcomes[[i]]<>" ("<>ToString[outcomeCounts[[i]]]<>" points)",{i,1,Length[outcomes]}]];
Print["Including ",Length[errors]," errors"];
(*Create plot*)simpleTooltips=Table[Table[If[res[[j,3]]==outcomes[[i]],Tooltip[res[[j,1;;2]],outcomes[[i]]<>": \[Beta]1="<>ToString[Round[res[[j,1]],0.01]]<>", \[Beta]2="<>ToString[Round[res[[j,2]],0.01]]],Nothing],{j,1,Length[res]}],{i,1,Length[outcomes]}];
finalPlot=ListPlot[simpleTooltips,GridLines->Automatic,PlotMarkers->Table[{Style["\[FilledCircle]",betterColors[[i]]],14},{i,1,Min[Length[betterColors],Length[outcomes]]}],PlotLegends->Table[outcomes[[i]]<>" ("<>ToString[outcomeCounts[[i]]]<>")",{i,1,Length[outcomes]}],AspectRatio->1,PlotLabel->"Parameter Space Analysis ("<>ToString[Length[res]]<>" points)",FrameLabel->{"Parameter "<>ToString[\[Beta]1Index],"Parameter "<>ToString[\[Beta]2Index]},LabelStyle->{FontSize->12,Black},ImageSize->450,PlotRange->{{\[Beta]1min,\[Beta]1max},{\[Beta]2min,\[Beta]2max}}];
(*Add R curves and overlays*)If[plot=!=Automatic,finalPlot=Show[plot,finalPlot,Sequence@@rCurves,Graphics[{EdgeForm[{Thick,Gray,Dashed}],FaceForm[None],Rectangle[{\[Beta]1min,\[Beta]2min},{\[Beta]1max,\[Beta]2max}],Blue,PointSize[0.015],Point[{{\[Beta]1Center,\[Beta]2Center}}]}],PlotRange->{{\[Beta]1min*0.9,\[Beta]1max*1.1},{\[Beta]2min*0.9,\[Beta]2max*1.1}}];,If[Length[rCurves]>0,finalPlot=Show[finalPlot,Sequence@@rCurves];];];
Print[Style["Summary:",FontWeight->Bold]];
Print[Grid[Table[{outcomes[[i]],outcomeCounts[[i]],ToString[Round[100.*outcomeCounts[[i]]/Length[res],1]]<>"%"},{i,1,Length[outcomes]}],Frame->All,Alignment->Left]];
{finalPlot,errors,res}];

(*Grid mode:scanPar[RHS,var,par,p0val,{1,2},30,plot,0.01,1/20,0.5,0.5,R01,R02,R21,R12]*)
(*Range mode:scanPar[RHS,var,par,p0val,{1,2},Automatic,plot,0.01,0.001,0.5,0.5,(*Fast parameter scanning using claP logic*)scanPar[RHS_,var_,par_,p0val_,plotInd_,gridRes_:Automatic,plot_:Automatic,hTol_:0.01,delta_:1/20,wRan_:1,hRan_:1,R01_:Automatic,R02_:Automatic,R21_:Automatic,R12_:Automatic]:=Module[{(*Parameter ranges and scanning*)\[Beta]1,\[Beta]2,\[Beta]1min,\[Beta]1max,\[Beta]2min,\[Beta]2max,stepSize,\[Beta]1vals,\[Beta]2vals,totalPoints,(*Results collection*)res,outcomes,outcomeCounts,(*Plotting and visualization*)finalPlot,plotData,betterColors,simpleTooltips,(*Mode detection*)useGridMode,\[Beta]1Index,\[Beta]2Index,\[Beta]1Center,\[Beta]2Center,(*Progress tracking*)progressVar,currentProgress,(*R curve variables*)R01F,R02F,R21F,R12F,rCurves,coPForR,(*Fast classification variables*)numpar,conPar,inP,sol,I1val,I2val,tol,currentType,errors,gridResults,nB1,nB2,jac,eigs,maxRealPart},(*Define better colors with old EE colors*)betterColors={RGBColor[0,0,1],(*Pure Blue for DFE*)RGBColor[0,1,0],(*Pure Green for E1*)RGBColor[0.8,0.5,0.2],(*Orange for E2*)RGBColor[0.6,0.3,0.8],(*Purple for EEstable*)RGBColor[1,0,0],(*Pure Red for EEunstable*)RGBColor[1.0,0.0,1.0](*Magenta for errors/NoSol*)};
(*Detect which mode we're using*)useGridMode=(gridRes=!=Automatic);
(*Extract parameter indices*)\[Beta]1Index=If[Length[plotInd]>=1,plotInd[[1]],1];
\[Beta]2Index=If[Length[plotInd]>=2,plotInd[[2]],2];
\[Beta]1Center=p0val[[\[Beta]1Index]];
\[Beta]2Center=p0val[[\[Beta]2Index]];
Print["Varying parameters at indices ",{\[Beta]1Index,\[Beta]2Index}," with center values: ",{\[Beta]1Center,\[Beta]2Center}];
If[useGridMode,(*GRID MODE*)Print["Using grid mode with ",gridRes,"\[Times]",gridRes," points..."];
\[Beta]1min=Max[\[Beta]1Center*(1-wRan),0.001];
\[Beta]1max=\[Beta]1Center*(1+wRan);
\[Beta]2min=Max[\[Beta]2Center*(1-hRan),0.001];
\[Beta]2max=\[Beta]2Center*(1+hRan);
\[Beta]1vals=Table[\[Beta]1min+(\[Beta]1max-\[Beta]1min)*k/(gridRes-1),{k,0,gridRes-1}];
\[Beta]2vals=Table[\[Beta]2min+(\[Beta]2max-\[Beta]2min)*k/(gridRes-1),{k,0,gridRes-1}];
stepSize="grid";,(*RANGE MODE*)Print["Using range mode with step size ",delta,"..."];
\[Beta]1min=Max[\[Beta]1Center*(1-wRan),0.001];
\[Beta]1max=\[Beta]1Center*(1+wRan);
\[Beta]2min=Max[\[Beta]2Center*(1-hRan),0.001];
\[Beta]2max=\[Beta]2Center*(1+hRan);
\[Beta]1vals=Table[\[Beta]1,{\[Beta]1,\[Beta]1min,\[Beta]1max,delta}];
\[Beta]2vals=Table[\[Beta]2,{\[Beta]2,\[Beta]2min,\[Beta]2max,delta}];
stepSize=delta;];
(*R curve preparation*)rCurves={};
If[R01=!=Automatic||R02=!=Automatic||R21=!=Automatic||R12=!=Automatic,coPForR=Thread[par->p0val];
coPForR=DeleteCases[coPForR,HoldPattern[par[[\[Beta]1Index]]->_]|HoldPattern[par[[\[Beta]2Index]]->_]];
If[R01=!=Automatic,R01F=R01/. coPForR;
AppendTo[rCurves,ContourPlot[R01F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Red,Thick,Opacity[0.8]],PlotPoints->30]]];
If[R02=!=Automatic,R02F=R02/. coPForR;
AppendTo[rCurves,ContourPlot[R02F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Blue,Thick,Opacity[0.8]],PlotPoints->30]]];
If[R21=!=Automatic,R21F=R21/. coPForR;
AppendTo[rCurves,ContourPlot[R21F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Green,Thick,Opacity[0.8]],PlotPoints->30]]];
If[R12=!=Automatic,R12F=R12/. coPForR;
AppendTo[rCurves,ContourPlot[R12F==1,{par[[\[Beta]1Index]],\[Beta]1min,\[Beta]1max},{par[[\[Beta]2Index]],\[Beta]2min,\[Beta]2max},ContourStyle->Directive[Orange,Thick,Opacity[0.8]],PlotPoints->30]]];];
totalPoints=Length[\[Beta]1vals]*Length[\[Beta]2vals];
Print["Total points to scan: ",totalPoints];
(*Initialize for fast scanning*)res={};
currentProgress=0;
progressVar=0;
tol=10^(-6);
inP=Table[{var[[j]],1/Length[var]},{j,Length[var]}];
Print[ProgressIndicator[Dynamic[progressVar]]];
(*MAIN SCANNING LOOP*)Do[Do[currentProgress++;
progressVar=N[currentProgress/totalPoints];
(*Create parameter values for this grid point*)numpar=p0val;
numpar[[\[Beta]1Index]]=\[Beta]1;
numpar[[\[Beta]2Index]]=\[Beta]2;
conPar=Thread[par->numpar];
(*Single FindRoot call with equal initial conditions*)sol=Quiet[FindRoot[RHS/. conPar,inP]];
If[Head[sol]===List,(*Extract I1,I2 and classify*){I1val,I2val}={I1,I2}/. sol;
currentType=Which[I1val<tol&&I2val<tol,"DFE",I1val>=tol&&I2val<tol,"E1",I1val<tol&&I2val>=tol,"E2",True,(*Both strains present-check EE stability*)(jac=D[RHS,{var}]/. conPar/. sol;
eigs=Chop[Eigenvalues[jac]//N];
maxRealPart=Max[Re[eigs]];
If[maxRealPart<0,"EEstable","EEunstable"])];
res=Append[res,{N[\[Beta]1],N[\[Beta]2],currentType}];,(*FindRoot failed*)res=Append[res,{N[\[Beta]1],N[\[Beta]2],"NoSol"}];];,{\[Beta]2,\[Beta]2vals}];,{\[Beta]1,\[Beta]1vals}];
Print["Scanning complete!"];
(*Simple error detection for grid mode*)errors={};
If[useGridMode&&gridRes>2,gridResults=Table[Table[Null,{j,1,gridRes}],{i,1,gridRes}];
Do[Do[pos=Position[res,{\[Beta]1vals[[i]],\[Beta]2vals[[j]],_}];
If[Length[pos]>0,gridResults[[i,j]]=res[[pos[[1,1]],3]]];,{j,1,gridRes}];,{i,1,gridRes}];
Do[Do[If[gridResults[[i,j]]=!=Null,neighbors={};
If[i>1,AppendTo[neighbors,gridResults[[i-1,j]]]];
If[i<gridRes,AppendTo[neighbors,gridResults[[i+1,j]]]];
If[j>1,AppendTo[neighbors,gridResults[[i,j-1]]]];
If[j<gridRes,AppendTo[neighbors,gridResults[[i,j+1]]]];
neighbors=DeleteCases[neighbors,Null];
If[Length[neighbors]>=3&&Length[DeleteDuplicates[Join[{gridResults[[i,j]]},neighbors]]]>=4,AppendTo[errors,{\[Beta]1vals[[i]],\[Beta]2vals[[j]],gridResults[[i,j]]}];
pos=Position[res,{\[Beta]1vals[[i]],\[Beta]2vals[[j]],gridResults[[i,j]]}];
If[Length[pos]>0,res[[pos[[1,1]],3]]="error"];];];,{j,2,gridRes-1}];,{i,2,gridRes-1}];];
(*Process results*)outcomes=DeleteDuplicates[Table[res[[i,3]],{i,1,Length[res]}]];
outcomeCounts=Table[Count[res,{_,_,outcomes[[i]]}],{i,1,Length[outcomes]}];
Print["Found equilibrium types: ",Table[outcomes[[i]]<>" ("<>ToString[outcomeCounts[[i]]]<>" points)",{i,1,Length[outcomes]}]];
Print["Including ",Length[errors]," errors"];
(*Create plot*)simpleTooltips=Table[Table[If[res[[j,3]]==outcomes[[i]],Tooltip[res[[j,1;;2]],outcomes[[i]]<>": \[Beta]1="<>ToString[Round[res[[j,1]],0.01]]<>", \[Beta]2="<>ToString[Round[res[[j,2]],0.01]]],Nothing],{j,1,Length[res]}],{i,1,Length[outcomes]}];
finalPlot=ListPlot[simpleTooltips,GridLines->Automatic,PlotMarkers->Table[{Style["\[FilledCircle]",betterColors[[i]]],14},{i,1,Min[Length[betterColors],Length[outcomes]]}],PlotLegends->Table[outcomes[[i]]<>" ("<>ToString[outcomeCounts[[i]]]<>")",{i,1,Length[outcomes]}],AspectRatio->1,PlotLabel->"Parameter Space Analysis ("<>ToString[Length[res]]<>" points)",FrameLabel->{"Parameter "<>ToString[\[Beta]1Index],"Parameter "<>ToString[\[Beta]2Index]},LabelStyle->{FontSize->12,Black},ImageSize->450,PlotRange->{{\[Beta]1min,\[Beta]1max},{\[Beta]2min,\[Beta]2max}}];
(*Add R curves and overlays*)If[plot=!=Automatic,finalPlot=Show[plot,finalPlot,Sequence@@rCurves,Graphics[{EdgeForm[{Thick,Gray,Dashed}],FaceForm[None],Rectangle[{\[Beta]1min,\[Beta]2min},{\[Beta]1max,\[Beta]2max}],Blue,PointSize[0.015],Point[{{\[Beta]1Center,\[Beta]2Center}}]}],PlotRange->{{\[Beta]1min*0.9,\[Beta]1max*1.1},{\[Beta]2min*0.9,\[Beta]2max*1.1}}];,If[Length[rCurves]>0,finalPlot=Show[finalPlot,Sequence@@rCurves];];];
Print[Style["Summary:",FontWeight->Bold]];
Print[Grid[Table[{outcomes[[i]],outcomeCounts[[i]],ToString[Round[100.*outcomeCounts[[i]]/Length[res],1]]<>"%"},{i,1,Length[outcomes]}],Frame->All,Alignment->Left]];
{finalPlot,errors,res}];

(*Grid mode:scanPar[RHS,var,par,p0val,{1,2},30,plot,0.01,1/20,1.0,1.0,R01,R02,R21,R12]*)
(*Range mode:scanPar[RHS,var,par,p0val,{1,2},Automatic,plot,0.01,0.05,1.0,1.0,R01,R02,R21,R12]*)R01,R02,R21,R12]*)
intEq[RHS_, var_, par_, p0val_, att_] := 
 Module[{eqs, sols, cleanSols, positiveSols, threshold = 10^(-10), 
   perturbedICs, attRules, coP}, 
  
  (* Create substitution rules *)
  coP = Thread[par -> p0val];
  eqs = Thread[RHS == 0] /. coP;
  
  (* Convert att list to replacement rules *)
  attRules = Thread[var -> att];
  
  (* Use attRules as seed instead of calling TS *)
  perturbedICs = Join[
   {var /. attRules}, 
   var /. # & /@ pertIC[attRules, var, 0.1, 15], 
   var /. # & /@ pertIC[attRules, var, 0.3, 8]   
   ];
   
  sols = Table[
   Quiet[FindRoot[eqs, 
    Evaluate[Thread[{var, perturbedICs[[i]]}]], 
    Method -> "Newton"]], {i, 1, Length[perturbedICs]}];
  
  sols = Select[sols, Head[#] === List && Length[#] == Length[var] &];
  cleanSols = sols /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
  cleanSols = DeleteDuplicatesBy[cleanSols, Round[Sort[var /. #], 10^(-6)] &];
  positiveSols = Select[cleanSols, And @@ (# > threshold & /@ (var /. #)) &];
  
  Print["Found ", Length[positiveSols], " equilibria"];
  positiveSols
  ];
(* Enhanced Continuation function with explicit parameter handling *)
cont[RHS_, var_, par_, p0val_, stepSize_,  plotInd_: {3,4}, bifInd_: 4, 
analyticalPlot_: None] := 
 Module[{posSols, bifParamIndex, currentVal, curve, newP0Val, paramVal, 
   complexEigs, eigs, eq, maxRealPart, complexPairs, bifPlot, plotC,
   forwardCurve, backwardCurve, i, maxIter = 10, successCount = 0, hopfPoints,
   testedPoints = {}, startingPoint, adaptiveStep, currentStep, overlayPossible,
   beta1Index, beta2Index, beta1Val, beta2Val, startPosInCurve},
  
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
   Print["\:2705 Will overlay on analytical plot (bifParam \[Element] plotInd)"];,
   If[analyticalPlot =!= None,
    Print["\:26a0\:fe0f Standalone plot only (bifParam \[NotElement] plotInd)"];
    ];
   ];
  
  (* Get starting equilibrium *)
  {posSols, complexEigs, eigs} = fpHopf[RHS, var, par, p0val];
  
  If[Length[posSols] == 0,
   Print["ERROR: No equilibrium found at starting p0val!"];
   Return[{}];
   ];
  
  (* Include starting equilibrium with \[Beta] coordinates *)
  {beta1Val, beta2Val} = p0val[[plotInd]];
  startingPoint = <|"Parameter" -> currentVal, "Values" -> posSols[[1]], 
    "MaxRealPart" -> If[Length[eigs] > 0, Max[Re[eigs]], 0], 
    "ComplexPairs" -> Length[complexEigs]/2, 
    "Eigenvalues" -> eigs, "ComplexEigs" -> complexEigs,
    "Beta1" -> beta1Val, "Beta2" -> beta2Val|>;
  
  curve = {};
  forwardCurve = {startingPoint};
  backwardCurve = {};
  
  (* Adaptive step size function *)
  adaptiveStep[stepNum_, baseStep_] := Which[
    stepNum <= 2, baseStep/5,
    stepNum <= 4, baseStep/3,
    stepNum <= 6, baseStep/2,
    True, baseStep
    ];
  
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
     Print["*** HOPF BIFURCATION FOUND at parameter = ", paramVal, " ***"];
     ];
    
    AppendTo[forwardCurve, <|"Parameter" -> paramVal, "Values" -> eq, 
      "MaxRealPart" -> maxRealPart, "ComplexPairs" -> complexPairs, 
      "Eigenvalues" -> eigs, "ComplexEigs" -> complexEigs,
      "Beta1" -> beta1Val, "Beta2" -> beta2Val|>];
    
    i++;
    If[i >= maxIter, Break[]];,
    Break[];
    ];
   ];
  
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
     Print["*** HOPF BIFURCATION FOUND at parameter = ", paramVal, " ***"];
     ];
    
    PrependTo[backwardCurve, <|"Parameter" -> paramVal, "Values" -> eq, 
      "MaxRealPart" -> maxRealPart, "ComplexPairs" -> complexPairs, 
      "Eigenvalues" -> eigs, "ComplexEigs" -> complexEigs,
      "Beta1" -> beta1Val, "Beta2" -> beta2Val|>];
    
    i++;
    If[i > maxIter, Break[]];,
    Break[];
    ];
   ];
  
  (* Combine curves *)
  curve = Join[backwardCurve, forwardCurve];
  
  Print["Found ", Length[curve], " points, tested ", Length[testedPoints], " parameters"];
  If[Length[curve] > 0,
   Print["Parameter range: ", curve[[1]]["Parameter"], " to ", curve[[-1]]["Parameter"]];
   ];
  
  (* Hopf analysis *)
  Print["\n=== HOPF BIFURCATION ANALYSIS ==="];
  hopfPoints = Select[curve, Length[#["ComplexEigs"]] > 0 && Abs[Max[Re[#["ComplexEigs"]]]] < 0.01 &];
  If[Length[hopfPoints] > 0,
   Print["HOPF BIFURCATIONS FOUND at parameters:"];
   Do[Print["  Parameter = ", pt["Parameter"]], {pt, hopfPoints}];,
   Print["NO HOPF BIFURCATIONS FOUND in parameter range"];
   ];
  
  (* Plotting *)
  If[Length[curve] > 1,
   (* S-shaped continuation curve *)
   curvePoints = {#["Parameter"], #["Values"][[1]]} & /@ curve;
   startPosInCurve = Length[backwardCurve] + 1;
   startPoint = curvePoints[[startPosInCurve]];
   
   bifPlot = ListLinePlot[
     curvePoints,
     PlotLegends -> {"S: " <> ToString[var[[1]]]},
     AxesLabel -> {"Parameter " <> ToString[bifParamIndex], ToString[var[[1]]]},
     PlotLabel -> "S-shaped Continuation Curve",
     ImageSize -> 600,
     Frame -> True,
     PlotStyle -> Blue,
     Epilog -> {
       Pink, PointSize[0.008],
       Point[{#, 0} & /@ testedPoints[[All, 2]]],
       Blue, PointSize[0.006],
       Point[curvePoints],
       Red, PointSize[0.012],
       Point[startPoint],
       Black, Arrowheads[0.015],
       Module[{leftArrow, rightArrow},
        leftArrow = If[startPosInCurve > 1, 
          Arrow[{startPoint, curvePoints[[startPosInCurve - 1]]}], {}];
        rightArrow = If[startPosInCurve < Length[curvePoints], 
          Arrow[{startPoint, curvePoints[[startPosInCurve + 1]]}], {}];
        {leftArrow, rightArrow}
        ],
       Text[Style["p\:2080", FontSize -> 10, FontWeight -> Bold, FontColor -> Red], 
            startPoint + {0.02, 0.05}]
       }
     ];
   Print[bifPlot];
   
   (* Overlay plot *)
   If[overlayPossible,
    Print["\|01f3a8 Creating overlay on analytical plot..."];
    
    overlayPoints = {#["Beta1"], #["Beta2"]} & /@ curve;
    overlayStart = overlayPoints[[startPosInCurve]];
    
    searchOverlayPoints = Table[
      Module[{testP0Val},
       testP0Val = p0val;
       testP0Val[[bifParamIndex]] = testedPoints[[j, 2]];
       testP0Val[[plotInd]]
       ],
      {j, 1, Length[testedPoints]}];
    
    plotC = Show[
     analyticalPlot,
     Graphics[{
       Thick, Purple,
       Line[overlayPoints],
       Red, PointSize[0.015],
       Point[overlayStart],
       Orange, PointSize[0.010],
       Point[searchOverlayPoints],
       Purple, PointSize[0.008],
       Point[overlayPoints],
       Black, Arrowheads[0.02],
       Module[{leftArrow, rightArrow},
        leftArrow = If[startPosInCurve > 1, 
          Arrow[{overlayStart, overlayPoints[[startPosInCurve - 1]]}], {}];
        rightArrow = If[startPosInCurve < Length[overlayPoints], 
          Arrow[{overlayStart, overlayPoints[[startPosInCurve + 1]]}], {}];
        {leftArrow, rightArrow}
        ],
       Text[Style["p\:2080", FontSize -> 10, FontWeight -> Bold, FontColor -> Red], 
            overlayStart + {0.02, 0.02}],
       Text[Style["Hopf Search Path", FontSize -> 10, FontWeight -> Bold, FontColor -> Purple], 
            {Mean[overlayPoints[[All, 1]]], Max[overlayPoints[[All, 2]]] + 0.05}]
       }],
     PlotLabel -> "Analytical Boundaries + Hopf Search Path"
     ];
    Print[plotC];
    ];
   ];
  
  curve
  ];

hopfD[curve_] := 
 Module[{hopfPoints, stabilityPlot, hopfPlot},
  If[Length[curve] < 3, 
   Print["hopfD needs \[GreaterEqual]3 points, got ", Length[curve]]; 
   Return[{}]];
  
  (* Check if curve has required fields *)
  If[!KeyExistsQ[curve[[1]], "MaxRealPart"],
   Print["Error: curve missing stability data. Run cont first."];
   Return[{}]];
  
  stabilityPlot = ListLinePlot[{#["Parameter"], #["MaxRealPart"]} & /@ curve,
    AxesLabel -> {"Parameter", "Max Re(\[Lambda])"},
    PlotLabel -> "Stability Analysis",
    PlotStyle -> Blue,
    GridLines -> {None, {0}},
    GridLinesStyle -> Directive[Red, Dashed],
    ImageSize -> 400];
  
  hopfPoints = Select[
    Table[
     Module[{prev, curr, next, prevReal, currReal, hasComplex},
      {prev, curr, next} = {curve[[i-1]], curve[[i]], If[i < Length[curve], curve[[i+1]], curve[[i]]]};
      {prevReal, currReal} = {prev["MaxRealPart"], curr["MaxRealPart"]};
      hasComplex = curr["ComplexPairs"] > 0;
      If[hasComplex && prevReal*currReal < 0,
       <|"parBif" -> curr["Parameter"], "Type" -> "Hopf"|>, Nothing]
      ], {i, 2, Length[curve] - 1}], # =!= Nothing &];
  
  If[Length[hopfPoints] > 0,
   hopfPlot = ListPlot[{#["parBif"], 0} & /@ hopfPoints, PlotStyle -> {Red, PointSize[0.015]}];
   Print[Show[stabilityPlot, hopfPlot, PlotLabel -> "Stability Analysis with Hopf Points"]];
   Print["hopfD found ", Length[hopfPoints], " Hopf points at parameters: ", #["parBif"] & /@ hopfPoints];,
   Print[stabilityPlot];
   Print["hopfD found no Hopf points"];];
  
  hopfPoints];

End[]
EndPackage[]
