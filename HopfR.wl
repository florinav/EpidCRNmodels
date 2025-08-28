(* ::Package:: *)

BeginPackage["HopfR`"]

(* Core Bifurcation Analysis Functions *)
intEq::usage = "equilibria = intEq[RHS, var, coP] finds interior equilibrium points - silent version."
intEqV::usage = "equilibria = intEqV[RHS, var, coP, E0, E1, E2] finds interior equilibrium points with verbose output and symbolic equilibria checking."
pertIC::usage = "initialConditions = pertIC[equilibrium, var, factor, n] creates n perturbed initial condition sets around equilibrium with perturbation factor."
cont::usage = "curve = cont[RHS, var, eqPoints, par, coP] performs continuation using first parameter and first equilibrium."
cont2::usage = "curve = cont2[RHS, var, eqPoints, par, coP, inp] performs continuation using parameter par[[inp]]."
cont3::usage = "curve = cont3[RHS, var, eqPoints, par, coP, inp, range] performs continuation with custom parameter range."
cont4::usage = "curve = cont4[RHS, var, eqPoints, par, coP, inp, range, ind] performs continuation with all options specified."
hopfD::usage = "hopfPoints = hopfD[curve] detects Hopf bifurcation points from continuation curve data."
TS::usage = "solution = TS[RHS, var, coP, R12, R21, tmax] performs time series analysis with oscillation detection."
phase::usage = "phasePlot = phase[RHS, var, coP, tmax] creates phase portrait with multiple trajectories around equilibrium."
pltBif::usage = "bifPlot = pltBif[curve, hopfPoints, paramName] creates enhanced bifurcation plot with Hopf points."
wF::usage = "wF[bifParam, coP, RHS, var, R12, R21] performs workflow analysis at bifurcation parameter."

plotBifurcationDiagram::usage = "plotBifurcationDiagram[curve_, hopfPoints_, bifParam_, opts___] creates a bifurcation diagram showing equilibria and Hopf bifurcations";
createFigure1::usage = "createFigure1[RHS_, var_, coP_, R12_, R21_, par_, opts___] reproduces Figure 1 style bifurcation analysis";
plotPhasePortrait::usage = "plotPhasePortrait[RHS_, var_, params_, xVar_, yVar_, opts___] creates phase portraits at specific parameter values";
quickBifurcationPlot::usage = "quickBifurcationPlot[curve_, hopfPoints_, opts___] creates a simple bifurcation plot from existing continuation data";




Begin["`Private`"]

(*===INTERIOR EQUILIBRIUM - SILENT VERSION===*)
intEq[RHS_, var_, coP_] := 
  Module[{eqs, sols, cleanSols, orthantSols, positiveSols, threshold = 10^(-10)}, 
   eqs = Thread[RHS == 0] /. coP;
   
   (* Improved precision: more attempts and multiple methods *)
   sols = Join[
     (* Try many random starting points with Newton *)
     Table[
       Quiet[FindRoot[eqs, Evaluate[Thread[{var, RandomReal[{0.01, 10}, Length[var]]}]], Method -> "Newton"]], 
       {40}],
     (* Also try with different method and range *)
     Table[
       Quiet[FindRoot[eqs, Evaluate[Thread[{var, RandomReal[{0.01, 5}, Length[var]]}]], Method -> "Broyden"]], 
       {20}]];
   
   sols = Select[sols, Head[#] === List && Length[#] == Length[var] &];
   cleanSols = sols /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
   cleanSols = DeleteDuplicatesBy[cleanSols, Round[Sort[var /. #], 10^(-8)] &];
   
   orthantSols = Select[cleanSols, And @@ (# >= 0 & /@ (var /. #)) &];
   positiveSols = Select[cleanSols, And @@ (# > threshold & /@ (var /. #)) &];
   
   (* PRINTS NOTHING, NEVER *)
   
   If[Length[positiveSols] > 0, positiveSols, orthantSols]];

(*===INTERIOR EQUILIBRIUM - VERBOSE VERSION===*)
intEqV[RHS_, var_, coP_, E0_, E1_, E2_] := 
  Module[{eqs, sols, cleanSols, orthantSols, positiveSols, boundaryFPs, threshold = 10^(-10), 
    symbolicEqs, identifiedFPs}, 
   eqs = Thread[RHS == 0] /. coP;
   
   (* Improved precision: more attempts, multiple methods, and seed with symbolic solutions *)
   symbolicEqs = {E0 /. coP, E1 /. coP, E2 /. coP};
   
   sols = Join[
     (* First try symbolic equilibria as starting points *)
     Table[
       Quiet[FindRoot[eqs, Evaluate[Thread[{var, var /. eq}]], Method -> "Newton"]], 
       {eq, symbolicEqs}],
     (* Then try many random starting points *)
     Table[
       Quiet[FindRoot[eqs, Evaluate[Thread[{var, RandomReal[{0.01, 10}, Length[var]]}]], Method -> "Newton"]], 
       {40}],
     (* Also try with different method *)
     Table[
       Quiet[FindRoot[eqs, Evaluate[Thread[{var, RandomReal[{0.01, 5}, Length[var]]}]], Method -> "Broyden"]], 
       {20}]];
   
   sols = Select[sols, Head[#] === List && Length[#] == Length[var] &];
   cleanSols = sols /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
   cleanSols = DeleteDuplicatesBy[cleanSols, Round[Sort[var /. #], 10^(-8)] &];
   
   orthantSols = Select[cleanSols, And @@ (# >= 0 & /@ (var /. #)) &];
   positiveSols = Select[cleanSols, And @@ (# > threshold & /@ (var /. #)) &];
   boundaryFPs = Select[orthantSols, !MemberQ[positiveSols, #] &];
   
   (* Identify which equilibrium corresponds to which symbolic one *)
   identifiedFPs = {};
   Module[{e0Num, e1Num, e2Num, tol = 0.01},
     e0Num = E0 /. coP; e1Num = E1 /. coP; e2Num = E2 /. coP;
     
     Do[
       Module[{fpVals = var /. fp},
         Which[
           Norm[fpVals - (var /. e0Num)] < tol, 
           AppendTo[identifiedFPs, {fp, "DFE (E0)"}],
           Norm[fpVals - (var /. e1Num)] < tol, 
           AppendTo[identifiedFPs, {fp, "E1 (strain 1 dominance)"}],
           Norm[fpVals - (var /. e2Num)] < tol, 
           AppendTo[identifiedFPs, {fp, "E2 (strain 2 dominance)"}],
           True, 
           AppendTo[identifiedFPs, {fp, "unknown"}]
         ]
       ], {fp, orthantSols}];
   ];
   
   (* EXTRA PRINTING LINE - interior and boundary fixed points separately, no total *)
   Print["  ", Length[positiveSols], " strictly positive fixed point: ", positiveSols];
   If[Length[positiveSols] > 0,
     Do[Module[{id = SelectFirst[identifiedFPs, #[[1]] == fp &]},
       If[id =!= Missing["NotFound"], Print["    -> ", id[[2]]]]
     ], {fp, positiveSols}]];
   
   Print["  ", Length[boundaryFPs], " boundary fixed point: ", boundaryFPs];
   If[Length[boundaryFPs] > 0,
     Do[Module[{id = SelectFirst[identifiedFPs, #[[1]] == fp &]},
       If[id =!= Missing["NotFound"], Print["    -> ", id[[2]]]]
     ], {fp, boundaryFPs}]];
   
   If[Length[orthantSols] < 3,
     Print["WARNING: Found only ", Length[orthantSols], " equilibria instead of expected 3"]];
   
   If[Length[positiveSols] > 0, positiveSols, orthantSols]];

(*===PERTURBED INITIAL CONDITIONS===*)
pertIC[equilibrium_, var_, factor_: 0.1, n_: 1] := 
  Module[{eqVals, pertVals}, 
   eqVals = var /. equilibrium;
   Table[
     pertVals = Max[#, 0.001] & /@ (eqVals + factor*RandomReal[{-1, 1}, Length[var]]*Abs[eqVals]);
     Thread[Through[var[0]] == pertVals], 
     {n}]];

(*===CONTINUATION FUNCTIONS===*)

(*cont: Basic continuation with first parameter and first equilibrium*)
cont[RHS_, var_, eqPoints_, par_, coP_] := 
  Module[{firstParamValue, paramRange},
   firstParamValue = par[[1]] /. coP;
   paramRange = If[NumericQ[firstParamValue] && firstParamValue > 1, 
                   {0.5*firstParamValue, 2*firstParamValue}, 
                   {0.1, 2}];
   Print["Using parameter range: ", paramRange, " for ", par[[1]]];
   cont4[RHS, var, eqPoints, par, coP, 1, paramRange, 1]
  ];

(*cont2: Continuation with specified parameter index*)
cont2[RHS_, var_, eqPoints_, par_, coP_, inp_] := 
  cont4[RHS, var, eqPoints, par, coP, inp, {0.1, 0.8}, 1];

(*cont3: Continuation with specified parameter index and range*)
cont3[RHS_, var_, eqPoints_, par_, coP_, inp_, range_] := 
  cont4[RHS, var, eqPoints, par, coP, inp, range, 1];

(*cont4: Full continuation with all options*)
cont4[RHS_, var_, eqPoints_, par_, coP_, inp_, range_, ind_] := 
  Module[{param, jac, paramVals, eqVals, eigVals, plot, testRHS, testJac}, 
   
   If[Length[eqPoints] < ind, 
     Print["Error: only ", Length[eqPoints], " equilibria available, need index ", ind];
     Return[{}]];
   
   param = par[[inp]];
   Print["Continuing with parameter: ", param, " from ", range[[1]], " to ", range[[2]]];
   Print["Starting equilibrium: ", eqPoints[[ind]]];
   
   testRHS = Quiet[N[RHS /. coP /. eqPoints[[ind]]]];
   If[!VectorQ[testRHS, NumericQ],
     Print["ERROR: RHS doesn't evaluate to numbers at starting point"];
     Print["RHS result: ", testRHS];
     Return[{}]];
   
   jac = D[RHS, {var}];
   
   testJac = Quiet[N[jac /. coP /. eqPoints[[ind]]]];
   If[!MatrixQ[testJac, NumericQ],
     Print["WARNING: Jacobian doesn't evaluate properly at starting point"];
     Print["Jacobian result: ", testJac]];
   
   paramVals = Range[range[[1]], range[[2]], (range[[2]] - range[[1]])/20];
   
   eqVals = Table[
     Module[{newParams, startVals, sol, rhsNumerical},
       newParams = coP /. (param -> _) -> (param -> pval);
       startVals = var /. eqPoints[[ind]];
       
       rhsNumerical = Quiet[N[RHS /. newParams /. Thread[var -> startVals]]];
       If[!VectorQ[rhsNumerical, NumericQ],
         Nothing, 
         sol = Quiet[FindRoot[Thread[RHS == 0] /. newParams, 
           Thread[{var, startVals}], 
           Method -> "Newton",
           MaxIterations -> 100]];
         If[Head[sol] === List && Length[sol] == Length[var] && 
            VectorQ[N[var /. sol], NumericQ], sol, Nothing]
       ]
     ], {pval, paramVals}];
   
   If[Length[eqVals] == 0, 
     Print["No equilibria found in continuation"];
     Print["Try different parameter (inp=", inp, ") or range ", range];
     Print["Available parameters: ", par];
     Return[{}]];
   
   Print["Found ", Length[eqVals], " equilibria along continuation curve"];
   
   eigVals = MapThread[
     Function[{eq, pval},
       Module[{newParams, jacEval, eigs},
         newParams = coP /. (param -> _) -> (param -> pval);
         jacEval = Quiet[N[jac /. eq /. newParams]];
         If[MatrixQ[jacEval, NumericQ],
           eigs = Quiet[Select[Eigenvalues[jacEval], NumericQ]];
           <|"Parameter" -> pval, "Equilibrium" -> eq, "Eigenvalues" -> eigs, 
             "MaxRealPart" -> If[Length[eigs] > 0, Max[Re[eigs]], Undefined], 
             "ComplexPairs" -> Count[eigs, _?((Im[#] != 0) &)]|>,
           Nothing
         ]
       ]], {eqVals, paramVals[[1 ;; Length[eqVals]]]}];
   
   eigVals = Select[eigVals, #["MaxRealPart"] =!= Undefined &];
   
   If[Length[eigVals] > 0,
     plot = ListLinePlot[{#["Parameter"], var[[1]] /. #["Equilibrium"]} & /@ eigVals, 
       AxesLabel -> {ToString[param], ToString[var[[1]]]}, 
       PlotLabel -> "Continuation Curve",
       ImageSize -> 500,
       PlotStyle -> Blue];
     Print[plot];,
     Print["No valid eigenvalue data for plotting"];];
   
   eigVals];

(*===HOPF DETECTION===*)
hopfD[curve_] := 
  Module[{hopfPoints, stabilityPlot, hopfPlot}, 
   If[Length[curve] < 3, 
     Print["hopfD needs \[GreaterEqual]3 points, got ", Length[curve]]; 
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

(*===TIME SERIES ANALYSIS===*)
TS[RHS_, var_, coP_, R12_, R21_, tmax_ : 50] := 
  Module[{intEqs, equilibrium, varT, odeSystem, sol, timePlot, 
    initialConds, infPos, plotVars, eqInfected, rhsT, plotTime}, 
   
   Print["R\:2081\:2082 = ", N[R12 /. coP, 4], ", R\:2082\:2081 = ", N[R21 /. coP, 4]];
   
   intEqs = intEq[RHS, var, coP];
   If[Length[intEqs] == 0, Return[$Failed]];
   equilibrium = intEqs[[1]];
   
   (* Create time-dependent variables *)
   varT = Through[var[t]];
   
   (* Get initial conditions *)
   initialConds = pertIC[equilibrium, var, 0.1, 1][[1]];
   
   (* Substitute parameters first, then variables *)
   rhsT = RHS /. coP /. Thread[var -> varT];
   
   (* Build ODE system carefully *)
   odeSystem = Join[
     Thread[D[varT, t] == rhsT],
     initialConds /. Thread[Through[var[0]] -> (varT /. t -> 0)]];
   
   sol = Quiet[NDSolve[odeSystem, varT, {t, 0, tmax}]];
   
   If[Length[sol] > 0 && Head[sol[[1]]] === List,
     infPos = If[Length[var] >= 5, {2, 5}, {2, Min[3, Length[var]]}];
     plotVars = varT[[infPos]];
     eqInfected = var[[infPos]] /. equilibrium;
     
     (* Plot only for first 1/5 of time horizon *)
     plotTime = tmax/25;
     
     timePlot = Plot[Evaluate[plotVars /. sol[[1]]], {t, 0, plotTime}, 
       PlotStyle -> {Blue, Red}, 
       PlotLegends -> {"Strain 1", "Strain 2"}, 
       AxesLabel -> {"Time", "Population"}, 
       PlotLabel -> "Two-Strain Dynamics", 
       ImageSize -> 600, 
       PlotRange -> All, 
       Epilog -> {Dashed, Blue, 
         Line[{{0, eqInfected[[1]]}, {plotTime, eqInfected[[1]]}}], 
         Dashed, Red, 
         Line[{{0, eqInfected[[2]]}, {plotTime, eqInfected[[2]]}}]}];
     Print[timePlot];
     
     Module[{finalVals, recentVals, variances, deviation}, 
       finalVals = varT /. sol[[1]] /. t -> tmax;
       deviation = Norm[finalVals - (var /. equilibrium)];
       
       recentVals = Table[
         Table[plotVars[[i]] /. sol[[1]] /. t -> tval, 
           {tval, Max[0, tmax - 20], tmax, 0.2}], 
         {i, Length[plotVars]}];
       variances = Variance /@ recentVals;
       
       If[Max[variances] > 10^(-8), 
         Print["OSCILLATIONS detected"], 
         If[deviation < 0.01, 
           Print["STABLE convergence"], 
           Print["Convergence to different state"]]];];
     (* Return Null to suppress InterpolatingFunction output *)
     Null, 
     Print["NDSolve failed"];
     $Failed]];

(*===PHASE PORTRAIT ANALYSIS===*)
phase[RHS_, var_, coP_, tmax_ : 50] := 
  Module[{intEqs, equilibrium, varT, initialCondsList, sol, phasePlot, 
    infPos, plotVars, eqPoint, rhsT}, 
   
   intEqs = intEq[RHS, var, coP];
   If[Length[intEqs] == 0, 
     equilibrium = Thread[var -> RandomReal[{0.01, 0.3}, Length[var]]];,
     equilibrium = intEqs[[1]];];
   
   varT = Through[var[t]];
   initialCondsList = pertIC[equilibrium, var, 0.15, 6];
   
   (* Substitute parameters first, then variables *)
   rhsT = RHS /. coP /. Thread[var -> varT];
   
   sol = Table[
     Module[{odeSystem}, 
       odeSystem = Join[
         Thread[D[varT, t] == rhsT],
         ic /. Thread[Through[var[0]] -> (varT /. t -> 0)]];
       Quiet[NDSolve[odeSystem, varT, {t, 0, tmax}]]], 
     {ic, initialCondsList}];
   
   sol = Select[sol, Length[#] > 0 &][[All, 1]];
   
   If[Length[sol] > 0,
     infPos = If[Length[var] >= 5, {2, 5}, {2, Min[3, Length[var]]}];
     plotVars = varT[[infPos]];
     eqPoint = var[[infPos]] /. equilibrium;
     
     phasePlot = Show[
       ParametricPlot[
         Evaluate[Table[plotVars /. sol[[i]], {i, Length[sol]}]], 
         {t, 0, tmax}, 
         AxesLabel -> {"Strain 1", "Strain 2"}, 
         PlotLabel -> "Phase Portrait", 
         ImageSize -> 500, 
         PlotStyle -> Table[Directive[Hue[i/Length[sol]], Thick], {i, Length[sol]}]],
       Graphics[{Red, PointSize[0.02], Point[eqPoint], 
         Text[Style["Equilibrium", 12, Bold], eqPoint + {0.01, 0.01}]}]];
     
     Print[phasePlot];
     Null,
     Print["Phase portrait failed"];
     $Failed]];

(*===IMPROVED BIFURCATION PLOT===*)
pltBif[curve_, hopfPoints_, paramName_] := 
  Module[{realPlot, hopfPlot, bifPlot}, 
   If[Length[curve] == 0, 
     Print[Graphics[Text["No data"], ImageSize -> 400]];
     Return[Null];];
   
   realPlot = ListLinePlot[
     Table[{curve[[i]]["Parameter"], curve[[i]]["MaxRealPart"]}, {i, Length[curve]}], 
     AxesLabel -> {ToString[paramName], "Max Re(\[Lambda])"}, 
     PlotStyle -> {Blue, Thick}, 
     GridLines -> {{}, {0}}, 
     GridLinesStyle -> {Gray, Dashed}, 
     ImageSize -> 500];
   
   bifPlot = If[Length[hopfPoints] > 0, 
     hopfPlot = ListPlot[
       Table[{hopfPoints[[i]]["parBif"], 0}, {i, Length[hopfPoints]}], 
       PlotStyle -> {Red, PointSize[0.02]}];
     Show[realPlot, hopfPlot], 
     realPlot];
   
   Print[bifPlot];
   Null];

(*===WORKFLOW FUNCTION===*)
wF[bifParam_, coP_, RHS_, var_, R12_, R21_] := 
  Module[{coPBif},
   Print["Testing dynamics at Hopf bifurcation parameter: ", bifParam];
   
   (* Create new parameter condition at bifurcation point *)
   coPBif = coP /. (coP[[1]][[1]] -> _) -> (coP[[1]][[1]] -> bifParam);
   
   (* Time series analysis at bifurcation *)
   TS[RHS, var, coPBif, R12, R21, 100];
   ];







(* Default options *)
Options[plotBifurcationDiagram] = {
  PlotRange -> Automatic,
  ImageSize -> {600, 400},
  Frame -> True,
  FrameLabel -> {"Parameter", "Population"},
  PlotStyle -> {Blue, Red, Green},
  HopfStyle -> {Red, PointSize[0.015]},
  GridLines -> Automatic,
  Background -> White
};

Options[createFigure1] = {
  BifurcationParameter -> \[Theta]1,
  ParameterRange -> {0, 1},
  NumPoints -> 200,
  PlotVariables -> {I1, I2, Y1, Y2},
  ShowHopf -> True,
  ShowPhasePortraits -> False
};

Options[plotPhasePortrait] = {
  TimeRange -> {0, 50},
  InitialConditions -> {{0.1, 0.1}},
  PlotRange -> {{0, 1}, {0, 1}},
  VectorField -> True
};

(* Main bifurcation diagram function *)
plotBifurcationDiagram[curve_, hopfPoints_, bifParam_, opts : OptionsPattern[]] := 
Module[{plotData, hopfData, stableBranches, unstableBranches, 
        plotOpts, hopfOpts},
  
  (* Extract plotting options *)
  plotOpts = FilterRules[{opts}, Options[Plot]];
  hopfOpts = OptionValue[HopfStyle];
  
  (* Process curve data - assuming curve contains {parameter, {equilibrium values}, stability} *)
  plotData = If[Length[curve] > 0,
    Table[
      {curve[[i, 1]], #} & /@ 
       Select[curve[[i, 2]], NumericQ[#] && # > 0 &], 
      {i, Length[curve]}
    ] // Flatten[#, 1] &,
    {}
  ];
  
  (* Separate stable and unstable branches based on eigenvalues *)
  stableBranches = Select[plotData, True &]; (* Simplified - add stability analysis *)
  unstableBranches = {}; (* Add unstable branch detection *)
  
  (* Extract Hopf bifurcation points *)
  hopfData = If[Length[hopfPoints] > 0,
    Table[{hopfPoints[[i, 1]], hopfPoints[[i, 2]]}, {i, Length[hopfPoints]}],
    {}
  ];
  
  (* Create the plot *)
  Show[
    (* Stable equilibria *)
    If[Length[stableBranches] > 0,
      ListPlot[stableBranches, 
        PlotStyle -> {OptionValue[PlotStyle][[1]], Thick},
        Joined -> True
      ],
      Graphics[{}]
    ],
    
    (* Unstable equilibria *)
    If[Length[unstableBranches] > 0,
      ListPlot[unstableBranches,
        PlotStyle -> {OptionValue[PlotStyle][[2]], Dashed, Thick},
        Joined -> True
      ],
      Graphics[{}]
    ],
    
    (* Hopf bifurcation points *)
    If[Length[hopfData] > 0,
      ListPlot[hopfData,
        PlotStyle -> hopfOpts
      ],
      Graphics[{}]
    ],
    
    (* Apply formatting options *)
    Frame -> OptionValue[Frame],
    FrameLabel -> OptionValue[FrameLabel],
    PlotRange -> OptionValue[PlotRange],
    ImageSize -> OptionValue[ImageSize],
    GridLines -> OptionValue[GridLines],
    Background -> OptionValue[Background],
    plotOpts
  ]
];

(* Comprehensive Figure 1 reproduction function *)
createFigure1[RHS_, var_, coP_, R12_, R21_, par_, opts : OptionsPattern[]] := 
Module[{bifParam, paramRange, numPoints, plotVars, equilibria, 
        bifurcationData, hopfData, finalPlot, parRules},
  
  (* Extract options *)
  bifParam = OptionValue[BifurcationParameter];
  paramRange = OptionValue[ParameterRange];
  numPoints = OptionValue[NumPoints];
  plotVars = OptionValue[PlotVariables];
  
  (* Convert par to replacement rules if it's not already *)
  parRules = If[MatchQ[par, {_Rule..}], 
    par, 
    Thread[par -> Table[0.1, Length[par]]] (* Default values if par is just symbols *)
  ];
  
  (* Parameter sweep for bifurcation analysis *)
  bifurcationData = Table[
    Module[{currentPar, currentEq, stability},
      (* Set current parameter value *)
      currentPar = Join[parRules, {bifParam -> paramVal}];
      
      (* Find equilibria at this parameter value *)
      currentEq = intEqV[RHS /. currentPar, var, coP];
      
      (* Return parameter value and equilibrium data *)
      {paramVal, currentEq, stability}
    ],
    {paramVal, paramRange[[1]], paramRange[[2]], 
     (paramRange[[2]] - paramRange[[1]])/numPoints}
  ];
  
  (* Detect Hopf bifurcations across parameter range *)
  hopfData = If[OptionValue[ShowHopf],
    Table[
      Module[{currentPar, jacobian, eigenvals},
        currentPar = Join[parRules, {bifParam -> paramVal}];
        (* Add Hopf detection logic here *)
        (* Return {paramVal, bifurcationValue} if Hopf exists *)
      ],
      {paramVal, paramRange[[1]], paramRange[[2]], 
       (paramRange[[2]] - paramRange[[1]])/(numPoints/10)}
    ] // Select[#, Length[#] > 0 &] &,
    {}
  ];
  
  (* Create the bifurcation diagram *)
  finalPlot = plotBifurcationDiagram[
    bifurcationData, 
    hopfData, 
    bifParam,
    FrameLabel -> {ToString[bifParam], "Population Density"},
    PlotRange -> {paramRange, Automatic}
  ];
  
  (* Add phase portraits if requested *)
  If[OptionValue[ShowPhasePortraits],
    (* Add phase portrait panels *)
    finalPlot,
    finalPlot
  ]
];

(* Phase portrait function *)
plotPhasePortrait[RHS_, var_, params_, xVar_, yVar_, opts : OptionsPattern[]] := 
Module[{odes, timeRange, initConds, vectorField, solutions},
  
  timeRange = OptionValue[TimeRange];
  initConds = OptionValue[InitialConditions];
  vectorField = OptionValue[VectorField];
  
  (* Create ODEs *)
  odes = Thread[D[var, t] == (RHS /. params)];
  
  (* Solve for trajectories *)
  solutions = Table[
    NDSolve[Join[odes, Thread[var /. t -> 0 == ic]], var, {t, 0, timeRange[[2]]}],
    {ic, initConds}
  ];
  
  (* Create phase portrait *)
  Show[
    (* Trajectories *)
    Table[
      ParametricPlot[
        Evaluate[{xVar, yVar} /. sol[[1]]], 
        {t, 0, timeRange[[2]]},
        PlotStyle -> Blue
      ],
      {sol, solutions}
    ],
    
    (* Vector field if requested *)
    If[vectorField,
      StreamPlot[
        Evaluate[{xVar, yVar} /. (RHS /. params)], 
        Evaluate[{xVar, OptionValue[PlotRange][[1]]}],
        Evaluate[{yVar, OptionValue[PlotRange][[2]]}],
        PlotStyle -> Gray
      ],
      Graphics[{}]
    ],
    
    Frame -> True,
    FrameLabel -> {ToString[xVar], ToString[yVar]},
    PlotRange -> OptionValue[PlotRange]
  ]
];

(* Utility function for parameter analysis *)
analyzeParameterSpace[RHS_, var_, par_, paramList_, coP_] := 
Module[{results, parRules},
  (* Convert par to replacement rules if needed *)
  parRules = If[MatchQ[par, {_Rule..}], 
    par, 
    Thread[par -> Table[0.1, Length[par]]]
  ];
  
  results = Table[
    Module[{currentPar, eq, stability, r0vals},
      currentPar = Join[parRules, Thread[paramList[[1]] -> paramVals]];
      eq = intEqV[RHS /. currentPar, var, coP];
      {paramVals, eq, stability}
    ],
    {paramVals, Tuples[paramList[[2]]]}
  ];
  results
];



(* Simple function to work directly with your existing continuation data *)
quickBifurcationPlot[curve_, hopfPoints_, opts : OptionsPattern[]] := 
Module[{plotData, hopfData},
  
  (* Extract parameter and equilibrium values from curve *)
  plotData = If[Length[curve] > 0 && Length[curve[[1]]] >= 2,
    Table[{curve[[i, 1]], curve[[i, 2]]}, {i, Length[curve]}],
    {}
  ];
  
  (* Extract Hopf points *)
  hopfData = If[Length[hopfPoints] > 0,
    hopfPoints,
    {}
  ];
  
  (* Create simple plot *)
  Show[
    If[Length[plotData] > 0,
      ListLinePlot[plotData, PlotStyle -> Blue],
      Graphics[{}]
    ],
    If[Length[hopfData] > 0,
      ListPlot[hopfData, PlotStyle -> {Red, PointSize[0.015]}],
      Graphics[{}]
    ],
    Frame -> True,
    FrameLabel -> {"Parameter", "Population"},
    ImageSize -> {600, 400},
    FilterRules[{opts}, Options[Graphics]]
  ]
];










End[]
EndPackage[]



