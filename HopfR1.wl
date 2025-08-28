(* ::Package:: *)

BeginPackage["HopfR1`"]

(* Core Bifurcation Analysis Functions *)
intEq::usage = "equilibria = intEq[RHS, var, coP] finds interior equilibrium points - silent version."
intEqV::usage = "equilibria = intEqV[RHS, var, coP] finds interior equilibrium points with verbose output and symbolic equilibria checking."
intEqG::usage = "equilibria = intEqG[RHS, var] finds interior equilibrium points for general models with verbose output and symbolic equilibria checking."
intEqAuto::usage =" intEqAuto[RHS, var] finds interior equilibrium with automatic parameter detection";
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
testHopf::usage ="testHopf[equilibrium_, params_] tests the occurance of Hopf."
plotBifurcationDiagram::usage = "plotBifurcationDiagram[curve_, hopfPoints_, bifParam_, opts___] creates a bifurcation diagram showing equilibria and Hopf bifurcations";
createFigure1::usage = "createFigure1[RHS_, var_, coP_, R12_, R21_, par_, opts___] reproduces Figure 1 style bifurcation analysis";
plotPhasePortrait::usage = "plotPhasePortrait[RHS_, var_, params_, xVar_, yVar_, opts___] creates phase portraits at specific parameter values";
quickBifurcationPlot::usage = "quickBifurcationPlot[curve_, hopfPoints_, opts___] creates a simple bifurcation plot from existing continuation data";
fpHopf::usage = "pos=fpHopf[RHS, var, conPar] yields all positive solutions for conPar, and examines
the eigenvalues of the first for possible closeness to Hopf";
scanPar::usage = "regions = scanPar[RHS, var, coP, par, ranges, att] scans parameter space to classify stability regions."


Begin["`Private`"]


intEq[RHS_, var_, coP_] := 
  Module[{eqs, sol, cleanSol, threshold = 10^(-10), attempts = 0, maxAttempts = 50},
   eqs = Thread[RHS == 0] /. coP;
   
   (* Try to find a single positive equilibrium with strategic starting points *)
   While[attempts < maxAttempts,
    attempts++;
    
    (* Use strategic starting points that are more likely to converge to positive equilibrium *)
    sol = Quiet[FindRoot[eqs, 
      Evaluate[Thread[{var, 
        (* Start with moderate positive values, scaled by variable importance *)
        RandomReal[{0.1, 2}, Length[var]]}]], 
      Method -> "Newton",
      MaxIterations -> 200,
      AccuracyGoal -> 12,
      PrecisionGoal -> 12]];
    
    (* Check if we got a valid solution *)
    If[Head[sol] === List && Length[sol] == Length[var],
     cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
     
     (* Check if this is a positive equilibrium *)
     If[And @@ (# > threshold & /@ (var /. cleanSol)),
      (* Found a positive equilibrium, return it *)
      Return[{cleanSol}];
      ];
     ];
    ];
   
   (* If no positive equilibrium found with Newton, try once with Solve *)
   sol = Quiet[Solve[eqs, var, Reals]];
   If[sol =!= {} && sol =!= $Failed,
    cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
    cleanSol = Select[cleanSol, And @@ (# > threshold & /@ (var /. #)) &];
    If[cleanSol =!= {}, Return[Take[cleanSol, 1]]];
    ];
   
   (* If still no positive solution, try with a more systematic approach *)
   (* Use the fact that for epidemiological models, equilibria often have *)
   (* certain scaling relationships *)
   sol = Quiet[FindRoot[eqs,
     Evaluate[Thread[{var, 
       Table[If[StringContainsQ[ToString[var[[i]]], "S"], 0.5, 0.1], {i, Length[var]}]}]],
     Method -> "Broyden",
     MaxIterations -> 500]];
   
   If[Head[sol] === List && Length[sol] == Length[var],
    cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
    If[And @@ (# > threshold & /@ (var /. cleanSol)),
     Return[{cleanSol}];
     ];
    ];
   
   (* If no positive equilibrium found, return empty *)
   Print["Warning: No positive equilibrium found"];
   {}];





           
           
  intEqG[RHS_, var_, paramValues_: {}] := 
 Module[{eqs, sol, cleanSol, threshold = 10^(-10), attempts = 0, 
   maxAttempts = 10, foundSolution = False, numericalRHS},
  
  (* Apply parameter values if provided *)
  If[paramValues =!= {},
   numericalRHS = RHS /. paramValues;
   Print["Debug: Using parameter values: ", paramValues];
   Print["Debug: Numerical RHS = ", numericalRHS];
   ,
   numericalRHS = RHS;
   Print["Debug: No parameter values provided, using symbolic RHS"];
   ];
  
  eqs = Thread[numericalRHS == 0];
  Print["Debug: equations = ", eqs];
  
  (* Check if equations are numerical *)
  If[!And @@ (NumericQ /@ (numericalRHS /. Thread[var -> RandomReal[{0.1, 1}, Length[var]]])),
   Print["Error: Equations contain symbolic parameters. Please provide numerical parameter values."];
   Print["Found parameters: ", Cases[eqs, _Symbol, Infinity] // DeleteDuplicates // Select[#, !MemberQ[var, #] &]];
   Return[{}];
   ];
  
  (* Try to find a positive equilibrium with strategic starting points *)
  While[attempts < maxAttempts && !foundSolution,
   attempts++;
   
   sol = Quiet[
     FindRoot[eqs, 
      Evaluate[Thread[{var, RandomReal[{0.1, 2}, Length[var]]}]], 
      Method -> "Newton", MaxIterations -> 100]];
   
   (* Check if we got a valid solution *)
   If[Head[sol] === List && Length[sol] == Length[var],
    cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
    
    (* Check if this is a positive equilibrium *)
    If[And @@ (# > threshold & /@ (var /. cleanSol)),
     Print["Found positive equilibrium: ", cleanSol];
     foundSolution = True;
     Return[{cleanSol}];
     ];
    ];
   ];
  
  (* If Newton failed, try Solve *)
  If[!foundSolution,
   sol = Quiet[Solve[eqs, var, Reals]];
   
   If[sol =!= {} && sol =!= $Failed,
    cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
    cleanSol = Select[cleanSol, And @@ (# > threshold & /@ (var /. #)) &];
    If[cleanSol =!= {}, 
     foundSolution = True;
     Print["Found positive equilibrium with Solve: ", cleanSol[[1]]];
     Return[Take[cleanSol, 1]];
     ];
    ];
   ];
  
  (* Try NSolve as backup *)
  If[!foundSolution,
   sol = Quiet[NSolve[eqs, var, Reals]];
   
   If[sol =!= {} && sol =!= $Failed,
    cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
    cleanSol = Select[cleanSol, And @@ (# > threshold & /@ (var /. #)) &];
    If[cleanSol =!= {}, 
     foundSolution = True;
     Print["Found positive equilibrium with NSolve: ", cleanSol[[1]]];
     Return[Take[cleanSol, 1]];
     ];
    ];
   ];
  
  (* Try specific initial conditions for epidemiological models *)
  If[!foundSolution,
   sol = Quiet[
     FindRoot[eqs, 
      Evaluate[Thread[{var, {1.0, 0.1, 0.1, 0.01}}]], 
      Method -> "Broyden", MaxIterations -> 200]];
   
   If[Head[sol] === List && Length[sol] == Length[var],
    cleanSol = sol /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
    If[And @@ (# > threshold & /@ (var /. cleanSol)),
     Print["Found positive equilibrium with specific IC: ", cleanSol];
     foundSolution = True;
     Return[{cleanSol}];
     ];
    ];
   ];
  
  (* If no positive equilibrium found *)
  Print["Warning: No positive interior equilibrium found after all attempts"];
  {}
  ]

(* Helper function to find interior equilibrium with automatic parameter detection *)
intEqAuto[RHS_, var_] := Module[{params, paramValues},
  (* Extract symbolic parameters *)
  params = Cases[RHS, _Symbol, Infinity] // DeleteDuplicates // Select[#, !MemberQ[var, #] &];
  
  Print["Detected parameters: ", params];
  Print["Please provide numerical values for these parameters."];
  Print["Example usage: intEq[RHS, var, {\[Beta]2 -> 0.1, \[Beta]3 -> 0.05, \[CapitalLambda] -> 1.0, \[Delta] -> 0.2, \[Mu]2 -> 0.1, \[Mu]3 -> 0.1, \[Mu]4 -> 0.1}]"];
  
  {}
]
           
           









(*===INTERIOR EQUILIBRIUM - VERBOSE VERSION===*)


intEqV[RHS_, var_, coP_] := 
 Module[{eqs, sols, cleanSols, orthantSols, positiveSols, 
   boundaryFPs, threshold = 10^-10, uniqueSols, equilibriaAssoc, 
   outputSummary},
  
  eqs = Thread[RHS == 0] /. coP;
  
  (* Comprehensive search targeting 3 boundary + 1 interior equilibria *)
  sols = Join[
    (* Search for Disease-Free Equilibrium (all infectious = 0) *)
    Table[
      Module[{startPoint = Table[0.001, Length[var]]},
        startPoint[[1]] = RandomReal[{2, 4}]; (* Large S *)
        Quiet[FindRoot[eqs, Evaluate[Thread[{var, startPoint}]], Method -> "Newton"]]
      ],
      {5}
    ],
    
    (* Search for strain 1 dominance (I1>0, I2=0) *)
    Table[
      Module[{startPoint = Table[0.001, Length[var]]},
        startPoint[[1]] = RandomReal[{1, 3}]; (* S *)
        startPoint[[2]] = RandomReal[{0.1, 0.8}]; (* I1 > 0 *)
        (* Keep I2-related compartments small *)
        If[Length[var] >= 5, startPoint[[5]] = 0.001]; (* I2 = 0 *)
        Quiet[FindRoot[eqs, Evaluate[Thread[{var, startPoint}]], Method -> "Newton"]]
      ],
      {8}
    ],
    
    (* Search for strain 2 dominance (I1=0, I2>0) *)
    Table[
      Module[{startPoint = Table[0.001, Length[var]]},
        startPoint[[1]] = RandomReal[{1, 3}]; (* S *)
        startPoint[[2]] = 0.001; (* I1 = 0 *)
        If[Length[var] >= 5, startPoint[[5]] = RandomReal[{0.1, 0.8}]]; (* I2 > 0 *)
        Quiet[FindRoot[eqs, Evaluate[Thread[{var, startPoint}]], Method -> "Newton"]]
      ],
      {8}
    ],
    
    (* Search for coexistence (both I1>0 and I2>0) *)
    Table[
      Module[{startPoint = Table[0.01, Length[var]]},
        startPoint[[1]] = RandomReal[{1, 2.5}]; (* S *)
        startPoint[[2]] = RandomReal[{0.05, 0.3}]; (* I1 > 0 *)
        If[Length[var] >= 5, startPoint[[5]] = RandomReal[{0.05, 0.3}]]; (* I2 > 0 *)
        Quiet[FindRoot[eqs, Evaluate[Thread[{var, startPoint}]], Method -> "Newton"]]
      ],
      {10}
    ],
    
    (* Additional boundary searches *)
    Table[
      Module[{startPoint = RandomReal[{0, 0.1}, Length[var]]},
        startPoint[[1]] = RandomReal[{1, 3}]; (* Ensure S > 0 *)
        Quiet[FindRoot[eqs, Evaluate[Thread[{var, startPoint}]], Method -> "Newton"]]
      ],
      {8}
    ],
    
    (* General searches with different methods *)
    Table[
      Quiet[FindRoot[eqs, 
        Evaluate[Thread[{var, RandomReal[{0.001, 0.5}, Length[var]]}]], 
        Method -> "Broyden"]],
      {5}
    ]
  ];
  
  (* Filter valid solutions *)
  sols = Select[sols, Head[#] === List && Length[#] == Length[var] &];
  
  (* Clean small values *)
  cleanSols = sols /. {x_?NumericQ :> If[Abs[x] < threshold, 0, x]};
  
  (* Aggressive duplicate removal *)
  uniqueSols = {};
  Do[
    Module[{fpVals = var /. sol, isDuplicate = False, isValid = True},
      
      (* Check if solution satisfies equations *)
      Module[{residuals = (RHS /. coP /. sol)},
        If[!VectorQ[residuals, NumericQ] || Max[Abs[residuals]] > 10^-5,
          isValid = False
        ]
      ];
      
      (* Check for duplicates *)
      If[isValid,
        Do[
          If[Norm[fpVals - (var /. existing)] < 10^-7,
            isDuplicate = True; Break[]
          ],
          {existing, uniqueSols}
        ];
        
        If[!isDuplicate, AppendTo[uniqueSols, sol]]
      ]
    ],
    {sol, cleanSols}
  ];
  
  cleanSols = uniqueSols;
  
  (* Separate into orthant and positive solutions *)
  orthantSols = Select[cleanSols, And @@ (# >= -threshold & /@ (var /. #)) &];
  positiveSols = Select[cleanSols, And @@ (# > threshold & /@ (var /. #)) &];
  boundaryFPs = Select[orthantSols, ! MemberQ[positiveSols, #] &];
  
  (* Systematic classification and naming *)
  equilibriaAssoc = Association[];
  
  (* Find Disease-Free Equilibrium (E0) *)
  Module[{dfeFound = False},
    Do[
      Module[{fpVals = var /. fp, infectiousVals},
        infectiousVals = Drop[fpVals, 1]; (* All except S *)
        If[And @@ (Abs[#] < threshold & /@ infectiousVals) && !dfeFound,
          equilibriaAssoc["E0"] = fp;
          dfeFound = True;
        ]
      ],
      {fp, boundaryFPs}
    ];
  ];
  
  (* Find strain dominance equilibria among remaining boundary points *)
  Module[{remainingBoundary, strain1Found = False, strain2Found = False},
    remainingBoundary = Select[boundaryFPs, 
      !KeyExistsQ[equilibriaAssoc, "E0"] || # =!= equilibriaAssoc["E0"] &];
    
    Do[
      Module[{fpVals = var /. fp, i1Val, i2Val},
        (* Assume I1 is var[[2]] and I2 is var[[5]] - adjust indices as needed *)
        i1Val = If[Length[var] >= 2, fpVals[[2]], 0];
        i2Val = If[Length[var] >= 5, fpVals[[5]], 0];
        
        Which[
          (* Strain 1 dominance: I1 > 0, I2 \[TildeTilde] 0 *)
          i1Val > threshold && Abs[i2Val] < threshold && !strain1Found,
          equilibriaAssoc["E1"] = fp; strain1Found = True,
          
          (* Strain 2 dominance: I1 \[TildeTilde] 0, I2 > 0 *)
          Abs[i1Val] < threshold && i2Val > threshold && !strain2Found,
          equilibriaAssoc["E2"] = fp; strain2Found = True,
          
          (* Other boundary case *)
          !KeyExistsQ[equilibriaAssoc, "E3"] && 
          (!strain1Found || !strain2Found || 
           (Abs[i1Val] < threshold && Abs[i2Val] < threshold)),
          equilibriaAssoc["E3"] = fp
        ]
      ],
      {fp, remainingBoundary}
    ];
  ];
  
  (* Find interior coexistence equilibrium *)
  If[Length[positiveSols] > 0,
    equilibriaAssoc["E4"] = First[positiveSols]
  ];
  
  (* If we don't have exactly 4, try to fill gaps *)
  Module[{allFound = orthantSols, counter = 0},
    Do[
      If[!MemberQ[Values[equilibriaAssoc], fp],
        counter++;
        If[counter <= 4 - Length[equilibriaAssoc],
          equilibriaAssoc["E" <> ToString[Length[equilibriaAssoc]]] = fp
        ]
      ],
      {fp, allFound}
    ];
  ];
  
  (* Print only the essential results *)
  Print["equ: ", equilibriaAssoc];
  Print["nbFP: ", Length[equilibriaAssoc]];
  Print["posEq: ", positiveSols];
  Print["var: ", var];
  
  (* Return values in the expected format for backward compatibility *)
  {equilibriaAssoc, Length[equilibriaAssoc], positiveSols, var}
];

(* Helper function for classification *)
ClassifyEquilibrium[fpVals_, var_, threshold_] := Module[{nVars, infectiousVals, zeroIndices, positiveIndices},
  nVars = Length[fpVals];
  infectiousVals = Drop[fpVals, 1];
  zeroIndices = Position[infectiousVals, x_ /; Abs[x] < threshold] // Flatten;
  positiveIndices = Position[infectiousVals, x_ /; x > threshold] // Flatten;
  
  Which[
    Length[zeroIndices] == Length[infectiousVals], "Disease-Free Equilibrium (DFE)",
    Length[positiveIndices] == 1, "Single strain dominance",
    Length[positiveIndices] == Length[infectiousVals], "Full coexistence",
    True, "Partial coexistence/boundary"
  ]
];


fpHopf[RHS_, var_, conPar_] := 
 Module[{eqns, allSols, posSols, eq, jac, eigs, complexEigs, realParts},
  
  (* Create equations properly *)
  eqns = Thread[(RHS /. conPar) == 0];
  
  (* Find all real solutions *)
  allSols = Quiet[NSolveValues[eqns, var, Reals]];
  Print[Length[allSols],"real solutions found"];
  
  (* Filter for positive solutions (coexistence equilibria) *)
  posSols = Select[allSols, AllTrue[#, # > 0 &] &];
  Print["Positive solutions (coexistence): ", Length[posSols]];
  
  If[Length[posSols] == 0, 
   Print["No coexistence equilibrium found"];
   Return["No equilibrium"]];
  
  (* Use first positive solution *)
  eq = posSols[[1]];
  Print["Using equilibrium: ", eq];
  
  (* Compute Jacobian and eigenvalues *)
  jac = D[RHS, {var}] /. conPar /. Thread[var->eq];
  eigs = Chop[Eigenvalues[jac] // N];
  complexEigs = Select[eigs, Im[#] != 0 &];
  
  Print["Eigenvalues: ", eigs];
  
  If[Length[complexEigs] >= 2,
   realParts = Re[complexEigs];
   If[Max[Abs[realParts]] < 0.1 && Max[Abs[Im[complexEigs]]] > 0.01,
    Print[Length[complexEigs],"complex eigs.","Result: Hopf - Complex eigenvalues with real 
parts near zero (potential oscillatory instability)"];
    "Hopf",
    Print["Result:Complex eigenvalues but not near Hopf conditions (stable spiral)"];
    "Complex"
    ],
   Print["Result: Real - Only real eigenvalues (stable/unstable node or saddle point)"];
   "Real"
   ];{posSols}
  ];

(* Usage:
{plot,  conPar} = ploR[par, coP, R0A, E0, E1, E2, R12, R21, bifInd];
result = testHopf[RHS, var, conPar];
*)





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






hopfD[curve_] := 
 Module[{hopfPoints, stabilityPlot, hopfPlot},
  If[Length[curve] < 3, 
   Print["hopfD needs \[GreaterEqual]3 points, got ", Length[curve]]; 
   Return[{}]];
  
  (* Check if curve has required fields *)
  If[!KeyExistsQ[curve[[1]], "MaxRealPart"],
   Print["Error: curve missing stability data. Run cont2 first."];
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




(* Hopf bifurcation test function - takes equilibrium and parameters *)
testHopf[equilibrium_, params_] := 
 Module[{jac, eigs, complexEigs, realParts},
  If[Length[equilibrium] == 0, Return["No equilibrium"]];
  
  jac = D[RHS, {var}] /. params /. equilibrium;
  eigs = Eigenvalues[jac] // N;
  complexEigs = Select[eigs, Im[#] != 0 &];
  
  If[Length[complexEigs] >= 2,
   realParts = Re[complexEigs];
   If[Max[Abs[realParts]] < 0.1 && Max[Abs[Im[complexEigs]]] > 0.01,
    "Hopf",
    "Complex"
    ],
   "Real"
   ]
  ]













(* Parameter scanning function for Hopf bifurcation detection *)
scanPar[RHS_, var_, par_, coP_, plotInd_, plot_, 
  wRan_: 0.5, hRan_: 0.5, gridRes_: 8, hTol_: 0.01] := 
 Module[{plotVar, beta1Symbol, beta2Symbol, beta1Center, beta2Center,
   beta1Range, beta2Range, gridPoints, hopfPoints, hopfRealParts, testCoP,
   fpResult, posSols, eigs, complexEigs, overlayPlot,
   beta1Min, beta1Max, beta2Min, beta2Max, totalTests, hopfCount},
  
  Print["Scanning ", gridRes^2, " points (wRan=", wRan, ", hRan=", hRan, ", hTol=", hTol, "] for Hopf bifurcations..."];
  
  (* Extract the two parameters to vary *)
  plotVar = par[[plotInd]];
  beta1Symbol = plotVar[[1]];
  beta2Symbol = plotVar[[2]];
  
  (* Get center values from coP *)
  beta1Center = beta1Symbol /. coP;
  beta2Center = beta2Symbol /. coP;
  
  (* Define ranges around center values *)
  beta1Min = beta1Center * (1 - wRan);
  beta1Max = beta1Center * (1 + wRan);
  beta2Min = beta2Center * (1 - hRan);
  beta2Max = beta2Center * (1 + hRan);
  
  (* Ensure positive values *)
  beta1Min = Max[beta1Min, 0.001];
  beta2Min = Max[beta2Min, 0.001];
  
  (* Create uniform grid *)
  beta1Range = Table[beta1Min + (beta1Max - beta1Min)*k/(gridRes - 1), {k, 0, gridRes - 1}];
  beta2Range = Table[beta2Min + (beta2Max - beta2Min)*k/(gridRes - 1), {k, 0, gridRes - 1}];
  
  gridPoints = Flatten[Table[{b1, b2}, {b1, beta1Range}, {b2, beta2Range}], 1];
  totalTests = Length[gridPoints];
  
  hopfPoints = {};
  hopfRealParts = {};
  hopfCount = 0;
  
  (* Test each grid point for Hopf bifurcations *)
  Do[
   (* Create parameter set for this grid point *)
   testCoP = Join[
     DeleteCases[coP, Rule[beta1Symbol, _] | Rule[beta2Symbol, _]],
     {beta1Symbol -> gridPoints[[j, 1]], beta2Symbol -> gridPoints[[j, 2]]}
     ];
   
   (* Test for equilibrium and eigenvalues *)
   fpResult = Quiet[fpHopf[RHS, var, testCoP]];
   
   (* Check if fpHopf succeeded *)
   If[Length[fpResult] == 3 && Length[fpResult[[1]]] > 0,
    {posSols, eigs, complexEigs} = fpResult;
    
    (* Check for Hopf bifurcation: complex eigenvalues with real part ~= 0 *)
    If[Length[complexEigs] > 0 && Abs[Max[Re[complexEigs]]] < hTol,
     hopfCount++;
     AppendTo[hopfPoints, {gridPoints[[j, 1]], gridPoints[[j, 2]]}];
     AppendTo[hopfRealParts, Max[Re[complexEigs]]];
     ];
    ];
   
   , {j, 1, totalTests}];
  
  If[Length[hopfPoints] > 0,
   Print["Found ", Length[hopfPoints], " Hopf points with real parts: ", hopfRealParts];,
   Print["Found ", Length[hopfPoints], " Hopf points"];
   ];
  
  (* Always create overlay plot with search rectangle *)
  overlayPlot = Show[
   plot,
   Graphics[{
     (* Search range rectangle *)
     EdgeForm[{Thick, Gray, Dashed}], FaceForm[None],
     Rectangle[{beta1Min, beta2Min}, {beta1Max, beta2Max}],
     
     (* Center point *)
     Blue, PointSize[0.012],
     Point[{{beta1Center, beta2Center}}],
     
     (* Hopf points if any found *)
     If[Length[hopfPoints] > 0, {Red, PointSize[0.012], Point[hopfPoints]}, {}],
     
     (* Labels *)
     Text[Style["Search Range" <> If[Length[hopfPoints] > 0, " (" <> ToString[Length[hopfPoints]] <> " Hopf)", " (No Hopf)"], 
                FontSize -> 8, FontWeight -> Bold, FontColor -> Gray], 
          {beta1Min + 0.02*(beta1Max - beta1Min), 
           beta2Max - 0.02*(beta2Max - beta2Min)}]
     }],
   PlotLabel -> "Parameter Scan: " <> ToString[Length[hopfPoints]] <> " Hopf Points Found",
   PlotRange -> {{beta1Min*0.9, beta1Max*1.1}, {beta2Min*0.9, beta2Max*1.1}}
   ];
  
  overlayPlot
  ];



End[]
EndPackage[]



