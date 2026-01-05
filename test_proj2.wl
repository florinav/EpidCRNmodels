(* proj2 - Project n-dimensional system onto 2D plane for visualization *)

proj2[RHS_, var_, projVars_, cN_:{}, opts___] := Module[{
  rhN, otherVars, varIndices, projIndices, otherIndices,
  otherSol, rhN2d, cFP, Xs, fpProj, eigVals, stabilities,
  pNull, pStream, pFP, r1, r2, xM},

  (* Substitute parameters *)
  rhN = If[Length[cN] > 0, RHS //. cN, RHS];

  (* Identify projection and other variables *)
  If[Length[projVars] != 2,
    Print["Error: projVars must have exactly 2 variables"];
    Return[{{}, {}, {}, Graphics[], Graphics[]}]
  ];

  otherVars = Complement[var, projVars];
  projIndices = Flatten[Position[var, #] & /@ projVars];
  otherIndices = Flatten[Position[var, #] & /@ otherVars];

  Print["Projecting ", Length[var], "D system onto (", projVars[[1]], ",", projVars[[2]], ")"];

  (* Try to solve nullclines of other variables *)
  If[Length[otherVars] > 0,
    Print["Solving nullclines for: ", otherVars];
    otherSol = Quiet[Solve[Thread[rhN[[otherIndices]] == 0], otherVars]];

    If[Length[otherSol] > 0,
      (* Use first solution *)
      rhN2d = rhN[[projIndices]] /. otherSol[[1]] // Simplify;
      Print["Projected 2D dynamics: ", rhN2d];,
      (* Cannot solve analytically - just use full system without projection *)
      Print["Warning: Cannot solve nullclines analytically, using full dynamics"];
      rhN2d = rhN[[projIndices]];
    ];,
    (* Already 2D *)
    rhN2d = rhN;
  ];

  (* Find fixed points of full n-D system *)
  cFP = Quiet[NSolve[Thread[rhN == 0], var, Reals]];
  Xs = DeleteDuplicates[
    Select[var /. cFP, AllTrue[#, NumericQ] && AllTrue[#, # >= 0 &] &],
    Norm[#1 - #2] < 10^-6 &];

  Print["fp (", Length[var], "D)=", Xs];

  (* Project onto 2D plane *)
  fpProj = If[Length[Xs] > 0, Xs[[All, projIndices]], {}];
  Print["fp proj (", projVars, ")=", fpProj];

  (* Compute eigenvalues and stability *)
  If[Length[Xs] > 0,
    Module[{jac},
      jac = Outer[D, rhN, var];
      eigVals = Table[Chop[Eigenvalues[jac /. Thread[var -> Xs[[i]]]]], {i, Length[Xs]}];
      stabilities = Table[
        Which[
          AllTrue[Re[eigVals[[i]]], # < 0 &], "Stable",
          AllTrue[Re[eigVals[[i]]], # > 0 &], "Unstable",
          True, "Saddle"
        ],
        {i, Length[Xs]}
      ];
      Print["eig=", eigVals, " stab=", stabilities];
    ];,
    eigVals = {}; stabilities = {};
  ];

  (* Plotting ranges *)
  If[Length[fpProj] > 0,
    xM = Max /@ Transpose[fpProj];
    r1 = {projVars[[1]], -0.05, xM[[1]] + 0.3};
    r2 = {projVars[[2]], -0.05, xM[[2]] + 0.3};,
    r1 = {projVars[[1]], 0, 1};
    r2 = {projVars[[2]], 0, 1};
  ];

  (* Nullcline plots *)
  pNull = ContourPlot[{rhN2d[[1]], rhN2d[[2]]},
    r1, r2,
    Contours -> {0},
    ContourStyle -> {Directive[Blue, Thick], Directive[Red, Thick]},
    ContourShading -> None,
    Frame -> True,
    FrameLabel -> {ToString[projVars[[1]]], ToString[projVars[[2]]]}
  ];

  (* Streamplot *)
  pStream = StreamPlot[{rhN2d[[1]], rhN2d[[2]]},
    r1, r2,
    StreamStyle -> Arrowheads[0.02],
    ColorFunction -> "Rainbow",
    StreamPoints -> Fine
  ];

  (* Fixed points *)
  If[Length[fpProj] > 0,
    pFP = Graphics[{
      EdgeForm[Directive[Black, Thick]],
      MapThread[
        {Switch[#2, "Stable", Green, "Unstable", Red, "Saddle", Orange, _, Black],
         Disk[#1, 0.015]} &,
        {fpProj, stabilities}
      ]
    }];,
    pFP = Graphics[];
  ];

  (* Return *)
  {fpProj, eigVals, stabilities, pFP, pStream}
];
