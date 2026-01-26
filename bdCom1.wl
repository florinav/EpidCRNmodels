(* ::Package:: *)

ClearAll[bdCom1];
bdCom1[RHS_, var_, prF_: Identity] := Module[{
    bdAnResult, RN, rts, spe, alp, bet, gam, rnR,
    par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infV, F, V,
    allInfectionVars, Esys, siphonVars, cEj, RHSEj, eqEj, varEj,
    timedOut, solSiph, rurSiph, nMsiph, mSiSorted
  },

  (* Convert ODE to reaction network - suppress output since already shown *)
  {RN, rts, spe, alp, bet, gam, rnR} = Block[{Print = (Null &)}, ODE2RN[RHS, var, Identity]];

  bdAnResult = bdAn[RN, rts, var];
  par = bdAnResult[[3]];
  cp = bdAnResult[[4]];
  mSi = bdAnResult[[5]];
  Jx = bdAnResult[[6]];
  Jy = bdAnResult[[7]];
  cDFE = bdAnResult[[8]];
  E0 = bdAnResult[[9]];
  F = bdAnResult[[10]];
  V = bdAnResult[[11]];
  K = bdAnResult[[12]];
  R0A = bdAnResult[[13]];
  infV = bdAnResult[[14]];

  allInfectionVars = ToExpression[Union[Flatten[mSi]]];

  (* First check for all equilibria using direct Solve *)
  allSol = {};
  Module[{eqAll, solAll},
    eqAll = Thread[RHS == 0];
    solAll = TimeConstrained[Quiet[Solve[eqAll, var]], 10, {}];

    (* Keep all solutions (filter later if needed) *)
    allSol = solAll;
  ];

  (* Reorder siphons by descending cardinality *)
  mSiSorted = SortBy[mSi, -Length[#]&];

  Esys = {};
  timedOut = False;
  solSiph = {};
  rurSiph = {};
  nMsiph = {};

  (* Build all system equations *)
  Do[
    siphonVars = ToExpression[mSiSorted[[j]]];
    cEj = Thread[siphonVars -> 0];
    RHSEj = RHS /. cEj;
    eqEj = Thread[RHSEj == 0];
    varEj = Complement[var, siphonVars];
    AppendTo[Esys, {eqEj, varEj}];
  , {j, Length[mSiSorted]}];

  (* Solve each siphon with timeout *)
  Do[
    {eqEj, varEj} = Esys[[j]];
    siphonVars = ToExpression[mSiSorted[[j]]];
    cEj = Thread[siphonVars -> 0];

    Module[{success},
      success = TimeConstrained[
        Module[{sols, compSols, rationalSols, irrationalSols, posSols, siphonIsAllInfVars},
          (* Step 1: Direct Solve *)
          sols = Solve[eqEj, varEj];

          (* Step 2: Complete with siphon constraints and filter parameters *)
          compSols = Table[
            Select[
              DeleteDuplicatesBy[Join[cEj, sol], First],
              MemberQ[var, First[#]]&
            ],
            {sol, sols}
          ];

          (* Step 3: Separate rational from irrational *)
          rationalSols = Select[compSols, FreeQ[#, Sqrt | Root | Power[_, Rational[_, _]]]&];
          irrationalSols = Select[compSols, !FreeQ[#, Sqrt | Root | Power[_, Rational[_, _]]]&];

          (* Step 4: Process rational solutions *)
          If[Length[rationalSols] > 0,
            siphonIsAllInfVars = ContainsAll[allInfectionVars, siphonVars] &&
                                 ContainsAll[siphonVars, allInfectionVars];
            posSols = Select[rationalSols,
              If[siphonIsAllInfVars,
                onlyNonneg[#, var],
                onlyNN[#, var, allInfectionVars]
              ] &];
            If[Length[posSols] > 0,
              AppendTo[solSiph, {mSiSorted[[j]], posSols}];
            ];
          ];

          (* Step 5: Process irrational solutions with RUR *)
          If[Length[irrationalSols] > 0,
            Module[{nonSiphonVars, RHSsub, simplifiedEqs, keepVar, keepVarIdx, eqsForFormulas,
                    elimVars, eliminated, polyLHS, collected, rationalFormulas},
              (* Use complementarity: divide by non-siphon variables *)
              nonSiphonVars = Complement[varEj, siphonVars];
              RHSsub = RHS /. cEj;

              (* Simplify equations by dividing out non-siphon factors *)
              simplifiedEqs = DeleteCases[Table[
                Module[{eq, num, factorList, simplified},
                  eq = Together[RHSsub[[i]]];
                  If[eq === 0,
                    Nothing,
                    num = Numerator[eq];
                    factorList = FactorList[num];
                    simplified = num;
                    Do[
                      If[MemberQ[nonSiphonVars, factor[[1]]],
                        simplified = Simplify[simplified / (factor[[1]]^factor[[2]])];
                      ];
                    , {factor, factorList}];
                    simplified / Denominator[eq]
                  ]
                ],
              {i, 1, Length[RHS]}], Nothing];

              If[Length[simplifiedEqs] > 0 && Length[nonSiphonVars] > 0,
                (* Keep x as the variable, drop its equation *)
                keepVar = First[var];
                If[MemberQ[nonSiphonVars, keepVar],
                  elimVars = DeleteCases[nonSiphonVars, keepVar];
                  (* Find index of keepVar in var and drop that equation *)
                  keepVarIdx = Position[var, keepVar];
                  eqsForFormulas = If[Length[keepVarIdx] > 0,
                    Delete[Thread[simplifiedEqs == 0], keepVarIdx[[1]]],
                    Thread[simplifiedEqs == 0]
                  ];

                  (* Eliminate other variables to get polynomial in keepVar *)
                  eliminated = Eliminate[eqsForFormulas, elimVars];
                  If[eliminated =!= False && eliminated =!= True,
                    polyLHS = If[Head[eliminated] === Equal, eliminated[[1]], eliminated];
                    collected = Collect[Expand[polyLHS], keepVar, FullSimplify];

                    (* Try to solve for eliminated variables rationally *)
                    rationalFormulas = {};
                    Module[{currentEqs, solForVar, remainingVars, foundSol},
                      remainingVars = elimVars;
                      currentEqs = eqsForFormulas;
                      While[Length[remainingVars] > 0 && Length[currentEqs] > 0,
                        foundSol = False;
                        Do[
                          If[!foundSol,
                            Do[
                              If[!foundSol,
                                solForVar = Solve[currentEqs[[eqIdx]], remainingVars[[varIdx]]];
                                If[Length[solForVar] > 0 && FreeQ[solForVar[[1]], Sqrt | Root | Power[_, Rational[_, _]]],
                                  rationalFormulas = Join[rationalFormulas, solForVar[[1]]];
                                  currentEqs = Delete[Simplify[currentEqs /. solForVar[[1]]], eqIdx];
                                  remainingVars = Delete[remainingVars, varIdx];
                                  foundSol = True;
                                ];
                              ];
                            , {varIdx, 1, Length[remainingVars]}];
                          ];
                        , {eqIdx, 1, Length[currentEqs]}];
                        If[!foundSol, Break[]];
                      ];
                    ];

                    (* Add RUR to rurSiph *)
                    AppendTo[rurSiph, {mSiSorted[[j]],
                      {"poly" -> collected, "formulas" -> rationalFormulas, "var" -> keepVar, "rur" -> True}}];
                  ];
                ];
              ];
            ];
          ];
          True
        ],
        15,
        $TimedOut
      ];

      If[success === $TimedOut,
        Print["Timeout on siphon ", j, ": ", mSiSorted[[j]]];
        timedOut = True;
      ];
    ];
  , {j, Length[mSiSorted]}];

  (* Print results *)
  Print["\n=== EQUILIBRIA ==="];

  (* Print DFE and endemic equilibria *)
  Module[{rationalAll, isDFE, dfe, endemic},
    rationalAll = Select[allSol, FreeQ[#, Sqrt | Root | Power[_, Rational[_, _]]]&];
    If[Length[rationalAll] > 0,
      dfe = Select[rationalAll, Function[sol, AllTrue[allInfectionVars, Function[v, (v /. sol) === 0]]]];
      endemic = Select[rationalAll, Function[sol, !AllTrue[allInfectionVars, Function[v, (v /. sol) === 0]]]];

      If[Length[dfe] > 0, Print["DFE: ", prF[dfe[[1]]]]];
      If[Length[endemic] > 0,
        Do[Print["  ", prF[Factor[endemic[[k]]]]], {k, 1, Length[endemic]}];
      ];
    ];
  ];

  (* RUR systems - only if no timeout *)
  If[!timedOut && Length[rurSiph] > 0,
    Do[
      {siph, polySys} = rurSiph[[i]];
      If[ListQ[polySys] && Length[polySys] > 0 && MatchQ[polySys, {___Rule}],
        Module[{poly, formulas, keepVar},
          poly = "poly" /. polySys;
          formulas = "formulas" /. polySys;
          keepVar = "var" /. polySys;
          If[poly =!= "poly",
            Print["RUR: ", prF[Collect[Factor[poly], keepVar]], " = 0"];
            If[formulas =!= "formulas" && Length[formulas] > 0,
              Print["  where: ", prF[formulas]];
            ];
          ];
        ];
      ];
    , {i, 1, Length[rurSiph]}];
  ];

  (* Summary *)
  Module[{nRational, nIrrational, isDFE, rationalAll},
    (* Count rational (excluding DFE) *)
    rationalAll = Select[allSol, FreeQ[#, Sqrt | Root | Power[_, Rational[_, _]]]&];
    nRational = Count[rationalAll, sol_ /; !AllTrue[allInfectionVars, Function[v, (v /. sol) === 0]]];

    (* Count irrational from RUR *)
    nIrrational = Total[Table[
      {siph, polySys} = rurSiph[[i]];
      If[ListQ[polySys] && Length[polySys] > 0 && MatchQ[polySys, {___Rule}],
        Module[{poly, keepVar},
          poly = "poly" /. polySys;
          keepVar = "var" /. polySys;
          If[poly =!= "poly", Exponent[poly, keepVar], 0]
        ],
        0
      ]
    , {i, 1, Length[rurSiph]}]];

    Print["Rational: ", nRational, ", Non-rational: ", If[nIrrational > 0, "up to " <> ToString[nIrrational], "0"]];
  ];

  {RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, Esys, solSiph, rurSiph, nMsiph}
];
