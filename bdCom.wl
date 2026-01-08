(* ::Package:: *)

ClearAll[bdCom];
bdCom[RHS_, var_, prF_: Identity] := Module[{
    bdAnResult, RN, rts, spe, alp, bet, gam, rnR,
    par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, EA, ngm, infV, F, V,
    allInfectionVars, Esys, EAsolved, siphonVars, cEj, RHSEj, eqEj, varEj, 
    timedOut,solSiph, rurSiph, nMsiph
  },

  (* Convert ODE to reaction network - suppress output since already shown *)
  {RN, rts, spe, alp, bet, gam, rnR} = Block[{Print = (Null &)}, ODE2RN[RHS, var, Identity]];

  bdAnResult = bdAn[RN, rts, var];
  (* RHS and var already provided as input, extract rest *)
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

  (* Reorder siphons by descending cardinality (larger siphons often simpler to solve) *)
  mSiSorted = SortBy[mSi, -Length[#]&];
  Print["Minimal siphons: ", mSiSorted];

  Esys = {};
  timedOut = False;
  solSiph = {};        (* {siphon, solutions} for siphons solved with replacement rules *)
  rurSiph = {};        (* {siphon, solutions} for siphons reduced to scalar equations *)
  timedOutIndices = {};   (* Siphons that timed out *)

  (* PASS 1: Build all Esys entries (fast - no solving) *)
  Do[
    siphonVars = ToExpression[mSiSorted[[j]]];
    cEj = Thread[siphonVars -> 0];
    RHSEj = RHS /. cEj;
    eqEj = Thread[RHSEj == 0];
    varEj = Complement[var, siphonVars];
    AppendTo[Esys, {eqEj, varEj}];
  , {j, Length[mSiSorted]}];


  (* PASS 2: Try to solve each system with timeout *)
  Do[
    {eqEj, varEj} = Esys[[j]];
    siphonVars = ToExpression[mSiSorted[[j]]];
    cEj = Thread[siphonVars -> 0];

    (* Try to solve and save rational solutions - with 15 sec timeout *)
    Module[{solutions, completeSolutions, hasIrrational, positiveSolutions,
            infVarsInSys, nonInfVarsInSys, scalarEq, solveResult},

      solveResult = TimeConstrained[
        Module[{sols, compSols, hasIrr, posSols, infVars, nonInfVars, scEq},
          sols = Solve[eqEj, varEj];
          compSols = Table[DeleteDuplicatesBy[Join[cEj, sol], First], {sol, sols}];
          (* Filter out any rules for parameters (not in var) - keep only rules for variables *)
          compSols = Map[Select[#, MemberQ[var, #[[1]]]&]&, compSols];
          hasIrr = !FreeQ[compSols, Sqrt | Root | Power[_, Rational[_, _]]];

          If[!hasIrr,
            (* Rational solutions - filter with onlyNN *)
            (* Special case: if siphon equals ALL infection vars, accept DFE *)
            siphonIsAllInfVars = ContainsAll[allInfectionVars, siphonVars] &&
                                 ContainsAll[siphonVars, allInfectionVars];
            posSols = Select[compSols,
              If[siphonIsAllInfVars,
                onlyNonneg[#, var],
                onlyNN[#, var, allInfectionVars]
              ] &];
            If[Length[posSols] > 0,
              {True, posSols},
              (* No positive solutions - check if any are DFE *)
              Module[{dfeSols},
                dfeSols = Select[compSols, Function[sol, AllTrue[allInfectionVars, (# /. sol) === 0 &]]];
                If[Length[dfeSols] > 0, {True, dfeSols}, {False, {}}]
              ]
            ]
          ,
            (* No rational solutions - use complementarity approach *)
            (* Variables not in siphon are positive, so divide by them *)
            nonSiphonVars = Complement[varEj, siphonVars];

            (* Simplify: substitute siphon=0, then divide by non-siphon factors *)
            RHSsub = RHS /. cEj;

            simplifiedEqs = DeleteCases[Table[
              Module[{eq, num, factorList, simplified},
                eq = Together[RHSsub[[i]]];

                (* Skip if equation is trivially 0 *)
                If[eq === 0,
                  Nothing,
                  (* Process non-trivial equation *)
                  num = Numerator[eq];
                  factorList = FactorList[num];
                  simplified = num;

                  (* Divide by non-siphon variable factors *)
                  Do[
                    If[MemberQ[nonSiphonVars, factor[[1]]],
                      simplified = Simplify[simplified / (factor[[1]]^factor[[2]])];
                    ];
                  , {factor, factorList}];

                  (* Return simplified/denominator *)
                  simplified / Denominator[eq]
                ]
              ],
            {i, 1, Length[RHS]}], Nothing];

            (* Solve simplified system *)
            If[Length[simplifiedEqs] > 0 && Length[nonSiphonVars] > 0,
              Module[{eqSimplified, compSolsNew},
                eqSimplified = Thread[simplifiedEqs == 0];
                sols = Solve[eqSimplified, nonSiphonVars];

                If[Length[sols] > 0,
                  compSolsNew = Table[DeleteDuplicatesBy[Join[cEj, sol], First], {sol, sols}];
                  hasIrr = !FreeQ[compSolsNew, Sqrt | Root | Power[_, Rational[_, _]]];

                  If[!hasIrr,

                    (* Rational solutions from complementarity *)
                    siphonIsAllInfVars = ContainsAll[allInfectionVars, siphonVars] &&
                                         ContainsAll[siphonVars, allInfectionVars];
                    posSols = Select[compSolsNew,
                      If[siphonIsAllInfVars,
                        onlyNonneg[#, var],
                        onlyNN[#, var, allInfectionVars]
                      ] &];
                    If[Length[posSols] > 0, {True, posSols}, {False, {}}]
                  ,
                    (* Still irrational - first try to extract linear solutions before RUR *)
                    Module[{linearRules, remainingEqs, remainingVars, foundLinear},
                      linearRules = {};
                      remainingEqs = eqSimplified;
                      remainingVars = nonSiphonVars;

                      (* Try to find variables that appear linearly in single equations *)
                      foundLinear = True;
                      While[foundLinear && Length[remainingEqs] > 0 && Length[remainingVars] > 0,
                        foundLinear = False;
                        (* Try each equation with each variable *)
                        Do[
                          If[!foundLinear,
                            Do[
                              If[!foundLinear && MemberQ[remainingVars, v],
                                Module[{eq, deg, sol},
                                  eq = remainingEqs[[eqIdx]];
                                  deg = Exponent[eq, v];
                                  If[deg == 1 && FreeQ[Coefficient[eq, v, 0], v] && FreeQ[Coefficient[eq, v, 1], v],
                                    (* Linear in v - solve it *)
                                    sol = Solve[eq, v];
                                    If[Length[sol] > 0 && FreeQ[sol[[1]], Sqrt | Root | Power[_, Rational[_, _]]],
                                      (* Found rational linear solution *)
                                      linearRules = Join[linearRules, sol[[1]]];
                                      remainingEqs = Delete[Simplify[remainingEqs /. sol[[1]]], eqIdx];
                                      remainingVars = DeleteCases[remainingVars, v];
                                      foundLinear = True;
                                    ];
                                  ];
                                ];
                              ];
                            , {v, remainingVars}];
                          ];
                        , {eqIdx, 1, Length[remainingEqs]}];
                      ];

                      (* If we found linear solutions, return them with remaining system *)
                      If[Length[linearRules] > 0,
                        (* Combine with siphon constraints *)
                        completeSols = DeleteDuplicatesBy[Join[cEj, linearRules], First];
                        (* Check if fully solved *)
                        If[Length[remainingVars] == 0,
                          (* Fully determined - return as replacement rules *)
                          {True, {completeSols}},
                          (* Partially solved - still need RUR for remaining vars *)
                          (* Continue to RUR with remaining system *)
                          Module[{keepVar, elimVars, eliminated, polyLHS, collected},
                            keepVar = First[remainingVars];
                            elimVars = Rest[remainingVars];
                            If[Length[elimVars] > 0,
                              eliminated = Eliminate[remainingEqs, elimVars];
                              If[eliminated =!= False && eliminated =!= True,
                                polyLHS = If[Head[eliminated] === Equal, eliminated[[1]], eliminated];
                                collected = Collect[Expand[polyLHS], keepVar, FullSimplify];
                                {True, {"poly" -> collected, "formulas" -> linearRules, "var" -> keepVar, "rur" -> True}},
                                {True, {"formulas" -> linearRules, "equations" -> remainingEqs}}
                              ],
                              (* Only one variable left - just polynomial *)
                              If[Length[remainingEqs] > 0,
                                polyLHS = remainingEqs[[1]];
                                collected = Collect[Expand[polyLHS], keepVar, FullSimplify];
                                {True, {"poly" -> collected, "formulas" -> linearRules, "var" -> keepVar, "rur" -> False}},
                                {True, {completeSols}}
                              ]
                            ]
                          ]
                        ],
                        (* No linear solutions found - proceed to original RUR logic *)
                        Module[{keepVar, otherVars, elimVars, eliminated, polyLHS, collected,
                                rationalSols, rationalFormulas},
                          keepVar = First[var];
                          If[MemberQ[nonSiphonVars, keepVar],
                            elimVars = DeleteCases[nonSiphonVars, keepVar];
                            eliminated = Eliminate[eqSimplified, elimVars];
                            If[eliminated =!= False && eliminated =!= True,
                              polyLHS = If[Head[eliminated] === Equal, eliminated[[1]], eliminated];
                              collected = Collect[Expand[polyLHS], keepVar, FullSimplify];
                              (* Build rational formulas by solving equations sequentially *)
                              (* Key: drop the equation for keepVar (x), solve remaining equations for other vars *)
                              Module[{eqsForFormulas, keepVarIdx},
                                (* Find index of keepVar in var to drop the corresponding equation *)
                                keepVarIdx = Position[var, keepVar];
                                eqsForFormulas = If[Length[keepVarIdx] > 0,
                                  Delete[eqSimplified, keepVarIdx[[1]]],
                                  eqSimplified
                                ];
                                rationalFormulas = {};
                                Module[{currentEqs, solForVar, remainingVars, foundSol, eqIdx, varIdx},
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
                                If[Length[rationalFormulas] >= Length[elimVars],
                                  {True, {"poly" -> collected, "formulas" -> rationalFormulas, "var" -> keepVar, "rur" -> True}},
                                  (* Fallback: try direct Solve on the reduced equation set *)
                                  Module[{rationalSols},
                                    rationalSols = Solve[eqsForFormulas, elimVars];
                                    If[Length[rationalSols] > 0,
                                      rationalFormulas = SelectFirst[rationalSols,
                                        FreeQ[#, Sqrt | Root | Power[_, Rational[_, _]]]&,
                                        {}
                                      ];
                                      If[rationalFormulas =!= {} && Length[rationalFormulas] > 0,
                                        {True, {"poly" -> collected, "formulas" -> rationalFormulas, "var" -> keepVar, "rur" -> True}},
                                        {True, {"poly" -> collected, "var" -> keepVar, "rur" -> False}}
                                      ],
                                      {True, {"poly" -> collected, "var" -> keepVar, "rur" -> False}}
                                    ]
                                  ]
                                ]
                              ],
                              {True, simplifiedEqs}
                            ],
                            keepVar = First[nonSiphonVars];
                            otherVars = Rest[nonSiphonVars];
                            eliminated = Eliminate[eqSimplified, otherVars];
                            If[eliminated =!= False && eliminated =!= True,
                              polyLHS = If[Head[eliminated] === Equal, eliminated[[1]], eliminated];
                              collected = Collect[Expand[polyLHS], keepVar, FullSimplify];
                              {True, {collected}},
                              {True, simplifiedEqs}
                            ]
                          ]
                        ]
                      ]
                    ]
                  ]
                ,
                  (* Could not solve - try elimination anyway *)
                  Module[{keepVar, otherVars, otherSols, eliminated, collected},
                    keepVar = First[nonSiphonVars];
                    otherVars = Rest[nonSiphonVars];
                    If[Length[otherVars] > 0,
                      otherSols = Solve[Thread[simplifiedEqs == 0], otherVars];
                      If[Length[otherSols] > 0,
                        eliminated = Simplify[(Thread[simplifiedEqs == 0])[[1]] /. otherSols[[1]]];
                        collected = Collect[Expand[eliminated], keepVar, FullSimplify];
                        {True, {"poly" -> collected, "formulas" -> otherSols[[1]], "var" -> keepVar}},
                        eliminated = Eliminate[Thread[simplifiedEqs == 0], otherVars];
                        If[eliminated =!= False && eliminated =!= True,
                          collected = Collect[eliminated[[1]], keepVar, FullSimplify];
                          {True, {collected}},
                          {True, simplifiedEqs}
                        ]
                      ],
                      {True, simplifiedEqs}
                    ]
                  ]
                ]
              ],
              {False, {}}
            ]
          ]
        ],
        15,
        $TimedOut
      ];

      If[solveResult === $TimedOut,
        AppendTo[timedOutIndices, j];
        timedOut = True;
        Break[];
      ];

      If[solveResult[[1]],
        (* Check if it's a scalar equation or list of rules *)
        If[MatchQ[solveResult[[2]], {{___Rule}..}],
          AppendTo[solSiph, {mSiSorted[[j]], solveResult[[2]]}];,
          AppendTo[rurSiph, {mSiSorted[[j]], solveResult[[2]]}];
        ];
      ];
    ];
  , {j, Length[mSiSorted]}];

  (* PASS 3: Try pairs of siphons that reduced to scalar equations *)
  nMsiph = {};
  scalarSiphons = rurSiph[[All, 1]];
  If[Length[scalarSiphons] >= 2 && !timedOut,
    Do[
      Do[
        siph1 = scalarSiphons[[i]];
        siph2 = scalarSiphons[[j]];

        (* Create combined siphon *)
        combinedSiphon = Union[siph1, siph2];
        siphonVars = ToExpression[combinedSiphon];
        cEj = Thread[siphonVars -> 0];

        (* Create system for combined siphon *)
        RHSEj = RHS /. cEj;
        eqEj = Thread[RHSEj == 0];
        varEj = Complement[var, siphonVars];

        (* Try to solve *)
        Module[{solveResult},
          solveResult = TimeConstrained[
            Module[{sols, compSols, hasIrr, posSols},
              sols = Solve[eqEj, varEj];
              compSols = Table[DeleteDuplicatesBy[Join[cEj, sol], First], {sol, sols}];
              hasIrr = !FreeQ[compSols, Sqrt | Root | Power[_, Rational[_, _]]];

              If[!hasIrr,
                siphonIsAllInfVars = ContainsAll[allInfectionVars, siphonVars] &&
                                     ContainsAll[siphonVars, allInfectionVars];
                posSols = Select[compSols,
                  If[siphonIsAllInfVars,
                    onlyNonneg[#, var],
                    onlyNN[#, var, allInfectionVars]
                  ] &];
                If[Length[posSols] > 0, {True, posSols}, {False, {}}],
                {False, {}}
              ]
            ],
            15,
            $TimedOut
          ];

          If[solveResult =!= $TimedOut && solveResult[[1]],
            If[MatchQ[solveResult[[2]], {{___Rule}..}],
              (* Found explicit solutions for the pair *)
              AppendTo[Esys, {eqEj, varEj}];
              AppendTo[nMsiph, {combinedSiphon, solveResult[[2]]}];
            ];
          ];
        ];
      , {j, i+1, Length[scalarSiphons]}];
    , {i, 1, Length[scalarSiphons]-1}];
  ];

  (* Print summary - identify categories and print *)
  Print["\n=== BOUNDARY EQUILIBRIA ==="];

  (* Individual siphons - solved with replacement rules *)
  Do[
    {siph, sols} = solSiph[[i]];
    (* Separate DFE and endemic equilibria *)
    endemicSols = Select[sols, Function[sol, !AllTrue[allInfectionVars, Function[v, (v /. sol) === 0]]]];
    If[Length[endemicSols] > 0,
      Print["Siphon face ", prF[siph], ": ", Length[endemicSols], " endemic equilibri", If[Length[endemicSols] == 1, "um", "a"]];
      Do[Print["  ", prF[endemicSols[[k]]]], {k, 1, Length[endemicSols]}];
    ];
  , {i, 1, Length[solSiph]}];

  (* Individual siphons - scalar equations *)
  Do[
    {siph, polySys} = rurSiph[[i]];
    If[ListQ[polySys] && Length[polySys] > 0,
      If[MatchQ[polySys, {___Rule}],
        (* New format with polynomial and formulas *)
        Module[{poly, formulas, keepVar, deg, sol, isRUR, equations},
          poly = "poly" /. polySys;
          formulas = "formulas" /. polySys;
          keepVar = "var" /. polySys;
          isRUR = "rur" /. polySys;
          equations = "equations" /. polySys;

          (* Check if we have formulas only (no polynomial yet) *)
          If[formulas =!= "formulas" && poly === "poly",
            (* Have formulas but no polynomial - show formulas and remaining equations *)
            Print["Siphon face ", prF[siph], " partial solution:"];
            Do[Print["  ", prF[formulas[[j]]]], {j, 1, Length[formulas]}];
            If[equations =!= "equations",
              Print["  with remaining equations: ", prF[equations]];
            ];,
            (* Have polynomial *)
            If[poly =!= "poly",
              (* Check if polynomial is in x (first variable in var) - then it's RUR *)
              If[keepVar === First[var] && poly =!= "poly",
                (* Polynomial in x - this is RUR *)
                Print["Siphon face ", prF[siph], " RUR:"];
                Print["  ", prF[Collect[Factor[poly], keepVar]], " = 0"];
                (* Show rational formulas if available *)
                If[formulas =!= "formulas",
                  Print["  where: ", prF[formulas]];
                ];,
                (* Not in x - check if linear or has formulas *)
                If[formulas =!= "formulas",
                  deg = Exponent[poly, keepVar];
                  If[deg == 1,
                    (* Linear - solve explicitly *)
                    sol = Solve[poly == 0, keepVar];
                    If[Length[sol] > 0,
                      Print["Siphon face ", prF[siph], ":"];
                      Print["  ", prF[keepVar -> (keepVar /. sol[[1]])]];
                      Print["  where:"];
                      Do[Print["    ", prF[formulas[[j]]]], {j, 1, Length[formulas]}];,
                      (* Couldn't solve - show polynomial *)
                      Print["Siphon face ", prF[siph], ":"];
                      Print["  ", prF[Collect[poly, keepVar]], " = 0"];
                      Print["  where:"];
                      Do[Print["    ", prF[formulas[[j]]]], {j, 1, Length[formulas]}];
                    ];,
                    (* Not linear - show polynomial *)
                    Print["Siphon face ", prF[siph], ":"];
                    Print["  ", prF[Collect[poly, keepVar]], " = 0"];
                    Print["  where:"];
                    Do[Print["    ", prF[formulas[[j]]]], {j, 1, Length[formulas]}];
                  ];,
                  (* No formulas - just polynomial *)
                  Print["Siphon face ", prF[siph], ":"];
                  Print["  ", prF[Collect[poly, keepVar]], " = 0"];
                ];
              ];,
              (* No poly - show error *)
              Print["Siphon face ", prF[siph],
              " reduces to polynomial system (", Length[polySys], " equations)"];
            ]
          ];
        ];,
        (* Old format *)
        If[Length[polySys] == 1,
          Print["Siphon face ", prF[siph], " reduces to: ", prF[Collect[polySys[[1]], Variables[polySys[[1]]]]], " = 0"];,
          Print["Siphon face ", prF[siph], " reduces to polynomial system (", Length[polySys], " equations)"];
        ];
      ];
    ];
  , {i, 1, Length[rurSiph]}];

  (* Pairs of siphons solved together *)
  Do[
    {combinedSiphon, sols} = nMsiph[[i]];
    numSols = Length[sols];
    Print["Siphon pair ", prF[combinedSiphon], ": ", numSols, If[numSols == 1, " solution", " solutions"]];
    Do[Print["  ", prF[sols[[k]]]], {k, 1, Min[numSols, 3]}];
  , {i, 1, Length[nMsiph]}];

  If[Length[timedOutIndices] > 0,
    Print["Timed out (", Length[timedOutIndices], " siphons after 15 sec):"];
    Do[
      idx = timedOutIndices[[i]];
      Print["  Siphon ", idx, " (", prF[mSiSorted[[idx]]], ")"];
    , {i, 1, Length[timedOutIndices]}];
    (* Add skipped siphons *)
    If[Length[timedOutIndices] > 0,
      firstTimeout = timedOutIndices[[1]];
      skipped = Range[firstTimeout + 1, Length[mSiSorted]];
      If[Length[skipped] > 0,
        Print["Skipped (", Length[skipped], " siphons after first timeout):"];
        Do[Print["  Siphon ", s, " (", mSiSorted[[s]], ")"], {s, skipped}];
      ];
    ];
  ];

  {RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, Esys, solSiph, rurSiph, nMsiph}
];

(* Tests for NGM 

RN = {"s" -> "i", "i" -> "s"};
rts = {be*s*i, ga*i};
{spe, alp, bet, gam, Rv, RHS, def} = extMat[RN];
var = ToExpression[spe];
RHS = gam . rts;
par = Par[RHS, var];
mod = {RHS, var, par};
infVars = {i};
ngm1 = NGM[mod, infVars];
Print["Test 1 - NGM without user F: K=", ngm1[[4]]];

Fuser = {{be*s}};
ngm2 = NGM[mod, infVars, Fuser];
Print["Test 2 - NGM with valid F: K=", ngm2[[4]]];

Finvalid = {{-be*s}};
ngm3 = NGM[mod, infVars, Finvalid];
Print["Test 3 - NGM with invalid F: K=", ngm3[[4]]];*)

