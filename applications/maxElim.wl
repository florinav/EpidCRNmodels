(* ::Package:: *)

(* ===================================================================== *)
(* maxElim - Find Maximal Eliminable Subset for Mass-Action ODEs        *)
(* Based on Saez-Feliu-Wiuf linear elimination theory                   *)
(* ===================================================================== *)

(* Identify species that appear linearly in their RHS *)
LinearSpecies[RN_List, rts_List] := Module[{RHS, spe, var, linear = {}, rhs},
  RHS = RHSFromCRN[RN, rts];
  spe = extSpe[RN];
  var = ToExpression /@ spe;
  Do[
    rhs = var[[i]] /. RHS;
    If[Exponent[rhs, var[[i]]] <= 1,
      AppendTo[linear, spe[[i]]]
    ],
  {i, Length[spe]}];
  DeleteDuplicates[linear]
];

(* Check if steady-state expressions are positive *)
PositiveSteadyStateQ[RN_List, rts_List, U_List] := Module[
  {rules, allPos = True, expr},
  rules = LinearEliminationRules[RN, rts, U];
  Do[
    expr = u /. rules;
    If[!FreeQ[expr, Alternatives[LinearSolve, PseudoInverse]],
      allPos = False; Break[]
    ],
  {u, ToExpression /@ U}];
  allPos
];

(* Find all maximal eliminable subsets *)
MaximalEliminableSubsets[RN_List, rts_List, checkPositive_:True] := Module[
  {linSpe, allMaxSets, positive},
  Print["Step 1: Identify linear species"];
  linSpe = LinearSpecies[RN, rts];
  Print["  Linear species: ", linSpe];
  Print["  (", Length[linSpe], " out of ", Length[extSpe[RN]], " total)\n"];

  Print["Step 2: Find all maximal eliminable subsets"];
  Print["  Checking FSW conditions (noninteracting, U-linear, spanning tree)\n"];

  allMaxSets = maxElim[RN, rts];
  Print["  Found ", Length[allMaxSets], " maximal eliminable set(s):"];
  Do[Print["    [", i, "] ", allMaxSets[[i]]], {i, Length[allMaxSets]}];

  If[Length[allMaxSets] == 0, Print["  No eliminable species found."]; Return[{}]];

  If[checkPositive,
    Print["\nStep 3: Check positivity of steady-state expressions"];
    Do[
      positive = PositiveSteadyStateQ[RN, rts, allMaxSets[[i]]];
      Print["  Set [", i, "] positive? ", positive],
    {i, Length[allMaxSets]}];
  ];

  Print["\nResult: ", Length[allMaxSets], " maximal eliminable set(s) found"];
  allMaxSets
];
