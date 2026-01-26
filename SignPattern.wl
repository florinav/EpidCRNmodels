(* ::Package:: *)

(* Check if a 2x2 submatrix has the bad pattern: one +1 and three -1 *)
isBadPattern[sub22_] := Module[{signs},
  signs = Sign[sub22];
  (* Bad pattern: exactly one +1 and three -1 entries *)
  Count[Flatten[signs], 1] == 1 && Count[Flatten[signs], -1] == 3
]

(* Find position of the +1 in a bad 2x2 submatrix *)
posOfPlus[sub22_] := Module[{signs, pos},
  signs = Sign[sub22];
  pos = Position[signs, 1];
  If[Length[pos] == 1, pos[[1]], Null]
]

(* Generate all pairs of row indices *)
allRowPairs[nRows_] := Flatten[Table[{i, j}, {i, 1, nRows - 1}, {j, i + 1, nRows}], 1]

(* Generate all pairs of column indices *)
allColPairs[nCols_] := Flatten[Table[{\[ScriptL], m}, {\[ScriptL], 1, nCols - 1}, {m, \[ScriptL] + 1, nCols}], 1]

(* Main function: find all bad 2x2 submatrices in S *)
findBadSub[S_] := Module[{nRows, nCols, rowPairs, colPairs, badList, sub22, plusPos},
  nRows = Length[S];
  nCols = Length[S[[1]]];
  rowPairs = allRowPairs[nRows];
  colPairs = allColPairs[nCols];

  badList = {};

  (* Check all 2x2 submatrices *)
  Do[
    Do[
      sub22 = S[[{rowPairs[[i, 1]], rowPairs[[i, 2]]},
                 {colPairs[[j, 1]], colPairs[[j, 2]]}]];

      If[isBadPattern[sub22],
        plusPos = posOfPlus[sub22];
        (* Store: {speciesRow1, speciesRow2, reactionCol1, reactionCol2,
                   positionOfPlus (1,1)=TL, (1,2)=TR, (2,1)=BL, (2,2)=BR} *)
        AppendTo[badList,
          {rowPairs[[i, 1]], rowPairs[[i, 2]],
           colPairs[[j, 1]], colPairs[[j, 2]],
           plusPos}]
      ],
    {j, 1, Length[colPairs]}],
  {i, 1, Length[rowPairs]}];

  badList
]

(* Step 2: Classify and fix bad submatrices *)

(* Given bad submatrix info, extract which species is B (produced),
   which reaction to split, and the production coefficient *)
classifyBad[S_, badInfo_] := Module[{row1, row2, col1, col2, plusPos,
    specA, specB, reacSplit, reacOther, prodCoef},

  {row1, row2, col1, col2, plusPos} = badInfo;

  (* The species with +1 entry is B (the one being produced) *)
  (* The reaction with +1 entry is the one to split *)
  Which[
    plusPos == {1, 1}, specB = row1; specA = row2; reacSplit = col1; reacOther = col2; prodCoef = S[[row1, col1]],
    plusPos == {1, 2}, specB = row1; specA = row2; reacSplit = col2; reacOther = col1; prodCoef = S[[row1, col2]],
    plusPos == {2, 1}, specB = row2; specA = row1; reacSplit = col1; reacOther = col2; prodCoef = S[[row2, col1]],
    plusPos == {2, 2}, specB = row2; specA = row1; reacSplit = col2; reacOther = col1; prodCoef = S[[row2, col2]],
    True, Return[$Failed]
  ];

  {specA, specB, reacSplit, reacOther, prodCoef}
]

(* Create new species name with "p" suffix, handling both symbols and indexed expressions *)
makePrime[spec_] := Module[{str},
  str = StringReplace[ToString[spec], {"[" -> "", "]" -> "", " " -> ""}];
  Symbol[str <> "p"]
]

(* Fix one bad submatrix by adding new species B' and new reaction B' -> p2*B *)
fixOneBad[S_, speciesNames_, rates_, badInfo_] := Module[{classif, specA, specB, reacSplit, reacOther, prodCoef,
    Snew, specNew, ratesNew, nRows, nCols, newSpecIdx, newReacIdx, Bprime, knew},

  classif = classifyBad[S, badInfo];
  If[classif === $Failed, Return[$Failed]];

  {specA, specB, reacSplit, reacOther, prodCoef} = classif;

  nRows = Length[S];
  nCols = Length[S[[1]]];

  (* Create new species name B' *)
  Bprime = makePrime[speciesNames[[specB]]];

  (* Create new rate constant with large value *)
  knew = Symbol["k" <> ToString[Length[rates] + 1]];

  (* Expand matrix: add one row (for B') and one column (for new reaction) *)
  Snew = ArrayPad[S, {{0, 1}, {0, 1}}];
  newSpecIdx = nRows + 1;
  newReacIdx = nCols + 1;

  (* Modify reaction reacSplit: B is no longer produced, B' is produced instead *)
  Snew[[specB, reacSplit]] = 0;  (* Remove B from reaction reacSplit *)
  Snew[[newSpecIdx, reacSplit]] = 1;  (* B' is produced instead *)

  (* Add new reaction B' -> prodCoef*B *)
  Snew[[newSpecIdx, newReacIdx]] = -1;  (* B' is consumed *)
  Snew[[specB, newReacIdx]] = prodCoef;  (* B is produced *)

  (* Update species and rates lists *)
  specNew = Append[speciesNames, Bprime];
  ratesNew = Append[rates, knew];

  {Snew, specNew, ratesNew}
]

(* Main function: iteratively fix all bad submatrices *)
signPatFix[S_, speciesNames_, rates_, maxIter_: 100] := Module[{Snew, specNew, ratesNew, badList, fixed, iter, diag},

  Snew = S;
  specNew = speciesNames;
  ratesNew = rates;
  diag = {};
  iter = 0;
  fixed = False;

  While[!fixed && iter < maxIter,
    badList = findBadSub[Snew];

    If[Length[badList] == 0,
      fixed = True;
      Break[]
    ];

    (* Fix first bad submatrix *)
    AppendTo[diag, {"iteration" -> iter + 1, "badFound" -> Length[badList], "badInfo" -> badList[[1]]}];
    {Snew, specNew, ratesNew} = fixOneBad[Snew, specNew, ratesNew, badList[[1]]];
    iter++;
  ];

  {Snew, specNew, ratesNew, diag}
]

(* High-level function: apply sign pattern fix to reaction network *)
Sfix[RN_, rts_, prF_: Identity] := Module[{spe, alp, bet, gam, Rv, RHSorig, defInfo,
    var, gamFix, varFix, rtsFix, diag, RvFix, RHSfix},

  (* Extract stoichiometric matrix and flux from RN *)
  {spe, alp, bet, gam, Rv, RHSorig, defInfo} = extMat[RN];
  var = ToExpression[spe];

  (* Apply sign pattern fix *)
  {gamFix, varFix, rtsFix, diag} = signPatFix[gam, var, rts];

  (* Build fixed flux vector *)
  If[Length[varFix] > Length[var],
    (* New species added, add new flux for new reaction *)
    RvFix = Join[Rv, {rtsFix[[-1]] * varFix[[-1]]}];
    RHSfix = gamFix . RvFix,
    (* No change, no bad submatrices *)
    RvFix = Rv;
    RHSfix = RHSorig
  ];

  {varFix, RHSfix, gamFix, rtsFix, diag}
]



