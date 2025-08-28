(* ::Package:: *)

(* EpidCRN Core Subpackage - Basic network analysis and utilities *)
BeginPackage["EpidCRN`Core`"];

(* ========================================================================= *)
(* CORE UTILITIES - Reaction parsing and species extraction *)
(* ========================================================================= *)


Begin["`Private`"];

(* ========================================================================= *)
(* CORE IMPLEMENTATION *)
(* ========================================================================= *)

(* Enhanced string to symbolic conversion *)
strToSymb[reactions_] := Module[{convertSide},
  convertSide[side_String] := Module[{result, terms, converted},
    If[side == "0" || side == "", Return[0]];
    terms = StringSplit[side, "+"];
    converted = Map[Module[{term, cleanTerm, coeff, species},
      term = StringTrim[#];
      Which[
        StringMatchQ[term, RegularExpression["^\\d+[A-Za-z0-9]+$"]],
        Module[{digits, letters},
          digits = StringTake[term, 1];
          letters = StringDrop[term, 1];
          ToExpression[digits]*letters
        ],
        StringMatchQ[term, RegularExpression["^[A-Za-z0-9]+$"]],
        term,
        True, term
      ]
    ] &, terms];
    If[Length[converted] == 1, converted[[1]], Apply[Plus, converted]]
  ];
  Map[Function[reaction,
    If[Head[reaction] === Rule,
      convertSide[First[reaction]] -> convertSide[Last[reaction]],
      {convertSide[First[reaction]], convertSide[Last[reaction]]}
    ]
  ], reactions]
];

(* Enhanced compToAsso function *)
compToAsso[side_] := Module[{coeffs}, 
  coeffs = <||>;
  If[side === 0, Return[coeffs]];
  
  If[Head[side] === Plus, 
   Do[
    If[Head[term] === Times, 
     Module[{parts, strings, nonStrings},
      parts = List @@ term;
      strings = Select[parts, StringQ];
      nonStrings = Select[parts, !StringQ[#] &];
      If[Length[strings] >= 1,
       Do[
        coeffs[str] = If[Length[nonStrings] == 0, 1, 
                        If[Length[nonStrings] == 1, nonStrings[[1]], Times @@ nonStrings]];
        , {str, strings}];,
       If[StringQ[term], coeffs[term] = 1]
       ]
      ], 
     If[StringQ[term], coeffs[term] = 1];
    ], {term, List @@ side}], 
   If[Head[side] === Times, 
    Module[{parts, strings, nonStrings},
     parts = List @@ side;
     strings = Select[parts, StringQ];
     nonStrings = Select[parts, !StringQ[#] &];
     If[Length[strings] >= 1,
      Do[
       coeffs[str] = If[Length[nonStrings] == 0, 1, 
                       If[Length[nonStrings] == 1, nonStrings[[1]], Times @@ nonStrings]];
       , {str, strings}];,
      If[StringQ[side], coeffs[side] = 1]
      ]
     ], 
    If[StringQ[side], coeffs[side] = 1]]];
  coeffs];

(* Main extMat function *)
extMat[reactions_] := Module[{
  spe, al, be, gamma, Rv, RHS, 
  numReactions, numSpecies, reactants, products, 
  var, rv, tk, complexes, linkageClasses, deficiency,
  defFormula, defTerms, defResult, Nc, l, s, leftSide, rightSide},
  
  (* Extract species *)
  spe = Module[{allSpecies, reactants, products}, 
    allSpecies = {};
    Do[
      If[Head[reactions[[i]]] === Rule,
        leftSide = First[reactions[[i]]];
        rightSide = Last[reactions[[i]]],
        leftSide = reactions[[i, 1]];
        rightSide = reactions[[i, 2]]
      ];
      reactants = compToAsso[leftSide];
      products = compToAsso[rightSide];
      allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
    , {i, Length[reactions]}];
    DeleteDuplicates[allSpecies]
  ];
  
  numSpecies = Length[spe];
  numReactions = Length[reactions];
  
  (* Build stoichiometric matrices *)
  al = Table[0, {numSpecies}, {numReactions}];
  be = Table[0, {numSpecies}, {numReactions}];
  
  Do[
    If[Head[reactions[[j]]] === Rule,
      leftSide = First[reactions[[j]]];
      rightSide = Last[reactions[[j]]],
      leftSide = reactions[[j, 1]];
      rightSide = reactions[[j, 2]]
    ];
    reactants = compToAsso[leftSide];
    products = compToAsso[rightSide];
    Do[
      If[KeyExistsQ[reactants, spe[[i]]], 
        al[[i, j]] = reactants[spe[[i]]]];
      If[KeyExistsQ[products, spe[[i]]], 
        be[[i, j]] = products[spe[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  gamma = be - al;
  var = ToExpression[spe];
  rv = expM[var, al // Transpose];
  tk = Array[Symbol["k" <> ToString[#]] &, numReactions];
  Rv = tk*rv;
  RHS = gamma . Rv;
  
  (* Deficiency calculation *)
  complexes = Module[{complexList},
    complexList = {};
    Do[
      If[Head[reactions[[i]]] === Rule,
        leftSide = First[reactions[[i]]];
        rightSide = Last[reactions[[i]]],
        leftSide = reactions[[i, 1]];
        rightSide = reactions[[i, 2]]
      ];
      If[!MemberQ[complexList, leftSide], AppendTo[complexList, leftSide]];
      If[!MemberQ[complexList, rightSide], AppendTo[complexList, rightSide]];
    , {i, Length[reactions]}];
    complexList
  ];
  
  Nc = Length[complexes];
  l = 1; (* Simplified - would need proper linkage class calculation *)
  s = If[gamma == {}, 0, MatrixRank[gamma]];
  
  defFormula = "\[Delta] = Nc - \[ScriptL] - s";
  defTerms = StringJoin[
    "Nc = ", ToString[Nc], " (complexes), ",
    "\[ScriptL] = ", ToString[l], " (linkage classes), ", 
    "s = ", ToString[s], " (stoich dimension)"
  ];
  defResult = StringJoin[
    "\[Delta] = ", ToString[Nc], " - ", ToString[l], " - ", ToString[s], 
    " = ", ToString[Nc - l - s]
  ];
  
  {spe, al, be, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
];

(* Other core functions *)
extSpe[reactions_] := Module[{allSpecies, reactants, products}, 
  allSpecies = {};
  Do[
   reactants = compToAsso[reactions[[i, 1]]];
   products = compToAsso[reactions[[i, 2]]];
   allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
   , {i, Length[reactions]}];
  DeleteDuplicates[allSpecies]];

convertReactionFormat[reactions_] := Module[{converted},
  If[Length[reactions] == 0, Return[{}]];
  If[Head[reactions[[1]]] === Rule,
    converted = Table[
      {First[reactions[[i]]], Last[reactions[[i]]]}, 
      {i, Length[reactions]}
    ];
    Return[converted],
    Return[reactions]
  ]
];

asoRea[RN_]:=Module[{parseSide,extractSpecies},
  extractSpecies[expr_]:=Module[{terms,species},
    If[expr===0,Return[{}]];
    terms=If[Head[expr]===Plus,List@@expr,{expr}];
    species={};
    Do[Which[
      StringQ[term]&&StringContainsQ[term," "],AppendTo[species,StringTrim[StringDrop[term,1]]],
      Head[term]===Times,AppendTo[species,Cases[term,_String][[1]]],
      StringQ[term],AppendTo[species,term]];,{term,terms}];
    DeleteDuplicates[species]];
  parseSide[expr_]:=extractSpecies[expr];
  Map[Function[r,Association["Substrates"->parseSide[r[[1]]],"Products"->parseSide[r[[2]]]]],RN]];

stoichiometricMatrices[reactions_] := Module[{
  species, numReactions, numSpecies, 
  alpha, beta, gamma, reactants, products, leftSide, rightSide},
  
  species = extSpe[reactions];
  numReactions = Length[reactions];
  numSpecies = Length[species];
  
  alpha = Table[0, {numSpecies}, {numReactions}];
  beta = Table[0, {numSpecies}, {numReactions}];
  
  Do[
    If[Head[reactions[[j]]] === Rule,
      leftSide = First[reactions[[j]]];
      rightSide = Last[reactions[[j]]],
      leftSide = reactions[[j, 1]];
      rightSide = reactions[[j, 2]]
    ];
    reactants = compToAsso[leftSide];
    products = compToAsso[rightSide];
    Do[
      If[KeyExistsQ[reactants, species[[i]]], 
        alpha[[i, j]] = reactants[species[[i]]]];
      If[KeyExistsQ[products, species[[i]]], 
        beta[[i, j]] = products[species[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  gamma = beta - alpha;
  {alpha, beta, gamma, species}
];

reaToRHS[reactions_] := Module[{alpha, beta, gamma, species, var, rv, tk, Rv, RHS, convertedReactions},
  convertedReactions = If[Length[reactions] > 0 && Head[reactions[[1]]] === Rule,
    Table[{First[reactions[[i]]], Last[reactions[[i]]]}, {i, Length[reactions]}],
    reactions
  ];
  {alpha, beta, gamma, species} = stoichiometricMatrices[convertedReactions];
  var = ToExpression[species];
  rv = expM[var, alpha // Transpose];
  tk = Array[Symbol["k" <> ToString[#]] &, alpha // Transpose // Length];
  Rv = tk*rv;
  RHS = gamma . Rv;
  {RHS, species, Rv}
];

expM=Inner[OperatorApplied[Power],#2,#1,Times]&;

End[];
EndPackage[];
(*Testing example-uncomment to test core package*)
(*RN={0->"S1","S1"+"I1"->2*"I1","I1"->0};
result=extMat[RN];
Print["Species: ",result[[1]]];
Print["Alpha matrix: ",result[[2]]//MatrixForm];
Print["Gamma matrix: ",result[[4]]//MatrixForm];*)
