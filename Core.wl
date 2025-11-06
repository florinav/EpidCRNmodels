(* ::Package:: *)

BeginPackage["EpidCRN1`"];
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

(* Backward compatibility alias *)
comp2Asso = compToAsso;

(* Main extMat function *)
extMat[reactions_] := Module[{
  spe, al, be, gamma, Rv, RHS, 
  numReactions, numSpecies, reactants, products, 
  var, rv, tk, complexes, linkageClasses, deficiency,
  defFormula, defTerms, defResult, Nc, l, s, leftSide, rightSide,
  compToAssoLocal, expMLocal},
  
  (* Local helper function to convert compound expressions to associations *)
  compToAssoLocal = Function[{expr},
    Module[{terms, result},
      If[expr === 0 || expr === Null,
        Association[], (* Empty association for zero or null *)
        terms = If[Head[expr] === Plus, List @@ expr, {expr}];
        result = Association[];
        Do[
          Which[
            Head[term] === Times && Length[term] >= 2 && NumericQ[First[term]],
            (* Handle cases like 2*A *)
            result[Last[term]] = First[term],
            Head[term] === Times,
            (* Handle cases like A*B (coefficient 1) *)
            result[term] = 1,
            True,
            (* Handle single species *)
            result[term] = 1
          ];
        , {term, terms}];
        result
      ]
    ]
  ];
  
  (* Local helper function for exponential matrix *)
  expMLocal = Function[{vars, matrix},
    Module[{},
      If[Length[vars] == 0 || Length[matrix] == 0,
        ConstantArray[1, Length[matrix]],
        (* Create the exponential terms *)
        Table[
          If[Length[vars] > 0 && Length[matrix[[j]]] > 0,
            Product[
              If[matrix[[j, i]] == 0, 1, vars[[i]]^matrix[[j, i]]], 
              {i, Length[vars]}
            ],
            1
          ],
          {j, Length[matrix]}
        ]
      ]
    ]
  ];
  
  (* Extract species *)
  spe = Module[{allSpecies = {}}, 
    Do[
      (* Handle different input formats *)
      {leftSide, rightSide} = Which[
        Head[reactions[[i]]] === Rule, 
        {reactions[[i, 1]], reactions[[i, 2]]},
        Head[reactions[[i]]] === List && Length[reactions[[i]]] >= 2,
        {reactions[[i, 1]], reactions[[i, 2]]},
        True,
        (Print["Warning: Unrecognized reaction format at position ", i]; 
         Continue[])
      ];
      
      reactants = compToAssoLocal[leftSide];
      products = compToAssoLocal[rightSide];
      allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
    , {i, Length[reactions]}];
    DeleteDuplicates[allSpecies]
  ];
  
  (* Check if species were found *)
  If[Length[spe] == 0,
    Print["Error: No species found in reactions"];
    Return[$Failed]
  ];
  
  numSpecies = Length[spe];
  numReactions = Length[reactions];
  
  (* Initialize stoichiometric matrices *)
  al = ConstantArray[0, {numSpecies, numReactions}];
  be = ConstantArray[0, {numSpecies, numReactions}];
  
  (* Build stoichiometric matrices *)
  Do[
    (* Handle different input formats *)
    {leftSide, rightSide} = Which[
      Head[reactions[[j]]] === Rule,
      {reactions[[j, 1]], reactions[[j, 2]]},
      Head[reactions[[j]]] === List && Length[reactions[[j]]] >= 2,
      {reactions[[j, 1]], reactions[[j, 2]]},
      True,
      (Print["Warning: Unrecognized reaction format at position ", j]; 
       Continue[])
    ];
    
    reactants = compToAssoLocal[leftSide];
    products = compToAssoLocal[rightSide];
    
    (* Fill in the matrices *)
    Do[
      If[KeyExistsQ[reactants, spe[[i]]], 
        al[[i, j]] = reactants[spe[[i]]]];
      If[KeyExistsQ[products, spe[[i]]], 
        be[[i, j]] = products[spe[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  (* Calculate derived matrices and terms *)
  gamma = be - al;
  
  (* Convert species to variables more safely *)
  var = Table[ToExpression[ToString[spe[[i]]]], {i, Length[spe]}];
  
  (* Calculate reaction rates *)
  rv = expMLocal[var, Transpose[al]];
  tk = Table[Symbol["k" <> ToString[j]], {j, numReactions}];
  Rv = tk*rv;
  RHS = gamma . Rv;
  
  (* Calculate complexes for deficiency *)
  complexes = Module[{complexList = {}},
    Do[
      {leftSide, rightSide} = Which[
        Head[reactions[[i]]] === Rule,
        {reactions[[i, 1]], reactions[[i, 2]]},
        Head[reactions[[i]]] === List && Length[reactions[[i]]] >= 2,
        {reactions[[i, 1]], reactions[[i, 2]]},
        True,
        Continue[]
      ];
      
      If[!MemberQ[complexList, leftSide], 
        AppendTo[complexList, leftSide]];
      If[!MemberQ[complexList, rightSide], 
        AppendTo[complexList, rightSide]];
    , {i, Length[reactions]}];
    complexList
  ];
  
  (* Deficiency calculation *)
  Nc = Length[complexes];
  l = 1; (* Simplified linkage class calculation *)
  s = Which[
    gamma === {}, 0,
    AllTrue[Flatten[gamma], # == 0 &], 0,
    True, MatrixRank[gamma]
  ];
  
  (* Format deficiency results *)
  defFormula = "\[Delta] = Nc - \[ScriptL] - s";
  defTerms = "Nc = " <> ToString[Nc] <> " (complexes), " <>
             "\[ScriptL] = " <> ToString[l] <> " (linkage classes), " <> 
             "s = " <> ToString[s] <> " (stoich dimension)";
  defResult = "\[Delta] = " <> ToString[Nc] <> " - " <> ToString[l] <> 
              " - " <> ToString[s] <> " = " <> ToString[Nc - l - s];
  
  (* Return results *)
  {spe, al, be, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
];










(* Other core functions *)
extSpe[reactions_] := Module[{allSpecies, reactants, products, leftSide, rightSide},
  allSpecies = {};
  Do[
   (* Handle both Rule and List formats *)
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
  DeleteDuplicates[allSpecies]];

arrow2pairReac[reactions_] := Module[{converted},
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
(*Basic parameter extraction*)
Par[RHS_,X_]:=Complement[Variables[RHS],X];
(*Par[{a*x+b*y,c*x+d*y},{x,y}] returns {a,b,c,d}*)
End[];
EndPackage[];
(*Testing example-uncomment to test core package*)
(*RN={0->"S1","S1"+"I1"->2*"I1","I1"->0};
result=extMat[RN];
Print["Species: ",result[[1]]];
Print["Alpha matrix: ",result[[2]]//MatrixForm];
Print["Gamma matrix: ",result[[4]]//MatrixForm];*)
