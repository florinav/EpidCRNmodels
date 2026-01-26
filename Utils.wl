(* ::Package:: *)

posM= 
Replace[#,{_?Negative->0,e_:>Replace[Expand[e],
{Times[_?Negative,_]->0,
t_Plus:>Replace[t,_?Negative|Times[_?Negative,_]->0,1]}]},{2}]&;
(*posL = Replace[#, {
  _?Negative -> 0,
  e_ :> Replace[Expand[e], {
    Times[_?Negative, _] -> 0,
    t_Plus :> Replace[t, _?Negative | Times[_?Negative, _] -> 0, 1]
  }]
}, {1}]&;*)


(*reCL[re_] :=DeleteCases[re, _Symbol > 0 | Subscript[_, __] > 0, Infinity];*)
remZ[li_]:=Select[li, # =!= 0 &];

remN[expr_] := Module[{withZeros},
  withZeros = expr //. {
    -Sqrt[x_] :> 0,
    Times[-1, Sqrt[x_], rest___] :> 0,
    Times[coef_?Negative, Sqrt[x_], rest___] :> 0
  };
  If[ListQ[withZeros], remZ[withZeros], withZeros]
];

redTim[con_, var_:{}, time_:300] := TimeConstrained[Reduce[con, var], time, "timeOut"];

whenP=Function[exprs,
  Module[{params,reducedConditions},
    params=Variables[exprs];
    reducedConditions=Reduce[And @@ (#>0 & /@ exprs)&&And @@ (#>0 & /@ params),params,Reals];
    Simplify[reducedConditions,Assumptions->And @@ (#>0 & /@ params)]
  ]
];

(*whenP[{\[Mu]3/\[Alpha]3,r (K \[Alpha]3-\[Mu]3)/(K \[Alpha]3^2)}]*)

Hur4M[mat_] := Module[{coefs, det, checks},
  coefs = CharacteristicPolynomial[mat, \[Lambda]] // CoefficientList[#, \[Lambda]] &;
  det = Det[mat];
  checks = Positive @@@ {coefs, {det}};
  {coefs, det, checks}
];
(* Hur4M[{{1, 2}, {3, 4}}] *)

makeLPM[mat_] := Module[{minors},
  minors = Table[Det[Take[Take[mat, i], All, i]], {i, Length[mat]}];
  minors
];
(* makeLPM[{{2, 1}, {1, 3}}] *)

GetVec[mat_] := Module[{vals, vecs},
  {vals, vecs} = Eigensystem[mat];
  Transpose[{vals, vecs}]
];
(* GetVec[{{0, 1}, {-2, 0}}] *)

Deg[poly_, var_] := Exponent[poly, var];
(* Deg[x^3 + 2 x^2 + 1, x] *)

Grobpol[RHS_, var_, par_, ind_, cn_ : {}] := Module[{dyn, X, eq, elim, pol, ratsub},
  dyn = RHS;
  X = var;
  eq = Thread[dyn == 0];
  elim = Complement[Range[Length[X]], ind];
  pol = Collect[GroebnerBasis[Numerator[Together[dyn /. cn]],
    Join[par, X[[ind]]], X[[elim]]], X[[ind]]];
  { pol}
]

(* Tests for onlyP (from Conv.wl)

onlyP[a + b]
onlyP[1/(a*b)]
onlyP[-a]
onlyP[a - b]
onlyP[0]
*)
