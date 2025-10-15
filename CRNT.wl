(* ::Package::*)

(* Extract complexes and edges from reaction network *)
Print["Load CRNT"]; 

getComE[RN_List] := Module[{complexes, edges},
  complexes = {};
  edges = {};
  Do[
    Module[{left, right},
      left = RN[[i, 1]];
      right = RN[[i, 2]];
      If[! MemberQ[complexes, left], AppendTo[complexes, left]];
      If[! MemberQ[complexes, right], AppendTo[complexes, right]];
      AppendTo[edges, {left, right}];
    ], 
    {i, Length[RN]}
  ];
  {complexes, edges}
];

(* 
Example usage:
getComE[{{a,b},{b,c},{c,a}}]
Returns: {{a, b, c}, {{a, b}, {b, c}, {c, a}}}
*)

(* Incidence matrix analysis for FHJ graphs *)
IaFHJ[vert_, edg_] := Module[{gg, oU, taF},
  gg[a_, b_] := Which[
    a === b[[1]], -1, 
    a === b[[2]], 1, 
    True, 0
  ];
  oU = Table[
    gg[vert[[i]], edg[[j]]], 
    {i, Length[vert]}, 
    {j, Length[edg]}
  ];
  taF = TableForm[
    oU, 
    TableHeadings -> {vert, edg}, 
    TableAlignments -> {Right, Top}
  ];
  {oU, taF}
];

(* 
Example usage:
IaFHJ[{a, b, c}, {{a, b}, {b, c}, {c, a}}]
Returns incidence matrix and formatted table for vertices {a,b,c} and edges {{a,b},{b,c},{c,a}}
*)

(* Ik matrix computation for FHJ analysis *)
IkFHJ[vert_, edg_, tk_] := Module[{tri, gg, oU},
  tri = MapThread[Append, {edg, tk}];
  gg[a_, b_] := Which[
    a === b[[1]], b[[3]], 
    a === b[[2]], 0, 
    True, 0
  ];
  oU = Table[
    gg[vert[[i]], tri[[j]]], 
    {i, Length[vert]}, 
    {j, Length[tri]}
  ] // Transpose
];

(* 
Example usage:
IkFHJ[{a, b, c}, {{a, b}, {b, c}, {c, a}}, {k1, k2, k3}]
Returns Ik matrix with rate constants k1, k2, k3 applied to edges
*)

(* Species-complex incidence matrix *)
SpeComInc[spec_, comp_] := Coefficient[#, spec] & /@ comp;

(* 
Example usage:
SpeComInc[{x, y}, {x + y, 2 x, y}]
Returns: {{1, 1}, {2, 0}, {0, 1}}
*)

(* Main Laplacian computation function *)
lapK[RN_, rates_] := Module[{complexes, edges, laplacian},
  {complexes, edges} = getComE[RN];
  laplacian = IkFHJ[complexes, edges, rates];
  laplacian
];

(* 
Example usage:
lapK[{{a, b}, {b, c}, {c, a}}, {k1, k2, k3}]
Returns the Laplacian matrix for the reaction network with given rate constants
*)
