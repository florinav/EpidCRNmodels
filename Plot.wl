(* ::Package:: *)

(* Plot.wl - CRN graph visualizations *)
(* Three graphs: crnGr (bipartite), intGr (interaction), qGr (quotient) *)
(* Port of https://github.com/tambysatya/crnplot to Mathematica/EpidCRN *)

(* merge directed pairs a->b, b->a into UndirectedEdge; keep self-loops *)
mergDir[edges_List] := Module[{de, done, result},
  de = DeleteDuplicates[edges];
  done = {}; result = {};
  Do[
    With[{a = e[[1]], b = e[[2]]},
      Which[
        MemberQ[done, {a, b}], Null,
        a === b,
          AppendTo[result, DirectedEdge[a, a]],
        MemberQ[de, DirectedEdge[b, a]],
          AppendTo[result, UndirectedEdge[a, b]];
          AppendTo[done, {a, b}]; AppendTo[done, {b, a}],
        True,
          AppendTo[result, DirectedEdge[a, b]]
      ]
    ],
    {e, de}
  ];
  result
];

(* directed edges for one reaction: self-loops on reactants, r->p, co-reactant pairs *)
rnEdges1[rcts_List, prds_List] := Join[
  DirectedEdge[#, #] & /@ rcts,
  Flatten[Outer[DirectedEdge, rcts, prds], 1],
  Flatten[Table[
    {DirectedEdge[rcts[[j]], rcts[[k]]], DirectedEdge[rcts[[k]], rcts[[j]]]},
    {j, Length[rcts]}, {k, j + 1, Length[rcts]}
  ]]
];

(* 1. bipartite species-reaction graph *)
crnGr[RN_List] := Module[{spe, rLab, edges, vSty},
  spe  = extSpe[RN];
  rLab = Table["R" <> ToString[i], {i, Length[RN]}];
  edges = Flatten[Table[
    Join[
      DirectedEdge[#, rLab[[i]]] & /@ Keys[comp2Asso[RN[[i, 1]]]],
      DirectedEdge[rLab[[i]], #] & /@ Keys[comp2Asso[RN[[i, 2]]]]
    ],
    {i, Length[RN]}
  ], 1];
  vSty = Join[
    Thread[spe  -> Directive[LightBlue,   EdgeForm[Gray]]],
    Thread[rLab -> Directive[LightOrange, EdgeForm[Gray]]]
  ];
  Graph[edges, VertexLabels -> "Name", VertexStyle -> vSty,
    ImageSize -> Medium, GraphLayout -> "LayeredDigraphEmbedding"]
];

(* 2. interaction graph: species only *)
intGr[RN_List] := Module[{spe, allDE},
  spe   = extSpe[RN];
  allDE = Flatten[Table[
    rnEdges1[Keys[comp2Asso[RN[[i, 1]]]], Keys[comp2Asso[RN[[i, 2]]]]],
    {i, Length[RN]}
  ], 1];
  Graph[mergDir[allDE], spe,
    VertexLabels -> "Name",
    VertexStyle -> Directive[LightBlue, EdgeForm[Gray]],
    ImageSize -> Medium]
];

(* 3. quotient graph: groups = {"G1"->{"s1","s2"}, "G2"->{"s3"}, ...} optional *)
(* without groups degenerates to intGr *)
qGr[RN_List, groups_List : {}] := Module[{spe, ungrp, grMap, allDE},
  spe   = extSpe[RN];
  ungrp = Complement[spe, If[groups === {}, {}, Flatten[groups[[All, 2]]]]];
  grMap = Association[Join[
    If[groups === {}, {},
      Flatten[Table[# -> groups[[i, 1]] & /@ groups[[i, 2]], {i, Length[groups]}]]
    ],
    Thread[ungrp -> ungrp]
  ]];
  allDE = Flatten[Table[
    rnEdges1[
      grMap /@ Keys[comp2Asso[RN[[i, 1]]]],
      grMap /@ Keys[comp2Asso[RN[[i, 2]]]]
    ],
    {i, Length[RN]}
  ], 1];
  Graph[mergDir[allDE],
    VertexLabels -> "Name",
    VertexStyle -> Directive[LightGreen, EdgeForm[Gray]],
    ImageSize -> Medium]
];

(* Test *)
(*
Module[{RNplt, g1, g2, g3},
  RNplt = {"S" + "I" -> 2*"I", "I" -> "R"};
  {g1, g2, g3} = {crnGr[RNplt], intGr[RNplt], qGr[RNplt]};
  Print[Row[{g1, g2, g3}, Spacer[10]]]
];
*)
