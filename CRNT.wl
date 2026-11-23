(* ::Package:: *)

(*Print["Load CRNT"]; *)
(* Extract complexes and edges from reaction network *)

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

(* whichInv: for each invasion block in R0A, test whether it can invade (R0 >= 1).
   R0A from NGMRN may be: scalar, list of scalars (one per block), or matrix K.
   - Scalar r:     returns {{r>=1, rval}}
   - List {r1,..}: returns {{r1>=1, v1}, {r2>=1, v2}, ...} one entry per block
   - Matrix K:     returns {{lam_i>=1, lam_i}} for each eigenvalue, sorted descending
   NOTE: R0>=1 is the invasion condition (local instability of the face equilibrium),
   NOT a theorem for Kdt membership in general.  Use vbCheck for the Kdt question.
   kv: parameter substitution rules.
   Returns: list of {canInvade, R0value} for each block. *)
whichInv[R0A_, kv_List] := Module[{r0val},
  r0val = R0A /. kv;
  Which[
    MatrixQ[r0val], Sort[Re[Eigenvalues[N[r0val]]], Greater],
    ListQ[r0val],   N[r0val],
    True,           {N[r0val]}
  ]
];

(* vbCheck: test vertex-balancing of the flux vector at a specific point.
   For (G,k) in Kdt there must exist an equilibrium xeq>0 with vertex-balanced flux.
   This checks a SPECIFIC (k, xeq) point; find xeq numerically first.
   rts: rate vector from ODE2RN, symbolic monomials k_j * x^alp_j.
   RN:  reaction network from ODE2RN, list of {sourceComplex, targetComplex}.
   kxv: substitution rules for both rate constants k and equilibrium xeq.
   Returns: {isVB, defectVector} where isVB=True iff max|defect| < 1e-8. *)
vbCheck[rts_, RN_, kxv_List] := Module[
  {com, edg, flux, outFlux, inFlux, defect},
  {com, edg} = getComE[RN];
  flux = N[rts /. kxv];
  outFlux = Table[
    Total[flux[[Flatten[Position[edg[[All, 1]], com[[i]]]]]]],
    {i, Length[com]}
  ];
  inFlux = Table[
    Total[flux[[Flatten[Position[edg[[All, 2]], com[[i]]]]]]],
    {i, Length[com]}
  ];
  defect = outFlux - inFlux;
  {Max[Abs[defect]] < 10^-8, defect}
];

(* kdtPoint: check Kdt membership at a specific numerical parameter point.
   Calls ODE2RN[RHS,var] internally to get rts and RN.
   RHS: concrete ODE right-hand side.
   var: species list (same order as ODE2RN[RHS,var]).
   kv:  all-numerical substitution rules for rate constants.
   Returns {isKdt, xStar, defect}.
   isKdt=False with xStar={} means no positive equilibrium exists at this kv. *)
kdtPoint[RHS_, var_, kv_List] := Module[
  {RN, rts, posSol, xStar, vb},
  Block[{$Output = {}}, {RN, rts} = ODE2RN[RHS, var][[{1, 2}]]];
  posSol = Select[
    NSolve[Thread[(RHS /. kv) == 0], var, Reals],
    AllTrue[var /. #, Positive] &
  ];
  If[posSol === {},
    {False, {}, {}},
    xStar = posSol[[1]];
    vb    = vbCheck[rts, RN, Join[kv, xStar]];
    {vb[[1]], xStar, vb[[2]]}
  ]
];

(* === BCH E-graph functions (Boros-Craciun-Henriksson-Jin-RojasLaLuz 2025) ===
   Input: reactions = list of directed edges {src,tgt} where src,tgt are
   stoichiometry vectors, e.g. {2,0}->{0,1} for 2X->Y.
   Rate constants: Subscript[\[Kappa],i].  Flux variables: Subscript[\[Beta],i]. *)

GetMassAction[reactions_] := Module[
  {n, G, rates, flux, vars, srcs, tgts, monomials, f},
  G         = Graph[reactions];
  n         = Length[reactions[[1, 1]]];
  rates     = Table[Subscript[\[Kappa], i], {i, 1, Length[reactions]}];
  flux      = Table[Subscript[\[Beta],  i], {i, 1, Length[reactions]}];
  vars      = Switch[n, 2, {Global`x, Global`y}, 3, {Global`x, Global`y, Global`z},
                4, {Global`x, Global`y, Global`z, Global`w},
                _, Table[Subscript[Global`x, i], {i, 1, n}]];
  srcs      = EdgeList[G][[All, 1]];
  tgts      = EdgeList[G][[All, 2]];
  monomials = rates * Apply[Times, vars^(srcs\[Transpose])];
  f         = Total[monomials * (tgts - srcs)];
  {f, rates, G, monomials, vars, flux}
];

(* Vertex-balance equilibrium equations: substitute \[Kappa]->\[Beta], set all vars=1 *)
GetEQ[reactions_] := Module[
  {f, vars, all1},
  {f, vars} = GetMassAction[reactions][[{1, 5}]];
  all1 = Table[vars[[i]] -> 1, {i, 1, Length[vars]}];
  (f /. {\[Kappa] -> \[Beta]} /. all1) == 0
];

(* Flux symbol vector \[Beta]_i *)
GetFlux[reactions_] := GetMassAction[reactions][[6]];

(* Dynamical equivalence (DE) and vertex-balancing (VB) conditions for rxnG~rxnH *)
GetDEandVB[rxnG_, rxnH_] := Module[
  {f1, vars1, f2, rates2, G2, flux2, DE, VB},
  {f1, vars1}       = GetMassAction[rxnG][[{1, 5}]];
  {f2, rates2, G2}  = GetMassAction[rxnH][[{1, 2, 3}]];
  flux2 = rates2 /. {\[Kappa] -> \[Gamma]};
  DE    = CoefficientList[
            (f1 /. {\[Kappa] -> \[Beta]}) - (f2 /. {\[Kappa] -> \[Gamma]}),
            vars1] == 0;
  VB    = IncidenceMatrix[G2] . flux2 == 0;
  {DE, VB}
];

(* Map flux \[Beta]_i to mass-action monomials \[Kappa]_i x^\[Alpha]_i *)
Get\[Beta]2\[Kappa][reactions_] := Module[
  {flux, monomials, b2k, vars, ks, i},
  flux              = GetFlux[reactions];
  {ks, monomials, vars} = GetMassAction[reactions][[{2, 4, 5}]];
  b2k = Table[flux[[i]] -> monomials[[i]], {i, 1, Length[reactions]}];
  {b2k, ks, vars}
];

(* maxEg: G^max = edges in G^comp for which G has a realization w.r.t. some weakly
   reversible H ⊆ G^comp.  For each edge e, tries all strongly connected subgraphs of
   G^comp containing e (largest first for speed), checks DE(G,H) with kappa>0, lambda>0. *)
maxEg[Eg_] := Module[{verts, Gcomp, f1, k1, vars1},
  verts = VertexList[Graph[Eg]];
  Gcomp = Select[Flatten[Outer[Rule, verts, verts, 1], 1], (#[[1]] =!= #[[2]]) &];
  {f1, k1, vars1} = GetMassAction[Eg][[{1, 2, 5}]];
  Select[Gcomp, mxInMax[f1, k1, vars1, Gcomp, #] &]
];

(* mxInMax: e in G^max iff exists kappa>0 (for G), lambda>=0 lambda_e>0
   satisfying VB(Gcomp) and DE(G, H_lambda).  By flow decomposition, support of
   lambda with VB is weakly reversible, so no subset enumeration needed.
   Reduces from exponential (Subsets) to a single linear feasibility check. *)
mxInMax[f1_, k1_, vars1_, Gcomp_, e_] := Module[
  {f2, k2, la, lam, eqs, vbEqs, eIdx},
  {f2, k2} = GetMassAction[Gcomp][[{1, 2}]];
  lam  = Array[la, Length[k2]];
  f2   = f2 /. Thread[k2 -> lam];
  eqs  = Flatten[Table[
    Thread[CoefficientRules[f1[[i]] - f2[[i]], vars1][[All, 2]] == 0],
    {i, Length[f1]}]];
  vbEqs = Thread[IncidenceMatrix[Graph[Gcomp]] . lam == 0];
  eIdx  = First[Position[Gcomp, e, 1]][[1]];
  FindInstance[
    Join[eqs, vbEqs, Thread[lam >= 0], Thread[k1 > 0], {lam[[eIdx]] > 0}],
    Join[k1, lam]
  ] =!= {}
];

(* Draw E-graph with stoichiometry vectors as vertex coordinates *)
DrawEGraph[rxn_, opts___] := Graph[rxn,
  VertexCoordinates -> Table[i -> i, {i, VertexList[rxn]}],
  opts, ImageSize -> Tiny
];

(* comp2RN: convert BCH E-graph to EpidCRN {RN, rts} format.
   Eg:  list of directed edges {src,tgt} between stoichiometry vectors,
        e.g. {{2,0}->{0,1}, {0,1}->{1,0}} for a 2-species network.
   RN:  list of {srcMonomial, tgtMonomial} with species from GetMassAction vars.
   rts: rate monomials Subscript[\[Kappa],i] * x^srci  (BCH rate-constant convention).
   Compatible with vbCheck, getComE, kdtPoint (pass Subscript[\[Kappa],i]->val in kv). *)
Eg2RN[Eg_] := Module[
  {res, G, srcs, tgts, vars, monomials, comSrc, comTgt},
  res      = GetMassAction[Eg];
  {G, monomials, vars} = res[[{3, 4, 5}]];
  srcs     = EdgeList[G][[All, 1]];
  tgts     = EdgeList[G][[All, 2]];
  comSrc   = Apply[Times, vars^#] & /@ srcs;
  comTgt   = Apply[Times, vars^#] & /@ tgts;
  Print["#Eg=", Length[Eg], " Eg2RN=",
    Prepend[Transpose[{MapThread[Rule, {comSrc, comTgt}], monomials}], {"RN1", "rts1"}] // MatrixForm];
  {MapThread[List, {comSrc, comTgt}], monomials}
];

(* RN2Eg: inverse of Eg2RN.  alp, bet = stoichiometry matrices from ODE2RN.
   Each row of alp/bet is the stoichiometry vector of the source/target complex. *)
RN2Eg[alp_, bet_] := MapThread[Rule, {Transpose[alp], Transpose[bet]}];

