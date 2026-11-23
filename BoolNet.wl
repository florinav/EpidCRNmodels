(* ::Package:: *)

(* ========================================================================= *)
(* BoolNet.wl - Boolean Network Analysis and EpidCRN Interface               *)
(* ========================================================================= *)
(* Converts Boolean regulatory networks to polynomial/Hill ODEs whose        *)
(* {var, RHS} output plugs directly into ODE2RN / NGMRN.                    *)
(* Sources: BoolDog (BoolCube/SQUAD) and BoolODE (two-layer Hill kinetics). *)

ClearAll[mkNet, readBnet, readBoolODE, evalNet, steadyBool, boolSTG,
         boolMultiRHS, boolHill2LRHS];

(* ========================================================================= *)
(* CORE DATA STRUCTURE                                                        *)
(* ========================================================================= *)

(* Build network from string gene names and Mathematica Boolean rule exprs.
   Regulators are inferred automatically from each rule. *)
mkNet[genes_List, rules_List] := Module[{regs},
  regs = Table[Select[genes, !FreeQ[rules[[i]], Symbol[#]] &], {i, Length[genes]}];
  <|"genes" -> genes,
    "rules" -> AssociationThread[genes -> rules],
    "regs"  -> AssociationThread[genes -> regs]|>
];

(* ========================================================================= *)
(* FILE I/O                                                                   *)
(* ========================================================================= *)

(* Convert Bnet-format rule string to Mathematica Boolean expression.
   Handles & -> &&, | -> ||, 0/1 literals, preserves existing && || !  *)
parseBnetRule["0"] := False;
parseBnetRule["1"] := True;
parseBnetRule[s_String] := ToExpression[
  StringReplace[StringTrim[s], {"&&" -> "&&", "||" -> "||",
                                 "&" -> "&&", "|" -> "||"}]
];

getRegs[ruleExpr_, genes_List] := Select[genes, !FreeQ[ruleExpr, Symbol[#]] &];

(* Read Bnet format: non-comment lines are "GeneName, BoolRule" *)
readBnet[file_String] := Module[{lines, pairs, genes, ruleExprs, regs},
  lines = Select[Import[file, "Lines"],
    StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &];
  pairs = StringTrim /@ StringSplit[#, ",", 2] & /@ lines;
  genes = First /@ pairs;
  ruleExprs = parseBnetRule /@ (Last /@ pairs);
  regs = Table[getRegs[ruleExprs[[i]], genes], {i, Length[genes]}];
  <|"genes" -> genes,
    "rules" -> AssociationThread[genes -> ruleExprs],
    "regs"  -> AssociationThread[genes -> regs]|>
];

(* Read BoolODE format: tab-separated with header line "Gene\tRule" *)
readBoolODE[file_String] := Module[{lines, data, genes, ruleExprs, regs},
  lines = Select[Import[file, "Lines"], StringLength[StringTrim[#]] > 0 &];
  data  = StringTrim /@ StringSplit[#, "\t"] & /@ Rest[lines];
  genes = First /@ data;
  ruleExprs = parseBnetRule /@ (Last /@ data);
  regs = Table[getRegs[ruleExprs[[i]], genes], {i, Length[genes]}];
  <|"genes" -> genes,
    "rules" -> AssociationThread[genes -> ruleExprs],
    "regs"  -> AssociationThread[genes -> regs]|>
];

(* ========================================================================= *)
(* BOOLEAN EVALUATION                                                         *)
(* ========================================================================= *)

(* Synchronous update: evaluate all rules at state (0/1 list, gene order) *)
evalNet[net_Association, state_List] := Module[{genes, rules, sub},
  genes = net["genes"];
  rules = net["rules"];
  sub = Thread[Symbol /@ genes -> (state /. {1 -> True, 0 -> False})];
  Table[Boole[rules[genes[[i]]] /. sub], {i, Length[genes]}]
];

(* Boolean fixed points: states s where synchronous update maps s -> s *)
steadyBool[net_Association] := Module[{n, allSt},
  n = Length[net["genes"]];
  allSt = Tuples[{0, 1}, n];
  Select[allSt, evalNet[net, #] === # &]
];

(* Synchronous state transition graph; vertices labeled as binary integers *)
boolSTG[net_Association] := Module[{n, allSt, lbl, edges},
  n = Length[net["genes"]];
  allSt = Tuples[{0, 1}, n];
  lbl = FromDigits[#, 2] &;
  edges = DirectedEdge[lbl[#], lbl[evalNet[net, #]]] & /@ allSt;
  Graph[edges, VertexLabels -> "Name"]
];

(* ========================================================================= *)
(* ODE TRANSFORMS — output {var, RHS} for EpidCRN pipeline                  *)
(* ========================================================================= *)

(* Multilinear (BoolCube) extension of Boolean rule over its regulators.
   Returns polynomial exact at {0,1} corners; regNames are gene-name strings. *)
multilinExt[rule_, regNames_List] := Module[{n, allV, regVars, boolVals, wts},
  n = Length[regNames];
  If[n == 0, Return[Boole[rule]]];
  allV     = Tuples[{0, 1}, n];
  regVars  = Symbol /@ regNames;
  boolVals = Table[
    Boole[rule /. Thread[Symbol /@ regNames -> (v /. {1 -> True, 0 -> False})]],
    {v, allV}
  ];
  wts = Table[
    Times @@ Table[If[allV[[j, i]] == 1, regVars[[i]], 1 - regVars[[i]]], {i, n}],
    {j, Length[allV]}
  ];
  Expand[boolVals . wts]
];

(* Multilinear ODE: dxi/dt = gam*(f_i^multi(x) - xi).
   Returns {var, RHS}; plugs directly into ODE2RN[RHS, var].
   Variables are Symbol/@genes; RHS is a polynomial list. *)
boolMultiRHS[net_Association, gam_: 1] := Module[{genes, rules, regs, vars, rhs},
  genes = net["genes"];
  rules = net["rules"];
  regs  = net["regs"];
  vars  = Symbol /@ genes;
  rhs   = Table[
    gam * (multilinExt[rules[genes[[i]]], regs[genes[[i]]]] - vars[[i]]),
    {i, Length[genes]}
  ];
  {vars, rhs}
];

(* Two-layer Hill ODE (BoolODE style): mRNA xG and protein pG per gene.
   dxG/dt = m*fHill(pRegs) - lx*xG,  dpG/dt = r*xG - dp*pG.
   fHill uses all 2^k regulator combinations; basal term alpha0 handled separately.
   pars: Association overriding defaults m=1, lx=0.5, r=2, dp=0.8, k=0.5, n=2.
   Returns {vars2L, RHS2L} with variables ordered {xA,pA,xB,pB,...}.
   Plugs into ODE2RN[RHS2L, vars2L] for CRN/NGM analysis. *)
boolHill2LRHS[net_Association, pars_Association : <||>] := Module[{
  genes, rules, regs, defs, m, lx, r, dp, k, nh,
  vars2L, rhs2L, g, rns, pvs, nreg, allV, nzV, bv, ht, alpha0, fH, xg, pg},
  genes = net["genes"];
  rules = net["rules"];
  regs  = net["regs"];
  defs  = Join[<|"m" -> 1., "lx" -> 0.5, "r" -> 2., "dp" -> 0.8, "k" -> 0.5, "n" -> 2|>, pars];
  {m, lx, r, dp, k, nh} = Lookup[defs, {"m", "lx", "r", "dp", "k", "n"}];
  vars2L = Flatten[Table[
    {Symbol["x" <> genes[[i]]], Symbol["p" <> genes[[i]]]},
    {i, Length[genes]}
  ]];
  rhs2L = Flatten[Table[
    g   = genes[[i]];
    xg  = Symbol["x" <> g];
    pg  = Symbol["p" <> g];
    rns = regs[g];
    pvs = Symbol["p" <> #] & /@ rns;
    nreg = Length[rns];
    fH = If[nreg == 0,
      N[Boole[rules[g]]],
      allV   = Tuples[{0, 1}, nreg];
      nzV    = Rest[allV];   (* exclude all-zeros combo *)
      bv     = Table[
        Boole[rules[g] /. Thread[Symbol /@ rns -> (v /. {1 -> True, 0 -> False})]],
        {v, nzV}
      ];
      ht = Table[
        Times @@ MapThread[If[#1 == 1, (#2/k)^nh, 1] &, {nzV[[j]], pvs}],
        {j, Length[nzV]}
      ];
      alpha0 = Boole[rules[g] /. Thread[Symbol /@ rns -> Table[False, nreg]]];
      (alpha0 + bv . ht) / (1 + Total[ht])
    ];
    {m * fH - lx * xg, r * xg - dp * pg},
    {i, Length[genes]}
  ]];
  {vars2L, rhs2L}
];

(* ========================================================================= *)
(* TEST: toggle switch A = !B, B = !A                                        *)
(* Expected: fps={{0,1},{1,0}}, rhs={1-A-B, 1-A-B}                          *)
(* ========================================================================= *)

(*
Module[{net, fps, vars, rhs, v2, r2},
  net = mkNet[{"A", "B"}, {!Symbol["B"], !Symbol["A"]}];
  fps = steadyBool[net];
  {vars, rhs} = boolMultiRHS[net, 1];
  {v2, r2}    = boolHill2LRHS[net];
  Print["gen=", net["genes"], " fps=", fps, " rhs=", rhs];
  Print["v2=", v2, " r2nz=", r2 =!= {}, " n2L=", Length[v2]]
];
*)
