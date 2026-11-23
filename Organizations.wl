(* ::Package:: *)

(* Organizations.wl -- Chemical Organization Theory (COT) and Elementary
   Reaction Closures (ERCs) for the EpidCRN package.

   RN representation: same as everywhere else in EpidCRN -- a list of
   lhs -> rhs rules (or {lhs,rhs} pairs), lhs/rhs a symbolic sum of
   coeff*"Species" terms or literal 0 for inflow/outflow. Reuses
   comp2Asso, extSpe, extMat from Core.wl. Does not touch minSiph or NGM.

   For a reaction rho, React(rho) = reactant species (Keys of
   comp2Asso[lhs]); rho is ACTIVATED on a species set X iff
   React(rho) subseteq X. Note: an inflow reaction (lhs=0) has
   React(rho)={}, which is a subset of every X -- so inflow reactions are
   always activated, and the closure of ANY set (including {}) already
   contains every freely-supplied species. This is a real structural fact,
   not a bug: open (demographic-inflow) networks have no "empty"
   organization. *)

actRea[RN_List, spe_List] := Select[RN, SubsetQ[spe, Keys[comp2Asso[#[[1]]]]] &];

speClo[RN_List, spe0_List] := Module[{cur = Union[spe0], active, prod, nxt},
  While[True,
    active = actRea[RN, cur];
    prod = Union @@ (Keys[comp2Asso[#[[2]]]] & /@ active);
    nxt = Union[cur, prod];
    If[nxt === cur, Break[]];
    cur = nxt;
  ];
  cur
];

isClo[RN_List, spe_List] := speClo[RN, spe] === Union[spe];

ercClo[RN_List, rxn_] := speClo[RN, Keys[comp2Asso[rxn[[1]]]]];

allErc[RN_List] := DeleteDuplicatesBy[{#, ercClo[RN, #]} & /@ RN, Last];

isReac[RN_List, spe_List] := isClo[RN, spe] && Length[actRea[RN, spe]] > 0;

(* SelfMaintainingQ: LP feasibility, v_rho>=1 for every activated reaction
   (strict positivity normalized to >=1, cf. task normalization), net
   production gamma.v >= 0 on the species rows of spe (Nc: rho epsilon
   scaling makes >=1 equivalent to >0 up to positive rescaling). *)
isSM[RN_List, spe_List] := Module[
  {spAll, alp, bet, gamma, active, colIdx, rowIdx, vars, cons, sol},
  {spAll, alp, bet, gamma} = extMat[RN][[1 ;; 4]];
  active = actRea[RN, spe];
  If[active === {}, Return[False]];
  colIdx = Flatten[Position[RN, #, {1}, 1] & /@ active];
  rowIdx = Flatten[Position[spAll, #, {1}, 1] & /@ spe];
  vars = Array[v, Length[active]];
  cons = Join[Thread[vars >= 1], Thread[(gamma[[rowIdx, colIdx]] . vars) >= 0]];
  sol = FindInstance[cons, vars, Reals];
  sol =!= {}
];

(* SemiSelfMaintainingQ: weaker LP, v_rho>=0 (not all zero), gamma.v>=0. *)
isSSM[RN_List, spe_List] := Module[
  {spAll, alp, bet, gamma, active, colIdx, rowIdx, vars, cons, sol},
  {spAll, alp, bet, gamma} = extMat[RN][[1 ;; 4]];
  active = actRea[RN, spe];
  If[active === {}, Return[False]];
  colIdx = Flatten[Position[RN, #, {1}, 1] & /@ active];
  rowIdx = Flatten[Position[spAll, #, {1}, 1] & /@ spe];
  vars = Array[v, Length[active]];
  cons = Join[Thread[vars >= 0], {Total[vars] > 0}, Thread[(gamma[[rowIdx, colIdx]] . vars) >= 0]];
  sol = FindInstance[cons, vars, Reals];
  sol =!= {}
];

isOrg[RN_List, spe_List] := isClo[RN, spe] && isSM[RN, spe];

Options[allOrg] = {Method -> "Exhaustive"};
allOrg[RN_List, OptionsPattern[]] := Module[{spAll, cands, ercs},
  spAll = extSpe[RN];
  cands = If[OptionValue[Method] === "ERCGenerated",
    ercs = allErc[RN][[All, 2]];
    DeleteDuplicates[speClo[RN, Union @@ #] & /@ Subsets[ercs, {1, Length[ercs]}]],
    Union /@ Subsets[spAll]
  ];
  Select[cands, isOrg[RN, #] &]
];

(* ERCSynergyQ[RN,ercList] -> {synQ,witnessRxns,joinClosure}. synQ=True iff
   some reaction is activated by the join of ercList's species but by none
   of the individual ERCs. *)
synQ[RN_List, ercs_List] := Module[{joinClo, actJoin, actEach, wit},
  joinClo = speClo[RN, Union @@ ercs];
  actJoin = actRea[RN, joinClo];
  actEach = Union @@ (actRea[RN, #] & /@ ercs);
  wit = Complement[actJoin, actEach];
  {wit =!= {}, wit, joinClo}
];

(* ERCComplementarity[RN,e1,e2] -> {cmpQ,supply12,supply21}. supply12 =
   products of R_e1 (all reactions activated by e1) intersected with
   req(e2) = reactant species of e2's own generating reaction(s), found by
   looking e2 up in allErc; symmetric for supply21. *)
cmpERC[RN_List, e1_List, e2_List] := Module[
  {gens1, gens2, req1, req2, prod1, prod2, supply12, supply21},
  gens1 = Cases[allErc[RN], {rxn_, x_} /; x === e1 :> rxn];
  gens2 = Cases[allErc[RN], {rxn_, x_} /; x === e2 :> rxn];
  req1 = Union @@ (Keys[comp2Asso[#[[1]]]] & /@ gens1);
  req2 = Union @@ (Keys[comp2Asso[#[[1]]]] & /@ gens2);
  prod1 = Union @@ (Keys[comp2Asso[#[[2]]]] & /@ actRea[RN, e1]);
  prod2 = Union @@ (Keys[comp2Asso[#[[2]]]] & /@ actRea[RN, e2]);
  supply12 = Intersection[prod1, req2];
  supply21 = Intersection[prod2, req1];
  {supply12 =!= {} || supply21 =!= {}, supply12, supply21}
];

(* Test: bare two-strain shared-susceptible network (Example A),
0->S,S->0,S+Ik->2Ik,Ik->0 (k=1,2). Expect: 3 distinct ERCs {S},{S,I1},{S,I2};
4 closed sets = 4 organizations (all of them, unlimited inflow); coinfection
network (Example C) adds I1+I2->I12, expect synQ True with witness that
reaction.

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

RNa = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> 0, "S" + "I2" -> 2*"I2", "I2" -> 0};
ercA = allErc[RNa][[All, 2]];
orgA = allOrg[RNa];
Print["ercA=", ercA, " orgA=", orgA, " nOrg=", Length[orgA]];

RNc = Join[RNa, {"I1" + "I2" -> "I12", "I12" -> 0}];
e1 = speClo[RNc, {"S", "I1"}];
e2 = speClo[RNc, {"S", "I2"}];
{synT, synW, synJ} = synQ[RNc, {e1, e2}];
Print["e1=", e1, " e2=", e2, " synT=", synT, " synW=", synW];
*)
