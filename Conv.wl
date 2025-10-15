(* Simple algebraic transformations *)
toSum = (# /. Times -> Plus) &;
(* toSum[a*b*c] returns a + b + c *)

toProd = (# /. Plus -> Times) &;
(* toProd[a + b + c] returns a*b*c *)

l2L = # &; (* Identity function placeholder *)
m2toM = # &; (* Identity function placeholder *)

(* List filtering utilities *)
remZ[li_] := Select[li, # =!= 0 &];
(* remZ[{1, 0, 3, 0, 5}] returns {1, 3, 5} *)

selZR[con_] := Select[con, MatchQ[#, Rule[_, 0]] &];
(* selZR[{a -> 1, b -> 0, c -> 3}] returns {b -> 0} *)

seZF[expr_] := Select[expr, FreeQ[#, 0] &];
(* seZF[{{1, 2}, {0, 3}, {4, 5}}] returns {{1, 2}, {4, 5}} *)

(* String conversion utilities *)
strEdg[edges_List] := Map[ToString, edges, {2}];
(* strEdg[{{A, B}, {B, C}}] returns {{"A", "B"}, {"B", "C"}} *)

rul2Str[rules_] := ToString[rules];
(* rul2Str[{a -> 1, b -> 2}] returns "{a -> 1, b -> 2}" *)

(* Mathematical analysis utilities *)
rtS[RHS_List] := DeleteDuplicates[Flatten[MonomialList /@ Expand[RHS]] /.
-1*x_ :> x];
(* rtS[{a*x + b*y, -a*x + c}] returns {a*x, b*y, c} *)

countMS[m_] := m // Together // NumeratorDenominator //
Map@CoefficientArrays //
   ReplaceAll[sa_SparseArray :> sa["NonzeroValues"]] // Flatten //
   Count[#, _?Negative] &;
(* countMS[{1, -2, 3, -4}] returns 2 *)

onlyP[m_] := m // Together // NumeratorDenominator //
Map@CoefficientArrays //
   ReplaceAll[sa_SparseArray :> sa["NonzeroValues"]] // Flatten //
   AllTrue[#, NonNegative] &;
(* onlyP[x + y + z] returns True, onlyP[x - y + z] returns False *)

