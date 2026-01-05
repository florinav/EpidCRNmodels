(* Format function to convert abbreviated notation to LaTeX subscripts *)

prF[expr_] := expr /. {
  (* Variables *)
  m0 -> Subscript[m, 0],
  m1 -> Subscript[m, 1],
  m2 -> Subscript[m, 2],
  mc -> Subscript[m, c],
  M1 -> Subscript[M, 1],
  M2 -> Subscript[M, 2],
  (* Parameters - Greek letters with subscripts *)
  Ph -> \[CapitalPhi],
  nu0 -> Subscript[\[Nu], 0],
  nu1 -> Subscript[\[Nu], 1],
  nu2 -> Subscript[\[Nu], 2],
  nuc -> Subscript[\[Nu], c],
  up1 -> Subscript[\[Upsilon], 1],
  up2 -> Subscript[\[Upsilon], 2],
  la1 -> Subscript[\[Lambda], 1],
  la2 -> Subscript[\[Lambda], 2],
  la1c -> Subscript[\[Lambda], "1,c"],
  la2c -> Subscript[\[Lambda], "2,c"],
  La1 -> Subscript[\[CapitalLambda], 1],
  La2 -> Subscript[\[CapitalLambda], 2]
};

(* Example: prF[m0' == Ph - nu0*m0 - m0*(la1 + la2)] displays with proper subscripts *)
