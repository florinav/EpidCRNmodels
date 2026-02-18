(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

(* SIR with vaccination (gas), waning (gar), treatment T[i], inhibition F via f[S,i].
   f[S,i] = be*S*i*F[i] is the total infection rate (opaque).
   T[i] is the treatment rate (opaque). Both handled by ODE2RNs. *)
var = {S, i, R}; RHS = 
{La + gar*R - f[S, i] - gas*S - mus*S + 
is*i,
  f[S, i] - (is+ir)*i - mui*i - T[i],
  gas*S +  T[i] - gar*R - mur*R+ ir*i};

prF[ex_] := ex /.
    {Derivative[1][T][i] -> Superscript[T, "\[Prime]"],
     Derivative[1, 0][f][S, i] -> Subscript[f, S],
     Derivative[0, 1][f][S, i] -> Subscript[f, i],
     La -> \[CapitalLambda], eta -> \[Eta], be -> \[Beta],
     mus -> Subscript[\[Mu], s], mui -> Subscript[\[Mu], i],
     mur -> Subscript[\[Mu], r], 
     ir -> Subscript[i, r],is -> Subscript[i, s],
     gai -> Subscript[\[Gamma], i], gar -> Subscript[\[Gamma], r],
     gas -> Subscript[\[Gamma], s]};
     con = {f[S_, i_] :> be*S*F[i]};
     conP = {f[S_, i_] :> be*S*i};
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNs[RHS, var, prF];
Print["chk=", Simplify[gam . rts - RHS]]
(*Print["alp=", alp // MatrixForm, " gam=", gam // MatrixForm]*)

(* Cauchy-Binet child selection analysis*)
{n, m} = {Length[var], Length[RN]};
RM = symRM[alp, n, m];
G = gam . RM;G//MatrixForm;
cpG = CharacteristicPolynomial[G, lam];
co = Rest@Reverse@CoefficientList[Expand[cpG], lam];
nCS = negCS[co, gam, n];upf = Select[nCS, upfQ];
{gH, edges} = hasseCS[nCS];

(*  ALGORITHM:
   1. symRM builds symbolic reactivity matrix RM (m x n, entries r[j,k] or 0).
      G = gam.RM is the symbolic Jacobian. Its char poly cpG has coefficients
      co[[k]] (k=1..n), for lam^(n-k), each a sum of k-minors of G (Cauchy-Binet).
   2. negCS extracts negative summands from each co[[k]] via allT+isN,
      then csInfo maps each to {k, term, sIdx, rIdx, Sk, det} where
      Sk = gam[[sIdx,rIdx]] is the INTEGER stoichiometric submatrix.
      Result: nCS list.
   3. upfQ  keeps entries with Sign[det]==(-1)^(k-1), det!=0
      (unstable positive feedback).
   4. hasseCS builds directed graph on nCS by r-set inclusion:
      edge a->b when rSet(a) is strict subset of rSet(b) and k(a)<k(b).
      TransitiveReductionGraph gives Hasse diagram.
   5. oscCS: for each UPF node u, finds descendants via VertexOutComponent,
      tests hurQ (Hurwitz stability of the integer Sk). If any descendant
      is stable, {u, stableDescendants} is an oscillatory core.
   6. Hopf via CS inheritance (Cor 5.21):
      csN = n-CS (full-size child selections from coeff of lam^0 = det G).
      stabN = Hurwitz-stable n-CS. If stabN nonempty AND upf nonempty -> Hopf.
   NOTE: gam here is integer (stoichiometric), so the numeric path works. *)
   Print["dim=", {n, m}, " #nCS=", Length[nCS], " #upf=", Length[upf], " #edges=", 
   Length[edges]]
upfRes = upfCS[nCS, gH, spe, rts];

(* Hopf via Cor 5.21 *)
csN = Select[nCS, #[[1]] == n &];
stabN = Select[csN, hurQ[#[[5]]] &];
fs = fixSgn[nCS, n];
Print["upfs exist=", Length[upf] > 0, " stable ", n, "x", n, " CS exist=",
      Length[stabN] > 0, ", Cor. 5.21 Hopf=", Length[stabN] > 0 && Length[upf] > 0]

(* Recipe I *)
riCores = oscRI[nCS, gH, spe, rts];
(* Recipe II *)
riiCores = oscRII[nCS, gH, spe, rts];

(* Symbolic Jacobian (with nonlinear f,T) and parametric Jacobian (con applied) *)
sJ = D[RHS, {var}];
sJc = Simplify[D[RHS /. con, {var}]];
Print["sJ=", sJ // MatrixForm//prF, " sJc=", sJc // MatrixForm//prF]


{co, res} = resHopf[sJ];FullSimplify/@co//prF


fu=FullSimplify/@res;fu//prF
cZF={gar->0,is->0};
fur=fu/.cZF//FullSimplify;fur//prF


fus=fur/.{gas->0,mur->mus}//FullSimplify;fus//prF
fus//Length
fus[[2]]//prF


conZ = {f[S_, i_] :> be*S*i, T[i_] :> eta*i/(1 + i b)};
RHS=RHS/.conZ;so=Solve[(RHS/.i->0)==0,{S,R}][[1]]//prF
sJ = Simplify[D[RHS /. con, {var}]];per={2,1,3};
Print["sJ=", sJ // MatrixForm//prF," sD=",sJ[[per,per]]/.i->0 // MatrixForm//prF]
{co, res} = resHopf[sJ];FullSimplify/@co//prF


fu=FullSimplify/@res//prF
fu//Length

