(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;

var = {S, i, R, W}; RHS = {La + gaw*W - be*S*i - mus*S,
  i*(be*S - gai - mui), gai*i + nu*i*W - (mur + gar)*R,
  gar*R - nu*i*W - gaw*W - muw*W};
  prF[ex_] := ex /.
  {x1 -> Subscript[x, 1], x2 -> Subscript[x, 2], x3 -> Subscript[x, 3],
  x0 -> Subscript[x, 0],
   S1 -> Subscript[S, 1], S2 -> Subscript[S, 2],
   B1 -> Subscript[B, 1], B2 -> Subscript[B, 2],
   et -> \[Eta], de -> \[Delta],x1H -> Subscript[OverHat[x], 1],
   x1T -> Subscript[OverTilde[x], 1],
   be -> \[Beta], la -> \[Lambda], La -> \[CapitalLambda], mu -> \[Mu],nu->\[Nu],
   mu1 -> Subscript[\[Mu], 1], mus -> Subscript[\[Mu], s],
   mur -> Subscript[\[Mu], r], muw -> Subscript[\[Mu], w],
   mui -> Subscript[\[Mu], i],
   ep1 -> Subscript[\[Epsilon], 1], ep2 -> Subscript[\[Epsilon], 2],
   gai -> Subscript[\[Gamma], i], gar -> Subscript[\[Gamma], r],gaw -> Subscript[\[Gamma], w],
   a1 -> Subscript[\[Alpha], 1], a2 -> Subscript[\[Alpha], 2], al -> \[Alpha],
   be1 -> Subscript[\[Beta], 1], be2 -> Subscript[\[Beta], 2],be3 -> Subscript[\[Beta], 3],
   om -> \[Omega]};
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var,prF];
RHS - gam . rts // FullSimplify  (*Should be {0,0,0,0}*)

{RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, sipVar} = NGMRN[RN, rts, var];
jac = Grad[RHS, var];
Print[" mSi=", mSi, " cE0=", Short[cE0]]
(* Cauchy-Binet child selection analysis
   ALGORITHM:
   1. symRM builds symbolic reactivity matrix RM (m x n, entries r[j,k] or 0).
      G = gam.RM is the symbolic Jacobian. Its char poly cpG has coefficients
      co[[k]] (k=1..n), each a sum of k-minors of G (Cauchy-Binet).
   2. negCS extracts negative summands from each co[[k]] via allT+isN,
      then csInfo maps each to {k, term, sIdx, rIdx, Sk, det} where
      Sk = gam[[sIdx,rIdx]] is the INTEGER stoichiometric submatrix.
      Result: nCS list (here 22 = 1+6+10+5 entries).
   3. upfQ keeps entries with Sign[det]==(-1)^(k-1), det!=0
      (unstable positive feedback). Here 13 of 22.
   4. hasseCS builds directed graph on nCS by r-set inclusion:
      edge a->b when rSet(a) is strict subset of rSet(b) and k(a)<k(b).
      TransitiveReductionGraph gives Hasse diagram. Here 71 edges.
   5. oscCS: for each UPF node u, finds descendants via VertexOutComponent,
      tests hurQ (Hurwitz stability of the integer Sk). If any descendant
      is stable, {u, stableDescendants} is an oscillatory core.
   HOW THE 4 OSCILLATORY CORES ARISE (all share rxn 3 = recovery):
   - osc#1: 1-CS {i},{3}. rSet={r[3,2]}. det(gam[{2},{3}])=1=(-1)^0: UPF.
     2 stable descendants (larger CS containing r[3,2] whose Sk is Hurwitz).
   - osc#2: 2-CS {i,R},{3,9}. rSet superset of osc#1. det=-1=(-1)^1: UPF.
     1 stable descendant.
   - osc#3: 2-CS {i,W},{3,10}. rSet superset of osc#1. det=-1=(-1)^1: UPF.
     1 stable descendant.
   - osc#4: 3-CS {i,R,W},{3,9,10}. rSet superset of osc#2,#3.
     det=1=(-1)^2: UPF. 1 stable descendant = full 4x4 system.
   Nested chain: osc#1 c osc#2,#3 c osc#4.
   Biological: the infection-recovery-waning cycle drives oscillations.
*)
{n, m} = {Length[var], Length[RN]};
RM = symRM[alp, n, m];
G = gam . RM;
cpG = CharacteristicPolynomial[G, lam];
co = Rest@Reverse@CoefficientList[Expand[cpG], lam];

nCS = negCS[co, gam, n];
upf = Select[nCS, upfQ];
{gH, edges} = hasseCS[nCS];
Print["dim=", {n, m}, " #nr of CS=", Length[nCS], " #nr of upf=", Length[upf], 
" #Hasse edges=", Length[edges]]
(* upfCS: UPF feedbacks with stable supersets *)
upfRes = upfCS[nCS, gH, spe, rts];

(* Hopf via Cor 5.21:
   csN: select n x n CS matrices from nCS (from co[[n]] = det G).
   stabN: among csN, select Hurwitz-stable ones. *)
csN = Select[nCS, #[[1]] == n &];
stabN = Select[csN, hurQ[#[[5]]] &];

fs = fixSgn[nCS, n];
Print["fixed sign by order=", fs,
      " #", n, "x", n, " CS=", Length[csN], " dets=", #[[6]] & /@ csN]
Print["upfs exist=", Length[upf] > 0, " stable ", n, "x", n, " CS exist=",
      Length[stabN] > 0, ", Cor. 5.21 Hopf=", Length[stabN] > 0 && Length[upf] > 0]

(* Recipe I: Oscillatory Cores of Class I (minimal stable CS containing a UPF) *)
riCores = oscRI[nCS, gH, spe, rts];
(* Recipe II: Oscillatory Cores of Class II (minimal UNF with stable submatrix) *)
riiCores = oscRII[nCS, gH, spe, rts];

(* Symbolic Jacobian) *)
sJ = D[RHS, {var}];
Print["sJ=", sJ // MatrixForm//prF]
{co, res} = resHopf[sJ];FullSimplify/@co//prF
(*fu=FullSimplify/@res aborted*)


var = {W, S2, i};
 RN = {
  0 -> "W",
  "W" -> "S2",
  "i" + "W" -> "i",
  "S2"+ "i" -> 2*"i",
  "W" -> 0,
  "S2" -> 0,
  "i" -> 0
};
rts = {
  La,
  ga*W,
  nu*i*W,
  be*i* S2,
  mu1*W,
  mu2*S2,
  mu3*i
};
 prF[ex_] := ex /. { S2 -> Subscript[S, 2],
  mu2 -> Subscript[\[Mu], 2], mu1 -> Subscript[\[Mu], 1], mu3 -> Subscript[\[Mu], 3],
   nu -> \[Nu], de -> \[Delta],ga->\[Gamma],
   be -> \[Beta], la -> \[Lambda], La -> \[CapitalLambda], mu -> \[Mu]
  };
  Print[rts//Length," reactions and transitions: ",Transpose[{RN,rts}]//MatrixForm//prF]

{RHS,  par, cp, mSi, Jx, Jy,Jxy, fx, cDFE, cE0, F, V, K, R0A, sipars} = 
NGMRN[RN, rts, var];
Print["RHS=", FullSimplify/@RHS//prF//MatrixForm," mSi=", 
mSi];
so=Solve[RHS==0,var];
Print[so//Length," sol="];
so//prF
Print["R0A=",R0A//prF," R0=",R0A/.so[[1]]//FullSimplify//prF];
(* Jacobian: differentiate first, then evaluate at endemic eq *)
jac = Grad[RHS, var];
eqE = {W -> La/(ga + nu*i + mu1),
       S2 -> mu3/be};
jEnd = jac /. eqE;
{co, res} = resHopf[jEnd]//FullSimplify;
Print["co=",  co // prF, " res=", res // prF]


so=selZF[Solve[res==0]]//FullSimplify;so//prF
co=CoefficientList[-Numerator[Together[RHS[[2]]/.eqE]],i]
cp=Thread[Variables[co]>0];FullSimplify[Reduce[And@@cp&&co[[1]]<0],cp]//prF



