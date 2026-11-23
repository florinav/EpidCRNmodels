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

pPrint[mat_, mSi_] := Module[
  {n, spos, cuts, unionCut, intraCuts, prMat, colSpec, rowStr, div, texStr},
  n = Length[mat];
  spos = Map[(First @ Flatten @ Position[Global`var,
                If[StringQ[#], ToExpression[#], #]]) &, mSi, {2}];
  cuts = Sort[Max /@ spos] + 1;
  unionCut = Last[cuts];   (* siphon-union / rest boundary *)
  intraCuts = Most[cuts];  (* intra-siphon boundaries *)
  prMat = Map[Global`prF, mat, {2}];
  (* LaTeX column spec: | between blocks, || at siphon/rest boundary *)
  colSpec = "|" <> StringJoin @ Table[
    "c" <> Which[
      k == n,                                   "|",
      MemberQ[intraCuts, k + 1],               "|",
      k + 1 == unionCut && unionCut <= n,      "||",
      True,                                     ""],
    {k, n}];
  (* LaTeX rows with \hline / \hline\hline separators *)
  rowStr = StringJoin @ Table[
    StringJoin[Riffle[Table[ToString[TeXForm[prMat[[i,j]]]], {j,n}], " & "]] <>
    " \\\\\n" <>
    Which[
      MemberQ[intraCuts, i + 1],             "\\hline\n",
      i + 1 == unionCut && unionCut <= n,    "\\hline\\hline\n",
      True,                                   ""],
    {i, n}];
  texStr = "\\begin{array}{" <> colSpec <> "}\n\\hline\n" <>
           rowStr <> "\\hline\n\\end{array}";
  Print[texStr];
  (* Visual Grid *)
  div = Join[{1 -> True, (n+1) -> True}, Thread[intraCuts -> True],
             If[unionCut <= n, {unionCut -> {True, True}}, {}]];
  Grid[prMat, Dividers -> {div, div}, Alignment -> Center, Spacings -> {0.8, 0.4}]];

zQ[v_, s_] := (v /. s) === 0;

SetAttributes[stabPrint, HoldFirst];
stabPrint[c_] := Module[{j, ch, lS, qS, hD, lF, qF, hd, dh, ine},
  Print["--- ", ToString[HoldForm[c]], " ---"];
  j = Global`jac /. c;
  ch = Numerator[Together[CharacteristicPolynomial[j, Global`u]]] // Factor;
  {lS, qS, hD, lF, qF} = staPP[ch, Global`par];
  Print["lSta=", FullSimplify /@ lS // Global`prF,
        "  qSta=", FullSimplify /@ qS // Global`prF, "  hDeg=", hD // Global`prF];
  If[hD =!= {},
    hd = hD[[1]]; dh = Exponent[hd, Global`u];
    Which[dh == 3, ine = Hur3[hd, Global`u]; Print["Hur3: ", FullSimplify /@ ine // Global`prF],
          dh == 4, ine = Hur4[hd, Global`u]; Print["Hur4: ", FullSimplify /@ ine // Global`prF],
          True, Print["Hur: deg=", dh, " not handled"]]]];

LyapPV[p_, svar_, zBlocks_, blkData_, cEp_, RHS_, var_] := Module[
  {m, VinvList, mRList, cList, sp, zp, np, zpstar, Gv, Hp, Lp,
   spos, zpPos, dLpdt, posOf},
  m = Length[blkData];
  Gv[x_] := x - 1 - Log[x];
  posOf = (First @ Flatten @ Position[var, #]) &;
  VinvList = Inverse /@ blkData[[All, 3]];
  (* mR_k = wInc_k^T V_k^{-1} wEnt_k; c_k = wInc_k^T V_k^{-1} / mR_k *)
  mRList = Table[blkData[[k, 1]] . VinvList[[k]] . blkData[[k, 2]], {k, m}];
  cList  = Table[blkData[[k, 1]] . VinvList[[k]] / mRList[[k]], {k, m}];
  sp     = svar /. cEp;
  zp     = zBlocks[[p]];
  np     = Length[zp];
  zpstar = zp /. cEp;
  spos   = posOf[svar];
  zpPos  = posOf /@ zp;
  (* H_p = s_p*G(s/s_p) + sum_j c_{p,j} * zpstar_j * G(zp_j / zpstar_j) *)
  Hp = sp*Gv[svar/sp] + Sum[cList[[p, j]]*zpstar[[j]]*Gv[zp[[j]]/zpstar[[j]]], {j, np}];
  (* L_p = H_p + sum_{k != p} c_k^T z_k *)
  Lp = Hp + Sum[If[k == p, 0, cList[[k]] . zBlocks[[k]]], {k, m}];
  (* dL_p/dt = (1-sp/s)*s' + sum_j c_{p,j}*(1-zpstar_j/zp_j)*zp_j' + sum_{k!=p} c_k . zk' *)
  dLpdt = (1 - sp/svar)*RHS[[spos]] +
          Sum[cList[[p, j]]*(1 - zpstar[[j]]/zp[[j]])*RHS[[zpPos[[j]]]], {j, np}] +
          Sum[If[k == p, 0, cList[[k]] . RHS[[(posOf /@ zBlocks[[k]])]]], {k, m}];
  {Lp, dLpdt, mRList, cList, zpstar}];

(* Tests for onlyP (from Conv.wl)

onlyP[a + b]
onlyP[1/(a*b)]
onlyP[-a]
onlyP[a - b]
onlyP[0]
*)

(* dulP: planar Dulac / second-compound test for a monomial metric.  In two
   dimensions J^[2]=tr J, so for p=x^a y^b the transformed second compound is
   pdot/p + tr J = div(p f)/p.  Strict negativity on the open quadrant is a
   Dulac certificate: no periodic orbit and no homoclinic loop.  Returns
   {TJ, trm} with trm the expanded summand list, so that a term-by-term sign
   test is available as a sufficient condition, and Reduce on TJ<0 gives the
   exact one. *)
dulP[RHS_, var_, exps_] := Module[{tr, pd, TJ},
  tr = Tr[D[RHS, {var}]];
  pd = exps[[1]] RHS[[1]]/var[[1]] + exps[[2]] RHS[[2]]/var[[2]];
  TJ = Expand[pd + tr];
  {TJ, If[Head[TJ] === Plus, List @@ TJ, {TJ}]}];

(* ACM: k-th additive compound matrix, ported from EpidCRNo.wl where it was the
   only definition, although EpidCRN.wl already declared its usage.  The
   auxiliary variable is localized; the original used a global t. *)
ACM[matrix_, k_] := Module[{t},
  D[Minors[IdentityMatrix[Length@matrix] + t*matrix, k], t] /. t -> 0];

(* Pdot: state-dependent derivative of a matrix metric along the flow,
   P_f = sum_i dP/dxi * fi.  D acts elementwise on a matrix P, so no
   Table/loop is needed. *)
Pdot[P_, f_, x_] := Sum[D[P, x[[i]]]*f[[i]], {i, Length[x]}];

(* TACM: the k-th transformed additive compound, B_k = P_f P^-1 + P J^[k] P^-1.
   Generalizes dulP's planar TJ (a 1x1, k=2 special case) to any order k and
   any invertible state-dependent metric P.  Reuses ACM for J^[k] and Pdot
   for P_f; does not touch either. *)
TACM[J_, k_, P_, f_, x_] := Module[{Pf, Pi, Jk},
  Pf = Pdot[P, f, x]; Pi = Inverse[P]; Jk = ACM[J, k];
  Pf.Pi + P.Jk.Pi];

(* matMu: the matrix measure (logarithmic norm) of B induced by the l1, l2 or
   linfinity vector norm, p in {1,2,Infinity}.  For a Metzler B, mu_1 reduces
   to the max column sum (no Abs needed, off-diagonal entries already >=0). *)
matMu[B_, p_] := With[{n = Length[B]},
  Which[
    p === 1, Max[Table[B[[j, j]] + Total[Abs[Delete[B[[All, j]], j]]], {j, n}]],
    p === Infinity, Max[Table[B[[i, i]] + Total[Abs[Delete[B[[i, All]], i]]], {i, n}]],
    p === 2, Max[Eigenvalues[(B + Transpose[B])/2]],
    True, $Failed]];

(* qLMI: fixed-rate feasibility test for a common quadratic contraction metric
   on a matrix polytope Jlist={J1,...,JN} (a single matrix is Jlist={J}).
   Builds a symmetric P from fresh Unique[] variables (n(n+1)/2 of them),
   normalizes P>=I to remove the trivial P=0 solution, and puts a fixed -I
   margin on every vertex LMI Jl^T.P+P.Jl<=-2 eta P - I.  Fixing eta rather
   than treating it as an unknown avoids the bilinear eta*P term, exactly as
   compu.tex recommends; qRate below bisects on eta using this as the oracle.
   Returns {True,P} (P a numeric symmetric matrix) if feasible, {False,Null}
   otherwise. *)
qLMI[Jlist_, eta_] := Module[{n, uv, P, vars, cons, sol},
  n = Length[First[Jlist]];
  uv = Table[Unique["qcp"], {n}, {n}];
  P = Table[If[j >= i, uv[[i, j]], uv[[j, i]]], {i, n}, {j, n}];
  vars = Union[Flatten[Table[uv[[i, j]], {i, n}, {j, i, n}]]];
  cons = Join[
    {VectorGreaterEqual[{P - IdentityMatrix[n], 0}, "SemidefiniteCone"]},
    Table[VectorGreaterEqual[
      {-(Transpose[Jl].P + P.Jl + 2 eta P) - IdentityMatrix[n], 0},
      "SemidefiniteCone"], {Jl, Jlist}]];
  sol = Quiet[SemidefiniteOptimization[0, cons, vars]];
  If[FreeQ[vars /. sol, Indeterminate], {True, P /. sol}, {False, Null}]];

(* qRate: largest common quadratic-contraction rate for the matrix polytope
   Jlist, by bisection on eta with qLMI as the feasibility oracle. Returns
   {eta,P} for the best certified rate found, to bisection tolerance tol
   (default 1e-4), capped by etaMax (default 10). *)
qRate[Jlist_, etaMax_: 10, tol_: 10^-4] := Module[{lo, hi, mid, feas, P, bestEta, bestP},
  {feas, P} = qLMI[Jlist, 0];
  If[! feas, Return[{0, Null}]];
  {bestEta, bestP} = {0, P};
  {feas, P} = qLMI[Jlist, etaMax];
  If[feas, Return[{etaMax, P}]];
  lo = 0; hi = etaMax;
  While[hi - lo > tol,
    mid = (lo + hi)/2;
    {feas, P} = qLMI[Jlist, mid];
    If[feas, lo = mid; bestEta = mid; bestP = P, hi = mid]];
  {bestEta, bestP}];

(* jPoly: matrix polytope vertices for J(x)=A0+Sum_j hj(x) Aj with each hj
   ranging over a known interval, given as hb={{lb1,ub1},...,{lbm,ubm}}. Only
   the mechanical vertex construction is automated (Tuples[hb] enumerates all
   2^m sign choices at the bounds); discovering the affine decomposition and
   its bounds from J,x symbolically -- what compu.tex's JacobianPolytope also
   proposes -- is left to the caller, as noted in the README plan. Feed the
   result straight to qRate for a polytope-robust quadratic certificate. *)
jPoly[A0_, As_, hb_] := (A0 + Total[MapThread[Times, {#, As}]]) & /@ Tuples[hb];

(* sosCon: builds (but does not solve) the SOS-matrix feasibility problem for
   a symmetric polynomial matrix M(x): v.M.v is SOS in (x,v) iff it equals
   z.Q.z for some PSD Gram matrix Q on the monomial basis z=v_i*(monomials of
   x up to degree d). Coefficient matching gives linear equalities on Q; kept
   separate from solving so several SOS-matrix conditions sharing the same
   outer unknowns (e.g. polyMet below) can be combined into one SDP call
   instead of solved independently. Returns {cons,vars,Q}. *)
sosCon[M_, x_, d_] := Module[{m, mono, v, z, nz, uv, Q, vars, expr, diff, ceqs, cons},
  m = Length[M];
  mono = (Times @@ (x^#)) & /@ Select[Tuples[Range[0, d], Length[x]], Total[#] <= d &];
  v = Table[Unique["sosv"], {m}];
  z = Flatten[Outer[Times, v, mono]];
  nz = Length[z];
  uv = Table[Unique["sosq"], {nz}, {nz}];
  Q = Table[If[j >= i, uv[[i, j]], uv[[j, i]]], {i, nz}, {j, nz}];
  vars = Union[Flatten[Table[uv[[i, j]], {i, nz}, {j, i, nz}]]];
  expr = Expand[v.M.v];
  diff = Expand[expr - z.Q.z];
  ceqs = (# == 0) & /@ (Last /@ CoefficientRules[diff, Join[x, v]]);
  cons = Join[{VectorGreaterEqual[{Q, 0}, "SemidefiniteCone"]}, ceqs];
  {cons, vars, Q}];

(* sosMat: standalone SOS-matrix certificate for a single symmetric polynomial
   matrix M(x) (a bare polynomial p(x) counts as the 1x1 case {{p}}), degree d
   monomial basis. Thin wrapper solving sosCon's single condition via
   SemidefiniteOptimization. Returns {True,Q} (numeric Gram matrix) if
   feasible, {False,Null} otherwise. Following Parrilo's SOS relaxation of
   polynomial matrix inequalities, cited in compu.tex. *)
sosMat[M_, x_, d_] := Module[{cons, vars, Q, sol, ok},
  {cons, vars, Q} = sosCon[M, x, d];
  sol = Quiet[SemidefiniteOptimization[0, cons, vars]];
  ok = FreeQ[vars /. sol, Indeterminate];
  If[ok, {True, Q /. sol}, {False, Null}]];

(* polyP: symmetric n x n matrix with entries that are independent degree-<=d
   polynomials in x, each with fresh Unique[] coefficients. Returns {P,vars};
   vars is exactly the free coefficient list, for feeding to a solver. *)
polyP[x_, n_, d_] := Module[{mono, nm, coefs, P, vars},
  mono = (Times @@ (x^#)) & /@ Select[Tuples[Range[0, d], Length[x]], Total[#] <= d &];
  nm = Length[mono];
  coefs = Table[Table[Unique["pa"], {nm}], {n}, {n}];
  P = Table[If[j >= i, coefs[[i, j]].mono, coefs[[j, i]].mono], {i, n}, {j, n}];
  vars = Flatten[Table[coefs[[i, j]], {i, n}, {j, i, n}], 2];
  {P, vars}];

(* polyMet: searches for a degree-<=d polynomial contraction metric P(x)>0
   certifying -(P_f+J^T.P+P.J+2*eta*P)>0 via SOS, both conditions with a -I
   margin exactly as qLMI uses for the constant case.  P is built by polyP,
   P_f by Pdot (reused unchanged), and the two SOS-matrix conditions
   (positivity of P, and the contraction inequality) are combined into one
   joint SDP via sosCon so they share the same metric coefficients -- this is
   what compu.tex's "PolynomialContractionMetric" asks for.  Fixing eta, as
   in qLMI, avoids the bilinear eta*P term.  At d=0, P is constant, P_f
   vanishes identically, and this reduces exactly to qLMI[{J},eta] (verified
   in the test below).  Returns {True,P} or {False,Null}. *)
polyMet[J_, f_, x_, d_, eta_] := Module[{n, P, avars, Pf, M1, M2, c1, v1, c2, v2,
   cons, vars, sol, ok},
  n = Length[J];
  {P, avars} = polyP[x, n, d];
  Pf = Pdot[P, f, x];
  M1 = P - IdentityMatrix[n];
  M2 = -(Pf + Transpose[J].P + P.J + 2 eta P) - IdentityMatrix[n];
  {c1, v1} = Take[sosCon[M1, x, d], 2];
  {c2, v2} = Take[sosCon[M2, x, d], 2];
  cons = Join[c1, c2];
  vars = Union[v1, v2, avars];
  sol = Quiet[SemidefiniteOptimization[0, cons, vars]];
  ok = FreeQ[vars /. sol, Indeterminate];
  If[ok, {True, P /. sol}, {False, Null}]];

(* optMon: exponent search for a monomial contraction metric p=x^a y^b via
   dulP, over box={{aLo,aHi},{bLo,bHi}}.  For fixed (x,y), dulP's TJ is
   AFFINE in (a,b) -- no a*b cross term, by inspection of the monomial
   compound formula -- so this is not the fragile nonconvex bi-level search
   compu.tex's OptimizeKMetric suggests in general; it is a semi-infinite
   LINEAR program, solved by a cutting-plane loop: maximize the margin delta
   over the exponent box subject to TJ<=-delta at a growing set of sample
   points (LinearOptimization), then verify the winning (a,b) EXACTLY via
   Reduce[ForAll[...]]; if verification fails, FindInstance extracts a real
   counterexample point, which is added to the sample set for the next LP.
   Terminates as soon as one exponent pair verifies (returns {True,{a,b}}),
   the LP itself becomes infeasible (returns {False,Null}, no certificate in
   this box), or maxIter is reached. Returns {found,{a,b}or Null,iters}. *)
optMon[RHS_, var_, box_, pts0_, maxIter_: 12] := Module[
  {a, b, d, TJsym, TJf, pts, cons, vars, sol, aStar, bStar, dStar,
   ok, iter, negReg, cex},
  a = Unique["oma"]; b = Unique["omb"]; d = Unique["omd"];
  TJsym = First[dulP[RHS, var, {a, b}]];
  TJf[x0_, y0_] := TJsym /. Thread[var -> {x0, y0}];
  pts = pts0; ok = False; iter = 0; aStar = bStar = Null;
  While[iter < maxIter && ! ok,
    cons = Join[
      (TJf[#[[1]], #[[2]]] <= -d) & /@ pts,
      {a >= box[[1, 1]], a <= box[[1, 2]], b >= box[[2, 1]], b <= box[[2, 2]]}];
    vars = {a, b, d};
    sol = Quiet[LinearOptimization[-d, cons, vars]];
    If[! FreeQ[vars /. sol, Indeterminate], Return[{False, Null, iter}]];
    {aStar, bStar, dStar} = vars /. sol;
    If[dStar <= 0, Return[{False, Null, iter}]];
    negReg = Reduce[var[[1]] > 0 && var[[2]] > 0 &&
       (TJsym /. {a -> aStar, b -> bStar}) >= 0, var, Reals];
    If[negReg === False,
      ok = True,
      cex = {var[[1]], var[[2]]} /. FindInstance[negReg, var, Reals][[1]];
      pts = Append[pts, cex]; iter++]];
  If[ok, {True, {aStar, bStar}, iter}, {False, Null, iter}]];

(* kTau: the k-shifted log norm tau_{p,k} of Dalin, Ofir, Bar-Shalom,
   Ovseevich, Bullo & Margaliot, "Verifying k-Contraction without Computing
   k-Compounds" (IEEE TAC 2024, arXiv:2209.01046), Definition 4 and
   Theorem 17 there:
     tau_{p,k}(A) = tr(A) + (n-k)*mu_q(T.(-A).T^-1),  1/p+1/q=1,
   which satisfies mu_{p,T^(k)}(A^[k]) <= tau_{p,k}(A) WITHOUT ever forming
   the k-compound A^[k] (T^(k) is the k-th MULTIPLICATIVE compound of the
   metric T, needed only to interpret the resulting norm, never computed
   here). So tau_{p,k}(A)<=-eta<0 certifies k-contraction directly from A.
   Reuses matMu (Phase 1) unchanged; T defaults to the identity. tau_{p,k}
   is convex in A (their Prop. 19(7)), so for a Jacobian ranging over a
   polytope or box, Max[kTau[#,k,p,T]&/@vertices] suffices, exactly as
   qLMI/qRate/jPoly already exploit convexity via vertex enumeration --
   this is compu.tex's CompoundFreeKContraction, Method->"DalinEtAl" (the
   compound-duality criterion); Method->"GeneralizedLyapunov"
   (Cecilia et al. 2026) is not yet implemented, so no Method-> dispatcher
   is built here -- a one-branch dispatcher would be premature. *)
kTau[A_, k_, p_, T_: Automatic] := Module[{n, q, Tm},
  n = Length[A];
  q = p /. {1 -> Infinity, Infinity -> 1, 2 -> 2};
  Tm = If[T === Automatic, IdentityMatrix[n], T];
  Tr[A] + (n - k)*matMu[Tm.(-A).Inverse[Tm], q]];

(* exKC: numeric-to-exact reconstruction of an SDP-solver candidate (a matrix,
   or a matrix with symbolic entries whose numeric COEFFICIENTS need
   rounding -- Rationalize touches only the numeric parts, symbols are left
   alone), per compu.tex's "P_hat -> P_rational" step: -2.999999997 should
   become -3, not stay a 9-digit float.  Elementwise Rationalize to within
   tol; the caller retries at several tol via exVerQ/exVerPoly below, since
   the RIGHT tolerance is whichever one still verifies exactly afterward. *)
exKC[val_, tol_: 10^-6] := Rationalize[val, tol];

(* posDefQ: exact positive-definiteness test for a numeric/rational constant
   matrix via Sylvester's criterion (all leading principal minors positive).
   No Reduce, no floating point once the matrix is exact -- genuinely exact
   linear algebra, as compu.tex specifies for the constant/rational case. *)
posDefQ[M_] := And @@ Table[Det[M[[1 ;; i, 1 ;; i]]] > 0, {i, Length[M]}];

(* verQ: exact verification of a constant quadratic k-contraction certificate
   Pex on the matrix polytope Jlist at rate eta -- Pex>0 and, for every
   vertex Jl, -(Jl^T.Pex+Pex.Jl+2*eta*Pex)>0, both via posDefQ. Companion to
   qLMI/qRate, which return the same shape of certificate but numeric. *)
verQ[Pex_, Jlist_, eta_] := posDefQ[Pex] &&
  And @@ (posDefQ[-(Transpose[#].Pex + Pex.# + 2*eta*Pex)] & /@ Jlist);

(* exVerQ: tries exKC at a sequence of tolerances (coarsest first, so the
   simplest exact candidate is tried before finer ones), keeps whichever
   ones verQ confirms exactly, and returns the simplest of those (smallest
   maximal denominator) -- compu.tex's ExactifyKCertificate+VerifyKCertificate
   combined for the constant case. Returns {"Exact",Pex} if some tolerance
   verified, {"Sampled",P} (the original numeric candidate, unverified) if
   none did.  "RigorousNumerical" (interval-arithmetic fallback) is not
   implemented -- only "Exact" and "Sampled" are genuine here. *)
exVerQ[P_, Jlist_, eta_, tols_: {10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-8}] :=
 Module[{cands, good},
  cands = exKC[P, #] & /@ tols;
  good = Select[cands, verQ[#, Jlist, eta] &];
  If[good =!= {},
    {"Exact", First[SortBy[good, Max[Denominator /@ Flatten[#]] &]]},
    {"Sampled", P}]];

(* verPoly: exact positive-definiteness of a symmetric polynomial matrix M(x)
   for all x satisfying cond (True = unrestricted), via Sylvester's leading-
   principal-minor test quantified over x with Reduce[ForAll[...]] -- the
   same exact-verification idiom optMon already uses, generalized from a
   scalar TJ to a matrix. x a one-element list is unwrapped to a bare symbol,
   since ForAll[{x},cond,expr] is not valid syntax (only ForAll[x,cond,expr]
   or ForAll[{x1,x2,...},cond,expr] are), a real syntax pitfall hit and fixed
   while building this. *)
verPoly[M_, x_, cond_: True] := Module[{minorsPos, xq},
  minorsPos = And @@ Table[Det[M[[1 ;; i, 1 ;; i]]] > 0, {i, Length[M]}];
  xq = If[Length[x] == 1, First[x], x];
  Reduce[ForAll[xq, cond, minorsPos], x, Reals] === True];

(* exVerPoly: the polyMet analogue of exVerQ.  Rebuilds M1=Pex-I and
   M2=-(Pdot[Pex,f,x]+J^T.Pex+Pex.J+2*eta*Pex)-I from each exactified
   candidate and checks both via verPoly; returns {"Exact",Pex} for the
   first tolerance (coarsest first) that verifies both globally, {"Sampled",
   Pnum} otherwise. Reuses Pdot (Phase 1) unchanged. *)
exVerPoly[J_, f_, x_, eta_, Pnum_, cond_: True,
   tols_: {10^-2, 10^-3, 10^-4, 10^-5, 10^-6}] := Module[
  {n, cands, good, test},
  n = Length[J];
  test[Pex_] := Module[{M1, M2},
    M1 = Pex - IdentityMatrix[n];
    M2 = -(Pdot[Pex, f, x] + Transpose[J].Pex + Pex.J + 2 eta Pex) - IdentityMatrix[n];
    verPoly[M1, x, cond] && verPoly[M2, x, cond]];
  cands = exKC[Pnum, #] & /@ tols;
  good = Select[cands, test];
  If[good =!= {}, {"Exact", First[good]}, {"Sampled", Pnum}]];

(* lyaTry: tries each mu in muCands (in order) as the shift in the Sylvester
   equation (A-mu I)^T.P+P.(A-mu I)=-I, returning the first {mu,P} for which
   LyapunovSolve gives a genuine unique numeric solution. A candidate mu can
   fail (LyapunovSolve::nuniq) when two eigenvalues of A-mu I sum to exactly
   zero -- e.g. the arithmetic midpoint of two real eigenvalues always causes
   this, a real degeneracy hit while building cecKC below -- so callers pass
   several candidates spanning the valid range, not just the midpoint.
   Internal helper for cecKC; not given a usage string. *)
lyaTry[A_, muCands_, n_] := Module[{P, ok},
  Do[
    Off[LyapunovSolve::nuniq];
    P = Quiet[LyapunovSolve[Transpose[N[A] - mu*IdentityMatrix[n]], -IdentityMatrix[n]]];
    On[LyapunovSolve::nuniq];
    ok = MatrixQ[P, NumericQ];
    If[ok, Return[{mu, P}, Module]],
    {mu, muCands}];
  {Null, Null}];

(* cecKC: Andreu/Cecilia, Zoboli, Astolfi, Serres & Andrieu, "Generalized
   Lyapunov Conditions for k-Contraction: Analysis and Feedback Design" (IEEE
   TAC 2026), Theorem k_contraction_LMI_cLTI specialized to l=2 -- the
   Method->"GeneralizedLyapunov" companion to kTau flagged as not yet
   implemented above. For xdot=A.x, finds constant symmetric P0 (In=(0,0,n),
   positive definite), P1 (In=(k-1,0,n-k+1)), and scalars mu0,mu1 with
   A^T.Pi+Pi.A < 2*mui*Pi (i=0,1) and mu1+(k-1)*mu0<0 -- NECESSARY AND
   SUFFICIENT for k-contraction of a constant A (their cLTI theorem), so
   strictly sharper than kTau (compound-free but sufficient-only) whenever
   the identity metric kTau defaults to is a bad fit for a non-normal A.
   Construction is DIRECT, not an SDP search, by their own Remark[algorithm]
   + Lemma[inertia_shifting] item (4) (the classical Ostrowski-Schneider
   Lyapunov inertia theorem): for Q=I and any mu not equal to any
   Re(eigenvalue of A), the Sylvester equation (A-mu I)^T.P+P.(A-mu I)=-I has
   a unique symmetric solution P, and expanding it gives A^T.P+P.A=2 mu P-I,
   i.e. exactly the required strict LMI with margin I; P's inertia is FORCED
   to (p,0,n-p) where p=#{eigenvalues of A with Re>mu} -- so mu0 chosen above
   the largest Re(eigenvalue) (p=0) and mu1 chosen strictly between the
   (k-1)th and kth largest Re(eigenvalue) (p=k-1) automatically yield P0,P1
   of the REQUIRED inertia, with no inertia-constrained SDP needed at all.
   Several mu0,mu1 candidates are tried via lyaTry to dodge the pairwise-
   eigenvalue-sum degeneracy noted there. Returns {ok,P0,P1,mu0,mu1}, ok=True
   iff a valid spectral gap at position (k-1,k) exists, both Lyapunov solves
   succeed, and mu1+(k-1)*mu0<0; requires 2<=k<=n (k=1 collapses In(P1) to
   In(P0), i.e. ordinary contraction -- use qLMI/qRate instead). The
   nonlinear/state-dependent case (their Theorem k_contraction_LMI proper,
   "for all x in S") needs a joint search across a matrix polytope and is not
   attempted here -- deferred, exactly as jPoly/qRate handle kTau's own
   analogous polytope extension, see README. *)
cecKC[A_, k_Integer] := Module[
  {n, eigs, reE, mu0Cands, mu1Cands, mu0, mu1, P0, P1, ok, gap},
  n = Length[A];
  If[k < 2 || k > n, Return[{False, Null, Null, Null, Null}]];
  eigs = Eigenvalues[N[A]];
  reE = Sort[Re[eigs], Greater];
  If[reE[[k - 1]] <= reE[[k]], Return[{False, Null, Null, Null, Null}]];
  gap = Max[1, Abs[reE[[1]]]];
  mu0Cands = reE[[1]] + # & /@ ({0.01, 0.1, 0.05, 0.2, 1, 2}*gap);
  mu1Cands = reE[[k]] + # (reE[[k - 1]] - reE[[k]]) & /@
    {0.3, 0.7, 0.2, 0.8, 0.4, 0.6, 0.1, 0.9, 0.25, 0.75};
  {mu0, P0} = lyaTry[A, mu0Cands, n];
  {mu1, P1} = lyaTry[A, mu1Cands, n];
  If[mu0 === Null || mu1 === Null, Return[{False, Null, Null, Null, Null}]];
  ok = mu1 + (k - 1)*mu0 < 0;
  {ok, P0, P1, mu0, mu1}];

(* inCnt: exact inertia triple {neg,zero,pos} of a symmetric matrix's
   eigenvalue signs, matching Cecilia et al.'s own In(P)=(pi_-,pi_0,pi_+)
   convention (negative count first, confirmed directly from their def.tex
   source). Genuinely exact when M has exact/rational entries -- Eigenvalues
   then returns algebraic numbers (Root objects for non-quadratic cases) with
   exactly-decidable signs, no numerics involved. *)
inCnt[M_] := Module[{ev = Eigenvalues[M]},
  {Count[ev, _?(# < 0 &)], Count[ev, _?(# == 0 &)], Count[ev, _?(# > 0 &)]}];

(* exCecKC: numeric-to-exact reconstruction and independent verification of a
   cecKC certificate, the same exKC-then-verify pattern as exVerQ. Exactifies
   P0,P1,mu0,mu1 elementwise at a sequence of tolerances (coarsest first) and
   keeps the first that verifies ALL FOUR of the theorem's own conditions
   independently: inCnt[P0]===(0,0,n), inCnt[P1]===(k-1,0,n-k+1) (bypassing
   cecKC's own inertia reasoning entirely, since inCnt just counts exact
   eigenvalue signs), posDefQ on both LMI residuals (reusing the exactify
   layer's own Sylvester-criterion code), and the scalar sum condition.
   Returns {"Exact",...} or {"Sampled",P0,P1,mu0,mu1} (the original numeric
   candidates, unverified) if no tolerance verifies. *)
exCecKC[A_, k_, P0_, P1_, mu0_, mu1_,
   tols_: {10^-2, 10^-3, 10^-4, 10^-5, 10^-6}] := Module[
  {n, test, cands0, cands1, candm0, candm1},
  n = Length[A];
  test[P0x_, P1x_, m0x_, m1x_] :=
    inCnt[P0x] === {0, 0, n} &&
    inCnt[P1x] === {k - 1, 0, n - k + 1} &&
    posDefQ[-(Transpose[A].P0x + P0x.A - 2*m0x*P0x)] &&
    posDefQ[-(Transpose[A].P1x + P1x.A - 2*m1x*P1x)] &&
    m1x + (k - 1)*m0x < 0;
  Do[
    {cands0, cands1, candm0, candm1} =
      {exKC[P0, tol], exKC[P1, tol], exKC[mu0, tol], exKC[mu1, tol]};
    If[test[cands0, cands1, candm0, candm1],
      Return[{"Exact", cands0, cands1, candm0, candm1}, Module]],
    {tol, tols}];
  {"Sampled", P0, P1, mu0, mu1}];

(* jntCec: the state-dependent/polytope extension of cecKC, deferred since
   Phase 3 began -- a single constant (P0,P1,mu0,mu1) valid at EVERY vertex
   of a matrix polytope Jlist (e.g. from jPoly), not just one frozen point,
   which is what Cecilia et al.'s own Theorem k_contraction_LMI actually
   asks for on a forward-invariant set (cecKC alone only ever certified a
   single matrix). Extensive search for a genuinely INDEFINITE P1 shared
   across several vertices via feasibility-only SDP (dropping the inertia
   constraint per Remark[algorithm], exactly as cecKC's own single-matrix
   construction does) consistently returns a POSITIVE DEFINITE point
   whenever anything is feasible at all, never an indefinite one -- the
   solver has no incentive to prefer one over the other, so this dead end
   was abandoned rather than forced. What DOES work: cecKC's own exact
   single-matrix certificate at ONE vertex, checked directly against every
   OTHER vertex's LMI residual -- valid by continuity whenever the
   polytope is not too large. This is not a new search at all, just reuse
   of cecKC exactly as built, tried at each vertex in turn until one
   certificate survives everywhere. Returns {True,P0,P1,mu0,mu1,i} (i =
   the vertex whose certificate worked) or {False,Null,Null,Null,Null,
   Null}. No claim of completeness: a polytope with no single vertex whose
   certificate extends everywhere may still admit some OTHER joint
   certificate this simple strategy does not find -- see the test below
   for a genuine case where it does not (and should not: two of that
   polytope's own vertices are, individually, not even 2-contractive). *)
jntCec[Jlist_, k_] := Module[{test, ok, P0, P1, mu0, mu1},
  test[P0x_, P1x_, m0x_, m1x_] := And @@ (Module[{m0, m1},
      m0 = Min[Re[Eigenvalues[-(Transpose[#].P0x + P0x.# - 2*m0x*P0x)]]];
      m1 = Min[Re[Eigenvalues[-(Transpose[#].P1x + P1x.# - 2*m1x*P1x)]]];
      m0 > 0 && m1 > 0] & /@ Jlist);
  Do[
    {ok, P0, P1, mu0, mu1} = cecKC[Jlist[[i]], k];
    If[ok && test[P0, P1, mu0, mu1],
      Return[{True, P0, P1, mu0, mu1, i}, Module]],
    {i, Length[Jlist]}];
  {False, Null, Null, Null, Null, Null}];

