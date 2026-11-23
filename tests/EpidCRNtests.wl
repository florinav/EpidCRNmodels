(* ::Package:: *)
(* EpidCRNtests.wl -- self-tests moved out of the package on 2026-11-20.
   The package no longer runs any test when it loads. Run this file only when a
   function it tests (Utils.wl: dulP, ACM, TACM, matMu, qRate, jPoly, sosMat,
   polyP, polyMet, optMon, kTau, exVerQ, exVerPoly, cecKC, jntCec; Buffering.wl:
   bufferBifurcationTest, bufLat, latHas, latBlk, rcMat, locChk, ccfRN, ccfIdl)
   has been changed:  Get["C:\\Users\\flori\\EpidCRNmodels\\EpidCRNtests.wl"]
   The Buffering tests call Julia (ccfRN) and take about a minute. *)
Get[FileNameJoin[{$HomeDirectory, "EpidCRNmodels", "EpidCRN.wl"}]];
Begin["EpidCRN`Private`"];
(* test: Bogdanov-Takens network with p=1/x loses the y dependence entirely *)
Module[{v = {dulx, duly}, R, TJ, trm},
  R = {-dk1 dulx + dk2 + dk3 dulx duly + dk4 dulx^2 - dk5 dulx^3,
       dk2 - dk3 dulx duly};
  {TJ, trm} = dulP[R, v, {-1, 0}];
  Print["dulP tst=", Simplify[TJ - (-dk2/dulx + dulx (dk4 - dk3 - 2 dk5 dulx))] === 0,
        " nTrm=", Length[trm], " freeY=", FreeQ[TJ, duly]]];

(* test: for 2x2, ACM[A,1]=A and the trace is recovered at k=2 *)
Module[{A = {{acma, acmb}, {acmc, acmd}}},
  Print["ACM tst=", ACM[A, 1] === A,
        " tr=", Simplify[Flatten[ACM[A, 2]] === {acma + acmd}],
        " n3=", Dimensions[ACM[Array[acmx, {3, 3}], 2]] === {3, 3}]];

(* test: on the same Bogdanov-Takens network as the dulP test above, TACM at
   k=2 with the matching 1x1 monomial metric must reproduce dulP's TJ
   exactly, since for a planar system J^[2]=trJ is already scalar. *)
Module[{v = {dulx, duly}, R, J, P, B, TJ},
  R = {-dk1 dulx + dk2 + dk3 dulx duly + dk4 dulx^2 - dk5 dulx^3,
       dk2 - dk3 dulx duly};
  J = D[R, {v}]; P = {{dulx^-1}};
  B = TACM[J, 2, P, R, v];
  TJ = dulP[R, v, {-1, 0}][[1]];
  Print["TACM tst=", Simplify[B[[1, 1]] - TJ] === 0,
        " dim=", Dimensions[B], " trM=", ACM[J, 2][[1, 1]] === Tr[J]]];

(* test: for a symmetric diagonal matrix all three measures agree with the
   largest diagonal entry; for a Metzler matrix mu_1 is the max column sum. *)
Module[{Bd = DiagonalMatrix[{-1, -3}], Bm = {{-1, 2}, {3, -4}}},
  Print["matMu tst=", matMu[Bd, 1] === -1 && matMu[Bd, 2] === -1 && matMu[Bd, Infinity] === -1,
        " metz1=", matMu[Bm, 1], " metz2=", N[matMu[Bm, 2]]]];

(* test: A=diag(-1,-2) has exact best common rate eta=1 at P=I (the slower
   eigenvalue), both alone and paired with a second vertex that only flips
   the off-diagonal sign, since the antisymmetric part cancels in Jl^T P+P Jl. *)
Module[{A1 = {{-1., 0.5}, {0., -2.}}, A2 = {{-1., -0.5}, {0., -2.}}, e1, e2, P1, P2},
  {e1, P1} = qRate[{{{-1., 0.}, {0., -2.}}}];
  {e2, P2} = qRate[{A1, A2}];
  Print["qRate tst=", Abs[e1 - 1] < 10^-3 && Abs[e2 - 1] < 10^-3,
        " eta1=", e1, " eta2=", e2]];

(* test: m=1 recovers the two endpoint matrices exactly; m=2 gives 2^2=4
   vertices and the first/last match hand computation at the corner bounds. *)
Module[{v1, v2},
  v1 = jPoly[{{0, 0}, {0, 0}}, {{{1, 0}, {0, 0}}}, {{2, 5}}];
  v2 = jPoly[{{0, 1}, {-1, 0}}, {{{1, 0}, {0, 0}}, {{0, 0}, {0, 1}}}, {{-1, 1}, {2, 3}}];
  Print["jPoly tst=", v1 === {{{2, 0}, {0, 0}}, {{5, 0}, {0, 0}}},
        " n1=", Length[v1], " n2=", Length[v2]]];

(* test: M=[[1,x],[x,x^2+1]] has v.M.v=(v1+x v2)^2+v2^2, manifestly SOS;
   M=[[x,0],[0,x]] has v.M.v=x(v1^2+v2^2), odd in x, so no SOS certificate
   exists at any even-in-v degree -- must fail. *)
Module[{svx = Symbol["svx"]},
  Print["sosMat tst=", First[sosMat[{{1, svx}, {svx, svx^2 + 1}}, {svx}, 1]] === True &&
                        First[sosMat[{{svx, 0}, {0, svx}}, {svx}, 1]] === False]];

(* test: symmetric by construction; n(n+1)/2 entries times nm monomials each
   gives the expected free-variable count. *)
Module[{P0, av0},
  {P0, av0} = polyP[{Symbol["svx"]}, 2, 1];
  Print["polyP tst=", P0 === Transpose[P0] && Length[av0] === 3*2]];

(* test: at d=0 (constant metric, P_f=0 identically) polyMet must agree with
   qLMI's own single-vertex feasibility exactly, at a rate below the known
   exact bound (eta=0.9<1, both should succeed) and above it (eta=1.5>1,
   both should fail), for A=diag(-1,-2) (qLMI's own validated test system). *)
Module[{Ac = {{-1., 0.}, {0., -2.}}, fc, okA, okA2, okQ, okQ2},
  fc = {-1. Symbol["px1"], -2. Symbol["px2"]};
  okA = First[polyMet[Ac, fc, {Symbol["px1"], Symbol["px2"]}, 0, 0.9]];
  okQ = First[qLMI[{Ac}, 0.9]];
  okA2 = First[polyMet[Ac, fc, {Symbol["px1"], Symbol["px2"]}, 0, 1.5]];
  okQ2 = First[qLMI[{Ac}, 1.5]];
  Print["polyMet tst=", okA === okQ === True && okA2 === okQ2 === False]];

(* test: on the LVA counterexample kappa^*=(1,1,4,1,1/4000,1/500), the family
   with b=-1 fixed correctly has no certifying a in [-3,0] (matches the
   proven witness TJ(1/4,1/20)>=9/400 for every such a, and the seed point IS
   that witness, so this branch never even calls Reduce); the free 2D search
   correctly finds a genuine certificate (p=x^-2y^-2 is known to work here,
   so at least one exponent pair must be found); and a hand-crafted
   adversarial single-exponent case, where RHS1/x=-x^2+6x-8 needs a<-1 at
   x=3 but a>-0.524 at x=1/2 -- contradictory -- forces at least one
   cutting-plane refinement (FindInstance+resolve) before correctly
   concluding no certificate exists. *)
Module[{k1, k2, k3, k4, k5, k6, xv, yv, RHSlva, RHSadv, ok1, ok2, okAdv, solAdv, itAdv},
  {k1, k2, k3, k4, k5, k6} = {1, 1, 4, 1, 1/4000, 1/500};
  xv = Symbol["olvax"]; yv = Symbol["olvay"];
  RHSlva = {k1 xv^2 - k2 xv^3 - k3 xv yv + k4 yv^2, k3 xv yv - k4 yv^2 - k5 yv + k6};
  ok1 = First[optMon[RHSlva, {xv, yv}, {{-3, 0}, {-1, -1}}, {{1/4, 1/20}}]];
  ok2 = First[optMon[RHSlva, {xv, yv}, {{-3, 0}, {-5, 0}}, {{1/4, 1/20}}]];
  RHSadv = {xv*(-xv^2 + 6 xv - 8), 0};
  {okAdv, solAdv, itAdv} = optMon[RHSadv, {xv, yv}, {{-3, 0}, {0, 0}}, {{3, 1}}];
  Print["optMon tst=", ok1 === False && ok2 === True && okAdv === False && itAdv >= 1,
        " b1infeas=", ok1, " fullFound=", ok2, " loopIters=", itAdv]];

(* test: A=I_n gives tau=k exactly for any n,k,p (their Example 10, where
   the bound is tight); and their own worked n=3 Hopfield network (Example
   12: w_ij=1 for all i,j, phi=tanh so phi_i' ranges over (0,1], r_i=r for
   all i -- J(x)=-(1/r)I+W.diag(phi') is affine in phi', so by convexity the
   worst case over phi' in [0,1]^3 is at a corner of the cube, checked by
   brute force since n=3 is trivial) certifies 2-contraction exactly on
   their own stated threshold r<1/2, despite the network provably having
   three equilibria (0 and +-1.2447*(1,1,1)), hence not being contractive
   w.r.t. any norm at all. *)
Module[{t10, W, JHop, corners, worstTau, v1, v2},
  t10 = kTau[IdentityMatrix[3], 2, 2] === 2;
  W = ConstantArray[1, {3, 3}];
  JHop[r_, phv_] := -(1/r)*IdentityMatrix[3] + W.DiagonalMatrix[phv];
  corners = Tuples[{0, 1}, 3];
  worstTau[r_] := Max[kTau[JHop[r, #], 2, Infinity] & /@ corners];
  v1 = worstTau[49/100]; v2 = worstTau[51/100];
  Print["kTau tst=", t10 && v1 < 0 && v2 > 0,
        " r049tau=", v1, " r051tau=", v2]];

(* test: a fresh qLMI[{A},eta] result on A=diag(-1,-2) at eta=0.9 (Phase 2a's
   own validated system) exactifies to a clean rational P, confirmed both by
   exVerQ itself and, independently, by direct Eigenvalues of the exact Pex
   substituted into the same LMI, bypassing verQ's own Sylvester code. *)
Module[{Ac, feas, Pnum, label, Pex, eigP, eigLMI},
  Ac = {{-1., 0.}, {0., -2.}};
  {feas, Pnum} = qLMI[{Ac}, 0.9];
  {label, Pex} = exVerQ[Pnum, {Ac}, 0.9];
  eigP = Eigenvalues[Pex];
  eigLMI = Eigenvalues[-(Transpose[Ac].Pex + Pex.Ac + 2*0.9*Pex)];
  Print["exVerQ tst=", label === "Exact" && And @@ Thread[eigP > 0] &&
        And @@ Thread[eigLMI > 0], " Pex=", Pex, " label=", label]];

(* test: on polyMet's own nonlinear validation system (f=-x-x^3, d=2, eta=0.5,
   Phase 2b's own test case), the numeric metric exactifies to a clean
   rational polynomial, confirmed by exVerPoly itself and, independently, by
   direct numeric evaluation of M1,M2 at seven sample points spanning
   [-5,5], bypassing verPoly's own Reduce-based code entirely. *)
Module[{xv3, fN, JN, okN, PN, label, Pex, M1f, M2f, samples, chk},
  xv3 = Symbol["exvx"];
  fN = {-xv3 - xv3^3}; JN = {{-1 - 3 xv3^2}};
  {okN, PN} = polyMet[JN, fN, {xv3}, 2, 0.5];
  {label, Pex} = exVerPoly[JN, fN, {xv3}, 0.5, PN];
  M1f[x0_] := (Pex - {{1}}) /. xv3 -> x0;
  M2f[x0_] := (-(Pdot[Pex, fN, {xv3}] + Transpose[JN].Pex + Pex.JN + 2*0.5*Pex)
      - {{1}}) /. xv3 -> x0;
  samples = {-5, -2, -1, 0, 1, 2, 5};
  chk = And @@ (First[First[M1f[#]]] > 0 && First[First[M2f[#]]] > 0 & /@ samples);
  Print["exVerPoly tst=", okN && label === "Exact" && chk, " Pex=", Pex]];

(* test: three cases. (1) A=diag(-1,-2), k=2, n=2 -- trace=-3<0 so genuinely
   2-contractive, and since k=n, In(P1)=(1,0,1) is honestly indefinite, not a
   degenerate corner case; exact eigenvalues of P0,P1 confirm (0,0,2) and
   (1,0,1) exactly. (2) a non-normal triangular shear system with KNOWN exact
   eigenvalues -3/10,-9/10,-3/2 (top-2 sum=-6/5<0, so 2-contractive) but
   large off-diagonal shear (200) that makes kTau's own default identity
   metric fail outright (tau>0) -- cecKC still certifies exactly via its own
   constructed (non-identity) metric, a genuine demonstration of the
   necessary-and-sufficient LTI theorem's advantage over the compound-free
   but sufficient-only kTau. (3) diag(1,-1/10,-1/10), top-2 eigenvalue sum
   =9/10>0, genuinely NOT 2-contractive -- must fail cleanly. *)
Module[{A1 = {{-1, 0}, {0, -2}}, Ash, Anot, ok1, P0a, P1a, m0a, m1a,
   lbl1, P0e, P1e, m0e, m1e, ok2, P0b, P1b, m0b, m1b, lbl2, P0f, P1f, m0f, m1f,
   ktauFail, ok3},
  {ok1, P0a, P1a, m0a, m1a} = cecKC[N[A1], 2];
  {lbl1, P0e, P1e, m0e, m1e} = exCecKC[A1, 2, P0a, P1a, m0a, m1a];
  Ash = {{-3, 200, 0}, {0, -9, 200}, {0, 0, -15}}/10;
  {ok2, P0b, P1b, m0b, m1b} = cecKC[N[Ash], 2];
  {lbl2, P0f, P1f, m0f, m1f} = exCecKC[Ash, 2, P0b, P1b, m0b, m1b];
  ktauFail = kTau[N[Ash], 2, 2] > 0;
  Anot = DiagonalMatrix[{1, -1/10, -1/10}];
  ok3 = First[cecKC[N[Anot], 2]];
  Print["cecKC tst=", ok1 && lbl1 === "Exact" && inCnt[P0e] === {0, 0, 2} &&
        inCnt[P1e] === {1, 0, 1} && ok2 && lbl2 === "Exact" && ktauFail &&
        ok3 === False, " lbl1=", lbl1, " lbl2=", lbl2, " ktauFail=", ktauFail]];

(* test: (1) the diag(-1,-2) system perturbed by a modest shear at a second
   vertex -- close enough that A1's own exact certificate should extend by
   continuity, confirmed both by jntCec succeeding and, independently, by
   direct eigenvalue evaluation of the SAME certificate's LMI residuals at
   A2, bypassing jntCec's own test[] code. (2) the four-dimensional network
   of Phase 3's own 4Pcert.tex, restricted to the box x1,x2,x3 in [2/5,8/5]
   (eight vertices) -- exactified via exCecKC and independently confirmed
   exact at all eight vertices via posDefQ, a genuine forward-invariant-
   region certificate, not a single point. (3) the same network's full box
   [0,2]^3 -- two of its own eight vertices individually fail cecKC (a
   corner where two coordinates are equal is an exact eigenvalue tie, see
   4Pcert.tex), so no vertex-certificate could possibly extend to the
   whole box; jntCec correctly returns False rather than a spurious
   result. *)
Module[{A1 = {{-1., 0.}, {0., -2.}}, A2, ok1, P0a, P1a, m0a, m1a, m0chk, m1chk,
   Jf4, pts8, Jlist8, ok2, P0b, P1b, m0b, m1b, srcV, lbl, P0e, P1e, m0e, m1e,
   JlistExact, allExact, pts0, Jlist0, ok3},
  A2 = {{-1.02, 0.03}, {0.02, -2.02}};
  {ok1, P0a, P1a, m0a, m1a} = cecKC[A1, 2];
  m0chk = Min[Eigenvalues[-(Transpose[A2].P0a + P0a.A2 - 2*m0a*P0a)]];
  m1chk = Min[Eigenvalues[-(Transpose[A2].P1a + P1a.A2 - 2*m1a*P1a)]];
  Jf4[x1v_, x2v_, x3v_] := {{-1 - 4*x1v, 0, 0, 2}, {0, -1 - 4*x2v, 0, 2},
     {0, 0, -1 - 4*x3v, 2}, {2*x1v, 2*x2v, 2*x3v, -4}};
  pts8 = Tuples[{2/5, 8/5}, 3];
  JlistExact = Jf4 @@@ pts8;
  Jlist8 = N /@ JlistExact;
  {ok2, P0b, P1b, m0b, m1b, srcV} = jntCec[Jlist8, 2];
  {lbl, P0e, P1e, m0e, m1e} = exCecKC[JlistExact[[srcV]], 2, P0b, P1b, m0b, m1b];
  allExact = And @@ (posDefQ[-(Transpose[#].P0e + P0e.# - 2*m0e*P0e)] &&
      posDefQ[-(Transpose[#].P1e + P1e.# - 2*m1e*P1e)] & /@ JlistExact);
  pts0 = Tuples[{0, 2}, 3];
  Jlist0 = N[Jf4 @@@ pts0];
  ok3 = First[jntCec[Jlist0, 2]];
  Print["jntCec tst=", ok1 && m0chk > 0 && m1chk > 0 && ok2 && lbl === "Exact" &&
        allExact && ok3 === False, " m0chk=", m0chk, " m1chk=", m1chk,
        " allExact=", allExact]];

End[];
Get[FileNameJoin[{$HomeDirectory, "EpidCRNmodels", "applications", "Buffering.wl"}]];
(* Test: two-strain downstream model of Okada.tex/COT.tex s:COT-buff.
Expect {P1},{P2},{P1,P2} to verdict DynamicTriangularBuffer with
hurw={-de1<0} etc (matching Corollary c:downstream-no-bifurcation), and
the whole-network buffer to verdict SteadyStateLocalizationOnly (u-block
empty). Runs inline (package already loaded); no header. *)
Module[{varT, RHST, RNT, rtsT, speT, alpT, betT, gamT, rnRT, parsT, resT, expT, vrdT},
  varT = {S, I1, I2, P1, P2};
  RHST = {
    La - mu*S - be1*S*I1 - be2*S*I2,
    (be1*S - v1)*I1,
    (be2*S - v2)*I2,
    si1*I1 - de1*P1,
    si2*I2 - de2*P2
  };
  {RNT, rtsT, speT, alpT, betT, gamT, rnRT} = ODE2RN[RHST, varT];
  parsT = {La, mu, be1, be2, v1, v2, si1, si2, de1, de2};
  resT = bufferBifurcationTest[RNT, RHST, varT, parsT];
  (* Silent on success.  This runs at every <<EpidCRN`, and printing the result
     of this example there made it look like output of the user's own model. *)
  expT = {{"P1"}, {"P2"}, {"P1", "P2"}, {"I1", "I2", "P1", "P2", "S"}};
  vrdT = {"DynamicTriangularBuffer", "DynamicTriangularBuffer",
          "DynamicTriangularBuffer", "SteadyStateLocalizationOnly"};
  If[resT[[All, 1]] =!= expT || resT[[All, -1]] =!= vrdT || Not[And @@ resT[[All, 5]]],
    Print["Buffering.wl SELFTEST FAILED m=", resT[[All, 1]],
      " tri=", resT[[All, 5]], " verdict=", resT[[All, -1]]]];
];

(* Test: five-species storage model of EE.tex (KTV tutorial).  Expect the chain
   {} < {e1,e2} < all species (9 reactions) < all species + inflow (10), three
   blocks with influence inflow -> hosts -> storage, two nonconstant factors of
   det L, and localization for {e1,e2}.  Silent on success. *)
Module[{varT, RHST, RNT, rtsT, speT, latT, covT, blkT, infT, fT, locT},
  varT = {s, i, e1, e2, p};
  RHST = {La - mu s - be s p, be s p - (ga + mu) i, q i - (k12 + c1) e1 + k21 e2,
    k12 e1 - (k21 + c2) e2, c1 e1 + c2 e2 - de p};
  {RNT, rtsT, speT} = ODE2RN[RHST, varT][[1 ;; 3]];
  latT = bufLat[RNT, speT];
  covT = latHas[latT];
  {blkT, infT} = latBlk[latT, covT];
  fT = rcMat[RNT, speT][[3]];
  locT = locChk[RNT, rtsT, varT, {"e1", "e2"}];
  If[Sort /@ latT[[All, 1]] =!= {{}, {"e1", "e2"}, Sort[speT], Sort[speT]} ||
     Length /@ latT[[All, 2]] =!= {0, 4, 9, 10} || Length[blkT] =!= 3 ||
     Length[fT] =!= 2 || !locT,
    Print["Buffering.wl bufLat SELFTEST FAILED lat=", latT, " blk=", blkT, " fL=", fT, " loc=", locT]];
];

(* Test: ccfRN/ccfIdl on the same five-species model must give the bufLat
   lattice (needs julia and LayeredMixedMatrices.jl; about 15 s).  Silent on success. *)
Module[{varT, RHST, RNT, speT, ccT, idT, bsT, latT},
  varT = {s, i, e1, e2, p};
  RHST = {La - mu s - be s p, be s p - (ga + mu) i, q i - (k12 + c1) e1 + k21 e2,
    k12 e1 - (k21 + c2) e2, c1 e1 + c2 e2 - de p};
  {RNT, speT} = ODE2RN[RHST, varT][[{1, 3}]];
  ccT = ccfRN[RNT, speT];
  idT = If[ListQ[ccT[[1]]], ccfIdl[Length[ccT[[1]]], ccT[[2]]], {}];
  bsT = {Sort[speT[[Union @@ Append[ccT[[1, #, 2]], {}]]]], Sort[Union @@ Append[ccT[[1, #, 1]], {}]]} & /@ idT;
  latT = {Sort[#[[1]]], Sort[#[[2]]]} & /@ bufLat[RNT, speT];
  If[Sort[bsT] =!= Sort[latT],
    Print["Buffering.wl ccfRN SELFTEST FAILED cc=", ccT, " idl=", idT, " lat=", latT]];
];


(* Test: hsGr and invSen on the two-resource model of Rappaport, Mirollo and
   Momeni at K=7 (paper Section 6): equilibria (7,0,0,0), (2,5/2,5/2,0),
   (2,2,8/3,1/3); graph {}->{n1}->{n1,n2}; r2(E1)=1/4; dr1(E0)/dK=1/2. *)
Module[{var, pop, RHS, par, hs, sen},
  var = {r1, r2, n1, n2}; pop = {n1, n2};
  RHS = {K - ga r1 - a n1 r1 + l b n2 r2, -ga r2 - b n2 r2 + l a n1 r1, n1 ((1 - l) a r1 - m1), n2 ((1 - l) b r2 - m2)};
  par = {K -> 7, ga -> 1, a -> 1, b -> 1, l -> 1/2, m1 -> 1, m2 -> 1};
  hs = hsGr[RHS, var, pop, par];
  sen = invSen[RHS, var, pop, par, {K}, hs];
  Print["hsGr tst=", hs[[1]] === {{}, {n1}, {n1, n2}} && hs[[4]] === {{1, 2}, {2, 3}} && hs[[5]] &&
     Abs[(n2 /. hs[[3, 2]]) - 1/4] < 10^-9, " invSen tst=", Abs[sen[[1, 3, 1]] - 1/2] < 10^-9, " nEq=", Length[hs[[2]]]]];
