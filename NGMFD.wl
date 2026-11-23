(* ::Package:: *)

(* NGMFD.wl -- NGM for Frequency-Dependent incidence models.

   MODEL CLASS: positive ODEs where some rate is normalized by a population
   total L(x) = c.x (c a 0/1 or nonneg row vector) whose support spans BOTH
   sides of a siphon partition, e.g. Nh = S+I1+...+R (susceptible AND
   infected host classes) in a frequency-dependent term beta*U*I/Nh. This is
   NOT an exact CRN conservation law (Nh isn't invariant: dNh/dt = LamN -
   mu*Nh - death terms), but a "quasi-conserved", bounded demographic total,
   as is standard in essentially every frequency-dependent/proportionate-
   mixing incidence model (as opposed to mass-action/density-dependent).

   WHY NGM[mod,infVars] (Core.wl) FAILS on this class: it builds the
   transversal removal matrix V in two SEPARATE zeroing stages -- first
   resident vars -> 0 (via X[[infc]]->0), then siphon vars -> 0 (via
   zeroOpaque[V,dfeRules]). Neither stage alone zeroes L(x), but since L(x)
   has support on both sides, the COMPOSITION of the two stages does, and a
   surviving L(x)^(-2) term (from differentiating a ratio S*I/L(x) w.r.t. a
   variable appearing only in L(x)) becomes literal 0^(-2): Power::infy.

   FIX: evaluate the transversal Jacobian Jx at the TRUE boundary
   equilibrium in ONE simultaneous substitution (siphon vars -> 0 AND
   resident vars -> their actual solved values on that face), never at the
   all-zero point. This never creates a 0 denominator for L(x), since the
   resident side of L(x) is evaluated at its genuine positive equilibrium
   value. Reuses posM (unchanged) for the F/V split, and minSiph (unchanged)
   to locate the siphon lattice and solve for boundary equilibria. *)

(* NGMFD[mod,sigVars,faceEq] -> {JxE,F,V,K,evals,R0}
   mod = {RHS,var,par}; sigVars = species set defining the siphon face;
   faceEq = substitution rules for the RESIDENT (non-sigVars) variables at
   the boundary equilibrium on that face (sigVars themselves are set to 0
   automatically). Get faceEq from minSiph's cE0 (DFE) or from solving
   RHS==0 on a smaller sub-siphon face (see NGMFDsub below). *)
NGMFD[mod_, sigVars_, faceEq_] := Module[
  {RHS, X, par, inf, Jx, JxE, F, V, Vinv, K, evals, R0},
  RHS = mod[[1]]; X = mod[[2]]; par = If[Length[mod] >= 3, mod[[3]], {}];
  inf = Flatten[Position[X, #] & /@ sigVars];
  Jx  = Grad[RHS[[inf]], X[[inf]]];
  JxE = Simplify[Jx /. faceEq /. Thread[X[[inf]] -> 0]];
  F   = posM[JxE];
  V   = F - JxE;
  Vinv = Quiet@Check[Inverse[V], $Failed];
  If[Vinv === $Failed || !FreeQ[Vinv, Indeterminate | DirectedInfinity | ComplexInfinity],
    K = 0; evals = {0}; R0 = 0,
    K = Simplify[F . Vinv];
    evals = Simplify /@ Eigenvalues[K];
    R0 = With[{realEvals = Select[evals, FreeQ[#, Complex] &]},
           If[realEvals === {}, 0, Max[realEvals]]]
  ];
  {JxE, F, V, K, evals, R0}
];

(* NGMFDsub[mod,sigVars] -> faceEq: solve RHS==0 on the boundary face
   sigVars==0 for the remaining (resident) variables, exactly like minSiph's
   own cE0 construction but for an ARBITRARY siphon face sigVars, not just
   the total DFE siphon -- use this to feed NGMFD an invasion number at a
   one-strain (or other) boundary equilibrium, not only the DFE. *)
NGMFDsub[mod_, sigVars_] := Module[
  {RHS, X, resid, sigRule, eqs, sol},
  RHS = mod[[1]]; X = mod[[2]];
  resid = Complement[X, sigVars];
  sigRule = Thread[sigVars -> 0];
  eqs = Thread[(RHS /. sigRule)[[Flatten[Position[X, #] & /@ resid]]] == 0];
  sol = Quiet[Solve[eqs, resid], Solve::svars];
  If[sol === {} || sol === $Failed, sigRule, Join[sigRule, First[sol]]]
];

(* NGMFDRN[RN,rts,var] -> {RHS,par,mSi,cDFE,cE0,R0A}
   Drop-in analog of NGMRN for frequency-dependent models: reuses minSiph
   (unchanged) for the siphon lattice and DFE, then computes a
   per-minimal-siphon reproduction number via NGMFD using minSiph's own
   cE0 as the (correct, nonzero-resident) evaluation point. *)
NGMFDRN[RN_, rts_, var_] := Module[
  {spe, alp, bet, gam, Rv, RHS, def, par, mS, mSi, cDFE, cE0, mod, R0A},
  {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
  RHS = gam . rts;
  par = Par[RHS, var];
  mS = minSiph[var, RN, RHS];
  mSi = mS[[1]]; cDFE = mS[[2]]; cE0 = mS[[3]];
  mod = {RHS, var, par};
  R0A = Table[
    Module[{sigSyms = ToExpression /@ mSi[[k]]},
      Last[NGMFD[mod, sigSyms, cE0]]
    ],
    {k, Length[mSi]}
  ];
  {RHS, par, mSi, cDFE, cE0, R0A}
];
