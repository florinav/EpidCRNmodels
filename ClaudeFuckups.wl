(* ::Package:: *)

(* ClaudeFuckups.wl — Experimental / prototype functions.
   Promote to Core.wl only after user approval.
   All new functions must be tested here first.
*)

(* =========================================================================
   Per-siphon NGM functions
   NGMsiph  — per-siphon K matrices and R_k(E0)
   R0siph   — invasion numbers R_k(E_j) at boundary equilibria
   invTable — formatted invasion table
   ========================================================================= *)

ClearAll[NGMsiph, R0siph, invTable];

(* -------------------------------------------------------------------------
   NGMsiph[RN_, rts_, var_]
   Computes one NGM per minimal siphon (at the DFE), unlike NGMRN which
   uses the union of all siphon variables and returns a single global K.

   Algorithm (per siphon Σ_k):
     Jsig_k  = D(RHS[Σ_k], Σ_k) |_{Σ_k=0}     (transversal Jacobian block)
     F_k     = posM[Jsig_k]                      (positive part, free vars)
     V_k     = posM[-Jsig_k]                     (non-positive negated, free vars)
     K_k     = F_k . Inverse[V_k]               (symbolic in free vars)
     R_k     = ρ(K_k)                            (symbolic spectral abscissa)

   Only Σ_k itself is zeroed. All other variables (S, U, other-strain vars)
   remain free. Use R0siph to evaluate at a specific equilibrium.

   Returns {mSi, JsigList, Klist, ChpList, R0list}
     mSi      — minimal siphons (list of string lists)
     JsigList — symbolic Jsig_k matrices (only Σ_k zeroed)
     Klist    — symbolic K_k matrices (functions of free variables)
     ChpList  — characteristic polynomials CharacteristicPolynomial[K_k, \[Lambda]]
     R0list   — symbolic R_k (spectral abscissae, functions of free variables)
   ------------------------------------------------------------------------- *)
NGMsiph[RN_, rts_, var_] :=
  Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, mS, mSi, cDFE, cE0,
          JsigList, Klist, ChpList, R0list, res},
    {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
    RHS  = gam . rts;
    par  = Par[RHS, var];
    cp   = Thread[par > 0];
    mS   = minSiph[var, RN, RHS];
    mSi  = mS[[1]];
    cDFE = mS[[2]];
    cE0  = mS[[3]];

    res = Table[
      Module[{invSyms, invIdx, JxSig, V1sig, F1sig, Jsig, Fsig, Vsig, Kk, Chp, Rk, lam, survVars},
        (* Force Global context: avoids EpidCRN`Private`lam$nnn in any output *)
        lam     = Symbol["Global`\[Lambda]"];
        invSyms = ToExpression /@ si;
        invIdx  = Flatten[Position[var, #] & /@ invSyms];
        survVars = Complement[var, invSyms];
        (* Two-step F-V split matching original NGM (EpidCRNo lines 1733-1736):
           Step 1: zero non-siphon (surv) vars → kills FOI terms, keeps pure transitions
           Step 2: Jx + V1 with siphon zeroed → transitions cancel, FOI survives *)
        JxSig   = Outer[D, RHS[[invIdx]], invSyms];
        V1sig   = -JxSig /. Thread[survVars -> 0];
        F1sig   = (JxSig + V1sig) /. Thread[invSyms -> 0];
        Jsig    = JxSig /. Thread[invSyms -> 0];
        Fsig    = posM[F1sig];
        Vsig    = Fsig - Jsig;
        Kk      = Simplify[Fsig . Inverse[Vsig]];
        Chp     = Factor[CharacteristicPolynomial[Kk, lam]];
        (* rhoFac: strips zero eigenvalues, no Root objects for degree ≥ 3 *)
        Rk      = rhoFac[Kk, lam];
        {Jsig, Kk, Chp, Rk}
      ],
      {si, mSi}
    ];

    JsigList = res[[All, 1]];
    Klist    = res[[All, 2]];
    ChpList  = res[[All, 3]];
    R0list   = res[[All, 4]];
    {mSi, JsigList, Klist, ChpList, R0list}
  ];

(* -------------------------------------------------------------------------
   R0siph[mSi_, JsigList_, var_, EA_]
   Evaluates per-siphon spectral abscissa at each boundary equilibrium.
   EA is the output of bdFp (indexed by siphon: EA[[j]] = sols on Σ_j=0 face).

   Entry [k, j] of the result = R_k(E_j): spectral abscissa of K_k
   evaluated at the boundary equilibrium on the Σ_j=0 face.
   Biologically: k = invader strain, j = resident face.
   Returns ns x ns matrix (ns = number of minimal siphons).
   ------------------------------------------------------------------------- *)
R0siph[mSi_, JsigList_, var_, EA_] :=
  Module[{ns = Length[mSi]},
    Table[
      Module[{ej, ejRules, invSyms, Jsig, Fsig, FsigSym, Vsig, Kk, Rk, survVars},
        ej = EA[[j]];
        (* EA[[j]] is a list of complete solutions (rule lists) *)
        ejRules = If[MatchQ[ej, {{__Rule}..}], First[ej],
                  If[MatchQ[ej, {__Rule}],     ej,
                     {}  (* no solution available *)
                  ]];
        If[ejRules === {},
          "no eq",
          invSyms  = ToExpression /@ mSi[[k]];
          Jsig     = JsigList[[k]];
          (* same survVar-based F-V split as NGMsiph: classify on symbolic Jsig, then substitute *)
          survVars = Complement[var, invSyms];
          FsigSym  = posM[Map[If[FreeQ[#, Alternatives @@ survVars], 0, #] &, Jsig, {2}]];
          Fsig     = FsigSym /. ejRules;
          Vsig     = (FsigSym - Jsig) /. ejRules;
          Kk       = Simplify[Fsig . Inverse[Vsig]];
          Rk      = FullSimplify[Max[Re[Eigenvalues[Kk]]]];
          Rk
        ]
      ],
      {k, ns}, {j, ns}
    ]
  ];

(* -------------------------------------------------------------------------
   invTable[mSi_, R0list_, JsigList_, var_, EA_]
   Returns a labeled association with all reproduction/invasion numbers.

   invM[k,j] = R_k(E_j) where E_j = face equilibrium with Sigma_j = 0.
   For k = j: R_k(E_j) = invasion number of strain k (standard R_{k,inv}).
   For k ≠ j: R_k(E_j) = R_k at the face where a different strain is absent
              (should equal 1 at endemic equilibria — useful as a check).

   Returns <|"R0" -> R0list, "invM" -> invM|> where:
     R0list[[k]]  = R_k(E0)  basic reproduction number at DFE
     invM[[k,k]]  = R_k(E_k) invasion number of strain k
   ------------------------------------------------------------------------- *)
invTable[mSi_, R0list_, JsigList_, var_, EA_] :=
  Module[{invM = R0siph[mSi, JsigList, var, EA]},
    <|"R0" -> R0list, "invM" -> invM|>
  ];

(* =========================================================================
   NGMsiph2 — variant of NGMsiph that factors the characteristic polynomial
   of K_k before computing the spectral abscissa.

   For each irreducible factor of degree d:
     d = 1 : eigenvalue extracted analytically (one root)
     d = 2 : eigenvalue extracted analytically (quadratic formula, two roots)
     d ≥ 3 : R_k is returned as {"fac", FactorList[charpol]} — no Root objects

   Helpers: eigsFac, rhoFac (private, not exported)
   ========================================================================= *)

ClearAll[NGMsiph2];

(* eigsFac[p, lam]: roots of an irreducible factor of degree 1 or 2 *)
eigsFac[p_, lam_] :=
  Module[{deg = Exponent[p, lam], co},
    co = CoefficientList[p, lam];  (* {c0, c1, ...} lowest degree first *)
    Which[
      deg == 1, {-co[[1]]/co[[2]]},
      deg == 2, Module[{a = co[[3]], b = co[[2]], c = co[[1]]},
                  {(-b + Sqrt[b^2 - 4*a*c])/(2*a),
                   (-b - Sqrt[b^2 - 4*a*c])/(2*a)}],
      True, {}  (* degree 0 constant — skip *)
    ]
  ];

(* rhoFac[Kk, lam]: spectral abscissa of Kk via char pol, stripping zero eigenvalues.

   Step 1: divide out lam^m (zero eigenvalues) using PolynomialQuotient,
           until the constant term is nonzero.  Does NOT rely on FactorList
           to detect the lam^k factor.
   Step 2: check degree of the remainder (the non-trivial part):
     deg <= 2 : extract eigenvalues analytically (Simplify, not FullSimplify,
                to prevent the quadratic-formula result being rewritten as Root)
     deg >= 3 : print a message; return Factor[remainder] — no Root objects *)
rhoFac[Kk_, lam_] :=
  Module[{chp, redp, maxDeg, fac, eigs},
    chp  = CharacteristicPolynomial[Kk, lam];
    (* strip lam^m: divide by lam while constant term is zero *)
    redp = chp;
    While[SameQ[Simplify[Coefficient[redp, lam, 0]], 0] && Exponent[redp, lam] > 0,
      redp = PolynomialQuotient[redp, lam, lam]
    ];
    maxDeg = Exponent[redp, lam];
    If[maxDeg >= 3,
      Print["NGMsiph2: degree-", maxDeg, " factor — spectral abscissa not computed"];
      Factor[redp],                   (* factored poly; no Root objects *)
      fac  = FactorList[redp];
      eigs = Flatten[
        eigsFac[#[[1]], lam] & /@ Select[fac[[2;;]], Exponent[#[[1]], lam] >= 1 &]];
      Simplify[Max[eigs]]
    ]
  ];

NGMsiph2[RN_, rts_, var_] :=
  Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, mS, mSi, cDFE, cE0,
          JsigList, Klist, R0list, res, lam},
    (* Force Global context: avoids EpidCRN`Private`lam$nnn in any output *)
    lam  = Symbol["Global`\[Lambda]"];
    {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
    RHS  = gam . rts;
    par  = Par[RHS, var];
    cp   = Thread[par > 0];
    mS   = minSiph[var, RN, RHS];
    mSi  = mS[[1]];
    cDFE = mS[[2]];
    cE0  = mS[[3]];

    res = Table[
      Module[{invSyms, invIdx, JxSig, V1sig, F1sig, Jsig, Fsig, Vsig, Kk, Rk, survVars},
        invSyms  = ToExpression /@ si;
        invIdx   = Flatten[Position[var, #] & /@ invSyms];
        survVars = Complement[var, invSyms];
        JxSig    = Outer[D, RHS[[invIdx]], invSyms];
        V1sig    = -JxSig /. Thread[survVars -> 0];
        F1sig    = (JxSig + V1sig) /. Thread[invSyms -> 0];
        Jsig     = JxSig /. Thread[invSyms -> 0];
        Fsig     = posM[F1sig];
        Vsig     = Fsig - Jsig;
        Kk       = Simplify[Fsig . Inverse[Vsig]];
        Rk      = rhoFac[Kk, lam];
        {Jsig, Kk, Rk}
      ],
      {si, mSi}
    ];

    JsigList = res[[All, 1]];
    Klist    = res[[All, 2]];
    R0list   = res[[All, 3]];
    {mSi, JsigList, Klist, R0list}
  ];

(* =========================================================================
   TEST — run after loading EpidCRN` in a model with 2 siphons.
   Activate when needed; keep commented to avoid auto-run on load.

   Example (two-strain SIS with s, i1, i2):
     RN = {...}; rts = {...}; var = {s, i1, i2};
     {mSiT, JsigT, KlistT, R0T} = NGMsiph[RN, rts, var];
     Print["NGMsiph R0list=", R0T];
   ========================================================================= *)
