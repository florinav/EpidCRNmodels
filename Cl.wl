(* ::Package:: *)

(* ========================================================================= *)
(* Cl.wl - Experimental/Untested Function Versions               *)
(* ========================================================================= *)
(* This subpackage contains new versions of functions that need testing     *)
(* before being promoted to the main trusted subpackages.                   *)
(* NO CONTEXT - subpackages don't need BeginPackage/EndPackage              *)
(* ========================================================================= *)


(* ========================================================================= *)
(* minSiph1: Copy of minSiph0, outputs cE0 when RHS present                 *)
(* ========================================================================= *)

minSiph1[vars_, reactions_, RHS_:Null] := Module[{
  species, alpha, beta, n, m, reactionsAsso, siphons, minimal, nonm,
  cDFEspecies, cDFE, cE0, nonSiphonVars, sol},

  species = toStr /@ vars;
  reactionsAsso = asoRea[reactions];

  n = Length[species];
  m = Length[reactions];

  alpha = ConstantArray[0, {n, m}];
  beta = ConstantArray[0, {n, m}];

  Do[
    Module[{subs, prods},
      subs = Lookup[reactionsAsso[[j]], "Substrates", {}];
      prods = Lookup[reactionsAsso[[j]], "Products", {}];
      subs = toStr /@ subs;
      prods = toStr /@ prods;

      Do[
        If[MemberQ[subs, species[[i]]], alpha[[i, j]] = 1];
        If[MemberQ[prods, species[[i]]], beta[[i, j]] = 1];
      , {i, n}];
    ];
  , {j, m}];

  siphons = Select[Subsets[species, {1, n}], isSiph[#, species, alpha, beta] &];

  minimal = Select[siphons,
    Function[s, AllTrue[siphons, Function[t, t === s || !SubsetQ[s, t]]]]
  ];
  nonm = Complement[siphons, minimal];

  (* Compute cDFE and E0 if RHS is provided *)
  cDFEspecies = DeleteDuplicates[Flatten[minimal]];
  cDFE = Thread[ToExpression[cDFEspecies] -> 0];

  If[RHS =!= Null,
    (* cE0: solve RHS==0 under cDFE constraint *)
    nonSiphonVars = Complement[vars, ToExpression[cDFEspecies]];
    If[Length[nonSiphonVars] == 0,
      (* All variables are siphon variables - nothing to solve *)
      cE0 = cDFE,
      (* Solve for non-siphon variables *)
      eqsAtDFE = RHS /. cDFE;
      (* Remove equations that became True (0==0) *)
      eqsNonTrivial = DeleteCases[eqsAtDFE, 0];
      If[Length[eqsNonTrivial] > Length[nonSiphonVars],
        (* cE0 cannot be computed: overdetermined system *)
        Print["E0 System with ", Length[nonSiphonVars], " vars and ", Length[eqsNonTrivial], " eqs does not have unique solution"];
        Print["vars=", nonSiphonVars, " eqs=", eqsNonTrivial];
        cE0 = cDFE,
        (* Solve with timeout *)
        sol = TimeConstrained[Solve[Thread[eqsNonTrivial == 0], nonSiphonVars], 2, $Failed];
        If[sol === {} || sol === $Failed,
          (* cE0 cannot be computed: Solve failed, timed out, or returned no solutions *)
          Print["E0 System with ", Length[nonSiphonVars], " vars and ", Length[eqsNonTrivial], " eqs does not have unique solution"];
          Print["vars=", nonSiphonVars, " eqs=", eqsNonTrivial];
          cE0 = cDFE,
          (* Take first solution *)
          cE0 = Join[cDFE, First[sol]]
        ]
      ]
    ];

    {minimal, cDFE, cE0, nonm},

    (* Legacy output if RHS not provided - still include cDFE *)
    {minimal, cDFE, nonm}
  ]
];


(* ========================================================================= *)
(* minSiph0: Exact copy from epi.wl minSiph, always outputs cDFE            *)
(* ========================================================================= *)

minSiph0[vars_, reactions_, RHS_:Null] := Module[{
  species, alpha, beta, n, m, reactionsAsso, siphons, minimal, nonm,
  cDFEspecies, cDFE, cE0, nonSiphonVars, sol},

  species = toStr /@ vars;
  reactionsAsso = asoRea[reactions];

  n = Length[species];
  m = Length[reactions];

  alpha = ConstantArray[0, {n, m}];
  beta = ConstantArray[0, {n, m}];

  Do[
    Module[{subs, prods},
      subs = Lookup[reactionsAsso[[j]], "Substrates", {}];
      prods = Lookup[reactionsAsso[[j]], "Products", {}];
      subs = toStr /@ subs;
      prods = toStr /@ prods;

      Do[
        If[MemberQ[subs, species[[i]]], alpha[[i, j]] = 1];
        If[MemberQ[prods, species[[i]]], beta[[i, j]] = 1];
      , {i, n}];
    ];
  , {j, m}];

  siphons = Select[Subsets[species, {1, n}], isSiph[#, species, alpha, beta] &];

  minimal = Select[siphons,
    Function[s, AllTrue[siphons, Function[t, t === s || !SubsetQ[s, t]]]]
  ];
  nonm = Complement[siphons, minimal];

  (* Compute cDFE and E0 if RHS is provided *)
  cDFEspecies = DeleteDuplicates[Flatten[minimal]];
  cDFE = Thread[ToExpression[cDFEspecies] -> 0];

  If[RHS =!= Null,
    (* cE0: solve RHS==0 under cDFE constraint *)
    nonSiphonVars = Complement[vars, ToExpression[cDFEspecies]];
    If[Length[nonSiphonVars] == 0,
      (* All variables are siphon variables - nothing to solve *)
      cE0 = cDFE,
      (* Solve for non-siphon variables *)
      eqsAtDFE = RHS /. cDFE;
      (* Remove equations that became True (0==0) *)
      eqsNonTrivial = DeleteCases[eqsAtDFE, 0];
      If[Length[eqsNonTrivial] > Length[nonSiphonVars],
        (* cE0 cannot be computed: overdetermined system *)
        Print["E0 System with ", Length[nonSiphonVars], " vars and ", Length[eqsNonTrivial], " eqs does not have unique solution"];
        Print["vars=", nonSiphonVars, " eqs=", eqsNonTrivial];
        cE0 = cDFE,
        (* Solve with timeout *)
        sol = TimeConstrained[Solve[Thread[eqsNonTrivial == 0], nonSiphonVars], 2, $Failed];
        If[sol === {} || sol === $Failed,
          (* cE0 cannot be computed: Solve failed, timed out, or returned no solutions *)
          Print["E0 System with ", Length[nonSiphonVars], " vars and ", Length[eqsNonTrivial], " eqs does not have unique solution"];
          Print["vars=", nonSiphonVars, " eqs=", eqsNonTrivial];
          cE0 = cDFE,
          (* Take first solution *)
          cE0 = Join[cDFE, First[sol]]
        ]
      ]
    ];

    {minimal, cDFE, cE0, nonm},

    (* Legacy output if RHS not provided - still include cDFE *)
    {minimal, cDFE, nonm}
  ]
];


(* clsEq1/clsEqb1 removed - use clsEq/clsEqb from epi.wl *)

(* ========================================================================= *)
(* FLaws: Four Laws of ABa models (Fall/Bonzi)                              *)
(* Inputs: A (n×n MH infected), vbe (m×n WAIFW), ba (n-vec entry probs),   *)
(*         AS (m×m MH susceptible), La (m-vec recruitment)                  *)
(* Output: {y0, mR, R0, Dw, ye, k, xe}                                     *)
(* ========================================================================= *)

FLaws[A_, vbe_, ba_, AS_, La_] := Module[
  {y0, Dw, mR, R0, kv, yeSym, kSol, ye, xe},
  y0    = LinearSolve[-AS, La];                           (* DFE: (-AS)^-1 La *)
  Dw    = LinearSolve[-A, ba];                            (* dwell times: (-A)^-1 ba *)
  mR    = FullSimplify[vbe . Dw];                         (* replacement vector *)
  R0    = FullSimplify[y0 . mR];                          (* Law 1: R0 = y0.mR *)
  yeSym = LinearSolve[kv*DiagonalMatrix[mR] - AS, La];   (* ye as function of kv *)
  Module[{ksols = Solve[FullSimplify[mR . yeSym] == 1, kv]},
    If[ksols === {}, Print["FLaws: Solve for k failed"]; Return[$Failed]];
    kSol = kv /. First[ksols]];                           (* Law 4: ye.mR=1 *)
  ye    = FullSimplify[yeSym /. kv -> kSol];              (* Law 3 *)
  xe    = FullSimplify[kSol * Dw];                        (* Law 2: xe = k*Dw *)
  Print["FLaws: y0=", y0, " mR=", mR, " R0=", R0];
  Print["FLaws: Dw=", Dw, " ye=", ye, " k=", kSol, " xe=", xe];
  {y0, mR, R0, Dw, ye, kSol, xe}
];

(* ========================================================================= *)
(* extABa: Extract ABa model matrices from a bilinear ODE                   *)
(* Input:  RHS, var (full system), yvar (susceptible), xvar (infected)      *)
(* Output: {La, AS, vbe, ba, A}                                             *)
(* Algorithm:                                                                *)
(*   vbe[i,j] = -Coeff(yRHS_i, y_i*x_j)    from y-equations               *)
(*   La = constant part of y-equations (at y=x=0)                          *)
(*   AS[i,j] = Coeff(yRHS_i - bilinear - La_i, y_j)                       *)
(*   A[k,l]  = Coeff(xRHS_k, x_l)           from x-equations               *)
(*   ba[k]   = Coeff(rem_k, y_i0*x_j0) / vbe[i0,j0]  (any nonzero entry) *)
(* ========================================================================= *)

extABa[RHS_, var_, yvar_, xvar_] := Module[
  {m, n, posY, posX, yRHS, xRHS, vbe, phi, La, AS, A, rem, ba, ij0},
  m    = Length[yvar];
  n    = Length[xvar];
  posY = Flatten[Position[var, #] & /@ yvar];
  posX = Flatten[Position[var, #] & /@ xvar];
  yRHS = RHS[[posY]];
  xRHS = RHS[[posX]];
  (* Step 1: extract vbe from bilinear terms in y-equations *)
  vbe  = Table[-Coefficient[yRHS[[i]], yvar[[i]]*xvar[[j]]], {i,m},{j,n}];
  If[AllTrue[Flatten[vbe], PossibleZeroQ],
    Message[extABa::novbe]; Return[$Failed]];
  (* Step 2: remove bilinear, extract La and AS *)
  phi  = Table[Expand[yRHS[[i]] + yvar[[i]]*(vbe[[i]] . xvar)], {i,m}];
  La   = phi /. Thread[Join[yvar,xvar] -> Table[0,{m+n}]];
  AS   = Table[Coefficient[phi[[i]] - La[[i]], yvar[[j]]], {i,m},{j,m}];
  Do[If[!PossibleZeroQ[Expand[phi[[i]] - La[[i]] - AS[[i]] . yvar]],
       Message[extABa::notABa, i, yvar[[i]], phi[[i]]-La[[i]]-AS[[i]].yvar]],
     {i,m}];
  (* Step 3: linear part A from x-equations — zero yvar first to avoid bilinear contamination *)
  A    = Table[Coefficient[(xRHS[[k]] /. Thread[yvar -> 0]), xvar[[l]]], {k,n},{l,n}];
  rem  = Table[Expand[xRHS[[k]] - A[[k]] . xvar], {k,n}];
  (* Step 4: recover ba from any nonzero vbe[i0,j0] *)
  ij0  = First[Cases[Flatten[Table[{i,j},{i,m},{j,n}],1],
                     {i_,j_}/;!PossibleZeroQ[vbe[[i,j]]]]];
  ba   = Table[
    FullSimplify[Coefficient[rem[[k]], yvar[[ij0[[1]]]]*xvar[[ij0[[2]]]]]
                / vbe[[ij0[[1]],ij0[[2]]]]],
    {k,n}];
  Print["extABa: La=", La, " ba=", ba, " vbe=", vbe];
  Print["extABa: AS=", AS, " A=", A];
  {La, AS, vbe, ba, A}
];

(* ========================================================================= *)
(* extABP: General DSI bilinear decomposition — extracts {La,AS,B,A,P}     *)
(* Input:  RHS, var (all variables), xvar (infected/x-class variables)      *)
(* yvars = complement of xvar in var (preserving order)                     *)
(* Output: {La, AS, B, A, P, yvar}                                         *)
(*   B (m×n WAIFW): B[i,j] = -Coeff(yRHS_i, y_i*x_j)                     *)
(*   A (n×n infected): A[k,l] = Coeff(xRHS_k, x_l)                        *)
(*   P (n×m insertion): P[k,i] = bilcoeff[k,i,j]/B[i,j] (first nonzero j) *)
(* rank(P)=1 iff ABa (ba*ones^t); rank(B)=1 iff rank-one WAIFW             *)
(* ========================================================================= *)

extABP[RHS_, var_, xvar_] := Module[
  {yvar, m, n, posY, posX, yRHS, xRHS, B, phi, La, AS, A, rem, bil, P,
   rowB, rowC, ratios, res, k, i, j},
  yvar = Select[var, !MemberQ[xvar, #]&];   (* preserve var order *)
  m    = Length[yvar];
  n    = Length[xvar];
  posY = Flatten[Position[var, #] & /@ yvar];
  posX = Flatten[Position[var, #] & /@ xvar];
  yRHS = Expand[RHS[[posY]]];
  xRHS = Expand[RHS[[posX]]];
  (* Step 1: B[i,j] from y-equations *)
  B   = Table[-Coefficient[yRHS[[i]], yvar[[i]]*xvar[[j]]], {i,m},{j,n}];
  phi = Table[Expand[yRHS[[i]] + yvar[[i]]*(B[[i]] . xvar)], {i,m}];
  Do[If[!PolynomialQ[phi[[i]], yvar] || Max[Exponent[phi[[i]], yvar]] > 1,
       Message[extABP::notbil, i]; Return[$Failed]], {i,m}];
  La  = phi /. Thread[yvar -> 0];
  AS  = Table[Coefficient[phi[[i]], yvar[[j]]], {i,m},{j,m}];
  (* Step 2: A[k,l] from x-equations — zero yvar first to avoid bilinear contamination *)
  A   = Table[Coefficient[(xRHS[[k]] /. Thread[yvar -> 0]), xvar[[l]]], {k,n},{l,n}];
  rem = Expand[xRHS - A . xvar];
  Do[res = Expand[rem[[k]] - Sum[Coefficient[rem[[k]], yvar[[i]]*xvar[[j]]]*yvar[[i]]*xvar[[j]], {i,m},{j,n}]];
     If[!PossibleZeroQ[res], Message[extABP::notbil, n+k]; Return[$Failed]],
     {k,n}];
  (* Step 3: bilcoeff[k,i,j] = P[k,i]*B[i,j] *)
  bil = Table[Coefficient[rem[[k]], yvar[[i]]*xvar[[j]]], {k,n},{i,m},{j,n}];
  (* Step 4: recover P *)
  P   = Table[0, {n},{m}];
  Do[
    rowB = B[[i]];
    rowC = bil[[k,i]];
    If[AllTrue[rowB, PossibleZeroQ],
      If[!AllTrue[rowC, PossibleZeroQ],
        Message[extABP::badP, k, i, rowC]; Return[$Failed]
        (* else P[k,i]=0 *)],
      ratios = DeleteDuplicates[
        Table[If[!PossibleZeroQ[rowB[[j]]], Simplify[rowC[[j]]/rowB[[j]]], Nothing], {j,n}]];
      If[Length[ratios] != 1,
        Message[extABP::badP, k, i, ratios]; Return[$Failed]];
      P[[k,i]] = First[ratios]],
    {k,n},{i,m}];
  Print["extABP: La=", La, " B=", B // MatrixForm, " P=", P // MatrixForm];
  Print["extABP: AS=", AS // MatrixForm, " A=", A // MatrixForm, " rank(P)=", MatrixRank[P], " rank(B)=", MatrixRank[B]];
  {La, AS, B, A, P, yvar}
];

(* ========================================================================= *)
(* invFQ: Invariant Facet Query - cleaner version of isInvariantFacet *)
(* ========================================================================= *)
(* MATHEMATICAL CONDITION:
   Let F = {x : x_i = 0 for i in facetVars} be a facet of the positive orthant.
   Let f(x) = RHS be the vector field (dx/dt = f(x)).

   F is FORWARD INVARIANT iff: for all i in facetVars, f_i|_F >= 0
   (i.e., when on the facet, the flow for facet variables is non-negative,
    so trajectories can't exit the positive orthant through F)

   WHY reaToRHS WORKS WITHOUT EXPLICIT RATE VALUES:
   reaToRHS[RN] constructs RHS using symbolic rates k1, k2, ... (assumed > 0).
   For mass action: rate_j = k_j * (product of substrates).
   On facet F (facet vars = 0):
   - Reactions CONSUMING facet var i: rate = 0 (because substrate i = 0)
   - Reactions PRODUCING facet var i: rate = k_j * (positive stuff) >= 0

   CONCLUSION: Invariance depends ONLY on reaction network STRUCTURE,
   not on specific rate values! Any positive rates give same answer.
*)
(* ========================================================================= *)

invFQ[facetVars_, RN_] := Module[{vf, vars, varSyms, facetIdx, facetRules, fi},
  (* Get RHS with symbolic rates from reaction network *)
  {vf, vars} = Take[reaToRHS[RN], 2];
  varSyms = ToExpression[vars];

  (* Find indices of facet variables *)
  facetIdx = Flatten[Position[vars, #] & /@ facetVars];
  If[facetIdx === {}, Return[True]]; (* Empty facet is trivially invariant *)

  (* Substitution rule: set facet variables to 0 *)
  facetRules = Thread[varSyms[[facetIdx]] -> 0];

  (* Check: for each facet var i, is f_i|_F >= 0? *)
  AllTrue[facetIdx, Function[i,
    fi = Simplify[vf[[i]] /. facetRules /. Power[0, _] -> 0];
    (* fi should be >= 0: either zero, or positive under positivity assumptions *)
    fi === 0 || PossibleZeroQ[fi] ||
    TrueQ[Simplify[fi >= 0, Assumptions -> And @@ (# >= 0 & /@ varSyms)]]
  ]]
];

(* TEST invFQ - uncomment to run
RNtest = {a -> b, b -> c, c -> a};  (* cycle: all facets invariant *)
Print["invFQ test cycle {a}=", invFQ[{"a"}, RNtest], " exp=True"];
Print["invFQ test cycle {b}=", invFQ[{"b"}, RNtest], " exp=True"];
RNtest2 = {a -> b + c};  (* {b,c} not invariant: f_b = k*a > 0 when a > 0 *)
Print["invFQ test {a}=", invFQ[{"a"}, RNtest2], " exp=True"];
Print["invFQ test {b}=", invFQ[{"b"}, RNtest2], " exp=False"];
*)

(* ========================================================================= *)
(* comp2Asso1: Handles parametric stoichiometric coefficients               *)
(* ========================================================================= *)
(* FIX: Original comp2Asso only handles Times with 2 elements.
   This version handles arbitrary products like be*et*"Vd" (3 elements).
   Extracts the string as species, everything else as coefficient. *)

comp2Asso1[expr_] := Module[{terms, result, coeff, species, stringParts, otherParts},
  If[expr === 0 || expr === Null,
    Association[],
    terms = If[Head[expr] === Plus, List @@ expr, {expr}];
    result = Association[];
    Do[
      Which[
        (* Handle Times with any number of elements *)
        Head[term] === Times,
        (* Separate strings (species) from non-strings (coefficients) *)
        stringParts = Select[List @@ term, StringQ];
        otherParts = Select[List @@ term, !StringQ[#] &];
        If[Length[stringParts] == 1,
          species = First[stringParts];
          coeff = If[Length[otherParts] == 0, 1,
                    If[Length[otherParts] == 1, First[otherParts], Times @@ otherParts]];
          result[species] = coeff,
          (* Multiple or no strings - fallback *)
          result[ToString[term]] = 1
        ],

        (* Handle string species alone (coefficient 1) *)
        StringQ[term],
        result[term] = 1,

        (* Handle symbol species alone (coefficient 1) *)
        AtomQ[term],
        result[ToString[term]] = 1,

        (* Handle other cases *)
        True,
        result[ToString[term]] = 1
      ];
    , {term, terms}];
    result
  ]
];

(* ========================================================================= *)
(* extMat1: Uses comp2Asso1 for parametric stoichiometry                    *)
(* ========================================================================= *)

extMat1[reactions_, speOrder_List:{}] := Module[{
  spe, alp, bet, gamma, Rv, RHS,
  numReactions, numSpecies, reactants, products,
  var, rv, tk, defFormula, defTerms, defResult, Nc, l, s, mSi},

  (* Extract species preserving case, use provided order if given *)
  If[Length[speOrder] > 0,
    spe = ToString[#] & /@ speOrder;
  ,
    spe = extSpe[reactions];
  ];

  If[Length[spe] == 0, Return[$Failed]];

  numSpecies = Length[spe];
  numReactions = Length[reactions];

  (* Initialize stoichiometric matrices - can contain symbolic entries *)
  alp = ConstantArray[0, {numSpecies, numReactions}];
  bet = ConstantArray[0, {numSpecies, numReactions}];

  (* Build stoichiometric matrices using comp2Asso1 *)
  Do[
    reactants = comp2Asso1[reactions[[j, 1]]];
    products = comp2Asso1[reactions[[j, 2]]];

    Do[
      If[KeyExistsQ[reactants, spe[[i]]], alp[[i, j]] = reactants[spe[[i]]]];
      If[KeyExistsQ[products, spe[[i]]], bet[[i, j]] = products[spe[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];

  gamma = bet - alp;
  var = ToExpression /@ spe;

  (* Calculate reaction rates *)
  rv = Table[
    Product[
      If[alp[[i, j]] == 0, 1,
        If[NumericQ[alp[[i, j]]], var[[i]]^alp[[i, j]], Power[var[[i]], alp[[i, j]]]]
      ],
      {i, numSpecies}
    ],
    {j, numReactions}
  ];

  tk = Table[Symbol["k" <> ToString[j]], {j, numReactions}];
  Rv = tk*rv;
  RHS = gamma . Rv;

  {Nc, l, s, mSi} = defi[reactions];

  defFormula = "\[Delta] = Nc - l - s";
  defTerms = "Nc=" <> ToString[Nc] <> ", l=" <> ToString[l] <> ", s=" <> ToString[s];
  defResult = If[mSi === "symbolic", "symbolic", ToString[mSi]];

  {spe, alp, bet, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
];

(* ========================================================================= *)
(* NGMRN1: Uses extMat1 for parametric stoichiometry                        *)
(* ========================================================================= *)

NGMRN1[RN_, rts_, var_] := Module[{spe, alp, bet, gam, Rv, RHS, def, par, cp, cv, ct, mS, mSi, inf, mod,
  K, R0A, cDFE, RDFE, eq0, var0, cE0, Jx, Jy, Jxy, fx, F, VV, eigenvals, ngm, infVarsStr, infV},

  (* Extract stoichiometric matrices with var order using extMat1 *)
  {spe, alp, bet, gam, Rv, RHS, def} = extMat1[RN, var];
  RHS = gam . rts;  (* rts = full rates (monomials * constants) *)
  par = Par[RHS, var];
  cp = Thread[par > 0];
  cv = Thread[var >= 0];
  ct = Join[cp, cv];
  mS = minSiph[var, RN, RHS];
  mSi = mS[[1]];
  cDFE = mS[[2]];
  cE0 = mS[[3]];

  (* infV: infection variable symbols (from minimal siphons) *)
  infVarsStr = Union[Flatten[mSi]];
  infV = ToExpression[infVarsStr];

  (* Call NGM with infection variables *)
  mod = {RHS, var, par};
  ngm = NGM[mod, infV];
  Jx = ngm[[1]];
  F = ngm[[2]];
  VV = ngm[[3]];
  K = ngm[[4]];
  Jy = ngm[[5]];
  Jxy = ngm[[6]];
  fx = ngm[[7]];

  If[K === 0,
    Module[{Vdfe, rnk, dim},
      Vdfe = VV /. Thread[infV -> 0];
      rnk = MatrixRank[Vdfe];
      dim = Length[Vdfe];
      Print["VV dim=", dim, " rank=", rnk, " singular"];
    ];
    R0A = {};
  ,
    eigenvals = Eigenvalues[K];
    R0A = Select[eigenvals, # =!= 0 &];
  ];

  V = VV;
  {RHS, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, F, V, K, R0A, infV}
];

(* ========================================================================= *)
(* TEST: comp2Asso1 with parametric coefficients                            *)
(* ========================================================================= *)
(*
Print["--- comp2Asso1 tests ---"];
Print["et*\"Vh\" = ", comp2Asso1[et*"Vh"]];
Print["be*et*\"Vd\" = ", comp2Asso1[be*et*"Vd"]];
Print["(et/ka)*\"Vh\" + (be+de)*(et/ka)*\"Vd\" = ",
      comp2Asso1[(et/ka)*"Vh" + (be+de)*(et/ka)*"Vd"]];
Print["\"Cs\" + \"Vh\" = ", comp2Asso1["Cs" + "Vh"]];
*)
