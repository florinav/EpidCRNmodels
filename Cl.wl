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
