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
