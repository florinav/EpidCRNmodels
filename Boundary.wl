(* ::Package:: *)

sta[pol_, par_] := Module[{
   factors, processedFactors, linearConditions, quadraticConditions, 
   higherFactors, linearFactors, quadraticFactors, deg, coeff, const, a, b, c, factor, var
   },
  
  (* Find the variable - it's in Variables[pol] but not in par *)
  var = First[Complement[Variables[pol], par]];
  
  (* Factor the polynomial completely *)
  factors = Factor[pol];
  
  (* Convert to list of factors, handling both single factors and products *)
  processedFactors = If[Head[factors] === Times, List @@ factors, {factors}];
  
  (* Remove constant factors (those not involving var) *)
  processedFactors = Select[processedFactors, !FreeQ[#, var] &];
  
  (* Initialize condition and factor lists *)
  linearConditions = {};
  quadraticConditions = {};
  higherFactors = {};
  linearFactors = {};
  quadraticFactors = {};
  
  (* Process each factor by degree *)
  Do[
   deg = Exponent[factor, var];
   
   Which[
    deg == 1,
    (* Linear factors: a*var + b *)
    coeff = Coefficient[factor, var, 1]; 
    const = Coefficient[factor, var, 0]; 
    AppendTo[linearConditions, const > 0];
    AppendTo[linearFactors, factor],
    
    deg == 2,
    (* Quadratic factors: a*var^2 + b*var + c *)
    a = Coefficient[factor, var, 2];
    b = Coefficient[factor, var, 1]; 
    c = Coefficient[factor, var, 0];
    AppendTo[quadraticConditions, b > 0];
    AppendTo[quadraticConditions, c > 0];
    AppendTo[quadraticFactors, factor],
    
    deg >= 3,
    (* Higher degree factors *)
    AppendTo[higherFactors, factor]
    ],
   
   {factor, processedFactors}
   ];
  
  (* Return as list {lSta, qSta, hDeg, linearFactors, quadraticFactors} *)
  {linearConditions, quadraticConditions, higherFactors, linearFactors, quadraticFactors}
 ]

ClearAll[bdFp];
bdFp[RHS_, var_, mSi_] := Module[{allInfectionVars, EA},
  (* Compute boundary equilibria for each siphon *)
  allInfectionVars = ToExpression[Union[Flatten[mSi]]];
  EA = {};

  Do[
    Module[{siphonVars, cEj, RHSEj, eqEj, varEj, solutions, completeSolutions, positiveSolutions},
      (* Set variables in siphon j to 0 *)
      siphonVars = ToExpression[mSi[[j]]];
      cEj = Thread[siphonVars -> 0];
      RHSEj = RHS /. cEj;
      eqEj = Thread[RHSEj == 0];
      varEj = Complement[var, siphonVars];

      (* Solve for non-siphon variables *)
      solutions = Solve[eqEj, varEj];
      (* Create complete solutions by adding siphon variables = 0 *)
      (* Ensure no duplicate rules by keeping first occurrence of each variable *)
      completeSolutions = Table[
        DeleteDuplicatesBy[Join[cEj, sol], First],
        {sol, solutions}
      ];

      (* Filter using onlyNN: remove DFE and negative solutions *)
      positiveSolutions = Select[completeSolutions,
        onlyNN[#, var, allInfectionVars] &
      ];

      (* Append either positive solutions or unsolved system *)
      If[Length[positiveSolutions] > 0,
        AppendTo[EA, positiveSolutions];,
        AppendTo[EA, {eqEj, varEj}];
      ];
    ];
  , {j, Length[mSi]}];

  EA
]

(* bdAn removed - use version from epi.wl instead *)

(* bdAnC - Boundary analysis for closed systems (N conserved) *)
bdAnC[RN_, rts_, var_] := Module[{
  spe, alp, bet, gam, Rv, RHS, RHSsum, varReduced, RHSreduced,
  RNreduced, rtsReduced, speReduced, alpReduced, betReduced, gamReduced,
  result, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ngm,
  mS, mod, ngmOrig, R0Aorig, DFEclosed
  },

  (* Extract stoichiometric matrices *)
  {spe, alp, bet, gam, Rv, RHS, def} = extMat[RN, var];
  RHS = gam . rts;

  (* Check if closed system: Sum[RHS] == 0 *)
  RHSsum = Total[RHS] // Simplify;

  If[RHSsum =!= 0,
    (* Not closed - use regular bdAn *)
    Print["System not closed (Sum[RHS] != 0), using regular bdAn"];
    Return[bdAn[RN, rts, var]]
  ];

  (* Closed system: First compute R0 from ORIGINAL system at s=1 *)
  par = Par[RHS, var];
  cp = Thread[par > 0];
  mS = minSiph[spe, RN];
  mSi = mS[[1]];
  infVars = Union[Flatten[mSi]];

  (* DFE for closed system: infection vars = 0, s = 1 *)
  DFEclosed = Join[Thread[ToExpression[infVars] -> 0], {var[[1]] -> 1}];

  (* Compute NGM at DFE with s=1 *)
  mod = {RHS, var, par};
  ngmOrig = NGM[mod, ToExpression[infVars]];
  K = ngmOrig[[4]] /. DFEclosed;
  R0Aorig = Select[Eigenvalues[K], # =!= 0 &];

  (* Now reduce system for other analyses *)
  varReduced = Drop[var, 1];
  RHSreduced = Drop[RHS, 1] /. {var[[1]] -> 1 - Total[varReduced]} // Simplify;

  (* Convert reduced RHS to reaction network *)
  {RNreduced, rtsReduced, speReduced, alpReduced, betReduced, gamReduced} =
    ODE2RN[RHSreduced, varReduced];

  (* Apply bdAn to reduced system for mSi, Jx, Jy, etc. *)
  result = bdAn[RNreduced, rtsReduced, varReduced];
  {RHSreduced, varReduced, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ngm} = result;

  (* Override R0A with the correctly computed value from original system *)
  R0A = R0Aorig;

  (* Return results (note: var and RHS are now for reduced system, but R0A is from original) *)
  {RHSreduced, varReduced, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ngm}
];


bd2[RN_, rts_] := 
  Module[{spe, al, be, gam, Rv, RHS, def, var, par, cp, cv, ct, mS, 
    mSi, inf, mod, K, eig, R0A, cDFE, RDFE, eq0, var0, E0, EA, cEj, 
    RHSEj, eqEj, varEj, E1t, E2t, Jx, Jy, eigenSystem, eigenvals, 
    eigenvecs, nonzeroIndices, relevantEigenvals, strainAssociation, 
    sortedPairs, mSiNGM, ngm}, 
   {spe, al, be, gam, Rv, RHS, def} = extMat[RN];
   var = ToExpression[spe];
   RHS = gam . rts;
   par = Par[RHS, var];
   cp = Thread[par > 0];
   cv = Thread[var >= 0];
   ct = Join[cp, cv];
   mS = minSiph[spe, asoRea[RN]];
   mSi = Map[Flatten[Position[spe, #] & /@ #] &, mS];
   inf = Union[Flatten[mSi]];
   
   (* Compute DFE *)
   cDFE = Flatten[Thread[ToExpression[#] -> 0] & /@ mS];
   RDFE = RHS /. cDFE;
   eq0 = Thread[RDFE == 0];
   var0 = Complement[var, var[[inf]]];
   E0 = Join[Solve[eq0, var0] // Flatten, Thread[var[[inf]] -> 0]];
   
   (* Compute NGM *)
   mod = {RHS, var, par};
   ngm = NGM[mod, inf];
   Jx = ngm[[1]] // FullSimplify;
   Jy = ngm[[5]] // FullSimplify;
   K = ngm[[4]] // FullSimplify;
   
   (* Get eigenvalues and organize by strain *)
   eigenSystem = Eigensystem[K];
   eigenvals = eigenSystem[[1]];
   eigenvecs = eigenSystem[[2]];
   nonzeroIndices = {};
   Do[If[eigenvals[[i]] =!= 0, AppendTo[nonzeroIndices, i]], {i, Length[eigenvals]}];
   
   If[Length[nonzeroIndices] > 0, 
    relevantEigenvals = eigenvals[[nonzeroIndices]];
    mSiNGM = Table[Flatten[Table[Position[inf, mSi[[i]][[j]]][[1, 1]], {j, Length[mSi[[i]]]}]], {i, Length[mSi]}];
    strainAssociation = Table[Module[{strain1Nonzeros, strain2Nonzeros, evec}, 
       evec = eigenvecs[[nonzeroIndices[[i]]]];
       strain1Nonzeros = Count[evec[[mSiNGM[[1]]]], Except[0]];
       strain2Nonzeros = Count[evec[[mSiNGM[[2]]]], Except[0]];
       If[strain1Nonzeros > strain2Nonzeros, 1, 
        If[strain2Nonzeros > strain1Nonzeros, 2, i]]], {i, Length[relevantEigenvals]}];
    sortedPairs = Sort[Transpose[{strainAssociation, relevantEigenvals}]];
    R0A = sortedPairs[[All, 2]];, 
    R0A = {};];
   
   (* Compute boundary equilibria for each strain *)
   EA = {};
   Do[elj = var[[Flatten[Delete[mSi, j]]]];
    cEj = Thread[elj -> 0];
    RHSEj = RHS /. cEj;
    eqEj = Thread[RHSEj == 0];
    varEj = Complement[var, elj];
    AppendTo[EA, {eqEj, varEj}], {j, mS // Length}];
   E1t = Solve[EA[[1]][[1]], EA[[1]][[2]]];
   E2t = Solve[EA[[2]][[1]], EA[[2]][[2]]];
   
   (* Filter rational solutions *)
   isRationalSolutionQ[sol_] := FreeQ[sol, Sqrt | Power[_, Except[_Integer]] | Root];
   E1tRational = Select[E1t, isRationalSolutionQ];
   E2tRational = Select[E2t, isRationalSolutionQ];
   {RHS, var, par, cp, mSi, Jx, Jy, E0, K, R0A, EA, E1tRational, E2tRational}];



(* Old NGM - commented out, replaced by new version with optional F parameter
NGM[mod_, infVars_] :=
  Module[{dyn, X, inf, infc, Jx, Jy, Jxy, Jyx, V1, F1, F, V, K, chpx, Kd},
   dyn = mod[[1]];
   X = mod[[2]];

   (* Convert infVars to positions in X *)
   inf = Flatten[Position[X, #] & /@ infVars];

   infc = Complement[Range[Length[X]], inf];

   (* Compute Jacobian blocks *)
   Jx = Grad[dyn[[inf]], X[[inf]]];
   Jy = If[Length[infc] > 0, Grad[dyn[[infc]], X[[infc]]], {}];
   Jxy = If[Length[infc] > 0, Grad[dyn[[inf]], X[[infc]]], {}];
   Jyx = If[Length[infc] > 0, Grad[dyn[[infc]], X[[inf]]], {}];

   chpx = CharacteristicPolynomial[Jx, #] &;

   (* NGM computation at DFE (all infection variables = 0) *)
   V1 = -Jx /. Thread[X[[infc]] -> 0];
   F1 = Jx + V1 /. Thread[X[[inf]] -> 0];
   F = posM[F1];
   V = F - Jx;

   K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // FullSimplify;
   Kd = (Inverse[V] . F) /. Thread[X[[inf]] -> 0] // FullSimplify;

   {Jx, F, V, K, Jy, Jxy, Jyx,  Kd}
  ]


NGM[mod_, infVars_, Fuser_: {}] :=
  Module[{dyn, X, par, inf, infc, Jx, Jy, Jxy, Jyx, V1, F1, F, V, K, Kd, Vinv, allPos, dfeRules},
   dyn = mod[[1]];
   X = mod[[2]];
   par = If[Length[mod] >= 3, mod[[3]], {}];
   inf = Flatten[Position[X, #] & /@ infVars];
   infc = Complement[Range[Length[X]], inf];
   Jx = Grad[dyn[[inf]], X[[inf]]];
   Jy = If[Length[infc] > 0, Grad[dyn[[infc]], X[[infc]]], {}];
   Jxy = If[Length[infc] > 0, Grad[dyn[[inf]], X[[infc]]], {}];
   Jyx = If[Length[infc] > 0, Grad[dyn[[infc]], X[[inf]]], {}];
   V1 = If[Length[infc] > 0, (-Jx) /. Thread[X[[infc]] -> 0], -Jx];
   F1 = (Jx + V1) /. Thread[X[[inf]] -> 0];
   If[Fuser === {},
     F = posM[F1];
     V = F - Jx;,
     V = Fuser - Jx;
     Vinv = Inverse[V] /. Thread[X[[inf]] -> 0];
     allPos = AllTrue[Flatten[Vinv], onlyP];
     If[allPos,
       F = Fuser;
       Print["NGM: user F accepted (Inverse[V] all positive)"],
       Print["NGM: user F rejected (Inverse[V] not all positive), using posM"];
       F = posM[F1];
       V = F - Jx;
     ];
   ];
   K = (F . Inverse[V]) /. Thread[X[[inf]] -> 0] // Simplify;
   Kd = (Inverse[V] . F) /. Thread[X[[inf]] -> 0] // Simplify;
   {Jx, F, V, K, Jy, Jxy, Jyx, Kd}
  ]
  *)


(* Helper function to compute infection indices directly from variable names *)
getInfectionIndices[variables_, siphonExpressions_] := 
  Module[{infectionVars, indices},
   (* Extract all infection variables from siphons *)
   infectionVars = Flatten[siphonExpressions];
   (* Find their positions in the variable list *)
   indices = Flatten[Position[variables, #] & /@ infectionVars];
   indices
  ];



DFE[mod_, inf_ : {}, cn_ : {}] := 
  Module[{dyn, X}, 
   dyn = mod[[1]] /. cn;
   X = mod[[2]];
   Quiet[Solve[Thread[dyn == 0] /. Thread[X[[inf]] -> 0], X]]];

mRts[RN_, ks_] := 
  Module[{rts, spe, al, be, var}, 
   {spe, al, be} = extMat[RN][[{1, 2, 3}]];
   var = ToExpression[spe];
   rts = Table[ks[[i]]*Product[var[[j]]^al[[j, i]], {j, Length[var]}], {i, Length[RN]}];
   rts];

JR0[pol_, u_] := 
  Module[{co, co1, cop, con, R0J}, 
   co = CoefficientList[pol, u];
   Print["the factor has degree ", Length[co] - 1];
   Print["its leading coefficient is ", co[[Length[co]]]];
   co1 = Expand[co[[1]]];
   Print["its constant coefficient is ", co1];
   cop = Replace[co1, _. _?Negative -> 0, {1}];
   con = cop - co1;
   Print["R0J is"];
   R0J = con/cop // FullSimplify;
   {R0J, co}];

extHD[poly_, var_] := 
  Module[{factored, factors, highDegree, linear}, 
   factored = Factor[poly];
   factors = If[Head[factored] === Times, List @@ factored, {factored}];
   Print[factors // Length, " factors: ", factors];
   highDegree = Collect[#, var, Simplify] & /@ 
     Select[factors, PolynomialQ[#, var] && Exponent[#, var] >= 2 &];
   linear = Collect[#, var, Simplify] & /@ 
     Select[factors, 
      PolynomialQ[#, var] && Exponent[#, var] == 1 && 
       MemberQ[List @@ Expand[# /. var -> 0], _?Negative, Infinity] &];
   Print["High degree factors (degree >= 2): ", highDegree];
   Print["Linear factors with possibly negative constant terms: ", linear];
   {highDegree, linear}];



(*bd1*)
bd1[RN_, rts_] := 
  Module[{spe, al, be, gam, Rv, RHS, def, var, par, cp, cv, ct, mS, 
    mSi, inf, mod, K, Jx, Jy, R0, R0A, E0, EA, E1, ngm, fps, 
    isRationalSolutionQ, isDFEQ}, 
   {spe, al, be, gam, Rv, RHS, def} = extMat[RN];
   var = ToExpression[spe];
   RHS = gam . rts // FullSimplify;
   par = Par[RHS, var];
   cp = Thread[par > 0];
   cv = Thread[var >= 0];
   ct = Join[cp, cv];
   mS = minSiph[spe, asoRea[RN]];

mSi=Map[Flatten[Position[var, #] & /@ #] &, mS]; 
inf=Union[Flatten[mSi]];
   
           
(* Compute DFE *)
   cDFE = Flatten[Thread[ToExpression[#] -> 0] & /@ mS];
   RDFE = RHS /. cDFE;
   eq0 = Thread[RDFE == 0];
   var0 = Complement[var, var[[inf]]];
   so0 = Solve[eq0, var0];
   E0 = Join[so0 // Flatten, Thread[var[[inf]] -> 0]];
   
  
   mSi=mS;(*simplest and safest*)
inf=Union[Flatten[mSi]];

   (* Compute NGM *)
   mod = {RHS, var, par};
   ngm = NGM[mod, inf];
   Jx = ngm[[1]] // FullSimplify;
   Jy = ngm[[5]] // FullSimplify;
   K = ngm[[4]] // FullSimplify;
   eig = Eigenvalues[K];
   R0A = Select[eig, (# =!= 0) &];
   
   (* Compute boundary equilibrium for single strain *)
   EA = {};
   RHSEj = RHS;
   eqEj = Thread[RHSEj == 0];
   varEj = var;
   AppendTo[EA, {eqEj, varEj}];
   
   (* Solve for fixed points *)
   fps = Solve[EA[[1]][[1]], EA[[1]][[2]]];
   
    mS = minSiph[spe, asoRea[RN]];

mSi=Map[Flatten[Position[var, #] & /@ #] &, mS]; 
inf=Union[Flatten[mSi]];

   (* Helper functions *)
   isRationalSolutionQ[sol_] := FreeQ[sol, Sqrt | Power[_, Except[_Integer]] | Root];
   isDFEQ[sol_] := Module[{infectionVars, vals}, 
     infectionVars = var[[inf]];
     vals = Simplify[infectionVars /. sol, cp];
     And @@ ((# === 0) &) /@ vals];
   
   (* Filter for non-DFE rational solutions *)
   E1NonDFE = Select[fps, (! isDFEQ[#]) &];
   E1 = Select[E1NonDFE, isRationalSolutionQ];
   {RHS, var, par, cp, mSi, Jx, Jy, E0, K, R0A, EA, E1}];
   
   


invN[E1_, E2_, R0A_, E0_, par_, cp_, fval_: {}, ins_: {}] := 
 Block[{R12, R21, coP, parSec, csi, findInstanceResult, nR0A},
  
  nR0A = Length[R0A];
  Print["R0A has ", nR0A, " elements"];
  
  R12 = Which[
    E1 === "nonRat" || E1 === F || !ListQ[E1] || Length[E1] == 0, "nonRat",
    nR0A >= 2, R0A[[2]] /. E1 // Factor,
    True, "nonRat"
    ];
  
  R21 = Which[
    E2 === "nonRat" || E2 === F || !ListQ[E2] || Length[E2] == 0, "nonRat",
    nR0A >= 1, R0A[[1]] /. E2 // Factor,
    True, "nonRat"
    ];
  
  Print["Invasion numbers: R12 = ", R12, ", R21 = ", R21];
  
  If[R12 === "nonRat" || R21 === "nonRat",
   coP = "unknown";
   Print["At least one equilibrium irrational - coexistence analysis not possible"],
   
   Print["Both equilibria rational - computing coexistence conditions"];
   
   If[Length[ins] == 0,
    coP = FindInstance[
      Join[cp, {R12 > 1, R21 > 1}, 
       Table[1 < (R0A[[k]] /. E0), {k, nR0A}],
       Thread[par != 1]], par];
    coP = If[coP === {}, "none found", Flatten[coP]],
    
    csi = Thread[par[[ins]] -> fval];
    parSec = Delete[par, List /@ ins];
    Print["Fixing parameters: ", csi];
    
    findInstanceResult = FindInstance[
      Join[Delete[cp, List /@ ins], 
       {(R12 /. csi) > 1, (R21 /. csi) > 1}, 
       Table[1 < (R0A[[k]] /. E0 /. csi), {k, nR0A}],
       Thread[parSec != 1]], parSec];
    
    Print["Thread result: ", Thread[parSec != 1]];
    Print["Manual result: ", Table[parSec[[i]] != 1, {i, Length[parSec]}]];
    
    coP = If[findInstanceResult === {}, "none found", 
      Join[Flatten[findInstanceResult], csi]];
    ];
   Print["Coexistence parameters: ", coP];
   ];
  
  {R12, R21, coP}
  ]


  (*invNr-computes invasion numbers and 
persistence condition*)
invN2[E1t_,E2t_,R0A_,E0_,par_,cp_,in1_,in2_,fval_:{},ins_:{}]:=
Module[{E1,E2,R12,R21,coP,parSec,csi,findInstanceResult},
(*Validate indices*)If[!(IntegerQ[in1]&&1<=in1<=Length[E1t]),Print["Error: in1 must be an integer between 1 and ",Length[E1t]];
Return[$Failed];];
If[!(IntegerQ[in2]&&1<=in2<=Length[E2t]),Print["Error: in2 must be an integer between 1 and ",Length[E2t]];
Return[$Failed];];
(*Select solutions based on user indices*)E1=E1t[[in1]]//Factor;
E2=E2t[[in2]]//Factor;
Print["Selected sol when i1=0 is (solution ",in1,"): ",E1];
Print["Selected sol when i2=0 is (solution ",in2,"): ",E2];
R12=R0A[[2]]/. E2//Factor;  (* strain 1 function at strain 2 equilibrium *)
R21=R0A[[1]]/. E1//Factor;  (* strain 2 function at strain 1 equilibrium *)
(*Print["invasion numbers are",{R12,R21}];
Handle parameter fixing based on ins*)If[ins==={}||Length[ins]==0,(*Case:no parameters to fix*)coP=FindInstance[Join[cp,{(R12)>1,(R21)>1,1<(R0A[[2]]/. E0)<(R0A[[1]]/. E0)}],par]//Flatten,(*Case:some parameters to fix*)csi=Thread[par[[ins]]->fval];
parSec=Delete[par,List/@ins];
Print["Fixing parameters at positions: ",ins," by csi",csi," leaves ",parSec];
(*Use FindInstance with constraints*)
findInstanceResult=
FindInstance[Join[Delete[cp,List/@ins],
{(R12/. csi)>1,(R21/. csi)>1,1<(R0A[[1]]/. E0/. csi)<(R0A[[2]]/. E0/. csi)}],parSec];
coP=Join[findInstanceResult//Flatten,csi];];
Print["under coP: ",coP," invasion nrs are",{R12,R21}/. coP//N," repr nrs are",
{(R0A[[1]]/. E0),(R0A[[2]]/. E0)}/. coP//N];
Print["END invNr OUTPUT"];
{E1,E2,R12,R21,coP}];



onlyNonneg[sol_, var_] := Module[{vals, defNeg},
  vals = var /. sol;
  (* Check if expression is DEFINITELY negative (all num coeffs neg AND all den coeffs pos, or vice versa *)
  defNeg[expr_] := Module[{num, den, numNeg, denNeg, numPos, denPos},
    {num, den} = {Numerator[#], Denominator[#]}& @ Together[expr];
    numNeg = onlyN[num];
    numPos = onlyN[-num];  (* Check if all coeffs positive *)
    denNeg = onlyN[den];
    denPos = onlyN[-den];
    (* Definitely negative if: (all num neg AND all den pos) OR (all num pos AND all den neg) *)
    (numNeg && denPos) || (numPos && denNeg)
  ];
  (* Only reject if definitely negative *)
  NoneTrue[vals, (# =!= 0 && defNeg[#]) &]
];

onlyNN[sol_, var_, infVars_] := Module[{vals, allInfZero},
  vals = FullSimplify[var /. sol];
  (* Check no non-zero component has all-negative coefficients *)
  If[!onlyNonneg[sol, var], Return[False]];
  (* Check not all infection variables are zero (DFE) *)
  allInfZero = And @@ ((# === 0) & /@ (infVars /. sol));
  !allInfZero
];
