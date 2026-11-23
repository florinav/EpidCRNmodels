(* ::Package:: *)

(* =========================================================================
   Lf.wl — Face-adapted Volterra-plus-linear Lyapunov functions
   EpidCRN sub-package.

   RECIPE (RahF.tex, Def. d:FaV):
   At boundary equilibrium E_Sigma on invariant face F_Sigma
   (species in Sigma vanish, complementary set SigmaBar positive):

       V_Sigma(z) = H_Sigma(z_SigmaBar) + Sum_{j in Sigma} b_j * z_j

   where
     H_Sigma  = Volterra function for the RESIDENT system on F_Sigma
              = Sum_{k in SigmaBar} a_k * z_k^* * G(z_k / z_k^*)
              G(u) = u - 1 - Log[u]  (non-negative, zero only at u=1)
     b_j > 0  = linear penalty weights for MISSING species

   The recipe is applied bottom-up over the siphon/invariant-face lattice:
   H_Sigma for a codim-p face uses the Volterra function of the
   (n-p)-dim resident subsystem (itself the Lyapunov function for the
   face of codimension p-1).

   PIPELINE (implemented below):
     volterraG    — kernel G(u)
     faceRHS      — restrict ODE to face F_Sigma
     faceBdEq     — find positive equilibrium E_Sigma on face
     volterraH    — build H_Sigma from resident vars + equilibrium + weights
     lfBuild      — assemble V_Sigma = H_Sigma + Sum b_j z_j
     lfDot        — compute dV/dt = Grad[V,var].RHS
     lfWeights    — find b_j by solving a LINEAR SYSTEM (see below)
     lfSignCheck  — verify dV/dt <= 0 via FullSimplify
     lfPrint      — pretty-print V_Sigma and dV/dt with verdict
     lfPipeline   — top-level driver over all invariant faces

   FINDING b_j (linear system approach):
   dV_Sigma is LINEAR in the b_j, so write
       dV = C0(z) + Sum_j b_j * Cj(z)
   where C0 = derivative of H_Sigma, Cj = RHS_j (rate of missing species j).
   For dV <= 0 we need the cross terms (monomials of z_Sigma appearing in C0)
   to be cancelled by the Sum b_j Cj term.

   Algorithm (lfWeights):
   1. Build V with symbolic bVars = {bf1, bf2, ...}.
   2. Compute dotV = lfDot[V, RHS, var]; expand.
   3. Extract coefficient array of dotV in sigmaVars (the missing species):
          {c0, c1mat} = CoefficientArrays[dotV, sigmaVars]
      c0   = part of dotV independent of sigmaVars (resident Volterra terms)
      c1mat[[j]] = coefficient of sigmaVars[[j]] in dotV   (linear in bVars)
   4. Each c1mat[[j]] is a polynomial in z_SigmaBar and linear in bVars.
      Expand c1mat[[j]] in monomials of z_SigmaBar; require each monomial
      coefficient to be <= 0.  The "cancellation" monomials give a linear
      system in bVars:
          A . bVars == rhs
      where A, rhs are extracted by CoefficientArrays[c1mat[[j]], bVars].
   5. Solve with LinearSolve (or Solve if underdetermined).
      Check solutions are positive; if not, the face cannot be handled by
      this recipe with the given H_Sigma weights (aVec), and report failure.

   NOTE on call order for lfDot:
   Always substitute eqSigma into V BEFORE calling lfDot, so that the
   equilibrium values z_k^* appear as parameters rather than variables
   (otherwise Grad differentiates with respect to them).

   DEPENDENCIES (all in EpidCRN):
     Par, allT    (Core.wl)
     invFace      (Siphons.wl)
     cSol         (Core.wl)

   STATUS:
     Implemented:  volterraG, faceRHS, volterraH, lfBuild, lfDot, faceBdEq
     Stub:         lfWeights, lfSignCheck, lfPrint, lfPipeline
   ========================================================================= *)

Begin["`Private`"]

(* -------------------------------------------------------------------------
   volterraG — Volterra kernel G(u) = u - 1 - Log[u]
   Non-negative for u > 0, zero only at u = 1.
   ------------------------------------------------------------------------- *)
volterraG[u_] := u - 1 - Log[u];

(* -------------------------------------------------------------------------
   faceRHS — restrict ODE system to invariant face F_Sigma
   Inputs:
     RHS       — full n-vector of rates
     var       — full variable list (length n)
     sigmaVars — list of variables in Sigma (those that vanish on the face)
   Returns {RHSres, varRes, rules, resIdx}:
     RHSres  — RHS with sigmaVars set to 0, restricted to resident equations
     varRes  — Complement[var, sigmaVars] in var-order
     rules   — Thread[sigmaVars -> 0]
     resIdx  — positions of varRes inside var
   ------------------------------------------------------------------------- *)
faceRHS[RHS_, var_, sigmaVars_] := Module[
  {rules, varRes, resIdx, RHSres},
  rules   = Thread[sigmaVars -> 0];
  varRes  = Select[var, !MemberQ[sigmaVars, #] &];          (* preserve var order *)
  resIdx  = Flatten[Position[var, #] & /@ varRes];
  RHSres  = (RHS /. rules)[[resIdx]];
  {RHSres, varRes, rules, resIdx}];

(* -------------------------------------------------------------------------
   faceBdEq — find positive equilibria of resident system on face F_Sigma
   Inputs:
     RHS, var   — full system
     sigmaVars  — missing species
     par        — parameter list (from Par[RHS, var])
   Returns list of rule sets {z_k -> z_k^*} for resident vars at E_Sigma,
   with sigmaVars -> 0 included. Returns {} if no positive solution found.
   ------------------------------------------------------------------------- *)
faceBdEq[RHS_, var_, sigmaVars_, par_] := Module[
  {RHSres, varRes, rules, resIdx, sols, posSols},
  {RHSres, varRes, rules, resIdx} = faceRHS[RHS, var, sigmaVars];
  sols = Quiet[Solve[Thread[RHSres == 0], varRes]];
  posSols = Select[sols, AllTrue[varRes /. #, # =!= 0 &] &&
                         AllTrue[varRes /. #, !NumericQ[#] || Positive[#] &] &];
  If[posSols === {}, {},
     Map[Join[#, rules] &, posSols]]];

(* -------------------------------------------------------------------------
   volterraH — build resident Volterra function H_Sigma
   Inputs:
     varRes  — resident variable list (z_SigmaBar)
     eqSigma — replacement rules for z_k^* (e.g. from faceBdEq); must already
               contain values for varRes (substitute BEFORE calling lfDot)
     aVec    — weight list a_k > 0, one per resident variable
               default: Table[1, Length[varRes]]
   Returns H_Sigma = Sum_k a_k * z_k^* * G(z_k / z_k^*)
   NOTE: eqSigma is substituted here so that z_k^* becomes a parameter
         rather than a variable — essential for correct Grad in lfDot.
   ------------------------------------------------------------------------- *)
volterraH[varRes_, eqSigma_, aVec_:Automatic] := Module[
  {weights, zStar, k},
  weights = If[aVec === Automatic, Table[1, {Length[varRes]}], aVec];
  Sum[Module[{zk = varRes[[k]], zkS = varRes[[k]] /. eqSigma},
        weights[[k]] * zkS * volterraG[zk / zkS]],
      {k, Length[varRes]}]];

(* -------------------------------------------------------------------------
   lfBuild — assemble V_Sigma = H_Sigma + Sum_{j in Sigma} b_j * z_j
   Inputs:
     var, sigmaVars, eqSigma — as above
     aVec   — Volterra weights for resident species (default: all 1)
     bVec   — linear penalty weights for missing species
              may be symbolic (e.g. {bf1, bf2}) for lfWeights
   Returns the full Lyapunov function expression.
   IMPORTANT: eqSigma must have been substituted so z_k^* are parameters.
   ------------------------------------------------------------------------- *)
lfBuild[var_, sigmaVars_, eqSigma_, aVec_:Automatic, bVec_:Automatic] :=
  Module[{varRes, H, bWeights},
    varRes   = Select[var, !MemberQ[sigmaVars, #] &];
    H        = volterraH[varRes, eqSigma, aVec];
    bWeights = If[bVec === Automatic, Table[1, {Length[sigmaVars]}], bVec];
    H + bWeights . sigmaVars];

(* -------------------------------------------------------------------------
   lfDot — time derivative dV/dt = Grad[V, var] . RHS
   Inputs:
     V   — Lyapunov function (with eqSigma already substituted; see NOTE above)
     RHS — full ODE right-hand sides
     var — full variable list
   Returns symbolic expression for dV/dt.
   ------------------------------------------------------------------------- *)
lfDot[V_, RHS_, var_] := Grad[V, var] . RHS;

(* =========================================================================
   STUB: lfWeights — find b_j by linear system
   See plan comment at top of file for the full algorithm.
   ========================================================================= *)
lfWeights[RHS_, var_, eqSigma_, sigmaVars_, par_, cp_,
          opts:OptionsPattern[{"Method" -> "LinearSystem"}]] :=
  Module[{bVars, varRes, V, dotV, c0, c1mat, eqs, sol, bNum},
    varRes = Select[var, !MemberQ[sigmaVars, #] &];
    (* 1. symbolic b_j unknowns *)
    bVars = Table[Unique["bf"], {Length[sigmaVars]}];
    (* 2. build V with symbolic b_j and substitute eqSigma *)
    V = lfBuild[var, sigmaVars, eqSigma, Automatic, bVars] /. eqSigma;
    (* 3. compute dV/dt *)
    dotV = Expand[lfDot[V, RHS, var]];
    (* 4. extract coefficients in sigmaVars (missing species) — linear in bVars *)
    {c0, c1mat} = {CoefficientArrays[dotV, sigmaVars][[1]],
                   Normal[CoefficientArrays[dotV, sigmaVars][[2]]]};
    (* c1mat[[j]] = coefficient of sigmaVars[[j]] in dotV, linear in bVars *)
    (* 5. for each j, extract monomials in varRes from c1mat[[j]]
          and build cancellation equations in bVars.
          *** TO BE COMPLETED: currently returns Identity weights *** *)
    Print["lfWeights: linear system extraction not yet implemented; ",
          "using b_j = 1 (Identity fallback)"];
    {Table[1, {Length[sigmaVars]}], Table[1, {Length[varRes]}], "Fallback"}];

(* =========================================================================
   STUB: lfSignCheck — verify dV/dt <= 0
   ========================================================================= *)
lfSignCheck[dotV_, var_, eqSigma_, par_, cp_, tmo_:30] :=
  Module[{assm, dotVS, result},
    assm = Join[cp, Thread[var > 0]];
    dotVS = TimeConstrained[
              FullSimplify[dotV, assm],
              tmo, dotV];
    result = TimeConstrained[
               FullSimplify[dotVS <= 0, assm],
               tmo, "Unknown"];
    Which[
      result === True,  {"NonPositive", dotVS},
      result === False, {"Positive",    dotVS},
      True,             {"Unknown",     dotVS}]];

(* =========================================================================
   STUB: lfPrint — pretty-print V_Sigma and dV/dt with sign verdict
   ========================================================================= *)
lfPrint[sigmaVars_, var_, eqSigma_, V_, dotV_, signResult_, prF_:Identity] :=
  Module[{col},
    col = Switch[First[signResult],
            "NonPositive", Darker[Green],
            "Positive",    Red,
            _,             Orange];
    Print["Face Sigma=", sigmaVars // prF,
          "  E_Sigma=", eqSigma // prF];
    Print["V_Sigma = ", V // prF];
    Print["dV/dt   = ", Last[signResult] // prF];
    Print[Style["Sign verdict: " <> First[signResult], col, Bold]]];

(* =========================================================================
   STUB: lfPipeline — top-level driver over all invariant faces
   ========================================================================= *)
Options[lfPipeline] = {"MaxCodim" -> 4, "WeightMethod" -> "LinearSystem",
                       "Timeout" -> 30, "PrettyPrint" -> Identity};
lfPipeline[RHS_, var_, reactions_, opts:OptionsPattern[]] :=
  Module[{par, cp, faces, results = {}, sigmaVars, eqs, eqSigma,
          weights, bVec, aVec, V, dotV, verdict},
    par   = Par[RHS, var];
    cp    = Thread[par > 0];
    faces = invFace[reactions, OptionValue["MaxCodim"]];
    Do[Module[{sigma = face},
      sigmaVars = ToExpression /@ sigma;
      eqs = faceBdEq[RHS, var, sigmaVars, par];
      If[eqs === {}, Print["Face ", sigma, ": no positive equilibrium, skipping"]; Return[]];
      eqSigma = First[eqs];
      {bVec, aVec, wOK} = lfWeights[RHS, var, eqSigma, sigmaVars, par, cp,
                                     "Method" -> OptionValue["WeightMethod"]];
      V     = lfBuild[var, sigmaVars, eqSigma, aVec, bVec] /. eqSigma;
      dotV  = lfDot[V, RHS, var];
      verdict = lfSignCheck[dotV, var, eqSigma, par, cp, OptionValue["Timeout"]];
      lfPrint[sigmaVars, var, eqSigma, V, dotV, verdict,
              OptionValue["PrettyPrint"]];
      AppendTo[results,
        <| "Face" -> sigmaVars, "Equilibrium" -> eqSigma,
           "V" -> V, "Vdot" -> Last[verdict],
           "SignResult" -> First[verdict], "Weights" -> {aVec, bVec} |>]],
    {face, SortBy[faces, Length]}];   (* increasing codimension *)
    results];

End[]
