(* ::Package:: *)

(* Amaral SEZR rumor spreading model *)
(* S: Susceptible/Ignorant *)
(* E: Exposed (not yet contagious) *)
(* Z: Spreaders/Infected (Zombies) *)
(* R: Removed/Stiflers *)

ClearAll["Global`*"];
SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];
<<EpidCRN`;

(* Pretty print function *)
prF[ex_] := ex /. {La -> \[Lambda], be -> \[Beta], ga -> \[Gamma],
   ka -> \[Kappa], mu -> \[Mu], ep -> \[Epsilon]};

(* ============================================ *)
(* SINGLE STRAIN MODEL *)
(* ============================================ *)

Print["\n=== SINGLE STRAIN SEZR MODEL ==="];

(* ODE system from paper *)
RHS = {
  La - be*s*z - mu*s,    (* dS/dt *)
  be*s*z - ep*e - mu*e,          (* dE/dt, using ep for the rate *)
  ga*e - ka*s*z- mu*z,          (* dZ/dt *)
  ka*s*z + (1-ga)*e - mu*r      (* dR/dt *)
};

var = {s, e, z, r};

Print["RHS=", RHS // prF, " var=", var // prF];

(* Convert ODE to reaction network *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var, prF];
Print["RN=", RN // prF];
Print["Stoichiometric matrix gam=", gam // MatrixForm];

(* Boundary analysis *)
{RHSb, varb, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars} =
bdAn[RN, rts, var];
Print["DFE=", E0 // prF];
Print["Infection variables=", infVars // prF];
Print["NGM K=", K // prF // MatrixForm];
Print["R0 functions=", R0A // prF];

(* Boundary equilibria *)
Print["\n--- Boundary Equilibria ---"];
Get["bdCom1.wl"];
{RHSb, varb, parb, cpb, mSib, Jxb, Jyb, cDFEb, E0b, Kb, R0Ab, Esys, solSiph, rurSiph, nMsiph} =
  bdCom1[RHS, var, prF];



(* ============================================ *)
(* TWO STRAIN MODEL *)
(* ============================================ *)

Print["\n\n=== TWO STRAIN SEZR MODEL ==="];

(* Two competing rumors:
   S: Susceptible to both
   E1, E2: Exposed to rumor 1 or 2
   Z1, Z2: Spreaders of rumor 1 or 2
   R: Removed from either rumor *)

RHS = {
  La - be1*s*z1 - be2*s*z2 - mu*s,           (* dS/dt *)
  be1*s*z1 - ep1*e1 - mu*e1,                 (* dE1/dt *)
  be2*s*z2 - ep2*e2 - mu*e2,                 (* dE2/dt *)
  ga1*e1 - ka1*s*z1 - mu*z1,                 (* dZ1/dt *)
  ga2*e2 - ka2*s*z2 - mu*z2,                 (* dZ2/dt *)
  ka1*s*z1 + ka2*s*z2 + (1-ga1)*e1 + (1-ga2)*e2 - mu*r  (* dR/dt *)
};

var = {s, e1, e2, z1, z2, r};

prF[ex_] := ex /. {La -> \[Lambda],
   be1 -> Subscript[\[Beta], 1], be2 -> Subscript[\[Beta], 2],
   ga1 -> Subscript[\[Gamma], 1], ga2 -> Subscript[\[Gamma], 2],
   ka1 -> Subscript[\[Kappa], 1], ka2 -> Subscript[\[Kappa], 2],
   ep1 -> Subscript[\[Epsilon], 1], ep2 -> Subscript[\[Epsilon], 2],
   mu -> \[Mu],
   e1 -> Subscript[e, 1], e2 -> Subscript[e, 2],
   z1 -> Subscript[z, 1], z2 -> Subscript[z, 2]};

Print["RHS=", RHS // prF];
Print["var=", var // prF];

(* Convert to reaction network *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var, prF];
Print["RN=", RN // prF];
Print["Stoichiometric matrix=", gam // MatrixForm];

(* Boundary analysis *)
{RHSb, varb, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars} =
  bdAn[RN, rts, var];
Print["DFE=", E0 // prF];
Print["Infection variables=", infVars // prF];
Print["NGM K=", K // prF // MatrixForm];
Print["R0 functions=", R0A // prF];

(* Boundary equilibria *)
Print["\n--- Two-Strain Boundary Equilibria ---"];
{RHSb, varb, parb, cpb, mSib, Jxb, Jyb, cDFEb, E0b, Kb, R0Ab,
  Esys, solSiph, rurSiph, nMsiph} = bdCom1[RHS, var, prF];


