(* ::Package:: *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;

(* Two-strain environmental ODE skeleton
   S  = susceptible hosts
   I1,I2 = infected hosts (strains 1,2)
   B1,B2 = environmental pathogen (strains 1,2)
   Parameters: lam=Lambda, mu, be1=beta1, be2=beta2,
               ga1=gamma1, ga2=gamma2, xi1=xi1, xi2=xi2,
               de1=delta1, de2=delta2
   NOTE: this is the ODE skeleton of the partly diffusive PDE system
         (B_j)_t = d_j Delta B_j + xi_j I_j - delta_j B_j.
   Under homogeneous Neumann BC the threshold R_j^PDE = R_j^ODE,
   so EpidCRN results transfer directly to the PDE setting. *)

var = {s, i1, i2, b1, b2};

RHS = {
  lam - mu*s - be1*s*b1 - be2*s*b2,
  be1*s*b1 - (mu + ga1)*i1,
  be2*s*b2 - (mu + ga2)*i2,
  xi1*i1 - de1*b1,
  xi2*i2 - de2*b2
};

prF[ex_] := ex /. {
  lam -> \[Lambda],
  mu  -> \[Mu],
  be1 -> Subscript[\[Beta], 1],  be2 -> Subscript[\[Beta], 2],
  ga1 -> Subscript[\[Gamma], 1], ga2 -> Subscript[\[Gamma], 2],
  xi1 -> Subscript[\[Xi], 1],    xi2 -> Subscript[\[Xi], 2],
  de1 -> Subscript[\[Delta], 1], de2 -> Subscript[\[Delta], 2],
  s   -> S,
  i1  -> Subscript[I, 1],  i2 -> Subscript[I, 2],
  b1  -> Subscript[B, 1],  b2 -> Subscript[B, 2]
};

(* \[HorizontalLine]\[HorizontalLine] CRN decomposition \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RNs[RHS, var, prF];
Print["verif RHS=", RHS - gam . rts // FullSimplify];

(* \[HorizontalLine]\[HorizontalLine] Minimal siphons \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
(* Expect two disjoint siphons: {i1,b1} and {i2,b2} *)
allSi = minSiph[var, RN];
Print["minSiphs=", allSi // prF];

(* \[HorizontalLine]\[HorizontalLine] NGM and R0 \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
{RHS2, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, FF, VV, K, R0A, infV} =
    NGMRN[RN, rts, var];
Print["E0=", cE0 // prF, "  infVars=", infV, "  R0=", R0A // prF];
Print["K=", K // prF // MatrixForm];
r1 = R0A[[1]];  r2 = R0A[[2]];
Print["R1=", r1 // prF, "  R2=", r2 // prF];

(* \[HorizontalLine]\[HorizontalLine] Invasion numbers on semitrivial faces \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
(* Strain-1 equilibrium: s1* = (mu+ga1)*de1 / (be1*xi1) = Lambda/(mu*R1)
   Linearize strain-2 subsystem around X1*:
     i2' = be2*s1**b2 - (mu+ga2)*i2
     b2' = xi2*i2 - de2*b2
   NGM gives R_{2|1} = be2*xi2*s1* / ((mu+ga2)*de2) = R2/R1            *)
s1star = lam/(mu*r1);
s2star = lam/(mu*r2);
r21 = r2/r1 // FullSimplify;   (* invasion of strain 2 into strain-1 face *)
r12 = r1/r2 // FullSimplify;   (* invasion of strain 1 into strain-2 face *)
Print["s1*=", s1star // prF,  "  s2*=", s2star // prF];
Print["R_{2|1}=", r21 // prF, "  R_{1|2}=", r12 // prF];
Print["R_{2|1}>1 iff R2>R1: ", Simplify[r21 - 1] // prF];

(* \[HorizontalLine]\[HorizontalLine] Child-selection and oscillation exclusion \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
{n, m} = {Length[var], Length[RN]};
RM = symRM[alp, n, m];
co = Rest @ Reverse @ CoefficientList[
         Expand[CharacteristicPolynomial[gam . RM, lam]], lam];
nCS = negCS[co, gam, n];
upf = Select[nCS, upfQ];
{gH, edges} = hasseCS[nCS];
Print["negCS: n=", n, " m=", m, "  #nCS=", Length[nCS],
      "  #UPF=", Length[upf], "  #Hasse-edges=", Length[edges]];
riCores  = oscRI[nCS, gH, spe, rts];
riiCores = oscRII[nCS, gH, spe, rts];
Print["oscRI #cores=", Length[riCores],
      "  oscRII #cores=", Length[riiCores]];



