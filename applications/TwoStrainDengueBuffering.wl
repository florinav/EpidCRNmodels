(* ::Package:: *)

(* ===== Cell 1: RN given directly, matching TwoStrainDengue.tex's own
   elementary-reaction split (NOT ODE2RN's auto-derivation: ODE2RN merges
   the two biologically separate events "V12 infects toward I1" and "V12
   infects toward I2" into one artifact reaction S+V12->I1+I2+S+V12 plus a
   compensating S+V12->V12, since both share the rate al*S*V12/M -- the
   paper instead treats these as two separate reactions S+V12->I1+V12 and
   S+V12->I2+V12, each consuming one S, matching its stated example for
   I1's production). M=U+V1+V2+V12, N=sum of the 10 host species are left
   as free symbols: quasi-constant normalizers, not true reactants. *)
ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;
Get[FileNameJoin[{$HomeDirectory, "EpidCRNmodels", "applications", "Buffering.wl"}]];

var = {S, I1, I2, R1, R2, S1, S2, I12, I21, R, U, V1, V2, V12};
spe = ToString /@ var;

RN = {
  0 -> "S", "S" -> 0,
  "S" + "V1" -> "I1" + "V1", "S" + "V2" -> "I2" + "V2",
  "S" + "V12" -> "I1" + "V12", "S" + "V12" -> "I2" + "V12",
  "I1" -> "R1", "I1" -> 0, "I2" -> "R2", "I2" -> 0,
  "R1" -> "S1", "R1" -> 0, "R2" -> "S2", "R2" -> 0,
  "S1" + "V2" -> "I12" + "V2", "S1" + "V12" -> "I12" + "V12", "S1" -> 0,
  "S2" + "V1" -> "I21" + "V1", "S2" + "V12" -> "I21" + "V12", "S2" -> 0,
  "I12" -> "R", "I12" -> 0, "I21" -> "R", "I21" -> 0, "R" -> 0,
  0 -> "U",
  "U" + "I1" -> "V1" + "I1", "U" + "I21" -> "V1" + "I21",
  "U" + "I2" -> "V2" + "I2", "U" + "I12" -> "V2" + "I12", "U" -> 0,
  "V1" + "I2" -> "V12" + "I2", "V1" + "I12" -> "V12" + "I12", "V1" -> 0,
  "V2" + "I1" -> "V12" + "I1", "V2" + "I21" -> "V12" + "I21", "V2" -> 0,
  "V12" -> 0
};

rts = {
  La, mu*S,
  al*S*V1/M, al*S*V2/M,
  al*S*V12/M, al*S*V12/M,
  ga*I1, mu*I1, ga*I2, mu*I2,
  nu*R1, mu*R1, nu*R2, mu*R2,
  si*S1*V2/M, si*S1*V12/M, mu*S1,
  si*S2*V1/M, si*S2*V12/M, mu*S2,
  ga*I12, (mu + de)*I12, ga*I21, (mu + de)*I21, mu*R,
  LaM,
  be*U*I1/N, be*U*I21/N,
  be*U*I2/N, be*U*I12/N, ka*U,
  be*V1*I2/N, be*V1*I12/N, ka*V1,
  be*V2*I1/N, be*V2*I21/N, ka*V2,
  ka*V12
};

RHS = {
  La - (al/M)*S*(V1 + V2 + 2*V12) - mu*S,
  (al/M)*S*(V1 + V12) - (ga + mu)*I1,
  (al/M)*S*(V2 + V12) - (ga + mu)*I2,
  ga*I1 - (nu + mu)*R1,
  ga*I2 - (nu + mu)*R2,
  nu*R1 - (si/M)*S1*(V2 + V12) - mu*S1,
  nu*R2 - (si/M)*S2*(V1 + V12) - mu*S2,
  (si/M)*S1*(V2 + V12) - (ga + mu + de)*I12,
  (si/M)*S2*(V1 + V12) - (ga + mu + de)*I21,
  ga*(I12 + I21) - mu*R,
  LaM - (be/N)*U*(I1 + I2 + I12 + I21) - ka*U,
  (be/N)*U*(I1 + I21) - (be/N)*V1*(I2 + I12) - ka*V1,
  (be/N)*U*(I2 + I12) - (be/N)*V2*(I1 + I21) - ka*V2,
  (be/N)*V1*(I2 + I12) + (be/N)*V2*(I1 + I21) - ka*V12
};

{spe2, alp, bet, gam} = extMat[RN, spe][[1 ;; 4]];
Print["nRxn=", Length[RN], " speMatch=", spe2 === spe, " check=", Simplify[gam . rts - RHS] == Table[0, {14}]];


(* ===== Cell 2: test buffS ===== *)
{mBuffS, allBuffS} = buffS[RN];
Print["mBuffS=", mBuffS[[All, 1]], " nMinimal=", Length[mBuffS], " nAll=", Length[allBuffS]];
Print["hasR=", MemberQ[mBuffS[[All, 1]], {"R"}], " hasV12=",
  AnyTrue[mBuffS[[All, 1]], MemberQ[#, "V12"] &], " allBuffS=", allBuffS[[All, 1]]];
