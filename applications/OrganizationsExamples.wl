(* ::Package:: *)

(* OrganizationsExamples.wl -- verified computations backing SubFF's
   Section 1 (COT.tex): organizations, ERCs, synergy, complementarity,
   buffering on Examples A-E and the dengue model. Every number quoted in
   COT.tex is transcribed from this file's Print output, not hand-derived.
   Run standalone: wolframscript -file OrganizationsExamples.wl *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;
Get[FileNameJoin[{$HomeDirectory, "EpidCRNmodels", "applications", "Buffering.wl"}]];

(* ---- Example A: bare two-strain shared-susceptible network ---- *)
RNa = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> 0, "S" + "I2" -> 2*"I2", "I2" -> 0};
ercA = allErc[RNa][[All, 2]];
orgA = allOrg[RNa];
{minA, nonmA} = minSiph[{S, I1, I2}, RNa];
buffA = allBuff[RNa];
Print["ercA=", ercA, " orgA=", orgA, " nOrgA=", Length[orgA]];
Print["minA=", minA, " nBuffA=", Length[buffA], " buffAm=", buffA[[All, 1]]];

(* ---- Example B: Okada.tex's S,I1,I2,P1,P2 network (buffering/siphons
   already in Okada.tex/OkadaTwoStrainBuffering.wl -- only ERC/org new) ---- *)
varB = {S, I1, I2, P1, P2};
RHSb = {La - mu*S - be1*S*I1 - be2*S*I2, (be1*S - v1)*I1, (be2*S - v2)*I2,
   si1*I1 - de1*P1, si2*I2 - de2*P2};
{RNb, rtsB, speB, alpB, betB, gamB, rnRB} = ODE2RN[RHSb, varB];
ercB = allErc[RNb][[All, 2]];
orgB = allOrg[RNb];
{minB, nonmB} = minSiph[varB, RNb];
{mBuffB, allBuffB} = buffS[RNb];
Print["ercB=", ercB, " orgB=", orgB, " minB=", minB];
Print["mBuffB=", mBuffB[[All, 1]], " disjointFromSiphon=", Intersection[mBuffB[[1, 1]], minB[[1]]]];

(* ---- Example C: coinfection synergy ---- *)
RNc = Join[RNa, {"I1" + "I2" -> "I12", "I12" -> 0}];
e1c = speClo[RNc, {"S", "I1"}];
e2c = speClo[RNc, {"S", "I2"}];
{synTc, synWc, synJc} = synQ[RNc, {e1c, e2c}];
Print["e1c=", e1c, " e2c=", e2c, " synTc=", synTc, " synWc=", synWc];

(* ---- Example D: sequential-infection complementarity ---- *)
RNd = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> "R1", "R1" -> "S1",
  "S1" + "I2" -> "I12" + "I2", "S" + "I2" -> 2*"I2", "I2" -> 0, "I12" -> 0};
ed1 = speClo[RNd, {"S", "I1"}];
gensD2 = Select[allErc[RNd], #[[1]] === ("S1" + "I2" -> "I12" + "I2") &];
ed2 = gensD2[[1, 2]];
{cmpTd, sup12d, sup21d} = cmpERC[RNd, ed1, ed2];
Print["ed1=", ed1, " ed2=", ed2, " cmpTd=", cmpTd, " sup12d=", sup12d];

(* ---- Example E: vector-host cycle, with vs. without demographic inflow ---- *)
RNe = {0 -> "SH", "SH" -> 0, 0 -> "SV", "SV" -> 0,
  "SH" + "IV" -> "IH" + "IV", "IH" -> 0,
  "SV" + "IH" -> "IV" + "IH", "IV" -> 0};
eHost = ercClo[RNe, "SH" + "IV" -> "IH" + "IV"];
eVec = ercClo[RNe, "SV" + "IH" -> "IV" + "IH"];
Print["withInflow: eHost=", eHost, " eVec=", eVec, " sameERC=", eHost === eVec];

RNe2 = {"SH" + "IV" -> "IH" + "IV", "IH" -> 0, "SV" + "IH" -> "IV" + "IH", "IV" -> 0};
eHost2 = ercClo[RNe2, "SH" + "IV" -> "IH" + "IV"];
eVec2 = ercClo[RNe2, "SV" + "IH" -> "IV" + "IH"];
{synTe, synWe, synJe} = synQ[RNe2, {eHost2, eVec2}];
Print["noInflow: eHost2=", eHost2, " eVec2=", eVec2, " joinE=", synJe, " synTe=", synTe];

(* ---- Dengue: reuse TwoStrainDengueBuffering.wl's exact 38-reaction RN
   (siphons/R0/buffering already in TwoStrainDengue.tex -- only ERC/synergy
   /complementarity new) ---- *)
RNden = {
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
distinctErcD = DeleteDuplicates[allErc[RNden][[All, 2]]];
e1den = ercClo[RNden, "S" + "V1" -> "I1" + "V1"];
e2den = ercClo[RNden, "S" + "V2" -> "I2" + "V2"];
{synTd, synWd, synJd} = synQ[RNden, {e1den, e2den}];
Print["nReactionsDen=", Length[RNden], " nDistinctErcD=", Length[distinctErcD], " synTd=", synTd];
Print["e1denHasV12=", MemberQ[e1den, "V12"], " e2denHasV12=", MemberQ[e2den, "V12"], " nWitnessD=", Length[synWd]];
