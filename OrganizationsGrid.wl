(* ::Package:: *)

(* OrganizationsGrid.wl -- every method (allErc, allOrg, minSiph, synQ,
   cmpERC, allBuff/buffS) applied to every example (A-E, OSN, dengue),
   backing SubFF's Summary.tex two summary tables: Table
   \ref{tab:summary-structural} (organizations, siphons, buffering -- all
   listed explicitly, orgSets/siphons always, buffSets when count<5) and
   Table \ref{tab:summary-erc} (ERC count, synergy, complementarity).
   Run standalone: wolframscript -file OrganizationsGrid.wl *)

ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;
Get[FileNameJoin[{$HomeDirectory,"Dropbox","EpidCRNmodels","Buffering.wl"}]];

RNa = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> 0, "S" + "I2" -> 2*"I2", "I2" -> 0};
nErcA = Length[DeleteDuplicates[allErc[RNa][[All, 2]]]];
orgA = allOrg[RNa];
{minA, nonmA} = minSiph[{S, I1, I2}, RNa];
buffA = allBuff[RNa];
eA1 = speClo[RNa, {"S", "I1"}]; eA2 = speClo[RNa, {"S", "I2"}];
{synA, synWA, joinA} = synQ[RNa, {eA1, eA2}];
{cmpA, sup12A, sup21A} = cmpERC[RNa, eA1, eA2];
Print["A: nErc=", nErcA, " orgSets=", orgA, " siphons=", minA, " nBuff=", Length[buffA] - 1];
Print["A: syn=", synA, " cmp=", cmpA, " sup12=", sup12A, " sup21=", sup21A];

varB = {S, I1, I2, P1, P2};
RHSb = {La - mu*S - be1*S*I1 - be2*S*I2, (be1*S - v1)*I1, (be2*S - v2)*I2, si1*I1 - de1*P1, si2*I2 - de2*P2};
{RNb, rtsB, speB, alpB, betB, gamB, rnRB} = ODE2RN[RHSb, varB];
nErcB = Length[DeleteDuplicates[allErc[RNb][[All, 2]]]];
orgB = allOrg[RNb];
{minB, nonmB} = minSiph[varB, RNb];
{mBuffB, allBuffB} = buffS[RNb];
eB1 = speClo[RNb, {"S", "I1"}]; eB2 = speClo[RNb, {"S", "I2"}];
{synB, synWB, joinB} = synQ[RNb, {eB1, eB2}];
{cmpB, sup12B, sup21B} = cmpERC[RNb, eB1, eB2];
Print["B: nErc=", nErcB, " orgSets=", orgB, " siphons=", minB, " nBuff=", Length[mBuffB]];
Print["B: syn=", synB, " cmp=", cmpB, " sup12=", sup12B];

RNc = Join[RNa, {"I1" + "I2" -> "I12", "I12" -> 0}];
nErcC = Length[DeleteDuplicates[allErc[RNc][[All, 2]]]];
orgC = allOrg[RNc];
{minC, nonmC} = minSiph[{S, I1, I2, I12}, RNc];
buffC = allBuff[RNc];
eC1 = speClo[RNc, {"S", "I1"}]; eC2 = speClo[RNc, {"S", "I2"}];
{synC, synWC, joinC} = synQ[RNc, {eC1, eC2}];
{cmpC, sup12C, sup21C} = cmpERC[RNc, eC1, eC2];
Print["C: nErc=", nErcC, " orgSets=", orgC, " siphons=", minC, " buffSets=", buffC];
Print["C: syn=", synC, " synW=", synWC, " cmp=", cmpC, " sup12=", sup12C];

RNd = {0 -> "S", "S" -> 0, "S" + "I1" -> 2*"I1", "I1" -> "R1", "R1" -> "S1",
  "S1" + "I2" -> "I12" + "I2", "S" + "I2" -> 2*"I2", "I2" -> 0, "I12" -> 0};
nErcD = Length[DeleteDuplicates[allErc[RNd][[All, 2]]]];
orgD = allOrg[RNd];
{minD, nonmD} = minSiph[{S, I1, I2, R1, S1, I12}, RNd];
buffD = allBuff[RNd];
ed1 = speClo[RNd, {"S", "I1"}];
gensD2 = Select[allErc[RNd], #[[1]] === ("S1" + "I2" -> "I12" + "I2") &];
ed2 = gensD2[[1, 2]];
{synD, synWD, joinD} = synQ[RNd, {ed1, ed2}];
{cmpD, sup12D, sup21D} = cmpERC[RNd, ed1, ed2];
Print["D: nErc=", nErcD, " orgSets=", orgD, " siphons=", minD, " nBuff=", Length[buffD] - 1];
Print["D: syn=", synD, " cmp=", cmpD, " sup12=", sup12D, " sup21=", sup21D, " buffSets=", buffD[[All, 1]]];

RNe = {0 -> "SH", "SH" -> 0, 0 -> "SV", "SV" -> 0,
  "SH" + "IV" -> "IH" + "IV", "IH" -> 0, "SV" + "IH" -> "IV" + "IH", "IV" -> 0};
nErcE = Length[DeleteDuplicates[allErc[RNe][[All, 2]]]];
orgE = allOrg[RNe];
{minE, nonmE} = minSiph[{SH, IH, SV, IV}, RNe];
buffE = allBuff[RNe];
eHost = ercClo[RNe, "SH" + "IV" -> "IH" + "IV"]; eVec = ercClo[RNe, "SV" + "IH" -> "IV" + "IH"];
Print["E(inflow): nErc=", nErcE, " orgSets=", orgE, " siphons=", minE, " nBuff=", Length[buffE] - 1, " sameERC=", eHost === eVec];

RNe2 = {"SH" + "IV" -> "IH" + "IV", "IH" -> 0, "SV" + "IH" -> "IV" + "IH", "IV" -> 0};
nErcE2 = Length[DeleteDuplicates[allErc[RNe2][[All, 2]]]];
orgE2 = allOrg[RNe2];
{minE2, nonmE2} = minSiph[{SH, IH, SV, IV}, RNe2];
buffE2 = allBuff[RNe2];
eHost2 = ercClo[RNe2, "SH" + "IV" -> "IH" + "IV"]; eVec2 = ercClo[RNe2, "SV" + "IH" -> "IV" + "IH"];
{synE2, synWE2, joinE2} = synQ[RNe2, {eHost2, eVec2}];
{cmpE2, sup12E2, sup21E2} = cmpERC[RNe2, eHost2, eVec2];
Print["E(noInflow): nErc=", nErcE2, " orgSets=", orgE2, " siphons=", minE2, " nBuff=", Length[buffE2] - 1];
Print["E(noInflow): syn=", synE2, " cmp=", cmpE2, " sup12=", sup12E2, " sup21=", sup21E2];

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
nErcDen = Length[DeleteDuplicates[allErc[RNden][[All, 2]]]];
orgDen = allOrg[RNden, Method -> "ERCGenerated"];
e1den = ercClo[RNden, "S" + "V1" -> "I1" + "V1"];
e2den = ercClo[RNden, "S" + "V2" -> "I2" + "V2"];
{synDen, synWDen, joinDen} = synQ[RNden, {e1den, e2den}];
{cmpDen, sup12Den, sup21Den} = cmpERC[RNden, e1den, e2den];
{mBuffDen, allBuffDen} = buffS[RNden];
Print["Dengue: nErc=", nErcDen, " orgSets(ERCGen)=", orgDen, " nBuffExhaustive=", Length[mBuffDen]];
Print["Dengue: syn=", synDen, " nWitness=", Length[synWDen], " cmp=", cmpDen, " sup12=", sup12Den];

(* OSN rumor-propagation model (Fakih-Halanay-Avram \cite{HA}),
   CholeraShuai13/mo3.tex: reactions from its table (rates omitted --
   only reactant/product species membership matters for these methods).
   omega>0: full network; omega=0: R->W removed, R a true sink, per the
   paper's own remark that R "does not feed back" when omega=0. *)
RNosnFull = {
  0 -> "Phi", "Phi" + "U" -> 2*"U", "U" + "W" -> 2*"W",
  "B1" + "U" -> "B1" + "U" + "S1", "B2" + "U" -> "B2" + "U" + "S2",
  "S1" -> "B1", "S2" -> "B2",
  "B1" -> "R", "B2" -> "R",
  "R" -> "W",
  "Phi" -> 0, "U" -> 0, "W" -> 0
};
speOSN = {Phi, U, S1, B1, S2, B2, R, W};

nErcF = Length[DeleteDuplicates[allErc[RNosnFull][[All, 2]]]];
orgF = allOrg[RNosnFull];
{minF, nonmF} = minSiph[speOSN, RNosnFull];
buffF = allBuff[RNosnFull];
e1F = speClo[RNosnFull, {"S1", "B1"}]; e2F = speClo[RNosnFull, {"S2", "B2"}];
{synF, synWF, joinF} = synQ[RNosnFull, {e1F, e2F}];
{cmpF, sup12F, sup21F} = cmpERC[RNosnFull, e1F, e2F];
Print["OSN(omega>0): nErc=", nErcF, " orgSets=", orgF, " siphons=", minF, " nBuff=", Length[buffF] - 1];
Print["OSN(omega>0): syn=", synF, " cmp=", cmpF, " buffSets=", buffF[[All, 1]]];

RNosn0 = DeleteCases[RNosnFull, "R" -> "W"];
nErc0 = Length[DeleteDuplicates[allErc[RNosn0][[All, 2]]]];
org0 = allOrg[RNosn0];
{min0, nonm0} = minSiph[speOSN, RNosn0];
buff0 = allBuff[RNosn0];
e10 = speClo[RNosn0, {"S1", "B1"}]; e20 = speClo[RNosn0, {"S2", "B2"}];
{syn0, synW0, join0} = synQ[RNosn0, {e10, e20}];
{cmp0, sup120, sup210} = cmpERC[RNosn0, e10, e20];
Print["OSN(omega=0): nErc=", nErc0, " orgSets=", org0, " siphons=", min0, " nBuff=", Length[buff0] - 1];
Print["OSN(omega=0): syn=", syn0, " cmp=", cmp0, " buffSets=", buff0[[All, 1]]];
