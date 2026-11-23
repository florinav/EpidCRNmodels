(* ::Package:: *)

(* ===== Cell 1: model (var,RHS) and derive reactions ===== *)
ClearAll["Global`*"];
SetDirectory[NotebookDirectory[]];
SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "EpidCRNmodels"}]];
<<EpidCRN`;
Get[FileNameJoin[{$HomeDirectory, "EpidCRNmodels", "applications", "Buffering.wl"}]];

var = {S, I1, I2, P1, P2};
RHS = {
  La - mu*S - be1*S*I1 - be2*S*I2,
  (be1*S - v1)*I1,
  (be2*S - v2)*I2,
  si1*I1 - de1*P1,
  si2*I2 - de2*P2
};
{RN, rts, spe, alp, bet, gam, rnR} = ODE2RN[RHS, var];
Print["spe=", spe, " nReactions=", Length[rts], " check=", Simplify[gam . rts - RHS] == Table[0, {5}]];


(* ===== Cell 2: test buffS ===== *)
{mBuffS, allBuffS} = buffS[RN];
Print["mBuffS=", mBuffS[[All, 1]], " allBuffS=", allBuffS[[All, 1]], " nFound=", Length[allBuffS]];
Print["hasP1=", MemberQ[mBuffS[[All, 1]], {"P1"}], " hasP2=", MemberQ[mBuffS[[All, 1]], {"P2"}],
  " noStrainBlock=", !MemberQ[mBuffS[[All, 1]], {"I1", "P1"}]];
