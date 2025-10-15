(* Diagnostic functions *)
WhereIs[funcName_String] := 
 Module[{usageSymbol, usage, match}, 
  usageSymbol = ToExpression[funcName <> "::usage"];
  usage = ToString[usageSymbol];
  If[usage === funcName <> "::usage", "Function not found", 
   match = StringCases[usage, 
     "(" ~~ source : WordCharacter .. ~~ ")" :> source, 1];
   If[match === {}, "Unknown", First[match]]]];
(* WhereIs["extMat"] *)

WhereIs[func_Symbol] := WhereIs[SymbolName[func]];
(* WhereIs[bd2] *)

ListFunctionsByPackage[] := 
 Module[{allFunctions, byPackage, validFunctions}, 
  allFunctions = 
   Select[Names["EpidCRN`*"], StringFreeQ[#, "`Private`"] &];
  validFunctions = 
   Select[allFunctions, 
    ToString[ToExpression[# <> "::usage"]] =!= (# <> "::usage") &];
  byPackage = GroupBy[validFunctions, WhereIs];
  KeyValueMap[Print[#1, ": ", StringRiffle[#2, ", "]] &, 
   KeySort[byPackage]];];
(* ListFunctionsByPackage[] *)

(* Validation functions *)
chRN[rn_List] := 
 Module[{validSpeciesQ, validTermQ, validComplexQ, validReactionQ},
  validSpeciesQ[sp_] := StringQ[sp];
  validTermQ[term_] := 
   StringQ[term] || (Head[term] === Times && Length[term] == 2 && 
      NumericQ[First[term]] && First[term] > 0 && 
      StringQ[Last[term]]);
  validComplexQ[expr_] := 
   expr === 0 || (Head[expr] === Plus && 
     AllTrue[List @@ expr, validTermQ]) || (Head[expr] =!= Plus && 
     validTermQ[expr]);
  validReactionQ[rxn_] := 
   Head[rxn] === Rule && validComplexQ[rxn[[1]]] && 
    validComplexQ[rxn[[2]]];
  AllTrue[rn, validReactionQ]];
(* chRN[{0 -> "S", "S" + "I" -> 2*"I", "I" -> 0}] *)

(* Format conversion functions *)
mat2Matl[matrix_List] := 
 Module[{matStr}, 
  matStr = 
   StringJoin[
    Riffle[StringJoin[Riffle[ToString /@ #, " "]] & /@ matrix, "; "]];
  StringJoin["[", matStr, "]"]];
(* mat2Matl[{{1, 2}, {3, 4}}] *)

matl2Mat[matrix_String] := 
 Module[{formattedMatrix}, 
  formattedMatrix = StringSplit[matrix, "\n"];
  formattedMatrix = StringReplace[formattedMatrix, Whitespace .. -> " "];
  formattedMatrix = StringReplace[#, " " -> ", "] & /@ formattedMatrix;
  formattedMatrix = "{" <> # <> "}" & /@ formattedMatrix;
  formattedMatrix = "{" <> StringRiffle[formattedMatrix, ",\n"] <> "}";
  ToExpression[formattedMatrix]];
(* matl2Mat["1 2\n3 4"] *)

matlr2Mat[str_String] := 
 Module[{formattedString, result}, 
  formattedString = 
   StringReplace[str, {"{" -> "", "}" -> "", "[" -> "", "]" -> ""}];
  formattedString = StringSplit[formattedString, " "];
  result = ToExpression[formattedString];
  DeleteCases[result, Null]];
(* matlr2Mat["[1 2 3 4]"] *)

(* Complex analysis functions *)
CreateMassActionRate[reactants_Association, kParam_] := 
 kParam*Product[
   reactants[spec]^coeff, {spec, Keys[reactants]}, {coeff, 
    Values[reactants]}];
(* CreateMassActionRate[<|"S" -> 1, "I" -> 1|>, k] *)

convNum[vertices_List] := 
 Module[{basis, processTerm, parseVertex}, 
  basis = Association[{"A" -> {1, 0}, "B" -> {0, 1}}];
  processTerm[term_] := 
   Module[{coef, letter}, 
    {coef, letter} = 
     StringCases[
       term, {a : DigitCharacter .. ~~ " " ~~ 
          l : ("A" | "B") :> {ToExpression[a], l}, 
        l : ("A" | "B") :> {1, l}}][[1]];
    coef*basis[letter]];
  parseVertex[vertex_String] := 
   Total[processTerm /@ StringSplit[vertex, " + "]];
  parseVertex /@ vertices];
(* convNum[{"A", "2 B", "A + B"}] *)

(* Advanced analysis functions *)
albe[RHS_, var_] := 
 Module[{rts, nvar, nrts, al, be, ga, i, j, rate, coeff}, 
  rts = rtS[RHS];
  nvar = Length[var];
  nrts = Length[rts];
  al = ConstantArray[0, {nvar, nrts}];
  be = ConstantArray[0, {nvar, nrts}];
  Do[Do[al[[i, j]] = Exponent[rts[[j]], var[[i]]], {i, nvar}], {j, 
    nrts}];
  Do[Do[rate = rts[[j]];
     coeff = Coefficient[RHS[[i]], rate];
     be[[i, j]] = al[[i, j]] + coeff, {j, nrts}], {i, nvar}];
  ga = be - al;
  {al, be, ga, rts}];
(* albe[RHS, vars] *)

RHS2RN[RHS_, var_] := 
 Module[{al, be, ga, spe, lhs, rhs, RN, rts}, 
  {al, be, ga, rts} = albe[RHS, var];
  spe = Symbol[ToUpperCase[ToString[#]]] & /@ var;
  lhs = spe . al;
  rhs = spe . be;
  lhs = lhs /. {x_ + 0 :> x, 0 + x_ :> x, 0*_ :> 0};
  rhs = rhs /. {x_ + 0 :> x, 0 + x_ :> x, 0*_ :> 0};
  RN = Thread[lhs -> rhs];
  {RN, rts, ga}];
(* RHS2RN[dynamics, variables] *)
