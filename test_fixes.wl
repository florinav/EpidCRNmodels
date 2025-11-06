(* Test script for Core.wl fixes *)

Get["Core.wl"];

Print["=== Testing extSpe and comp2Asso ==="];

(* Test 1: extSpe with Rule format *)
Print["\n Test 1: extSpe with Rule format"];
RN1 = {0 -> "S1", "S1" + "I1" -> 2*"I1", "I1" -> 0};
spe1 = extSpe[RN1];
Print["Species extracted: ", spe1];
Print["All strings? ", AllTrue[spe1, StringQ]];

(* Test 2: extSpe with List format *)
Print["\nTest 2: extSpe with List format"];
RN2 = {{0, "S1"}, {"S1" + "I1", 2*"I1"}, {"I1", 0}};
spe2 = extSpe[RN2];
Print["Species extracted: ", spe2];
Print["All strings? ", AllTrue[spe2, StringQ]];

(* Test 3: comp2Asso backward compatibility *)
Print["\nTest 3: comp2Asso backward compatibility"];
testExpr = "v1" + 2*"v2" + 3*"v3";
result1 = compToAsso[testExpr];
result2 = comp2Asso[testExpr];
Print["compToAsso result: ", result1];
Print["comp2Asso result: ", result2];
Print["Results equal? ", result1 === result2];

(* Test 4: Converting species to variables *)
Print["\nTest 4: Converting species to variables"];
RN3 = {"v1" -> "x1", "v2" -> "x2", "v3" -> "x3", "x1" + "x2" -> "z", "z" -> 0};
spe3 = extSpe[RN3];
Print["Species (strings): ", spe3];
var3 = ToExpression[spe3];
Print["Variables (symbols): ", var3];
Print["All symbols? ", AllTrue[var3, SymbolQ]];

(* Test 5: Full workflow like user's code *)
Print["\nTest 5: Full workflow"];
RN4 = {"v1" -> "v1" + "x1", "v2" -> "v2" + "x2", "v3" -> "v3" + "x3",
       "x1" + "v1" -> "v1", "x2" + "v2" -> "v2", "x3" + "v3" -> "v3",
       "v1" + "v2" + "v3" -> "v1" + "v2" + "v3" + "z", "z" -> 0};
result = extMat[RN4];
spe4 = result[[1]];
RHS = result[[6]];
Print["Species: ", spe4];
Print["RHS: ", RHS];
var4 = ToExpression[spe4];
Print["Variables: ", var4];

(* Simulate user's code *)
varV = Select[var4, StringMatchQ[ToString[#], "v" ~~ __] &];
Print["varV (v-variables): ", varV];
cInf = Thread[varV -> 0];
Print["cInf (constants): ", cInf];
eqD = Thread[(RHS /. cInf) == 0];
Print["eqD (equations at DFE): ", eqD];

(* Now solve *)
Print["\nSolving system..."];
cDFE = Solve[eqD, var4];
Print["Solution: ", cDFE];

Print["\n=== All tests completed ==="];
