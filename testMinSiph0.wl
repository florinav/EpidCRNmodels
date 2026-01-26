(* Test minimal minSiph0 *)
<<EpidCRN`;

(* Define simple version *)
minSiph0Test[vars_, reactions_, RHS_] := {{{1}}, {2}, {3}, {4}};

(* Test it *)
result = minSiph0Test[{s, i}, {}, {}];
Print["result=", result];

(* Now check if minSiph0 from ClaudeFuckups is defined *)
Print["minSiph0 DownValues=", DownValues[minSiph0]];
Print["minSiph0 defined?=", ValueQ[minSiph0[{s,i},{},{}]]];
