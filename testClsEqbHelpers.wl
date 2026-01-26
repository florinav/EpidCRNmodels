(* ::Package:: *)

(* Test auxiliary functions for clsEqb *)

ClearAll["Global`*"];
SetDirectory["C:\\Users\\flori\\Dropbox\\EpidCRNmodels"];

(* OSN test data *)
mSi = {{x2}, {B1, S1}, {B2, S2}};
currentSiphon = {B2, S2};

(* Test solution from Siphon 4 *)
sol = {x1 -> (et*La)/((al + et)*mu), x2 -> mu/et,
       x3 -> (al*et*La - al*mu*mun - et*mu*mun)/(et*(al + et)*mu),
       B1 -> 0, S1 -> 0, R -> 0};

Print["mSi=", mSi];
Print["currentSiphon=", currentSiphon];
Print["sol=", sol];

(* Test 1: Check if other siphon is zero *)
hasOtherSiphonZero[s_, currentSip_, allSiphons_] := Module[{result},
  result = False;
  Do[
    If[sip =!= currentSip,
      Module[{vals},
        vals = sip /. s;
        Print["  Checking sip ", sip, " vs current ", currentSip, ": vals=", vals];
        If[AllTrue[vals, # === 0 &],
          Print["  --> Found other siphon zero: ", sip];
          result = True;
        ];
      ]
    ],
    {sip, allSiphons}
  ];
  result
];

Print["\nTest hasOtherSiphonZero:"];
result = hasOtherSiphonZero[sol, currentSiphon, mSi];
Print["Result: ", result];
Print["Expected: True (because {B1,S1} is zero)"];
