(* ::Package:: *)

Print["chemYAk is loading..."]

(* mr3: Function returning the rank of matrix m.
For compatibility with earlier versions of Mathematica *)

mr3[m_] := Length[Select[RowReduce[m], ! (Dot[#, #] === 0) &]];

(* uUnion performs an unsorted union *)

uUnion[x_] := Reap[Sow[1, x], _, #1&][[2]]


outOf1Make2[vec_] := {Map[ Max[0, #] &, vec], - Map[ Min[0, #] &, vec]}
  
(* matrixY takes the Stoichiometric matrix as intput and yields a 
list of three matrices. 
The first is matrix Y, the second is a list of the input complexes for all 
the reactions (with repetitions), and the third is a list of the output 
complexes for all the reactions (with repetitions). *)

matrixY[Smat_] := 
 Module[{transposeS = Transpose[Smat], InputMatrix, OutputMatrix, ActualY},
  InputMatrix = {};
  OutputMatrix = {};
  YY = {};
  YY = Flatten[Map[outOf1Make2[#] &, transposeS], 1];
  InputMatrix = Map[outOf1Make2[#][[2]] &, transposeS];
  OutputMatrix = Map[outOf1Make2[#][[1]] &, transposeS];
  ActualY = Transpose[uUnion[YY]];
  {ActualY, InputMatrix, OutputMatrix}]
  
(* matrixG takes in matrix Y, the input complexes matrix, 
and the output complexes matrix, and yields matrix G. *)

matrixG[{YYmat_, Inputmat_, Outputmat_}] := 
 Module[{ActualY, InputMatrix, OutputMatrix, ActualK, autvec},
  ActualY = YYmat;
  InputMatrix = Inputmat;
  OutputMatrix = Outputmat;
  ActualK = {};
  For[i = 1, i <= Length[InputMatrix], i++,
   (* We will have Y.autvec = i'th column of S. InputMatrix[[
   i]] is the input complex of the i'th reaction. 
   If the j'th column of matrix Y is InputMatrix[[i]], 
   then the j'th entry of autvec will be -1.
   *)
   autvec = 
    Table[If[j == 
       Flatten[Position[Transpose[ActualY], InputMatrix[[i]]]][[1]], -1,
      (* OutputMatrix[[i]] is the output complex of the i'th reaction. 
      If the j'th column of matrix Y is OutputMatrix[[i]], 
      then the j'th entry of autvec will be 1. Otherwise, it will be 0.
      *)
      If[j == Flatten[Position[Transpose[ActualY], OutputMatrix[[i]]]][[1]], 
       1, 0]],
     {j, Length[ActualY[[1]]]}];
   (* ActualK will be the matrix consisting of all the autvec's generated at \
each iteration in our For loop.
   *)
   ActualK = Append[ActualK, autvec];
   ];
  Transpose[ActualK]]
  
makeMonomial[Ymat_, symb_] := 
  Module[{tmp, lista, pow, mat}, 
   mat = Refine[Sign[Ymat], Map[# > 0 &, Variables[-Ymat]]];
   pow = Table[0, {Length[mat[[1]]]}]; 
   tmp = Table[
     If[mat[[i, j]] >= 0, 1, 
      If[pow[[j]] > 0, symb[i] (-Ymat[[i, j]]), pow[[j]]++;  
       symb[i]^(-Ymat[[i, j]])]], {i, 1, Length[mat]}, {j, 1, 
      Length[mat[[1]]]}];
   lista = Map[Apply[Times, #] &, Transpose[tmp]]; lista];
   
(* matrixK takes matrix G as input and yields matrix K. *)

matrixK[GGmat_] := Module[{GG, ActualKK, tvec, position, n, i}, Clear[k];
  	GG = {};
  	ActualKK = Transpose[GGmat];
  	n = 1;
  (* K is a #(reactions)x #(complexes) matrix. For each complex... 
  *)
  	For[i = 1, i <= Length[Transpose[ActualKK]], i++,
   (* If complex i participates in reaction j, then the j'th entry of tvec 
   is the rate constant k[i,j]. It's 0 otherwise.
   *)
       tvec = Table[If[Transpose[ActualKK][[i]][[j]] == -1, position = 
       	Flatten[Position[ActualKK[[j]], 1]][[1]]; 
      		k[i, position], 0], {j, 1, Length[ActualKK]}];
       GG = Append[GG, tvec];
       ];
   Transpose[GG]
  ]
  
(* The desired decomposition Y Ak Psi(x). *)

decomposition[Smatt_] := 
 Module[ {Ymatt, Gmatt, Kmatt, Akmatt, psimatt}, Ymatt = matrixY[Smatt]; Gmatt = matrixG[Ymatt]; 
  Kmatt = matrixK[Gmatt]; 
  Akmatt = Gmatt.Kmatt; psimatt = makeMonomial[-Ymatt[[1]], x]; {Ymatt[[1]], Akmatt, psimatt}]

(* linkageClasses takes the Laplacian of a graph (e.g. the A_k matrix \
in the Y A_k Psi(x)) as an input. It returns a vector with \
two components. The first component gives a list of all the linkage \
classes. For example, {1,3,4} would be a linkage class where the \
first, third and fourth complex participate in reactions with each \
other. By definition, the linkage classes are disjoint. 

The second component in our output is a list of 3-tuples for each \
different complex in our chemical network. A 3-tuple will list all \
the reactions a certain complex particpates in. The first vector in \
the 3-tuple gives the index i of a particular complex y_i. The second \
vector in the 3-tuple gives indices for of each complex that y_i \
particates in a reaction with, and where y_i is an input for that \
reaction.
The third vector in the 3-tuple gives indices for of each complex \
that y_i particates in a reaction with, and where y_i is an output \
for that reaction. For example, the 3-tuple {{5},{1,2},{2,4,6}} means \
we have the reactions y_5->y_1, y_5->y_2, y_2->y_5, y_4->y_5, and \
y_6->y_2. *)

linkageClasses[Akkmat_] := 
 Module[{diag, A, newA, n, compgraph, classes, class, testvec, 
   componentIn, componentOut, i, j},
  A = Akkmat;
  diag = Tr[A, List];
  newA = A - DiagonalMatrix[diag];
  n = Length[A[[1]]];
  compgraph = {};
  classes = {};
  For[i = 1, i <= n, i++,
   testvec = Table[If[j == i, 1, 0], {j, n}];
   componentIn = 
    Flatten[Position[Map[ToString, newA.testvec], _?(# != "0" &)]];
   componentOut = 
    Flatten[Position[Map[ToString, testvec.newA], _?(# != "0" &)]];
   compgraph = Append[compgraph, {{i}, componentIn, componentOut}];
   class = Union[componentIn, componentOut, {i}];
   classes = Append[classes, class];
   ];
  For[i = 1, i <= Length[classes], i++,
   For[j = i + 1, j <= Length[classes], j++,
     If[Intersection[classes[[i]], classes[[j]]] != {}, 
       classes[[i]] = Union[classes[[i]], classes[[j]]]; 
       classes = Delete[classes, j]; j = i];
     ];];
  {classes, compgraph}
  ]
  
(* The next formula gives the topological deficiency of a chemical \
reaction network. Its input is the stoichiometric matrix. *)

topDeficiency[SSmat_] := 
 Module[ {ymat, akmat, psimatt}, {ymat, akmat, psimatt} = decomposition[SSmat];
  Length[ akmat ] - Length[  linkageClasses[ akmat][[1]] ] - mr3[SSmat]]
  
(* The next formula gives the matrix deficiency of a chemical \
reaction network. Its input is the stoichiometric matrix. Two versions are given.
A probabilistic version that runs on moderately big examples and a symbolic one. *)

matDeficiency[Smmat_] := 
 Module[{tmpak, ymat, akmat, psimat}, {ymat, akmat, psimat} = decomposition[Smmat];
  Max[Table[ tmpak=akmat/.Map[# -> RandomInteger[{1, 3}] &, Variables[akmat]];
  Length[NullSpace[ymat.tmpak]] - Length[NullSpace[tmpak]],{5}]] ]

matDeficiencySymbolic[Smaat_] := 
 Module[ {yymat, akkmat, ppsimat}, {yymat, akkmat, ppsimat} = decomposition[Smaat];
  Length[ NullSpace[ yymat.akkmat ]] - Length[ NullSpace[akkmat]] ]

Print["chemYAk has loaded"]
