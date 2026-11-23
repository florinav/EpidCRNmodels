(* ::Package:: *)

Print["chemJac is loading..."]

(* mr3: Function returning the rank of matrix m.  For compatibility with earlier versions of Mathematica *)

mr3[m_] := Length[Select[RowReduce[m], ! (Dot[#, #] === 0) &]];

(* makeMonomial: This command takes the stoichiometric matrix (Smatrix) as an argument to create the flux vector v, using mass action kinetics.  This is a column vector indexed by the reactions of the chemical network (i.e., the columns of Smatrix) whose entries will be the products of the concentrations of the reactants of the corresponding reaction, with multiplicity, multiplied by a rate constant.  The chemical species are given the names a[1], a[2], ... while the rate constants are given the names k[1], k[2], ... *)
   
makeMonomial[Smatrix_] := Module[{tmp, lista, pow,
   mat}, refine[Smat_] := Smat /. Map[# -> 1 &, Variables[Smat]]; mat = 
    refine[Smatrix];
    pow = Table[0, {Length[mat[[1]]]}]; tmp = Table[
          If[mat[[i, j]] >= 0, 1, If[pow[[j]] > 0, a[
      i]^(-Smatrix[[i, j]]), pow[[j]]++; k[j] a[
    i]^(-Smatrix[[i, j]])]], {i, 1, Length[mat]}, {j, 1, Length[mat[[1]]]}];
    lista = Map[Apply[Times, #] &, Transpose[tmp]]; lista]

(* makeMonomialNMA: Command to make the flux vector v in the more general situation of non-mass-action kinetics.  The entry of the flux vector for the jth reaction will be of the form k[j][a[i_1], a[i_2], ..., a[i_n]], where a[i_1], ..., a[i_n] are the reactants.  In other words, we assume the reaction rate is an indeterminate function of the reactant concentrations. *)

makeMonomialNMA[Smatrix_] :=
  Module[{reactantlist, tmp}, tmp = Smatrix; tmp = Map[If[# < 0, 1, 0] &, 
      tmp, {2}]; tmp = Transpose[tmp]; reactantlist = Table[
      tmp[[i, j]]*a[j], {i, 1, Dimensions[tmp][[1]]}, {
        j, 1, Dimensions[tmp][[2]]}];
    reactantlist = DeleteCases[reactantlist, 0, {2}];
    Table[Apply[k[i], reactantlist[[i]]], {i, 1, Length[reactantlist]}]]

(* variables: This function takes a vector of monomials (fluxes) as input, and returns a vector of all the variables (species) appearing in the vector of monomials.  WARNING: if the stoichiometric matrix is not fully reversible, then not all species may appear as variables in the monomial vector, so this command may return too few variables. *)

variables[monom_] :=
     Module[{list},  list = Select[Variables[monom], Head[#] == Head[a[
        1]] || (Head[#] == Power && Head[Level[#, 1][[1]]] == Head[a[1]]) &]; 
                
      Table[ If[ Length[list[[
              i]]] == 1, list[[i]], ( Level[list[[i]], 
                  1])[[1]] ], {i, 1, Length[list]}]] // Union

(* svars: The command 'svars' can be used to generate the species variables whenever the flux vector has been generated using the makeMonomial or makeMonomialNMA commands.  It takes a stoichiometric matrix as input and simply returns a list {a[1], a[2], ..., a[n]} where n is the number of rows in the stoichiometric matrix (a.k.a. the number of species).  This can be preferable to using the above command since it avoids the problem of returning too few variables. *) 

svars[Smatrix_] := Table[a[i], {i, 1, Length[Smatrix]}];

(* makeReversible: The command 'makeReversible' takes a stoichiometric matrix that only has one direction entered for each reaction and makes it fully reversible (so the number of columns will be doubled) *)

makeReversible[Smat_] := Module[{tmp}, tmp = Transpose[S]; Flatten[ 
        Table[ {tmp[[i]], -tmp[[i]]}, {i, 1, Length[tmp]}], 1] // Transpose]

(* jac: This command returns the Jacobian matrix obtained from the list of functions Fns and the list of variables Vars. *)

jac[Fns_, Vars_] := Table[D[Fns[[i]], Vars[[j]]], {
    i, 1, Length[Fns]}, {j, 1, Length[Vars]}]

(* cfDet: cfDet takes a square Jacobian matrix Jac and returns the associated Craciun-Feinberg determinant, defined as Det[Jac-I] where I is the identity  matrix. *) 

cfDet[Jac_] := Det[Jac - IdentityMatrix[Length[Jac]]] // Expand;

(* coeffsList, coeffs: The functions 'coeffsList' and 'coeffs' count the number of terms in a determinant expansion, find all the numerical coefficients of all the terms in the determinant expansion, and count the number of times each coefficient appears in the expansion. *)

coeffsList[expansion_] := DeleteCases[Flatten[CoefficientList[expansion, Variables[expansion]]], 0]

coeffs[expansion_] := Module[{tmp}, tmp = coeffsList[expansion];
    Print["The number of terms in the det expansion is ", Length[tmp], ","];
    Print["and (a,b) says that the number of terms with coef a is b:"];
    Map[{#, Count[tmp, #]} &, Union[tmp]]]

(* coreDet: This command computes the core determinant of the Jacobian, which is natural for analyzing chemical networks with no outflows (see Helton-Klep-Gomez).  The inputs are the Jacobian as computed by the 'jac' command, and the stoichiometric matrix. *)  

coreDet[Jac_, Smatrix_] := Module[{d, cdet}, d = Length[Jac] - mr3[Smatrix];
    cdet = (Det[Jac - t IdentityMatrix[Length[Jac]]])/t^d // Expand;
    cdet /. t -> 0]


Print["chemJac has loaded"]
