(* ::Package:: *)

BeginPackage["EpidCRN`"];

(* Global variables used throughout the package *)
Global`ome;

(* ========================================================================= *)
(* CORE UTILITIES - Reaction parsing and species extraction *)
(* ========================================================================= *)

extMat::usage = "{spe, al, be, gamma, Rv, RHS, def} = extMat[reactions] extracts comprehensive stoichiometric and network information. Returns species list, alpha matrix (reactants), beta matrix (products), gamma matrix (net stoichiometric), reaction rate vector, RHS of mass action ODEs, and deficiency as {formula, terms, result}. Combines functionality of extSpe, stoichiometricMatrices, reaToRHS, and diagnostics.";
checkPersistence::usage = "checkPersistence[RN] determines persistence status of reaction network RN. Returns {status, analysis} where status is 'Persistent', 'Unknown', or 'Non-persistent'.";
persistenceReport::usage = "persistenceReport[RN] provides comprehensive persistence analysis of reaction network RN with detailed output.";
isCatalytic::usage = "isCatalytic[RN] checks if the reaction network contains catalytic sets.";

findCatalyticSets::usage = "findCatalyticSets[RN] identifies all catalytic sets in the reaction network.";
convertReactionFormat::usage = "convertReactionFormat[reactions] converts arrow format reactions {A->B, C->D} to pair format {{A,B}, {C,D}}. Returns input unchanged if already in pair format.";
compToAsso::usage = "compToAsso[side] parses a reaction side (left or right) and returns an association of species names to stoichiometric coefficients. Example: compToAsso[k*\"i\" + 2*\"s\"] returns <|\"i\"->k, \"s\"->2|>";
extSpe::usage = "extSpe[reactions] extracts all species names from a reaction network. Returns a list of unique species strings. Example: species = extSpe[reactions]";
asoRea::usage = "asoRea[RN] transforms classic reaction network format into association format with \"Substrates\" and \"Products\" keys. Example: associations = asoRea[reactions]";
strToSymb::usage = "strToSymb[reactions] converts string format reactions to symbolic format for internal processing.";

(* ========================================================================= *)
(* STOICHIOMETRIC ANALYSIS *)
(* ========================================================================= *)

stoichiometricMatrices::usage = "{alpha, beta, gamma, species} = stoichiometricMatrices[reactions] creates stoichiometric matrices for a reaction network. Returns {alpha, beta, gamma, species} where alpha is reactant matrix, beta is product matrix, gamma=beta-alpha is net matrix, and species is the species list";
reaToRHS::usage = "{RHS, species, Rv} = reaToRHS[reactions] generates the right-hand side of the ODE system for a reaction network using mass action kinetics. Returns {RHS, species, gamma, Rv} where RHS is the vector field, species is the species list, gamma is the net stoichiometric matrix, and Rv is the reaction rate vector";
expM::usage = "expM[var,expo] gives the vector var at power in matrix expo using Inner[OperatorApplied[Power],#2,#1,Times]&";

(* =======================================================*)
(*DRAINABILITY ANALYSIS-Based on Deshpande& Gopalkrishnan (2014)*)
(* =========================================================================*)
autocatalysisReport::usage = "autocatalysisReport[reactions] provides comprehensive analysis of autocatalytic behavior and persistence properties. Options: Verbose -> True/False. Returns detailed association with network info, persistence analysis, catalysis analysis, and theoretical insights based on Deshpande & Gopalkrishnan (2014).";

isDrainable::usage = "isDrainable[reactions, speciesSet] checks if speciesSet is drainable for the given reaction network. A set is drainable if there exists a reaction pathway that decreases all species in the set. Returns True if such a pathway exists, False otherwise.";

isSelfReplicable::usage = "isSelfReplicable[reactions, speciesSet] checks if speciesSet is self-replicable for the given reaction network. A set is self-replicable if there exists a reaction pathway that increases all species in the set. Returns True if such a pathway exists, False otherwise.";

isCritical::usage = "isCritical[reactions, speciesSet] checks if speciesSet is critical (no positive conservation law has support contained in the set). Returns True if the set is critical, False otherwise.";

siphonAnalysis::usage = "siphonAnalysis[reactions] provides comprehensive siphon classification for a reaction network. Returns an association with each siphon classified as drainable, self-replicable, critical, and categorized by type. This implements the core theoretical insights from Deshpande & Gopalkrishnan (2014).";

persistenceAnalysis::usage = "persistenceAnalysis[reactions] analyzes network persistence based on drainable siphons. Returns detailed information about persistence, drainable siphons, and extinction threats. Networks without drainable siphons are persistent.";

catalysisAnalysis::usage = "catalysisAnalysis[reactions] identifies autocatalytic behavior in reaction networks by finding self-replicable critical siphons. Returns information about catalytic sets and autocatalytic potential.";

findMinimalCriticalSiphons::usage = "findMinimalCriticalSiphons[reactions] finds all minimal critical siphons in a reaction network. These are the siphons that pose potential threats to persistence according to the theory.";
(* ========================================================================= *)

(* ENDOTACTIC ANALYSIS - Geometric and structural properties *)
(* ========================================================================= *)

endo::usage = "endo[reactions] analyzes a reaction network for endotactic properties. Options: \"ShowPlot\" -> True/False/Automatic (default: Automatic), \"Verbose\" -> True/False (default: False).";
diagnostics::usage = "diagnostics[reactions] provides detailed network analysis including deficiency.";
getComplexes::usage = "getComplexes[reactions] extracts all complexes from a list of reactions.";
isEndotactic::usage = "isEndotactic[reactions, speciesList] checks if network is endotactic.";
isStronglyEndotactic::usage = "isStronglyEndotactic[reactions, speciesList] checks if network is strongly endotactic.";

(* ========================================================================= *)
(* VISUALIZATION *)
(* ========================================================================= *)

EucFHJ::usage = "EucFHJ[reactions] draws the Euclidean Feinberg-Horn-Jackson graph for a two-species reaction network, showing complexes as red points, reactions as blue arrows with red heads, and the Newton polytope in gray. Example: EucFHJ[{{0,\"A\"},{\"A\"+\"B\",2*\"B\"}}]";
NewtPol::usage = "NewtPol[complexCoords] generates the Newton polytope (convex hull) graphics primitive from an association of complex names to 2D coordinates. Returns {LightGray, EdgeForm[{Black,Thin}], ConvexHullMesh[points]} or {} if insufficient points.";
FHJ::usage="FHJ[comp_List,edges_List,rates_List, ver_:{},groups_List:{}] generates the Feinberg-Horn-Jackson graph.";

(* ========================================================================= *)
(* MULTISCALE ANALYSIS *)
(* ========================================================================= *)

ACK::usage = "ACK[RN, continuousIndices] yields a list of four outputs: discrete network reactions, continuous network reactions, continuous species ODE, and quasi-equilibrium Poisson parameter in complex balanced case (or not Poisson), by Anderson-Cappelletti-Kurtz multiscale theory.";

(* ========================================================================= *)
(* EPIDEMIOLOGICAL AND DYNAMICAL ANALYSIS *)
(* ========================================================================= *)

(* Next Generation Matrix analysis *)
NGM::usage = "{Jx,F,V,K,Jy,Jxy,Jyx,chp,Kd} = 
NGM[mod_,inf_] yields {Jx,F,V,K,Jy,Jxy,Jyx,chp(u),Kd}; 
they are the infectious Jacobian, the new infections, transitions, 
and next generation matrices, char. pol,  other blocks of Jacobian, and the alt K"; 
NGMs::usage = "NGMs[mod_,inf_] simpler version of NGM[mod_,inf_], treats incorrectly denominators and exponents"; 
JR0::usage = "{R0,co} = JR0[pol] computes basic reproduction number and coefficients";
DFE::usage = "{diseaseFreeeEquilibrium} = DFE[mod_,inf_] yields the DFE of the model";
extP::usage ="{extinctionProb} = extP[mod_,inf_] yields the Bacaer equation for approximate extinction probability";

(* Stability analysis *)
Stab::usage = "{stabilityResult} = Stab[mod_,cfp_,cn_:{}] analyzes stability at fixed point";
ACM::usage = "A2 = ACM[A,k] yields additive compound matrix";
CofP::usage = "co = CofP[list] yields coefficients of a polynomial as required by Routh-Hurwitz theory, ie normalized so the free coefficient is 1 (see for example R\[ODoubleDot]st, Tekeli, Thm 4A)";
CofRH::usage = "co = CofRH[mat] yields coefficients of CharacteristicPolynomial, as required by Routh-Hurwitz theory, ie normalized so the free coefficient is 1 (see for example R\[ODoubleDot]st, Tekeli, Thm 4A): Drop[Reverse[CoefficientList[(-1)^(Length@A) CharacteristicPolynomial[A,x],x]],1]";
Hur2::usage = "ine = Hur2[co] yields stability conditions";
Hur3M::usage = "{co,h3,ine} = Hur3M[A] yields Hurwitz analysis where ine=Append[inec,h3>0]";
Hur4M::usage = "{co,h4,ine} = Hur4M[A] yields fourth-order Hurwitz analysis";
Hur5M::usage = "{co,h5,ine,H5} = Hur5M[jac] yields fifth-order Hurwitz analysis";
H4::usage = "H4[co] gives the 4th Hurwitz det, needed in Routh-Hurwitz theory (see for example R\[ODoubleDot]st, Tekeli, Thm 4A). H4[CofRH[M]] gives the 4th Hurwitz det of the matrix M, and could be used in Hur4M[mat]";

(* Bifurcation and phase analysis *)
Bifp::usage = "{bifurcationPlot} = Bifp[mod_,cN_,indX_,bifv_,pl0_:0,pL_:10,y0_:-1, yM_:10,cR0_:0] gives the bifurcation plot of the dynamics wrt to one parameter";
fix::usage = "{fixedPoints} = fix[mod_,cn_:{}] finds fixed points of the model with conditions cn";
phasePl2::usage = "{phasePlot} = phasePl2[mod_,plc_:{},cn_:{}] plots a 2dim phase-plot of mod, for the components not excluded in plc";

(* ========================================================================= *)
(* STRUCTURAL ANALYSIS - Conservation, siphons, invariant facets *)
(* ========================================================================= *)

cons::usage = "{conservationVectors} = cons[mat,cp_:{}] parametrizes positively the left kernel of mat, using also conditions cp; cp is not necessary if mat is numeric";
minSiph::usage = "{minimalSiphons} = minSiph[species,reactions] finds minimal siphons in a reaction network";
isSiph::usage = "isSiph[species,reactions,siphon] checks if a given set forms a siphon in the reaction network";
invFacet::usage = "invF = invFacet[reactions,maxCodim] finds invariant facets of the positive orthant for reaction networks up to specified codimension";
isInvariantFacet::usage = "isInvariantFacet[facetSet,reactions] checks if a given set of species forms an invariant facet";

(* ========================================================================= *)
(* SIMULATION AND NUMERICAL TOOLS *)
(* ========================================================================= *)

mSim::usage = "{simulationResult} = mSim[mod,cN, cInit,T,excluded] performs model simulation";
Sta::usage = "Sta[] numeric stability analysis";

(* ========================================================================= *)
(* ADVANCED ANALYSIS TOOLS *)
(* ========================================================================= *)

JTD::usage = "{jtdResult} = JTD[mod,cn_:{}] performs JTD analysis";
JTDP::usage = "{jtdpResult} = JTDP[mod,\[Zeta]_:\[Zeta],cn_:{}] performs JTDP analysis with parameter \[Zeta]";
verHir::usage ="{reductionResult} = verHir[RHS,var,intRows] checks whether a network can be reduced using the species indexed by intRows";
RUR::usage = "{ratsub,pol,ln} = RUR[mod,ind,cn_:{}] attempts to reduce the fixed point system to one with variables specified by the list ind; only singleton ind is allowed currently; outputs are ratsub,pol,and ln=pol//Length";
GBH::usage = "{gbhResult} = GBH[pol_,var_,sc_,cn_:{}] performs Gr\[ODoubleDot]bner basis analysis";
L1Planar::usage = "{l1Result} = L1Planar[fg,eq:{}] performs L1 planar analysis; eq is condition";
DerEq::usage = "{derivativeEq} = DerEq[fg,eq:{}] computes derivative equations; eq is condition";

(* ========================================================================= *)
(* UTILITY FUNCTIONS *)
(* ========================================================================= *)

red::usage = "recl = red[re,cond] erases from the output of a Reduce all the conditions in cond";
reCL::usage = "recl = reCL[re,cond] erases from the output of a Reduce all the conditions in cond";
seZF::usage = "seZF[so_] removes in a list of lists those with a 0";
onePR::usage = "onePR[cof_,cp_:{}] outputs conditions that the first and last coefs of a list have different signs";
expon::usage= "Exponent[p,Variable[p]] computes the maximum power of an expanded form p";
posM::usage="posM[matrix] keeps all syntactically positive terms";
FposEx::usage="FposEx[matrix] extracts first syntactically positive term in a nonnegative matrix";
perR::usage="perR[M_, i_, j_] = ReplacePart[M, {i -> M[[j]], j -> M[[i]]}] swaps rows i and j";
perC::usage="perC[matrix_, cycle_List] performs a cyclic permutation on the rows (or columns) of the input matrix based on the indices in the list cycle_List. The rows (or columns) specified by cycle are rearranged according to the right rotation of cycle, and the modified matrix is returned";
IaFHJ::usage = "{oU,taF} = IaFHJ[vert,edg] analyzes Feinberg-Horn-Jackson graph structure";
IkFHJ::usage = "Ik = IkFHJ[vert,edg,tk] computes Ik matrix for FHJ analysis";
SpeComInc::usage = "SpeComInc[comp,spec] creates species-complex incidence matrix";
makeLPM::usage = "makeLPM[mat_] := Table[Det@mat[[1 ;; i, 1 ;; i]], {i, 1, Length@mat}] yields the leading principal minors";
onlyP::usage ="onlyP[m_] checks whether all the coefficients of the numerator of a rational expression m are nonnegative";
CreateMassActionRate::usage = "CreateMassActionRate[reactants,kParam] creates mass action rate expression from reactant association and rate parameter";
Par::usage = "Par[dyn,var] extracts parameters from dynamics";

(* Additional utility functions *)
H6::usage = "H6[co] computes the 6th Hurwitz determinant for stability analysis";
Idx::usage = "Idx[args] indexing function for complex manipulations";
sym2Str::usage = "sym2Str[expr] converts symbolic expressions to string format";
str2Sym::usage = "str2Sym[expr] converts string expressions to symbolic format";
rul2Str::usage = "rul2Str[rules] converts replacement rules to string format";
toSum::usage = "toSum[expr] converts expression to sum format";
toProd::usage = "toProd[expr] converts expression to product format";
strEdg::usage = "strEdg[edges] processes edge structures for graph operations";
countMS::usage = "countMS[matrix] counts negative coefficients in a matrix";
mat2Matl::usage = "mat2Matl[matrix] converts Mathematica matrix to MATLAB format string";
matl2Mat::usage = "{matrix} = matl2Mat[string] converts MATLAB format string to Mathematica matrix";
matlr2Mat::usage = "{list} = matlr2Mat[string] converts MATLAB row vector string to Mathematica list";
l2L::usage = "l2L[list] converts lowercase list format to uppercase format";
m2toM::usage = "m2toM[expr] converts m2 expressions to M format";
Stodola::usage = "{stodolaResult} = Stodola[args] implements Stodola method for eigenvalue problems";
DerL::usage = "DerL[expr] computes derivatives in L format";
convNum::usage = "convNum[expr] converts numerical expressions";
Hirono::usage = "{hironoResult} = Hirono[args] implements Hirono method";
L13::usage = "{l13Coeff} = L13[args] computes L13 coefficient";
L23::usage = "{l23Coeff} = L23[args] computes L23 coefficient";
Res1F::usage = "{residueForm} = Res1F[args] computes first residue form";
Deg::usage = "Deg[poly] computes degree of polynomial";
GetVec::usage = "{vectorResult} = GetVec[A,om] extracts vectors, used in L13,L23";

Begin["`Private`"];

(* ========================================================================= *)
(* CORE UTILITIES IMPLEMENTATION - COMBINED VERSION *)
(* ========================================================================= *)

(* Enhanced string to symbolic conversion *)
strToSymb[reactions_] := Module[{convertSide},
  
  (* Convert one side of a reaction *)
  convertSide[side_String] := Module[{result, terms, converted},
    If[side == "0" || side == "", Return[0]];
    
    (* Split by + to get individual terms *)
    terms = StringSplit[side, "+"];
    
    (* Process each term *)
    converted = Map[Module[{term, cleanTerm, coeff, species},
      term = StringTrim[#];
      
      Which[
        (* Case: "3A" or "2B" - number followed by letters *)
        StringMatchQ[term, RegularExpression["^\\d+[A-Za-z0-9]+$"]],
        Module[{digits, letters},
          digits = StringTake[term, 1];
          letters = StringDrop[term, 1];
          ToExpression[digits]*letters
        ],
        
        (* Case: just letters like "A" or "B" *)
        StringMatchQ[term, RegularExpression["^[A-Za-z0-9]+$"]],
        term,
        
        (* Default: return as is *)
        True, term
      ]
    ] &, terms];
    
    (* Combine terms with Plus *)
    If[Length[converted] == 1, 
      converted[[1]], 
      Apply[Plus, converted]
    ]
  ];
  
  (* Apply to entire reaction network *)
  Map[Function[reaction,
    If[Head[reaction] === Rule,
      convertSide[First[reaction]] -> convertSide[Last[reaction]],
      {convertSide[First[reaction]], convertSide[Last[reaction]]}
    ]
  ], reactions]
];

(* Add symbToStr conversion for visualization *)
symbToStr[complex_] := Module[{},
  Which[
    complex === 0, "0",
    Head[complex] === Plus, 
      StringJoin[Riffle[Map[
        Which[
          Head[#] === Times,
          Module[{parts, strings, nonStrings},
            parts = List @@ #;
            strings = Select[parts, StringQ];
            nonStrings = Select[parts, !StringQ[#] &];
            If[Length[nonStrings] == 1 && nonStrings[[1]] == 1,
              strings[[1]],
              ToString[nonStrings[[1]]] <> "*" <> strings[[1]]
            ]
          ],
          StringQ[#], #,
          True, ToString[#]
        ] &, List @@ complex], " + "]],
    Head[complex] === Times,
      Module[{parts, strings, nonStrings},
        parts = List @@ complex;
        strings = Select[parts, StringQ];
        nonStrings = Select[parts, !StringQ[#] &];
        If[Length[nonStrings] == 1 && nonStrings[[1]] == 1,
          strings[[1]],
          ToString[nonStrings[[1]]] <> "*" <> strings[[1]]
        ]
      ],
    StringQ[complex], complex,
    True, ToString[complex]
  ]
];

(* Enhanced compToAsso function *)
compToAsso[side_] := Module[{coeffs}, 
  coeffs = <||>;
  If[side === 0, Return[coeffs]];
  
  If[Head[side] === Plus, 
   (* Multiple species - handle each term *)
   Do[
    If[Head[term] === Times, 
     (* Manually separate strings from non-strings *)
     Module[{parts, strings, nonStrings},
      parts = List @@ term;
      strings = Select[parts, StringQ];
      nonStrings = Select[parts, !StringQ[#] &];
      
      If[Length[strings] >= 1,
       (* For each string, assign coefficient as product of non-strings *)
       Do[
        coeffs[str] = If[Length[nonStrings] == 0, 1, 
                        If[Length[nonStrings] == 1, nonStrings[[1]], Times @@ nonStrings]];
        , {str, strings}];,
       (* No strings - treat whole term as species *)
       If[StringQ[term], coeffs[term] = 1]
       ]
      ], 
     (* Not Times - single species *)
     If[StringQ[term], coeffs[term] = 1];
    ], {term, List @@ side}], 
   (* Single term *)
   If[Head[side] === Times, 
    (* Manually separate strings from non-strings *)
    Module[{parts, strings, nonStrings},
     parts = List @@ side;
     strings = Select[parts, StringQ];
     nonStrings = Select[parts, !StringQ[#] &];
     
     If[Length[strings] >= 1,
      (* For each string, assign coefficient as product of non-strings *)
      Do[
       coeffs[str] = If[Length[nonStrings] == 0, 1, 
                       If[Length[nonStrings] == 1, nonStrings[[1]], Times @@ nonStrings]];
       , {str, strings}];,
      (* No strings - treat whole side as species *)
      If[StringQ[side], coeffs[side] = 1]
      ]
     ], 
    (* Single species *)
    If[StringQ[side], coeffs[side] = 1]]];
  coeffs];

(* {spe, al, be, gamma, Rv, RHS, def} =extMat[RN] *)
extMat[reactions_] := Module[{
  spe, al, be, gamma, Rv, RHS, 
  numReactions, numSpecies, reactants, products, 
  var, rv, tk, complexes, linkageClasses, deficiency,
  defFormula, defTerms, defResult, Nc, l, s},
  
  (* 1. Extract species directly from Rule format *)
  spe = Module[{allSpecies, reactants, products}, 
    allSpecies = {};
    
    Do[
     (* Extract left and right sides of the rule directly *)
     reactants = compToAsso[reactions[[i, 1]]];
     products = compToAsso[reactions[[i, 2]]];
     allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
     , {i, Length[reactions]}];
    DeleteDuplicates[allSpecies]
  ];
  
  numSpecies = Length[spe];
  numReactions = Length[reactions];
  
  (* 2. Build stoichiometric matrices *)
  al = Table[0, {numSpecies}, {numReactions}];  (* alpha - reactants *)
  be = Table[0, {numSpecies}, {numReactions}];  (* beta - products *)
  
  Do[
    (* Extract left and right sides directly *)
    reactants = compToAsso[reactions[[j, 1]]];
    products = compToAsso[reactions[[j, 2]]];
    
    Do[
      If[KeyExistsQ[reactants, spe[[i]]], 
        al[[i, j]] = reactants[spe[[i]]]];
      If[KeyExistsQ[products, spe[[i]]], 
        be[[i, j]] = products[spe[[i]]]];
    , {i, numSpecies}];
  , {j, numReactions}];
  
  (* 3. Net stoichiometric matrix *)
  gamma = be - al;
  
  (* 4. Reaction rate vector *)
  var = ToExpression[spe];
  rv = expM[var, al // Transpose];
  tk = Array[Symbol["k" <> ToString[#]] &, numReactions];
  Rv = tk*rv;
  
  (* 5. Right-hand side *)
  RHS = gamma . Rv;
  
  (* 6. Deficiency calculation *)
  complexes = Module[{complexList},
    complexList = {};
    Do[
      Module[{reactant, product},
        (* Extract left and right sides directly *)
        reactant = reactions[[i, 1]];
        product = reactions[[i, 2]];
        If[!MemberQ[complexList, reactant], AppendTo[complexList, reactant]];
        If[!MemberQ[complexList, product], AppendTo[complexList, product]];
      ], {i, Length[reactions]}
    ];
    complexList
  ];
  
  Nc = Length[complexes];
  
  (* Number of linkage classes *)
  l = Module[{adjMatrix, components},
    If[Length[complexes] == 0, Return[0]];
    
    (* Build adjacency matrix *)
    adjMatrix = ConstantArray[0, {Length[complexes], Length[complexes]}];
    Do[
      Module[{reactant, product, reactantPos, productPos},
        reactant = reactions[[i, 1]];
        product = reactions[[i, 2]];
        reactantPos = FirstPosition[complexes, reactant][[1]];
        productPos = FirstPosition[complexes, product][[1]];
        adjMatrix[[reactantPos, productPos]] = 1;
        adjMatrix[[productPos, reactantPos]] = 1; (* Make undirected *)
      ], {i, Length[reactions]}
    ];
    
    (* Find connected components *)
    components = ConnectedComponents[AdjacencyGraph[adjMatrix]];
    Length[components]
  ];
  
  (* Dimension of stoichiometric subspace *)
  s = If[gamma == {}, 0, MatrixRank[gamma]];
  
  (* Deficiency components *)
  defFormula = "\[Delta] = Nc - \[ScriptL] - s";
  defTerms = StringJoin[
    "Nc = ", ToString[Nc], " (complexes), ",
    "\[ScriptL] = ", ToString[l], " (linkage classes), ", 
    "s = ", ToString[s], " (stoich dimension)"
  ];
  defResult = StringJoin[
    "\[Delta] = ", ToString[Nc], " - ", ToString[l], " - ", ToString[s], 
    " = ", ToString[Nc - l - s]
  ];
  
  (* Return all results *)
  {spe, al, be, gamma, Rv, RHS, {defFormula, defTerms, defResult}}
];
(* Enhanced extSpe function *)
extSpe[reactions_] := Module[{allSpecies, reactants, products}, 
  allSpecies = {};
  
  Do[
   reactants = compToAsso[reactions[[i, 1]]];
   products = compToAsso[reactions[[i, 2]]];
   allSpecies = Join[allSpecies, Keys[reactants], Keys[products]];
   , {i, Length[reactions]}];
  DeleteDuplicates[allSpecies]];

(* Standard conversion function *)
convertReactionFormat[reactions_] := Module[{converted},
  If[Length[reactions] == 0, Return[{}]];
  
  If[Head[reactions[[1]]] === Rule,
    converted = Table[
      {First[reactions[[i]]], Last[reactions[[i]]]}, 
      {i, Length[reactions]}
    ];
    Return[converted],
    Return[reactions]
  ]
];

(* Other core functions *)
asoRea[RN_]:=Module[{parseSide,extractSpecies},
  extractSpecies[expr_]:=Module[{terms,species},
    If[expr===0,Return[{}]];
    terms=If[Head[expr]===Plus,List@@expr,{expr}];
    species={};
    Do[Which[
      StringQ[term]&&StringContainsQ[term," "],AppendTo[species,StringTrim[StringDrop[term,1]]],
      Head[term]===Times,AppendTo[species,Cases[term,_String][[1]]],
      StringQ[term],AppendTo[species,term]];,{term,terms}];
    DeleteDuplicates[species]];
  parseSide[expr_]:=extractSpecies[expr];
  Map[Function[r,Association["Substrates"->parseSide[r[[1]]],"Products"->parseSide[r[[2]]]]],RN]];

(* Stoichiometric matrices *)
stoichiometricMatrices[reactions_] := Module[{
  convertedReactions, species, numReactions, numSpecies, 
  alpha, beta, gamma, reactants, products},
  
  convertedReactions = EpidCRN`convertReactionFormat[reactions];
  species = extSpe[convertedReactions];
  numReactions = Length[convertedReactions];
  numSpecies = Length[species];
  
  alpha = Table[0, {numSpecies}, {numReactions}];
  beta = Table[0, {numSpecies}, {numReactions}];
  
  Do[
   reactants = compToAsso[convertedReactions[[j, 1]]];
   products = compToAsso[convertedReactions[[j, 2]]];
   
   Do[
    If[KeyExistsQ[reactants, species[[i]]], 
     alpha[[i, j]] = reactants[species[[i]]]];
    If[KeyExistsQ[products, species[[i]]], 
     beta[[i, j]] = products[species[[i]]]];
    , {i, numSpecies}];
   , {j, numReactions}];
  
  gamma = beta - alpha;
  {alpha, beta, gamma, species}
];

(* RHS generation *)
reaToRHS[reactions_] := Module[{alpha, beta, gamma, species, var, rv, tk, Rv, RHS, convertedReactions},
  convertedReactions = If[Length[reactions] > 0 && Head[reactions[[1]]] === Rule,
    Table[{First[reactions[[i]]], Last[reactions[[i]]]}, {i, Length[reactions]}],
    reactions
  ];
  
  {alpha, beta, gamma, species} = stoichiometricMatrices[convertedReactions];
  var = ToExpression[species];
  rv = expM[var, alpha // Transpose];
  tk = Array[Symbol["k" <> ToString[#]] &, alpha // Transpose // Length];
  Rv = tk*rv;
  RHS = gamma . Rv;
  {RHS, species, Rv}
];

(* Keep these core utility functions *)
convertReactionFormat[reactions_] := Module[{converted},
  If[Length[reactions] == 0, Return[{}]];
  
  If[Head[reactions[[1]]] === Rule,
    converted = Table[
      {First[reactions[[i]]], Last[reactions[[i]]]}, 
      {i, Length[reactions]}
    ];
    Return[converted],
    Return[reactions]
  ]
];

asoRea[RN_]:=Module[{parseSide,extractSpecies},
  extractSpecies[expr_]:=Module[{terms,species},
    If[expr===0,Return[{}]];
    terms=If[Head[expr]===Plus,List@@expr,{expr}];
    species={};
    Do[Which[
      StringQ[term]&&StringContainsQ[term," "],AppendTo[species,StringTrim[StringDrop[term,1]]],
      Head[term]===Times,AppendTo[species,Cases[term,_String][[1]]],
      StringQ[term],AppendTo[species,term]];,{term,terms}];
    DeleteDuplicates[species]];
  parseSide[expr_]:=extractSpecies[expr];
  Map[Function[r,Association["Substrates"->parseSide[r[[1]]],"Products"->parseSide[r[[2]]]]],RN]];

expM=Inner[OperatorApplied[Power],#2,#1,Times]&;

(* ========================================================================= *)
(* CORE DRAINABILITY AND SELF-REPLICABILITY FUNCTIONS *)
(* ========================================================================= *)

isDrainable[RN_, speciesSet_List] := Module[{
  RND, spe, gamma, speciesIndices, coeffs, constraints, solution},
  
  RND = extMat[RN];
  spe = RND[[1]];
  gamma = RND[[4]];
  
  speciesIndices = Flatten[Position[spe, #] & /@ speciesSet];
  
  If[Length[speciesIndices] != Length[speciesSet],
    Print["Warning: Some species not found in network"];
    Return[False]
  ];
  
  coeffs = Array[c, Dimensions[gamma][[2]]];
  linearComb = coeffs . Transpose[gamma];
  constraints = Join[
    Thread[coeffs >= 0],
    Thread[linearComb[[speciesIndices]] < 0],
    {Total[coeffs] > 0}
  ];
  
  solution = FindInstance[constraints, coeffs, Reals];
  Length[solution] > 0
];

(* Function to test if a species set is self-replicable *)
isSelfReplicable[RN_, speciesSet_List] := Module[{
  RND, spe, gamma, speciesIndices, coeffs, constraints, solution},
  
  RND = extMat[RN];
  spe = RND[[1]];
  gamma = RND[[4]];
  
  speciesIndices = Flatten[Position[spe, #] & /@ speciesSet];
  
  If[Length[speciesIndices] != Length[speciesSet],
    Print["Warning: Some species not found in network"];
    Return[False]
  ];
  
  coeffs = Array[c, Dimensions[gamma][[2]]];
  linearComb = coeffs . Transpose[gamma];
  constraints = Join[
    Thread[coeffs >= 0],
    Thread[linearComb[[speciesIndices]] > 0],
    {Total[coeffs] > 0}
  ];
  
  solution = FindInstance[constraints, coeffs, Reals];
  Length[solution] > 0
];

(* Function to check if a set is critical *)
isCritical[RN_, speciesSet_List] := Module[{
  RND, spe, conservationLaws, speciesIndices},
  
  RND = extMat[RN];
  spe = RND[[1]];
  speciesIndices = Flatten[Position[spe, #] & /@ speciesSet];
  conservationLaws = cons[RND[[4]], {}];
  
  If[Length[conservationLaws] == 0,
    True,
    !AnyTrue[conservationLaws,
      Function[law, 
        AllTrue[speciesIndices, law[[#]] >= 0 &] && 
        AnyTrue[speciesIndices, law[[#]] > 0 &] &&
        AllTrue[Complement[Range[Length[spe]], speciesIndices], law[[#]] == 0 &]
      ]
    ]
  ]
];

(* Function to find all minimal critical siphons *)
findMinimalCriticalSiphons[RN_] := Module[{
  RND, spe, allSiphons, criticalSiphons},
  
  RND = extMat[RN];
  spe = RND[[1]];
  allSiphons = minSiph[spe, asoRea[RN]];
  criticalSiphons = Select[allSiphons, isCritical[RN, #] &];
  criticalSiphons
];

(* Function to classify all siphons *)
classifySiphons[RN_] := Module[{
  RND, spe, allSiphons, classification},
  
  RND = extMat[RN];
  spe = RND[[1]];
  allSiphons = minSiph[spe, asoRea[RN]];
  
  classification = Association[];
  
  Do[
    Module[{isDrain, isSelfRep, isCrit},
      isDrain = isDrainable[RN, siphon];
      isSelfRep = isSelfReplicable[RN, siphon];
      isCrit = isCritical[RN, siphon];
      
      classification[siphon] = <|
        "Drainable" -> isDrain,
        "SelfReplicable" -> isSelfRep,
        "Critical" -> isCrit,
        "Type" -> Which[
          isDrain && isSelfRep, "Both drainable and self-replicable",
          isDrain && !isSelfRep, "Drainable only",
          !isDrain && isSelfRep, "Self-replicable only",
          True, "Neither drainable nor self-replicable"
        ]
      |>;
    ], {siphon, allSiphons}
  ];
  
  classification
];

(* Function to check persistence - returns status string instead of boolean *)
checkPersistence[RN_] := Module[{
  RND, spe, allSiphons, drainableSiphons, analysis, 
  persistenceStatus, competing},
  
  RND = extMat[RN];
  spe = RND[[1]];
  allSiphons = minSiph[spe, asoRea[RN]];
  drainableSiphons = Select[allSiphons, isDrainable[RN, #] &];
  competing = Select[drainableSiphons, isSelfReplicable[RN, #] &];
  
  persistenceStatus = Which[
    Length[drainableSiphons] == 0, 
    "Persistent",
    Length[competing] == Length[drainableSiphons], 
    "Unknown",
    True, 
    "Non-persistent"
  ];
  
  analysis = <|
    "AllSiphons" -> allSiphons,
    "DrainableSiphons" -> drainableSiphons,
    "CompetingSiphons" -> competing,
    "PersistenceStatus" -> persistenceStatus
  |>;
  
  {persistenceStatus, analysis}
];

(* Function to generate persistence report - verbose by default, minimal output *)
persistenceReport[RN_] := Module[{
  RND, spe, persistenceCheck, siphonClassification,
  minimalCriticalSiphons, result},
  
  RND = extMat[RN];
  spe = RND[[1]];
  
  persistenceCheck = checkPersistence[RN];
  siphonClassification = classifySiphons[RN];
  minimalCriticalSiphons = findMinimalCriticalSiphons[RN];
  
  result = <|
    "Species" -> spe,
    "NumberOfSpecies" -> Length[spe],
    "NumberOfReactions" -> Length[RN],
    "Deficiency" -> RND[[7]],
    "PersistenceStatus" -> persistenceCheck[[1]],
    "PersistenceAnalysis" -> persistenceCheck[[2]],
    "SiphonClassification" -> siphonClassification,
    "MinimalCriticalSiphons" -> minimalCriticalSiphons
  |>;
  
  Print["=== PERSISTENCE ANALYSIS ==="];
  Print["Persistence: ", persistenceCheck[[1]]];
  
  If[persistenceCheck[[1]] == "Unknown",
    Print["WARNING: Competing forces - drainable siphons are also self-replicable"];
    Print["This represents the central challenge of the persistence conjecture"];
  ];
  
  Print[""];
  Print["Siphon Analysis:"];
  Print["  Drainable siphons: ", persistenceCheck[[2]]["DrainableSiphons"]];
  Print["  Competing siphons: ", persistenceCheck[[2]]["CompetingSiphons"]];
  
  Print[""];
  Print["Siphon Classification:"];
  Do[
    Print["  ", siphon, " -> ", siphonClassification[siphon]["Type"]],
    {siphon, Keys[siphonClassification]}
  ];
  
  result
];
    
(* ========================================================================= *)
(* COMPREHENSIVE SIPHON ANALYSIS *)
(* ========================================================================= *)

siphonAnalysis[reactions_, OptionsPattern[{"Verbose" -> False}]] := Module[{
  species, allSiphons, classification, verbose},
  
  verbose = OptionValue["Verbose"];
  
  (* Get basic network information *)
  species = extSpe[reactions];
  
  If[verbose, Print["Analyzing network with ", Length[species], " species: ", species]];
  
  (* Find all minimal siphons *)
  allSiphons = minSiph[species, asoRea[reactions]];
  
  If[verbose, Print["Found ", Length[allSiphons], " siphons: ", allSiphons]];
  
  (* Classify each siphon with clean variable handling *)
  classification = Association[];
  
  Do[
    Module[{siphon, drainableResult, selfReplicableResult, criticalResult, siphonType, significance},
      siphon = allSiphons[[i]];
      
      (* Calculate properties separately to avoid scoping issues *)
      drainableResult = isDrainable[reactions, siphon];
      selfReplicableResult = isSelfReplicable[reactions, siphon];
      criticalResult = isCritical[reactions, siphon];
      
      (* Determine siphon type *)
      siphonType = Which[
        drainableResult && selfReplicableResult, "Competing (both drainable and self-replicable)",
        drainableResult && !selfReplicableResult, "Extinction risk (drainable only)",
        !drainableResult && selfReplicableResult, "Autocatalytic growth (self-replicable only)",
        True, "Neutral (neither drainable nor self-replicable)"
      ];
      
      (* Assess significance *)
      significance = Which[
        !criticalResult, "Non-critical (protected by conservation)",
        drainableResult && !selfReplicableResult, "High extinction risk",
        drainableResult && selfReplicableResult, "Competition between growth and extinction",
        !drainableResult && selfReplicableResult, "Pure autocatalytic potential",
        True, "No direct persistence threat"
      ];
      
      classification[siphon] = <|
        "IsDrainable" -> drainableResult,
        "IsSelfReplicable" -> selfReplicableResult,
        "IsCritical" -> criticalResult,
        "Type" -> siphonType,
        "Significance" -> significance,
        "IsMinimal" -> True
      |>;
      
      If[verbose,
        Print["Siphon ", siphon, ": ", siphonType];
        Print["  Drainable: ", drainableResult, ", Self-replicable: ", selfReplicableResult, ", Critical: ", criticalResult];
        Print["  Significance: ", significance];
      ];
    ]
  , {i, Length[allSiphons]}];
  
  classification
];
(* ========================================================================= *)
(* PERSISTENCE ANALYSIS *)
(* ========================================================================= *)

persistenceAnalysis[reactions_, OptionsPattern[{"Verbose" -> False}]] := Module[{
  siphonClass, drainableSiphons, criticalDrainableSiphons, 
  isPersistent, threats, analysis, verbose},
  
  verbose = OptionValue["Verbose"];
  
  If[verbose, Print["=== PERSISTENCE ANALYSIS ==="]];
  
  (* Get comprehensive siphon classification *)
  siphonClass = siphonAnalysis[reactions, "Verbose" -> False];
  
  (* Extract drainable siphons *)
  drainableSiphons = Select[Keys[siphonClass], 
    siphonClass[#]["IsDrainable"] &];
  
  (* Extract critical drainable siphons (main threats) *)
  criticalDrainableSiphons = Select[drainableSiphons,
    siphonClass[#]["IsCritical"] &];
  
  (* Determine persistence based on Deshpande & Gopalkrishnan theorem *)
  isPersistent = Length[drainableSiphons] == 0;
  
  (* Assess specific threats *)
  threats = Table[
    siphon -> siphonClass[siphon]["Significance"], 
    {siphon, criticalDrainableSiphons}
  ];
  
  analysis = <|
    "Persistent" -> isPersistent,
    "AllSiphons" -> Keys[siphonClass],
    "DrainableSiphons" -> drainableSiphons,
    "CriticalDrainableSiphons" -> criticalDrainableSiphons,
    "NumberOfDrainableSiphons" -> Length[drainableSiphons],
    "ExtinctionThreats" -> threats,
    "TheoreticalBasis" -> "Networks without drainable siphons are persistent (Deshpande & Gopalkrishnan, 2014)",
    "SiphonClassification" -> siphonClass
  |>;
  
  If[verbose,
    Print["Persistent: ", isPersistent];
    Print["Total siphons: ", Length[Keys[siphonClass]]];
    Print["Drainable siphons: ", Length[drainableSiphons]];
    Print["Critical drainable siphons: ", Length[criticalDrainableSiphons]];
    If[Length[criticalDrainableSiphons] > 0,
      Print["Main extinction threats: ", criticalDrainableSiphons];
    ];
  ];
  
  analysis
];

(* ========================================================================= *)
(* CATALYSIS ANALYSIS *)
(* ========================================================================= *)

catalysisAnalysis[reactions_, OptionsPattern[{"Verbose" -> False}]] := Module[{
  siphonClass, selfReplicableSiphons, catalyticSets, strictlyCatalyticSets,
  autocatalyticPotential, analysis, verbose},
  
  verbose = OptionValue["Verbose"];
  
  If[verbose, Print["=== CATALYSIS ANALYSIS ==="]];
  
  (* Get comprehensive siphon classification *)
  siphonClass = siphonAnalysis[reactions, "Verbose" -> False];
  
  (* Extract self-replicable siphons *)
  selfReplicableSiphons = Select[Keys[siphonClass],
    siphonClass[#]["IsSelfReplicable"] &];
  
  (* According to Deshpande & Gopalkrishnan: 
     self-replicable critical siphons correspond to catalytic sets *)
  catalyticSets = Select[selfReplicableSiphons,
    siphonClass[#]["IsCritical"] &];
  
  (* Strictly catalytic sets (self-replicable but not drainable) *)
  strictlyCatalyticSets = Select[catalyticSets,
    !siphonClass[#]["IsDrainable"] &];
  
  (* Assess overall autocatalytic potential *)
  autocatalyticPotential = Which[
    Length[catalyticSets] == 0, "No autocatalytic behavior detected",
    Length[strictlyCatalyticSets] > 0, "Strong autocatalytic potential (pure growth)",
    True, "Autocatalytic potential with competition (growth vs extinction)"
  ];
  
  analysis = <|
    "HasCatalyticSets" -> Length[catalyticSets] > 0,
    "CatalyticSets" -> catalyticSets,
    "StrictlyCatalyticSets" -> strictlyCatalyticSets,
    "SelfReplicableSiphons" -> selfReplicableSiphons,
    "AutocatalyticPotential" -> autocatalyticPotential,
    "TheoreticalBasis" -> "Self-replicable critical siphons correspond to catalytic sets (Deshpande & Gopalkrishnan, 2014)",
    "SiphonClassification" -> siphonClass
  |>;
  
  If[verbose,
    Print["Catalytic sets found: ", Length[catalyticSets]];
    Print["Strictly catalytic sets: ", Length[strictlyCatalyticSets]];
    Print["Autocatalytic potential: ", autocatalyticPotential];
    If[Length[catalyticSets] > 0,
      Print["Catalytic sets: ", catalyticSets];
    ];
  ];
  
  analysis
];

(* ========================================================================= *)
(* MINIMAL CRITICAL SIPHONS *)
(* ========================================================================= *)

findMinimalCriticalSiphons[reactions_] := Module[{
  species, allSiphons, criticalSiphons},
  
  species = extSpe[reactions];
  allSiphons = minSiph[species, asoRea[reactions]];
  
  (* Filter for critical siphons *)
  criticalSiphons = Select[allSiphons, isCritical[reactions, #] &];
  
  (* All siphons from minSiph are already minimal *)
  criticalSiphons
];

(* ========================================================================= *)
(* COMPREHENSIVE NETWORK ANALYSIS *)
(* ========================================================================= *)

autocatalysisReport[reactions_, OptionsPattern[{"Verbose" -> True}]] := Module[{
  matResults, persistenceResults, catalysisResults, 
  minimalCriticalSiphons, summary, verbose},
  
  verbose = OptionValue["Verbose"];
  
  If[verbose, Print["=== COMPREHENSIVE AUTOCATALYSIS ANALYSIS ==="]];
  
  (* Get basic network information *)
  matResults = extMat[reactions];
  
  (* Perform all analyses *)
  persistenceResults = persistenceAnalysis[reactions, "Verbose" -> False];
  catalysisResults = catalysisAnalysis[reactions, "Verbose" -> False];
  minimalCriticalSiphons = findMinimalCriticalSiphons[reactions];
  
  (* Create summary *)
  summary = <|
    "NetworkInfo" -> <|
      "Species" -> matResults[[1]],
      "NumberOfSpecies" -> Length[matResults[[1]]],
      "NumberOfReactions" -> Length[reactions],
      "Deficiency" -> matResults[[7]]
    |>,
    "Persistence" -> persistenceResults,
    "Catalysis" -> catalysisResults,
    "MinimalCriticalSiphons" -> minimalCriticalSiphons
  |>;
  
  If[verbose,
    Print["Network: ", Length[matResults[[1]]], " species, ", Length[reactions], " reactions"];
    Print["Persistent: ", persistenceResults["Persistent"]];
    Print["Has catalytic sets: ", catalysisResults["HasCatalyticSets"]];
    Print["Minimal critical siphons: ", Length[minimalCriticalSiphons]];
  ];
  
  summary
];


isCatalytic[RN_] := Module[{catalyticSets},
  catalyticSets = findCatalyticSets[RN];
  Length[catalyticSets] > 0
];


findCatalyticSets[RN_] := Module[{
  RND, spe, allSiphons, selfReplicableCriticalSiphons, 
  reactionCatalysts, allCatalyticSets},
  
  RND = extMat[RN];
  spe = RND[[1]];
  allSiphons = Quiet[minSiph[spe, asoRea[RN]]];
  
  If[allSiphons === $Failed || !ListQ[allSiphons],
    allSiphons = {}
  ];
  
  (* Find self-replicable critical siphons *)
  selfReplicableCriticalSiphons = Select[allSiphons,
    isSelfReplicable[RN, #] && isCritical[RN, #] &
  ];
  
  (* Find traditional catalysts (species unchanged in reactions) *)
  reactionCatalysts = {};
  Do[
    Module[{reactants, products, catalysts},
      reactants = compToAsso[RN[[i, 1]]];
      products = compToAsso[RN[[i, 2]]];
      catalysts = Select[Keys[reactants], 
        KeyExistsQ[products, #] && reactants[#] == products[#] &];
      If[Length[catalysts] > 0,
        reactionCatalysts = Union[reactionCatalysts, {catalysts}]
      ];
    ], {i, Length[RN]}
  ];
  
  Union[selfReplicableCriticalSiphons, reactionCatalysts]
];
(* ========================================================================= *)
(* ENHANCED ENDOTACTIC ANALYSIS IMPLEMENTATION *)
(* ========================================================================= *)

(* Convert complex to vector in species space - FIXED *)
complexToVector[complex_, speciesList_] := Module[{parsed, vector},
  parsed = compToAsso[complex];
  vector = Table[
    If[KeyExistsQ[parsed, speciesList[[i]]], 
      parsed[speciesList[[i]]], 
      0
    ], {i, Length[speciesList]}];
  vector
];

(* Compute reaction vector (product - reactant) - FIXED *)
reactionVector[reaction_, speciesList_] := Module[{reactantVec, productVec},
  reactantVec = complexToVector[reaction[[1]], speciesList];
  productVec = complexToVector[reaction[[2]], speciesList];
  productVec - reactantVec
];

(* Compute stoichiometric subspace *)
stoichiometricSubspace[reactions_, speciesList_] := Module[{reactionVectors},
  reactionVectors = Table[reactionVector[reactions[[i]], speciesList], {i, Length[reactions]}];
  reactionVectors = DeleteDuplicates[Select[reactionVectors, # != Table[0, {Length[speciesList]}] &]];
  If[Length[reactionVectors] == 0, {}, reactionVectors]
];

(* Enhanced getMaximalElements with better error handling *)
getMaximalElements[vectors_, w_] := Module[{dotProducts, maxValue, validVectors, result},
  If[Length[vectors] == 0, Return[{}]];
  
  (* Ensure vectors is a proper list of lists *)
  validVectors = Select[vectors, VectorQ[#] && Length[#] == Length[w] &];
  
  If[Length[validVectors] == 0, Return[{}]];
  
  (* Compute dot products safely *)
  dotProducts = Map[w . # &, validVectors];
  
  (* Check if we have any finite dot products *)
  If[Length[dotProducts] == 0 || !AllTrue[dotProducts, NumericQ], Return[{}]];
  
  maxValue = Max[dotProducts];
  
  (* Use Position for robust selection *)
  result = validVectors[[Flatten[Position[dotProducts, maxValue]]]];
  
  result
];

(* Other endotactic helper functions *)
getReactantVectors[reactions_, speciesList_] := Module[{reactants},
  reactants = Table[complexToVector[reactions[[i, 1]], speciesList], {i, Length[reactions]}];
  DeleteDuplicates[reactants]
];

isWEssential[reaction_, w_, speciesList_] := Module[{reacVec},
  reacVec = reactionVector[reaction, speciesList];
  w . reacVec != 0
];

getWEssentialReactions[reactions_, w_, speciesList_] := 
  Select[reactions, isWEssential[#, w, speciesList] &];

getWSupport[reactions_, w_, speciesList_] := Module[{essentialReactions, essentialReactants},
  essentialReactions = getWEssentialReactions[reactions, w, speciesList];
  If[Length[essentialReactions] == 0, Return[{}]];
  essentialReactants = Table[complexToVector[essentialReactions[[i, 1]], speciesList], {i, Length[essentialReactions]}];
  getMaximalElements[essentialReactants, w]
];

isWEndotactic[reactions_, w_, speciesList_] := Module[{essentialReactions, wSupport},
  essentialReactions = getWEssentialReactions[reactions, w, speciesList];
  wSupport = getWSupport[reactions, w, speciesList];
  
  AllTrue[essentialReactions, 
    Module[{reactantVec, reacVec, inSupport},
      reactantVec = complexToVector[#[[1]], speciesList];
      reacVec = reactionVector[#, speciesList];
      inSupport = MemberQ[wSupport, reactantVec];
      If[inSupport, w . reacVec < 0, True]
    ] &
  ]
];

(* Generate test directions *)
generateTestDirections[speciesList_] := Module[{n, directions},
  n = Length[speciesList];
  directions = {};
  
  (* Standard basis vectors and their negatives *)
  Do[
    AppendTo[directions, UnitVector[n, i]];
    AppendTo[directions, -UnitVector[n, i]];
  , {i, n}];
  
  (* Some random directions *)
  Do[
    AppendTo[directions, Normalize[RandomReal[{-1, 1}, n]]];
  , {20}];
  
  (* Some systematic combinations *)
  If[n >= 2,
    Do[
      Do[
        AppendTo[directions, Normalize[UnitVector[n, i] + UnitVector[n, j]]];
        AppendTo[directions, Normalize[UnitVector[n, i] - UnitVector[n, j]]];
      , {j, i + 1, n}];
    , {i, 1, n - 1}];
  ];
  
  DeleteDuplicates[directions, Norm[#1 - #2] < 10^(-10) &]
];

(* Main endotactic checks *)
isEndotactic[reactions_, speciesList_] := Module[{testDirections},
  testDirections = generateTestDirections[speciesList];
  AllTrue[testDirections, isWEndotactic[reactions, #, speciesList] &]
];

isStronglyEndotactic[reactions_, speciesList_] := Module[{stoichSubspace, testDirections},
  If[!isEndotactic[reactions, speciesList], Return[False]];
  
  stoichSubspace = stoichiometricSubspace[reactions, speciesList];
  testDirections = Select[generateTestDirections[speciesList], 
    !AllTrue[stoichSubspace, # . # == 0 &] &];
  
  AllTrue[testDirections,
    Module[{w, reactantVectors, maximalReactants, hasGoodReaction},
      w = #;
      reactantVectors = getReactantVectors[reactions, speciesList];
      maximalReactants = getMaximalElements[reactantVectors, w];
      
      hasGoodReaction = False;
      Do[
        Module[{maxReactant, candidateReactions},
          maxReactant = maximalReactants[[i]];
          candidateReactions = Select[reactions, complexToVector[#[[1]], speciesList] == maxReactant &];
          Do[
            Module[{reacVec},
              reacVec = reactionVector[candidateReactions[[j]], speciesList];
              If[w . reacVec < 0, hasGoodReaction = True; Break[]];
            ], {j, Length[candidateReactions]}];
          If[hasGoodReaction, Break[]];
        ], {i, Length[maximalReactants]}];
      hasGoodReaction
    ] &
  ]
];

(* Updated endo function with 1-endotactic report and direction search *)
Options[endo] = {"ShowPlot" -> Automatic, "Verbose" -> False};

endo[reactions_, OptionsPattern[]] := Module[{
  matResults, speciesList, isEndo, isStrongEndo, result, showPlot, verbose, processedReactions,
  axisDirections, axisEndotacticResults, failingSpecies, oneEndotacticReport,
  testDirections, endotacticDirections, endotacticDirectionResults},
  
  (* Use extMat to get species list *)
  matResults = extMat[reactions];
  speciesList = matResults[[1]];
  
  (* Convert Rule format to pair format directly - no other conversions *)
  processedReactions = Table[{reactions[[i, 1]], reactions[[i, 2]]}, {i, Length[reactions]}];
  
  showPlot = OptionValue["ShowPlot"];
  verbose = OptionValue["Verbose"];
  
  If[verbose,
    Print["Analyzing network with species: ", speciesList];
    Print["Reactions: ", processedReactions];
  ];
  
  (* 1-Endotactic analysis: check each axis direction *)
  axisDirections = Table[UnitVector[Length[speciesList], i], {i, Length[speciesList]}];
  axisEndotacticResults = Table[
    isWEndotactic[processedReactions, axisDirections[[i]], speciesList],
    {i, Length[speciesList]}
  ];
  
  (* Find failing species for 1-endotactic property *)
  failingSpecies = {};
  Do[
    If[!axisEndotacticResults[[i]], 
      AppendTo[failingSpecies, speciesList[[i]]]
    ], {i, Length[speciesList]}
  ];
  
  (* Create 1-endotactic report *)
  oneEndotacticReport = Table[
    speciesList[[i]] -> axisEndotacticResults[[i]], 
    {i, Length[speciesList]}
  ];
  
  (* Search for endotactic directions *)
  testDirections = Module[{n, directions},
    n = Length[speciesList];
    directions = {};
    
    (* Standard basis vectors and their negatives *)
    Do[
      AppendTo[directions, UnitVector[n, i]];
      AppendTo[directions, -UnitVector[n, i]];
    , {i, n}];
    
    (* Diagonal and anti-diagonal directions for 2D *)
    If[n == 2,
      AppendTo[directions, Normalize[{1, 1}]];   (* (1,1) direction *)
      AppendTo[directions, Normalize[{1, -1}]];  (* (1,-1) direction *)
      AppendTo[directions, Normalize[{-1, 1}]];  (* (-1,1) direction *)
      AppendTo[directions, Normalize[{-1, -1}]]; (* (-1,-1) direction *)
    ];
    
    (* Some additional systematic combinations for higher dimensions *)
    If[n >= 2,
      Do[
        Do[
          AppendTo[directions, Normalize[UnitVector[n, i] + UnitVector[n, j]]];
          AppendTo[directions, Normalize[UnitVector[n, i] - UnitVector[n, j]]];
        , {j, i + 1, n}];
      , {i, 1, n - 1}];
    ];
    
    DeleteDuplicates[directions, Norm[#1 - #2] < 10^(-10) &]
  ];
  
  (* Test which directions are endotactic *)
  endotacticDirectionResults = Table[
    {testDirections[[i]], isWEndotactic[processedReactions, testDirections[[i]], speciesList]},
    {i, Length[testDirections]}
  ];
  
  (* Extract only the endotactic directions *)
  endotacticDirections = Select[endotacticDirectionResults, #[[2]] == True &][[All, 1]];
  
  (* Overall endotactic checks *)
  isEndo = isEndotactic[processedReactions, speciesList];
  isStrongEndo = If[isEndo, isStronglyEndotactic[processedReactions, speciesList], False];
  
  result = <|
    "Species" -> speciesList,
    "Reactions" -> processedReactions,
    "Endotactic" -> isEndo,
    "StronglyEndotactic" -> isStrongEndo,
    "OneEndotacticReport" -> oneEndotacticReport,
    "FailingSpeciesDirections" -> failingSpecies,
    "EndotacticDirections" -> endotacticDirections
  |>;
  
  If[verbose,
    Print["Results:"];
    Print["  Endotactic: ", isEndo];
    Print["  Strongly Endotactic: ", isStrongEndo];
    Print["  1-Endotactic by species:"];
    Do[
      Print["    ", speciesList[[i]], ": ", axisEndotacticResults[[i]]],
      {i, Length[speciesList]}
    ];
    If[Length[failingSpecies] > 0,
      Print["  Failing species directions: ", failingSpecies],
      Print["  All species directions are endotactic"]
    ];
    Print["  Endotactic directions found: ", Length[endotacticDirections]];
    If[Length[endotacticDirections] > 0,
      Do[
        Print["    ", N[endotacticDirections[[i]], 3]],
        {i, Min[5, Length[endotacticDirections]]} (* Show first 5 *)
      ];
      If[Length[endotacticDirections] > 5,
        Print["    ... and ", Length[endotacticDirections] - 5, " more"]
      ]
    ];
  ];
  
  (* Show plot for 2D cases *)
  If[(showPlot === True || (showPlot === Automatic && Length[speciesList] == 2)) && Length[speciesList] == 2,
    If[verbose, Print[""]];
    Print[EucFHJ[reactions]];,
    If[showPlot === True && Length[speciesList] != 2,
      Print["Note: Visualization only available for 2-species networks. This network has ", Length[speciesList], " species."];
    ]
  ];
  
  result
];
(* ========================================================================= *)
(* ENHANCED VISUALIZATION *)
(* ========================================================================= *)
(* EucFHJ - FIXED to properly handle Rule format *)
EucFHJ[reactions_] := EucFHJ[reactions, {}]

EucFHJ[reactions_, opts___] := Module[{
  matResults, species, complexes, reactantComplexes, productOnlyComplexes, 
  reactionPairs, coords, vertices, edges, labels, reactionPolygon, result},
  
  (* Use extMat to get species *)
  matResults = extMat[reactions];
  species = matResults[[1]];
  
  If[Length[species] > 2, species = Take[species, 2]];
  If[Length[species] < 2, species = Join[species, Table["dummy" <> ToString[i], {i, 2 - Length[species]}]]];
  
  (* FIXED: Properly extract complexes from Rule format *)
  complexes = {};
  reactantComplexes = {};
  reactionPairs = {};
  
  Do[
    Module[{reactant, product},
      (* Handle Rule format properly *)
      If[Head[reactions[[i]]] === Rule,
        reactant = reactions[[i]] /. Rule -> List // First;
        product = reactions[[i]] /. Rule -> List // Last,
        (* Pair format *)
        reactant = reactions[[i, 1]];
        product = reactions[[i, 2]]
      ];
      
      If[!MemberQ[complexes, reactant], AppendTo[complexes, reactant]];
      If[!MemberQ[complexes, product], AppendTo[complexes, product]];
      If[!MemberQ[reactantComplexes, reactant], AppendTo[reactantComplexes, reactant]];
      
      AppendTo[reactionPairs, {reactant, product}];
    ], {i, Length[reactions]}
  ];
  
  productOnlyComplexes = Complement[complexes, reactantComplexes];
  
  (* Create coordinate mapping *)
  coords = <||>;
  Do[
    Module[{prodAssoc, coord},
      prodAssoc = compToAsso[complex];
      coord = Table[
        If[KeyExistsQ[prodAssoc, species[[j]]], 
          prodAssoc[species[[j]]], 
          0
        ], {j, 2}];
      coords[complex] = coord;
    ], {complex, complexes}
  ];
  
  (* Create reaction polygon *)
  reactionPolygon = Module[{reactantPoints, convexHull},
    reactantPoints = Values[KeyTake[coords, reactantComplexes]];
    If[Length[reactantPoints] >= 3,
      convexHull = ConvexHullMesh[reactantPoints];
      {Lighter[Blue, 0.8], EdgeForm[{Blue, Thick}], convexHull},
      {}
    ]
  ];
  
  (* Create vertices *)
  vertices = Module[{reactantVertices, productVertices},
    reactantVertices = Table[
      {Red, EdgeForm[{Red, Thick}], Disk[coords[complex], 0.03]}, 
      {complex, reactantComplexes}
    ];
    productVertices = Table[
      {PointSize[0.015], Red, Point[coords[complex]]}, 
      {complex, productOnlyComplexes}
    ];
    Join[reactantVertices, productVertices]
  ];
  
  (* Create labels and edges *)
  labels = Table[
    Text[Style[ToString[complex], 12, Black], coords[complex], {-1.5, -1.5}], 
    {complex, complexes}
  ];
  
  edges = Table[
    Module[{start, end},
      start = coords[pair[[1]]];
      end = coords[pair[[2]]];
      If[start == end,
        {Blue, Red, Arrowheads[Large], Arrow[{start, start + {0.2, 0.2}}]},
        {Blue, Red, Arrowheads[Large], Arrow[{start, end}]}
      ]
    ], {pair, reactionPairs}
  ];
  
  (* Final graphics *)
  result = Graphics[{reactionPolygon, edges, vertices, labels}, 
    PlotRange -> All, AspectRatio -> 1, Axes -> True, 
    AxesLabel -> {ToString[species[[1]]], ToString[species[[2]]]}, 
    PlotLabel -> "Euclidean Feinberg-Horn-Jackson Graph", 
    ImageSize -> 400, GridLines -> Automatic, GridLinesStyle -> LightGray, opts];
  result
];

NewtPol[complexCoords_Association, opts___] := Module[{points, convexHull, polytope},
  points = Values[complexCoords];
  If[Length[points] >= 3 && Length[Dimensions[points]] == 2,
    convexHull = ConvexHullMesh[points];
    polytope = {LightGray, EdgeForm[{Black, Thin}], convexHull},
    polytope = {}
  ];
  polytope
];

(* ========================================================================= *)
(* ALL OTHER FUNCTIONS FROM ORIGINAL PACKAGE *)
(* ========================================================================= *)

ACM[matrix_,k_]:=D[Minors[IdentityMatrix[Length@matrix]+t*matrix,k],t]/. t->0;

mat2Matl[matrix_List]:=Module[{matStr},
matStr=StringJoin[Riffle[StringJoin[Riffle[ToString/@#," "]]&/@matrix,"; "]];
StringJoin["[",matStr,"]"]];

perR[M_, i_, j_] := ReplacePart[M, {i -> M[[j]], j -> M[[i]]}];
perC[matrix_,cycle_List]:=Module[{tempMatrix=matrix},tempMatrix[[cycle]]=tempMatrix[[RotateRight[cycle]]];tempMatrix];

Par[RHS_,X_]:=Complement[Variables[RHS],X];

m2toM[a_List]:=ReplaceAll[str_String:>Total[StringSplit[str,"+"]]][Rule@@@StringSplit[First@(List@@@a),"-->"]]
SpeComInc[comp_,spec_]:=Coefficient[#,spec]&/@comp;
countMS[m_]:=m//Together//NumeratorDenominator//Map@CoefficientArrays//ReplaceAll[sa_SparseArray:>sa["NonzeroValues"]]//Flatten//Count[#, _?Negative] &;

onlyP[m_]:=m//Together//NumeratorDenominator//Map@CoefficientArrays//ReplaceAll[sa_SparseArray:>sa["NonzeroValues"]]//Flatten//AllTrue[#,NonNegative]&;
   
CofRH[A_?MatrixQ]:=Module[{x},Drop[Reverse[CoefficientList[(-1)^(Length@A) CharacteristicPolynomial[A,x],x]],1]];

CofP[co_?ListQ]:=Drop[Reverse[(-1)^(Length@co) *co],1];
red[re_,cond_:{}]:=re/. (# -> True & /@ cond);
reCL[re_] :=DeleteCases[re, _Symbol > 0 | Subscript[_, __] > 0, Infinity];

seZF[expr_] := Select[expr, FreeQ[#, 0] &];

onePR[cof_,cp_:{}]:=Append[cp,(cof[[#]]//First) (cof[[#]]//Last)<0]&/@Range[cof//Length];
makeLPM[mat_] := Table[Det@mat[[1 ;; i, 1 ;; i]], {i, 1, Length@mat}];
Hur2[co_]:=Module[{co3, ine}, co3=co[[3]];ine= {co[[1]] co3>0,co[[2]] co3>0};ine];

cons[gamma_, {}] := Module[{
  leftKernel, positiveConservationLaws, nullSpace, dims},
  
  dims = Dimensions[gamma];
  
  (* Find the left nullspace of gamma (vectors w such that w.gamma = 0) *)
  nullSpace = Quiet[NullSpace[Transpose[gamma]]];
  
  If[nullSpace === {} || nullSpace === $Failed,
    (* No conservation laws *)
    {},
    (* Filter for positive conservation laws *)
    positiveConservationLaws = Select[nullSpace,
      AllTrue[#, # >= 0 &] && AnyTrue[#, # > 0 &] &
    ];
    
    (* If no positive ones found, try to make them positive *)
    If[Length[positiveConservationLaws] == 0,
      positiveConservationLaws = Select[nullSpace,
        AllTrue[-#, # >= 0 &] && AnyTrue[-#, # > 0 &] &
      ];
      positiveConservationLaws = -# & /@ positiveConservationLaws;
    ];
    
    positiveConservationLaws
  ]
];
   
matl2Mat[matrix_String]:=Module[{formattedMatrix},
formattedMatrix=StringSplit[matrix,"\n"];
formattedMatrix=StringReplace[formattedMatrix,Whitespace..->" "];
formattedMatrix=StringReplace[#," "->", "]&/@formattedMatrix;
formattedMatrix="{"<>#<>"}"&/@formattedMatrix;
formattedMatrix="{"<>StringRiffle[formattedMatrix,",\n"]<>"}";
ToExpression[formattedMatrix]];

matlr2Mat[str_String]:=Module[{formattedString,result},
formattedString=StringReplace[str,{"{"->"","}"->"","["->"","]"->""}];
formattedString=StringSplit[formattedString," "];
result=ToExpression[formattedString];
DeleteCases[result,Null]];

expon:=Exponent[#,Variables[#]]&;

isSiph[species_List, reactions_List, siphon_List] := Module[{ns, sm, siphonSet, isSiphonQ, subIdx, prodIdx, substrates, products},
  ns = Length[species];sm = AssociationThread[species -> Range[ns]];siphonSet = siphon;isSiphonQ = True;
  Do[substrates = reaction["Substrates"];products = reaction["Products"];
    subIdx = If[substrates === {} || substrates === {""}, {}, Select[Lookup[sm, substrates, Nothing], IntegerQ]];
    prodIdx = If[products === {} || products === {""}, {}, Select[Lookup[sm, products, Nothing], IntegerQ]];
    Do[If[MemberQ[siphonSet, p],If[Length[subIdx] == 0,isSiphonQ = False; Break[],If[!AnyTrue[subIdx, MemberQ[siphonSet, #] &],isSiphonQ = False; Break[]]]]; , {p, prodIdx}];
    If[!isSiphonQ, Break[]];, {reaction, reactions}];
  isSiphonQ];

minSiph[species_List,reactions_List]:=Module[{ns,sm,specs,constraints,solutions,siphons,minimal},
ns=Length[species];sm=AssociationThread[species->Range[ns]];specs=Array[Symbol["s"<>ToString[#]]&,ns];constraints={Or@@specs};
Do[Module[{subIdx,prodIdx,substrates,products},substrates=reactions[[i]]["Substrates"];products=reactions[[i]]["Products"];
subIdx={};If[substrates=!={}&&substrates=!={""},Do[If[KeyExistsQ[sm,sub],AppendTo[subIdx,sm[sub]];];,{sub,substrates}];];
prodIdx={};If[products=!={}&&products=!={""},Do[If[KeyExistsQ[sm,prod],AppendTo[prodIdx,sm[prod]];];,{prod,products}];];
Do[Module[{constraint},If[Length[subIdx]==0,constraint=Not[specs[[p]]];,constraint=Implies[specs[[p]],If[Length[subIdx]==1,specs[[subIdx[[1]]]],Or@@specs[[subIdx]]]];];
AppendTo[constraints,constraint];];,{p,prodIdx}];];,{i,Length[reactions]}];
solutions=FindInstance[constraints,specs,Integers,50];If[solutions==={},Return[{}];];
siphons=Map[Flatten@Position[specs/. #,True]&,solutions];siphons=DeleteDuplicates[siphons];siphons=Select[siphons,Length[#]>0&];
minimal={};Do[If[Not[AnyTrue[siphons,Function[other,other=!=siphon&&SubsetQ[siphon,other]]]],AppendTo[minimal,siphon]];,{siphon,siphons}];
Map[species[[#]]&/@#&,minimal]];

isInvariantFacet[facetSet_,reactions_]:=Module[{vf,vars,facetIndices,facetRules,isInvariant,derivative,varSymbols},
{vf,vars}=Take[reaToRHS[reactions],2];varSymbols=ToExpression[vars];facetIndices=Flatten[Position[vars,#]&/@facetSet];
facetRules=Table[varSymbols[[facetIndices[[j]]]]->0,{j,Length[facetIndices]}];isInvariant=True;
Do[derivative=Simplify[vf[[facetIndices[[j]]]]/. facetRules];derivative=derivative/.Power[0,_]:>0;derivative=Simplify[derivative];
If[!(PossibleZeroQ[derivative]||TrueQ[Simplify[derivative<=0]]||TrueQ[Simplify[derivative<=0,Assumptions->And@@(#>=0&/@varSymbols)]]||MatchQ[derivative,_?NonPositive]||MatchQ[derivative,_?Negative]||derivative===0),isInvariant=False;Break[]];
,{j,Length[facetIndices]}];isInvariant];

invFacet[reactions_,maxCodim_]:=Module[{species,n,invariantFacets,subsets},species=extSpe[reactions];n=Length[species];invariantFacets={};
Do[subsets=Subsets[species,{k}];Do[candidateSet=subsets[[i]];If[isInvariantFacet[candidateSet,reactions],AppendTo[invariantFacets,candidateSet]];,{i,Length[subsets]}];,{k,1,Min[maxCodim,n]}];
invariantFacets];



ACK[RN_,continuousSpeciesIndices_]:=Module[{species,continuous,discrete,discreteRN,continuousRN,reaction,left,right,discreteLeft,discreteRight,continuousLeft,continuousRight,poissonParam},
species=DeleteDuplicates[Flatten[Cases[RN,_String,Infinity]]];continuous=species[[continuousSpeciesIndices]];discrete=Complement[species,continuous];
removeDiscreteTerms[expr_]:=Module[{terms},If[expr===0,Return[0]];If[Head[expr]===Plus,terms=List@@expr;terms=Select[terms,FreeQ[#,Alternatives@@discrete]&];
If[Length[terms]==0,Return[0]];If[Length[terms]==1,Return[terms[[1]]]];Return[Plus@@terms];];If[FreeQ[expr,Alternatives@@discrete],Return[expr],Return[0]];];
removeContinuousTerms[expr_]:=Module[{terms},If[expr===0,Return[0]];If[Head[expr]===Plus,terms=List@@expr;terms=Select[terms,FreeQ[#,Alternatives@@continuous]&];
If[Length[terms]==0,Return[0]];If[Length[terms]==1,Return[terms[[1]]]];Return[Plus@@terms];];If[FreeQ[expr,Alternatives@@continuous],Return[expr],Return[0]];];
discreteRN={};continuousRN={};Do[reaction=RN[[i]];left=reaction[[1]];right=reaction[[2]];
discreteLeft=removeContinuousTerms[left];discreteRight=removeContinuousTerms[right];
If[discreteLeft=!=discreteRight&&!(discreteLeft===0&&discreteRight===0),AppendTo[discreteRN,discreteLeft->discreteRight];];
continuousLeft=removeDiscreteTerms[left];continuousRight=removeDiscreteTerms[right];
If[continuousLeft=!=continuousRight&&!(continuousLeft===0&&continuousRight===0),AppendTo[continuousRN,continuousLeft->continuousRight];];,{i,Length[RN]}];
discreteRN=DeleteDuplicates[discreteRN];continuousRN=DeleteDuplicates[continuousRN];
<|"DiscreteNetwork"->discreteRN,"ContinuousNetwork"->continuousRN,"PoissonParameter"->"Solve detailed balance equations"|>];

Bifp[mod_,cN_,indX_,bifv_,pl0_:0,pL_:10,y0_:-1, yM_:10,cR0_:0]:=Module[{dyn, X,fp,pl,epi,plf},dyn=mod[[1]]/.cN;X=mod[[2]];
fp=Quiet[Solve[Thread[(dyn)==0],X]//N];epi={Text["\!\(\*SubscriptBox[\(c\), \(R0\)]\)",Offset[{10,10},{ cR0,0}]],{PointSize[Large],Style[Point[{ cR0,0}],Purple]}};
pl=Plot[Evaluate@(X[[indX]]/.fp),{bifv,pl0,pL}, PlotStyle->{Blue,Green,Red,Brown}];
plf=Show[{pl},Epilog->epi,PlotRange->{{pl0,pL},{y0,yM}},AxesLabel->{bifv,"Fixed points"}];{fp,plf}];

Idx[set_,n_PositiveInteger]:=Module[{seq},seq=(Table[Count[set,i],{i,n}]/.List->Sequence);seq];

FHJ[comp_List,edges_List,rates_List, ver_:{},groups_List:{}]:=Module[{colorList,shapeList,vertexColors,options,vertexShapes,defaultColor=Yellow},
colorList={Green,Red,Yellow,Purple,Orange};shapeList={"Square","Circle","ConcaveHexagon","Triangle","Hexagon","Pentagon","Star"};
vertexColors=Join[Flatten[MapIndexed[Thread[#1->colorList[[#2[[1]]]]]&,groups]],#->defaultColor&/@Complement[comp,Flatten[groups]]];
vertexShapes=Flatten[MapIndexed[Thread[#1->shapeList[[#2[[1]]]]]&,groups]];
options={VertexShapeFunction->vertexShapes,VertexStyle->vertexColors,VertexSize->ver,VertexLabels->{_ -> Placed[Automatic, Center]},EdgeStyle -> {{Black, Thick}},
PerformanceGoal->"Quality",EdgeLabels->Thread[edges->rates],EdgeLabelStyle->Directive[Black,Bold,Background->White]};LayeredGraphPlot[edges,Right,options]];
   
IaFHJ[vert_,edg_]:=Module[{gg,oU,taF},gg[a_,b_]:=Which[a===b[[1]],-1,a===b[[2]],1,True,0];oU=Outer[gg,vert,edg];
taF=TableForm[oU,TableHeadings->{vert,edg},TableAlignments->{Right,Top}];{oU,taF}];

IkFHJ[vert_,edg_,tk_]:=Module[{tri,gg,oU},tri=MapThread[Append,{edg,tk}];gg[a_,b_]:=Which[a===b[[1]],b[[3]],a===b[[2]],0,True,0];oU=Outer[gg,vert,tri]//Transpose];

convNum[vertices_List] := Module[{basis, processTerm, parseVertex},
  basis = Association[{"A" -> {1, 0}, "B" -> {0, 1}}];
  processTerm[term_] := Module[{coef, letter},{coef, letter} = StringCases[term, {a : DigitCharacter .. ~~ " " ~~ l : ("A" | "B") :> {ToExpression[a], l}, l : ("A" | "B") :> {1, l}}][[1]];coef * basis[letter]];
  parseVertex[vertex_String] := Total[processTerm /@ StringSplit[vertex, " + "]];parseVertex /@ vertices];

sym2Str=Replace[Thread[#1->#2],x_Symbol:>ToString[x],All]&;
str2Sym= #//. s_String :>ToExpression[s]&; 
varLS= (#//. s_String :>ToExpression[s])//Variables&;
rul2Str=# /. r_Rule :> ToString /@ r &;

Hur3M[A_]:=Module[{co,h3,inec,ineSys,\[Omega]},co=CoefficientList[(-1)^Length[A] CharacteristicPolynomial[A,\[Omega]],\[Omega]];
h3=co[[ 2]]* co[[ 3]]-co[[ 1]] *co[[ 4]];inec={co[[ 1]]>0,co[[ 2]]>0};ineSys=Append[inec,h3>0];{co,h3,ineSys}];

Hur4M[mat_]:=Module[{lm,ch,cot,co,H4,h4,ine},lm=mat//Length;ch=((-1)^lm * CharacteristicPolynomial[mat,\[Lambda]]//Factor);
cot=CoefficientList[ch,\[Lambda]];co=Reverse[Drop[cot,-1]];H4={{co[[1]],1,0,0},{co[[3]],co[[2]],co[[1]],1},{0,co[[4]],co[[3]],co[[2]]},{0,0,0,co[[4]]}};
h4=Det[H4];ine=Thread[co>0];{co,h4,ine}];

H4[co_]:={{co[[1]],1,0,0},{co[[3]],co[[2]],co[[1]],1},{0,co[[4]],co[[3]],co[[2]]},{0,0,0,co[[4]]}};

Hur5M[jac_]:=Module[{lm,ch,cot,co,H5,h5,ine},lm=jac//Length;ch=((-1)^lm * CharacteristicPolynomial[jac,\[Lambda]]//Factor);
cot=CoefficientList[ch,\[Lambda]];co=Reverse[Drop[cot,-1]];H5={{co[[1]],1,0,0,0},{co[[3]],co[[2]],co[[1]],1,0},
{co[[5]],co[[4]],co[[3]],co[[2]],co[[1]]},{0,0,co[[5]],co[[4]],co[[3]]},{0,0,0,0,co[[5]]}};h5=Det[H5];
ine=Append[Thread[co>0],co[[1]] co[[2]]>co[[3]]];{co,h5,ine,H5}];

H6[co_]:=Module[{hm},hm={{co[[1]],1,0,0,0,0},{co[[3]],co[[2]],co[[1]],1,0,0},{co[[5]],co[[4]],co[[3]],co[[2]],co[[1]],1},
{0,co[[6]],co[[5]],co[[4]],co[[3]],co[[2]]},{0,0,0,co[[6]],co[[5]],co[[4]]},{0,0,0,0,0,co[[6]]}}];

JTD[mod_,cn_:{}]:=Module[{dyn,X,jac,tr,det},dyn=mod[[1]];X=mod[[2]];jac=Grad[dyn,X]/.cn;tr=Tr[jac];det=Det[jac];{jac,tr,det}];

JTDP[mod_,\[Zeta]_:\[Zeta],cn_:{}]:=Module[{dyn,X,jac,tr,det,chp,cof},dyn=mod[[1]];X=mod[[2]];jac=Grad[dyn,X]/.cn;tr=Tr[jac];det=Det[jac];
chp=CharacteristicPolynomial[jac,\[Zeta]];cof=CoefficientList[chp,\[Zeta]];{jac,tr,det,cof,chp}];

Res1F[mod_,csr_,pol_,in_,cn_:{}]:=Module[{jac,det,res,chp,cof},jac=JTDP[mod][[1]]/.csr/.cn;det=Numerator[Together[Det[jac]]];res=Resultant[det,pol,in]//Factor];

DFE[mod_,inf_:{},cn_:{}]:=Module[{dyn,X},dyn=mod[[1]]/.cn;X=mod[[2]];Quiet[Solve[Thread[dyn==0]/.Thread[X[[inf]]->0],X]]];

fix[mod_,cn_:{}]:=Module[{dyn,X,fp,Xp},dyn=mod[[1]]//.cn;X=mod[[2]];fp=X/.Quiet[Solve[Thread[(dyn)==0],X]];
If[cn!={},Xp=Cases[_?(AllTrue[NonNegative]@#&)]@fp;fp=SortBy[Xp,{ #[[1]]&,#[[2]]&}]];fp];
   
phasePl2[mod_,cn_:{},plc_:{},in_:1]:=Module[{dyn,X,pl,fp,jac,jacE,Xp,Xs,sp,Gp,cP,xM,yM,r1,r2},
dyn=mod[[1]]//.cn;X=mod[[2]];pl=Complement[Range[Length[X]],plc];fp=X/.Quiet[NSolve[Thread[(dyn)==0],X]];
jac=Grad[dyn,X]; jacE=jac/.{Thread[X->fp[[1]]]};Xp=Cases[_?(AllTrue[NonNegative]@#&)]@fp;xM=Max/@Transpose[Xp];
Xs=SortBy[Xp,{ #[[1]]&,#[[2]]&}];r1={X[[pl[[1]]]],-xM[[pl[[1]]]]-.5,xM[[pl[[1]]]]+.5};r2={X[[pl[[2]]]],-xM[[pl[[2]]]]-.5,xM[[pl[[2]]]]+.5};
sp=StreamPlot[{dyn[[pl[[1]]]],dyn[[pl[[2]]]]},r1,r2,StreamStyle->Arrowheads[Medium],ColorFunction->"Rainbow",StreamPoints->Fine,
Frame->True,FrameLabel->{"x[t]","y[t]"},PlotLabel->Style["Phase portrait",Large],LabelStyle->18];
Gp=Graphics[{PointSize[0.03],{Red,Black,Cyan},Point[Xp]}];cP=ContourPlot[{dyn[[1]],dyn[[2]]},r1,r2,FrameLabel->{"x[t]","y[t]"},ContourStyle->{Blue,Red},
LabelStyle->Directive[Black,Medium]];{Xs, jacE,Show[sp,cP,Gp]}];
   
posM= Replace[#,{_?Negative->0,e_:>Replace[Expand[e],{Times[_?Negative,_]->0,t_Plus:>Replace[t,_?Negative|Times[_?Negative,_]->0,1]}]},{2}]&;
 
FposEx=With[{pos=First@SparseArray[#]["NonzeroPositions"]},SparseArray[{pos->Extract[#,pos]},Dimensions@#]]&;
 
NGMs[mod_,inf_:{}]:=Module[{dyn,X,infc,M,V,F,F1,V1,K,chp,Jy,Kd},dyn=mod[[1]];X=mod[[2]];infc=Complement[Range[Length[X]],inf];
Jy=Grad[dyn[[inf]],X[[inf]]];chp=CharacteristicPolynomial[Jy,u];V1=-Jy/.Thread[X[[infc]]->0];F1=Jy+V1/.Thread[X[[inf]]->0];
F=ReplaceAll[F1, _. _?Negative -> 0];V=F-Jy;K=(F . Inverse[V])/.Thread[X[[inf]]->0]//FullSimplify;Kd=( Inverse[V] . F)/.Thread[X[[inf]]->0]//FullSimplify;
{Jy,V1,F1,F,V,K,Kd,chp}];
 
NGM[mod_,inf_:{}]:=Module[{dyn,X,infc,M,V,F,F1,V1,K,chp,Jy,Jx,Jxy,Jyx,Kd},dyn=mod[[1]];X=mod[[2]];infc=Complement[Range[Length[X]],inf];
Jx=Grad[dyn[[inf]],X[[inf]]];Jy=Grad[dyn[[infc]],X[[infc]]];Jxy=Grad[dyn[[inf]],X[[infc]]];Jyx=Grad[dyn[[infc]],X[[inf]]];
chp=CharacteristicPolynomial[Jx,#]&;V1=-Jx/.Thread[X[[infc]]->0];F1=Jx+V1/.Thread[X[[inf]]->0];F=posM[F1];V=F-Jx;
K=(F . Inverse[V])/.Thread[X[[inf]]->0]//FullSimplify;
Kd=( Inverse[V] . F)/.Thread[X[[inf]]->0]//FullSimplify;
{Jx,F,V,K,Jy,Jxy,Jyx,chp,Kd}];

JR0[pol_,u_]:=Module[{co,co1,cop,con,R0J},co=CoefficientList[pol,u];Print["the  factor  has degree ",Length[co]-1];
Print["its leading  coefficient  is ",co[[Length[co]]]];co1=Expand[co[[1]] ];Print["its  constant coefficient  is ",co1];
cop=Replace[co1, _. _?Negative -> 0, {1}];con=cop-co1;Print["R0J is"];R0J=con/cop//FullSimplify;{R0J,co}];

Hirono[S_, intRows_, intCols_] :=Module[{S11, S12, S21, S22, S11plus, Sred},
  S11 = S[[intRows, intCols]];S12 = S[[intRows, Complement[Range[Dimensions[S][[2]]], intCols]]];
  S21 = S[[Complement[Range[Dimensions[S][[1]]], intRows], intCols]];S22 = S[[Complement[Range[Dimensions[S][[1]]], intRows], Complement[Range[Dimensions[S][[2]]], intCols]]];
  S11plus = PseudoInverse[S11];Sred = Simplify[S22 - S21 . S11plus . S12];Sred];

verHir[RHS_,var_,intRows_]:=Module[{extRows,sub,rhsRed,fpRed},extRows=Complement[Range[Length[var]],intRows];
sub=Solve[RHS[[intRows]]==0,var[[intRows]]][[1]];rhsRed=RHS[[extRows]]/. sub//Simplify;fpRed=Solve[rhsRed==0,var[[extRows]]][[1]];{sub,rhsRed,fpRed}];

extP[mod_,inf_]:=Module[{X,Xi,qv,ov,ngm,fv,eq},X=mod[[2]];Xi=X[[inf]];qv=Array[q,Length[Xi]];ov=Table[1,{j,Length[Xi]}];
ngm=NGM[mod,inf];F=ngm[[4]];V=ngm[[5]];fv=ov . F;eq=(qv . F)*qv-qv*fv+(ov-qv) . V];

RUR[mod_, ind_:{1}, cn_ : {}] := Module[{RHS, var, par, elim,ratsub,pol,rat1},RHS = mod[[1]]/.cn; var = mod[[2]]; par = mod[[3]]; 
elim = Complement[Range[Length[var]], ind];ratsub = seZF[Solve[Delete[Thread[RHS == 0], ind], var[[elim]]]];pol =Numerator[Together[RHS//.ratsub]];
RHS[[ind]]/.ratsub;rat1=Append[(ratsub/.var[[ind]]->1),var[[ind]]->1];{ratsub, pol,rat1}];
      
GBH[pol_,var_,sc_,cn_:{}]:=Module[{li,pa},li={pol,sc};pa=Complement[Variables[li],{var}];
GroebnerBasis[{Numerator[Together[pol]],Numerator[Together[sc]]}/.cn,pa,{var},MonomialOrder->EliminationOrder]];

Stodola[pol_,var_] :=Equal@@Sign[CoefficientList[pol,var]];
mSim[mod_,cN_, cInit_,T_:100,exc_:{}]:=Module[{dyn, X,vart,diff,diffN,initcond,eqN,ndesoln,ind}, 
dyn=mod[[1]];X=mod[[2]];vart=Through[X[t]];diff= D[vart,t] \[Minus] (dyn/.Thread[X->vart]);diffN=diff//.cN;
initcond = (vart/.t->0)- cInit;eqN=Thread[Flatten[{diffN, initcond}] == 0];ndesoln = Chop[NDSolveValue[eqN,vart,{t, 0, T}]];
ind=Complement[Range[Length[X]],exc];pl=Plot[ndesoln[[ind]],{t,0,T},AxesLabel->{"t"," "}];pl];

Stab[mod_,cfp_,cn_:{}]:=Module[{dyn,X,par,jac,jacfp,eig},dyn=mod[[1]];X=mod[[2]];par=mod[[3]];jac=Grad[dyn,X];jacfp=jac//.cfp;eig=Eigenvalues[jacfp/.cn]]

Sta[jac_,X_,Xv_]:=Map[Max[Re[Eigenvalues[jac/.Thread[X->#]]]]&,Xv];

L1Planar[fg_,var_,equilcon_:{}]:=Module[{J,xyshift,Tm,Tinvuv,FG,derivatives,a,b,i,j,L1,x,y,F,G,normForm},
{x,y}=var;J=Simplify[D[fg,{var}]/. equilcon];xyshift={x->x+(x/. equilcon),y->y+(y/. equilcon)};
Tm={{1,0},{-a/Global`ome,-b/Global`ome}};Tinvuv=Inverse[Tm] . {Global`u,Global`v};
FG=(Tm . fg/. xyshift)/. {x->Tinvuv[[1]],y->Tinvuv[[2]]}/. {a->J[[1,1]],b->J[[1,2]]};derivatives={};
For[i=0,i<=3,i++,For[j=0,j<=3-i,j++,derivatives=Join[derivatives,{Subscript[F,i,j]->(D[FG[[1]],{Global`u,i},{Global`v,j}]/. {Global`u->0,Global`v->0}),
Subscript[G,i,j]->(D[FG[[2]],{Global`u,i},{Global`v,j}]/. {Global`u->0,Global`v->0})}]]];
L1=Subscript[F,3,0]+Subscript[F,1,2]+Subscript[G,0,3]+Subscript[G,2,1]+1/Global`ome*(Subscript[F,1,1]*(Subscript[F,2,0]+Subscript[F,0,2])-Subscript[G,1,1]*(Subscript[G,2,0]+Subscript[G,0,2])+Subscript[F,0,2]*Subscript[G,0,2]-Subscript[F,2,0]*Subscript[G,2,0]);
L1=L1/. derivatives;F1=Sum[Subscript[F,i,j] Global`u^i Global`v^j,{i,0,3},{j,0,3-i}];G1=Sum[Subscript[G,i,j] Global`u^i Global`v^j,{i,0,3},{j,0,3-i}];
normForm={F1,G1}/. derivatives;{L1,FG,normForm}];

DerEq[f_,var_,equilcon_]:=Module[{derivatives,order,deriv,i,j,k,A,B,CC,DD,EE,x,y,z,par,cp},n=Length[f];A=D[f,{var}]/.equilcon;
{x,y,z}=var;par=Par[f,var];cp=Thread[par>0];derivatives={};order=2;
For[i=0,i<=order,i++,For[j=0,j<=order-i,j++,For[k=0,k<=order-i-j,k++,deriv=Simplify[D[f,{x,i},{y,j},{z,k}]/.equilcon];
derivatives=Join[derivatives,Simplify[{Subscript[F, i,j,k]->deriv[[1]],Subscript[G, i,j,k]->deriv[[2]],Subscript[H, i,j,k]->deriv[[3]]}]];]]];
B[x_,y_]:=Sum[{Subscript[F, Idx[{k,l},n]],Subscript[G, Idx[{k,l},n]],Subscript[H, Idx[{k,l},n]]}x[[k]]y[[l]]/.derivatives,{k,n},{l,n}];
CC[x_,y_,z_]:={0,0,0};DD[x_,y_,z_,s_]:={0,0,0};EE[x_,y_,z_,s_,t_]:={0,0,0};{A,B,CC}];

GetVec[A_,om_]:=Module[{n,mtx,pconj,q,qconj,normalize},n=Length[A];mtx=A-om I IdentityMatrix[n];q=NullSpace[mtx[[Range[1,n-1]]]][[1]];
mtx=A\[Transpose]-om I IdentityMatrix[n];pconj=NullSpace[mtx[[Range[1,n-1]]]][[1]];normalize=FullSimplify[pconj . q];pconj=pconj/normalize;
qconj=FullSimplify[ComplexExpand[q\[Conjugate]]];{pconj,q,qconj}];

L13[A_,B_,CC_,cp_]:=Module[{n,pconj,q,qconj,v1,v2,v3,c1,numer,denom,a,b,c,d,L1\[Kappa]\[Omega]},n=Length[A];
{pconj,q,qconj}= GetVec[A,ome];v1=CC[q,q,qconj];v2=B[q,Inverse[-A] . B[q,qconj]];v3=B[qconj,Inverse[2I ome IdentityMatrix[n]-A] . B[q,q]];
c1=pconj . (1/2 v1+v2+1/2 v3);numer=Numerator[c1];denom=Denominator[c1];a=Simplify[ComplexExpand[Re[numer]],cp];
b=Simplify[ComplexExpand[Im[numer]],cp];c=Simplify[ComplexExpand[Re[denom]],cp];d=Simplify[ComplexExpand[Im[denom]],cp];
L1\[Kappa]\[Omega]=Simplify[(a c+b d)/(c^2+d^2)]];

L23[A_,B_,CC_:{0,0,0},DD_:{0,0,0},EE_:{0,0,0}]:=Module[{n,Id,omega,invA,inv2,inv3,pconj,q,qconj,h,prec,c,invbig},
n=Length[A];Id=IdentityMatrix[n];omega=Sqrt[Det[A]/Tr[A]];invA=Inverse[A];inv2=Simplify[Inverse[2 omega I Id-A]];
inv3=Simplify[Inverse[3 omega I Id-A]];{pconj,q,qconj}= GetVec[A,omega];q=FullSimplify[q/. {\[Omega]->omega}];pconj=FullSimplify[pconj/. {\[Omega]->omega}];
Subscript[h, 2,0]=FullSimplify[inv2 . B[q,q]];Subscript[h, 1,1]=FullSimplify[-invA . B[q,q\[Conjugate]]];
prec=FullSimplify[CC[q,q,q\[Conjugate]]+2 B[q,Subscript[h, 1,1]]+B[q\[Conjugate],Subscript[h, 2,0]]];Subscript[c, 1]=FullSimplify[1/2 (pconj . prec)];
invbig=FullSimplify[Inverse[Join[Join[omega I Id-A,{q}\[Transpose],2],{Join[pconj,{0}]}]]];
Subscript[h, 2,1]=FullSimplify[invbig . Join[FullSimplify[prec-2 Subscript[c, 1] q],{0}]][[1;;n]];
Subscript[h, 3,0]=FullSimplify[inv3 . (CC[q,q,q]+3 B[q,Subscript[h, 2,0]])];
Subscript[h, 3,1]=FullSimplify[inv2 . (DD[q,q,q,q\[Conjugate]]+3 CC[q,q,Subscript[h, 1,1]]+3 CC[q,q\[Conjugate],Subscript[h, 2,0]]+3 B[Subscript[h, 2,0],Subscript[h, 1,1]]+B[q\[Conjugate],Subscript[h, 3,0]]+3 B[q,Subscript[h, 2,1]]-6 Subscript[c, 1] Subscript[h, 2,0])];
Subscript[h, 2,2]=FullSimplify[-invA . (DD[q,q,q\[Conjugate],q\[Conjugate]]+4 CC[q,q\[Conjugate],Subscript[h, 1,1]]+CC[q\[Conjugate],q\[Conjugate],Subscript[h, 2,0]]+CC[q,q,Subscript[h, 2,0]\[Conjugate]]+2 B[Subscript[h, 1,1],Subscript[h, 1,1]]+2 B[q,Subscript[h, 2,1]\[Conjugate]]+2 B[q\[Conjugate],Subscript[h, 2,1]]+B[Subscript[h, 2,0]\[Conjugate],Subscript[h, 2,0]])];
Subscript[c, 2]=FullSimplify[1/12 (pconj . (EE[q,q,q,q\[Conjugate],q\[Conjugate]]+DD[q,q,q,Subscript[h, 2,0]\[Conjugate]]+3 DD[q,q\[Conjugate],q\[Conjugate],Subscript[h, 2,0]]+6 DD[q,q,q\[Conjugate],Subscript[h, 1,1]]+CC[q\[Conjugate],q\[Conjugate],Subscript[h, 3,0]]+3 CC[q,q,Subscript[h, 2,1]\[Conjugate]]+6 CC[q,q\[Conjugate],Subscript[h, 2,1]]+3 CC[q,Subscript[h, 2,0]\[Conjugate],Subscript[h, 2,0]]+6 CC[q,Subscript[h, 1,1],Subscript[h, 1,1]]+6 CC[q\[Conjugate],Subscript[h, 2,0],Subscript[h, 1,1]]+2 B[q\[Conjugate],Subscript[h, 3,1]]+3 B[q,Subscript[h, 2,2]]+B[Subscript[h, 2,0]\[Conjugate],Subscript[h, 3,0]]+3 B[Subscript[h, 2,1]\[Conjugate],Subscript[h, 2,0]]+6 B[Subscript[h, 1,1],Subscript[h, 2,1]]))];
ComplexExpand[Re[Subscript[c, 2]]]];

toSum= (# /. Times -> Plus) &;
toProd=(# /.  Plus->Times ) &;

CreateMassActionRate[reactants_Association, kParam_] := kParam * Product[reactants[spec]^coeff, {spec, Keys[reactants]}, {coeff, Values[reactants]}];

strEdg[edges_List] := Map[ToString, edges, {2}];

Deg[poly_, var_] := Exponent[poly, var];

End[];
EndPackage[];
$ContextPath=DeleteDuplicates[Append[$ContextPath,"model`Private`"]];
