(* ::Package:: *)

(* EpidCRN1 Main Package - Loader for all subpackages *)
BeginPackage["EpidCRN`"];
(* Global variables *)
Global`ome;

(* Core1 functions *)
extMat::usage = "(Core){spe, al, be, gamma, Rv, RHS, def} = extMat[reactions] Returns species list, alpha matrix (reactants), beta matrix (products), gamma matrix (net stoichiometric), reaction rate vector, RHS of mass action ODEs, and deficiency as {formula, terms, result}.";
compToAsso::usage = "(Core1)compToAsso[side] parses a reaction side (left or right) and returns an association of species names to stoichiometric coefficients. Example: compToAsso[k*\"i\" + 2*\"s\"] returns <|\"i\"->k, \"s\"->2|>";
extSpe::usage = "(Core1)extSpe[reactions] extracts all species names from a reaction network. Returns a list of unique species strings.";
asoRea::usage = "(Core1)asoRea[RN] transforms classic reaction network format into association format with \"Substrates\" and \"Products\" keys.";
arrow2pairReac::usage = "(Core1)arrow2pairReac[reactions] converts arrow format reactions {A->B, C->D} to pair format {{A,B}, {C,D}}. Returns input unchanged if already in pair format.";
strToSymb::usage = "(Core1)strToSymb[reactions] converts string format reactions to symbolic format for internal processing.";
symbToStr::usage = "(Core1)symbToStr[complex] converts symbolic expressions to string format for visualization and internal processing.";
convertReactionFormat::usage = "(Core1)convertReactionFormat[reactions] converts between different reaction network formats";
str2Sym::usage = "(Core1)str2Sym[string] converts string expressions to symbolic format";
sym2Str::usage = "(Core1)sym2Str[symbolic] converts symbolic expressions to string format";
stoichiometricMatrices::usage = "(Core1){alpha, beta, gamma, species} = stoichiometricMatrices[reactions] creates stoichiometric matrices for a reaction network.";
reaToRHS::usage = "(Core1){RHS, species, Rv} = reaToRHS[reactions] generates the right-hand side of the ODE system for a reaction network using mass action kinetics.";
expM::usage = "(Core1)expM[var,expo] gives the vector var at power in matrix expo using
Inner[OperatorApplied[Power],#2,#1,Times]&";
Par::usage="(Core1)Par[RHS,var] extracts parameters from dynamics";
allT::usage="(Core1)allT[expr] extracts all terms from expanded polynomial as list";
isN::usage="(Core1)isN[term] returns True if term has negative sign, False otherwise";
ODE2RN::usage="(Core){RN,rts,spe,alp,bet,gam}=ODE2RN[RHS,var] converts ODE system to reaction network. RHS is list of expressions, var is list of variables. Returns reaction network RN, rates vector rts, species list spe, reactant matrix alp (alpha), product matrix bet (beta), and net stoichiometric matrix gam (gamma). Algorithm: identifies monomials, uses negative terms to define rates, constructs stoichiometric matrices where sources use alpha=0, then beta=gamma+alpha. Inverse operation to extMat";
ODE2WY::usage="(Core){W,Y}=ODE2WY[RHS,var] converts ODE system to W,Y matrix format for wr0RealODE. W is coefficient matrix, Y is exponent matrix. Satisfies RHS==W.x^Y.";
ab2RN::usage="(Core)RN=ab2RN[alp,bet,spe] converts alpha/beta matrices to reaction network format. alp is reactant matrix, bet is product matrix, spe is species list.";

(* CRNT functions *)
IaFHJ::usage="(CRNT)IaFHJ[vert,edg] analyzes Feinberg-Horn-Jackson graph structure, returns {incidenceMatrix, tableForm}";
IkFHJ::usage="(CRNT)IkFHJ[vert,edg,tk] computes Ik matrix for FHJ analysis";
SpeComInc::usage="(CRNT)SpeComInc[comp,spec] creates species-complex incidence matrix";
getComE::usage="(CRNT)getComE[RN] extracts complexes and edges from reaction network, returns {complexes, edges}";
lapK::usage="(CRNT)lapK[RN, rates] computes Laplacian matrix of reaction network";

(* Boundary1 analysis functions *)
NGM::usage = "NGM[mod, infVars, F] computes Next Generation Matrix. mod={RHS,var,par}. If F provided, validates Inverse[V] all positive under par>0. If not or invalid, uses posM. Returns {Jx,F,V,K,Jy,Jxy,Jyx,Kd}.";

getInfectionIndices::usage = "getInfectionIndices[variables, siphonExpressions] extracts infection compartment indices from minimal siphons. variables is the list of system variables, siphonExpressions is the list of minimal siphon expressions. Returns a list of integer indices corresponding to infection compartments. Used to automatically compute the inf parameter for NGM analysis from siphon structure.";
DFE::usage="(Boundary1)DFE[mod,inf] yields the disease-free equilibrium of the model";
mRts::usage="(Boundary1)mRts[RN,ks] creates mass action rates with names ks";
JR0::usage="(Boundary1){R0,co} = JR0[pol,u] computes basic reproduction number from polynomial";
extHD::usage="(Boundary1){hd,li} = extHD[pol,u] factors polynomial, extracts high degree and linear factors";
bd1::usage="(Boundary1)bd1[RN,rts] boundary analysis for single strain: expects DFE and one endemic equilibrium. Returns {RHS,var,par,cp,mSi,Jx,Jy,E0,ngm,R0A,EA,E1}";
bd2::usage="(Boundary1){RHS,var,par,cp,mSi,Jx,Jy,E0,ngm,R0A,EA,E1t,E2t}=bd2[RN,rts] 
boundary analysis for two strains (assuming DFE, 
two boundary points, and one coexistence equilibrium).";
invN2::usage = "(Boundary1)invNr[E1t,E2t,R0A,E0,par,cp,in1,in2,fval,ins] computes invasion numbers 
and persistence conditions for two-strain models. Takes lists of boundary equilibria E1t, E2t, basic reproduction numbers R0A, disease-free equilibrium E0, parameters par with constraints cp, solution indices in1,in2 to select specific equilibria, optional parameter values fval and indices ins for parameter fixing. Returns {E1,E2,R12,R21,coP} where E1,E2 are selected equilibria, R12,R21 are invasion numbers, 
and coP are parameter conditions ensuring both invasion numbers >1 for coexistence.";
invN::usage = "invN[E1t, E2t, R0A, E0, par, cp, in1, in2, fval, ins] computes invasion numbers and persistence conditions for multi-strain epidemic models. E1t and E2t are lists of boundary fixed points for strains 1 and 2. R0A contains basic reproduction numbers. E0 is the disease-free equilibrium. par and cp are model parameters and constraints. in1 and in2 are indices selecting equilibria from E1t and E2t (use -1 for irrational equilibria). Optional: fval gives values for fixed parameters, ins specifies which parameters to fix. Returns {E1, E2, R12, R21, coP} where R12, R21 are invasion numbers and coP contains coexistence parameter values. For irrational equilibria, returns 'nonRat' 
for the equilibrium and invasion number, 'unknown' for coP.";bdCo::usage="(Boundary1){RHS, var, par, cp, mSi, Jx, Jy, E0, ngm, R0A, EA, E1t, E2t}=
bdCo[RN,rts] boundary analysis for two strains with nondisjoint minSiph";
bdAn::usage="(Boundary1){RHS,var,par,cp,mSi,Jx,Jy,cDFE,E0,K,R0A,infVars,alMat,beMat,gam,ngm}=
bdAn[RN,rts] boundary analysis for possibly nonrational bdFp. Returns RHS (right-hand side), var (variables), par (parameters), cp (parameter constraints), mSi (minimal siphons), Jx (Jacobian x-block), Jy (Jacobian y-block), cDFE (disease-free equilibrium conditions), E0 (DFE solution), K (next generation matrix), R0A (basic reproduction numbers), infVars (infection variables), alMat (reactant matrix alpha), beMat (product matrix beta), gam (stoichiometric matrix gamma), ngm (full NGM output).";
bdFp::usage = "bdFp[RHS, var, mSi] computes boundary fixed points on siphon facets for epidemic models. RHS is the right-hand side vector of the ODE system, var is the list of all variables as symbols, and mSi is the list of minimal siphons as variable names (strings). Returns a list of pairs {{rationalSols, scalarEq}, ...}, one for each siphon facet, where rationalSols contains explicit rational solutions as replacement rules and scalarEq is either None (if all solutions are rational), a factored polynomial equation 
(for non-rational solutions), or \"froze\" (if computation timed out).";
sta::usage = "sta[pol] analyzes polynomial stability by 
factoring pol and examining linear and quadratic factors separately. Assumes all parameters are positive. For linear factors (au + b), stability requires b > 0. For quadratic factors (au\.b2 + bu + c), stability requires both b > 0 and c > 0 (Routh-Hurwitz criteria). Returns {lSta, qSta, hDeg} where lSta contains linear stability conditions, qSta contains quadratic stability conditions, and hDeg contains higher-degree factors (degree \[GreaterEqual] 3) requiring manual stability analysis. Polynomial variable should be u.";

(* Bifurcation analysis functions *)
fpHopf::usage="(Bifurcation)fpHopf[RHS,var,par,p0val] finds positive fixed points and analyzes for Hopf bifurcations. Returns {posSols,complexEigs,angle,eigs} where angle is ArcTan[Re/Im]*180/Pi of complex eigenvalues";
simpleOptHopf::usage="(Bifurcation)simpleOptHopf[RHS,var,par,coP,optInd,numTries] simple optimization to maximize angle (find Hopf bifurcations) by varying parameters at optInd positions. Returns {bestAngle,bestValues,finalP0Val}";
optHopf::usage="(Bifurcation)optHopf[RHS,var,par,coP,optInd,timeLimit,method,accGoal,precGoal,maxIter] sophisticated optimization using NMaximize to find Hopf bifurcations. Returns {bestAngle,bestValues,finalP0Val}";
cont::usage="(Bifurcation)cont[RHS,var,par,p0val,stepSize,plotInd,bifInd,analyticalPlot] continuation analysis tracking equilibria along parameter bifInd, plotting in parameter space plotInd. Returns curve data with stability analysis";
hopfD::usage="(Bifurcation)hopfD[curve] analyzes continuation curve for Hopf bifurcation points by detecting stability changes. Returns list of Hopf points with parameters";

pertIC::usage="(Bifurcation)pertIC[equilibrium,var,factor,minq,n] generates n perturbed initial conditions around equilibrium with perturbation factor";
TS::usage="(Bifurcation)TS[RHS,var,par,p0val,tmax] time series simulation using NDSolve with BDF method. Returns solution functions";
intEq::usage="(Bifurcation)intEq[RHS,var,par,p0val,att] finds equilibria using perturbed initial conditions around attractor att. Returns list of positive solutions";
Bifp::usage="(Bifurcation)Bifp[RHS,var,par] performs bifurcation point analysis";
mSim::usage="(Bifurcation)mSim[RHS,var,par,p0val,opts] performs multiple simulations 
for bifurcation analysis";
scan::usage = "(Bifurcation)scan[RHS,var,par,coP,plotInd,mSi,mode,steTol,staTol,choTol,
R01,R02,R12,R21] performs  parameter space scanning for multi-strain epidemic models 
with siphon-aware analysis. Takes ODE system RHS, variables var, parameters par, 
constraints coP (as rules), plot parameter indices plotInd, minimal siphons mSi for 
invasion analysis, scanning mode, step tolerance steTol for numerical integration, 
stability tolerance staTol for eigenvalue analysis, choice tolerance choTol for 
equilibrium classification. 
Reproduction numbers R01,R02 for strains and invasion numbers R12,R21 
for strain interaction analysis. Returns {fPl,unC,res} where fPl is final plot with 
equilibrium classification, unC contains uncategorized points requiring further analysis, 
res is raw scanning results with detailed siphon and invasion number computations. 
Designed for models where minimal siphons and invasion dynamics are central to the analysis. 
Uses specialized tolerance controls for multi-scale epidemic dynamics.";

scanPar::usage = "(Bifurcation)scanPar[RHS,var,par,p0val,plotInd,gridRes,plot,hTol,delta,
wRan,hRan,R01,R02,R21,R12] performs comprehensive parameter space scanning with 
flexible grid resolution and range controls. Takes ODE system RHS, variables var, 
parameters par, reference parameter values p0val (as list), plot parameter indices plotInd,
 grid resolution gridRes (or Automatic for range mode), base plot to overlay on, 
Hopf tolerance hTol for bifurcation detection, step size delta for range mode, width range wRan and height range hRan as fractions of reference values. Reproduction numbers R01,R02 for individual strains and invasion numbers R21,R12 (note reversed order from scan). Returns {finalPlot,errors,res} where finalPlot is publication-ready plot with color-coded equilibrium types and counts, errors contains diagnostic information, res is scanning results with equilibrium stability classification. Supports both grid mode (fixed NxN resolution) and range mode (adaptive stepping). Optimized for general dynamical systems analysis with automatic equilibrium detection and stability analysis via eigenvalue computation. Does not require siphon analysis, 
making it suitable for broader classes of models beyond epidemic systems.";

(* Siphon and persistence functions *)
minSiph::usage="(Siphons)minSiph[species,reactions] finds minimal siphons in a reaction network";
isSiph::usage="(Siphons)isSiph[species,reactions,siphon] checks if a given set forms a siphon in the reaction network";
isDrainable::usage="(Siphons)isDrainable[reactions, speciesSet] checks if speciesSet is drainable for the given reaction network. A set is drainable if there exists a reaction pathway that decreases all species in the set";
isSelfReplicable::usage="(Siphons)isSelfReplicable[reactions, speciesSet] checks if speciesSet is self-replicable for the given reaction network. A set is self-replicable if there exists a reaction pathway that increases all species in the set";
isCritical::usage="(Siphons)isCritical[reactions, speciesSet] checks if speciesSet is critical (no positive conservation law has support contained in the set)";
siphonAnalysis::usage="(Siphons)siphonAnalysis[reactions] provides comprehensive siphon classification for a reaction network. Returns an association with each siphon classified as drainable, self-replicable, critical, and categorized by type";
persistenceAnalysis::usage="(Siphons)persistenceAnalysis[reactions] analyzes network persistence based on drainable siphons. Returns detailed information about persistence, drainable siphons, and extinction threats. Networks without drainable siphons are persistent";
catalysisAnalysis::usage="(Siphons)catalysisAnalysis[reactions] identifies autocatalytic behavior in reaction networks by finding self-replicable critical siphons. Returns information about catalytic sets and autocatalytic potential";
autocatalysisReport::usage="(Siphons)autocatalysisReport[reactions] provides comprehensive analysis of autocatalytic behavior and persistence properties. Returns detailed association with network info, persistence analysis, catalysis analysis, and theoretical insights";
checkPersistence::usage="(Siphons)checkPersistence[RN] determines persistence status of reaction network RN. Returns {status, analysis} where status is 'Persistent', 'Unknown', or 'Non-persistent'";
persistenceReport::usage="(Siphons)persistenceReport[RN] provides comprehensive persistence analysis of reaction network RN with detailed output";
isCatalytic::usage="(Siphons)isCatalytic[RN] checks if the reaction network contains catalytic sets";
findCatalyticSets::usage="(Siphons)findCatalyticSets[RN] identifies all catalytic sets in the reaction network";
findMinimalCriticalSiphons::usage="(Siphons)findMinimalCriticalSiphons[reactions] finds all minimal critical siphons in a reaction network. These are the siphons that pose potential threats to persistence according to the theory";
invFace::usage="(Siphons)invF = invFace[reactions,maxCodim] finds invariant facets of the positive orthant for reaction networks up to specified codimension";
isInvariantFacet::usage="(Siphons)isInvariantFacet[facetSet,reactions] checks if a given set of species forms an invariant facet";

(* Conv functions - Simple one-liner conversions *)
toSum::usage="(Conv)toSum[expr] converts expression to sum format";
toProd::usage="(Conv)toProd[expr] converts expression to product format";
l2L::usage="(Conv)l2L[list] converts lowercase list format to uppercase format";
m2toM::usage="(Conv)m2toM[expr] converts m2 expressions to M format";
remZ::usage="(Conv)remZ[li] removes zeroes from list";
isNNe::usage="(Utils)isNNe[expr] checks if expression is nonnegative (no negative terms and nonzero)";
selZR::usage="(Conv)selZR[con] selects zero rules from conditions";
seZF::usage="(Conv)seZF[so] removes in a list of lists those with a 0";
strEdg::usage="(Conv)strEdg[edges] processes edge structures for graph operations";
countMS::usage="(Conv)countMS[matrix] counts negative coefficients in a matrix";
onlyP::usage="(Conv)onlyP[m] checks whether all the coefficients of the numerator of a rational expression m are nonnegative";
rtS::usage="(Conv)rtS[RHS] extracts reaction terms from RHS";

(* Extra functions - Complex utilities *)
WhereIs::usage="(Extra)WhereIs[functionName] identifies which subpackage contains a function. Works with strings or symbols.";
ListFunctionsByPackage::usage="(Extra)ListFunctionsByPackage[] organizes all EpidCRN1 functions by their source subpackage.";
chRN::usage="(Extra)chRN[reactionList] validates reaction network format. Returns True if all reactions are properly formatted with string species names and positive coefficients, False otherwise.";
mat2Matl::usage="(Extra)mat2Matl[matrix] converts Mathematica matrix to MATLAB format string";
matl2Mat::usage="(Extra)matl2Mat[string] converts MATLAB format string to Mathematica matrix";
matlr2Mat::usage="(Extra)matlr2Mat[string] converts MATLAB row vector string to Mathematica list";
CreateMassActionRate::usage="(Extra)CreateMassActionRate[reactants,kParam] creates mass action rate expression from reactant association and rate parameter";
convNum::usage="(Extra)convNum[expr] converts numerical expressions";
albe::usage="(Extra)albe[RHS,var] extracts stoichiometric matrices al,be,ga from RHS";
RHS2RN::usage="(Extra)RHS2RN[RHS,var] extracts reactions representation of an ODE";

(* Utils functions - Mathematical and stability analysis *)
CofRH::usage="(Utils)CofRH[A] yields coefficients of CharacteristicPolynomial, as required by Routh-Hurwitz theory, normalized so the free coefficient is 1";
CofP::usage="(Utils)CofP[list] yields coefficients of a polynomial as required by Routh-Hurwitz theory, normalized so the free coefficient is 1";
Hur2::usage="(Utils)Hur2[co] yields stability conditions for second-order system";
Hur3M::usage="(Utils){co,h3,ine} = Hur3M[A] yields third-order Hurwitz analysis";
Hur4M::usage="(Utils){co,h4,ine} = Hur4M[A] yields fourth-order Hurwitz analysis";
Hur5M::usage="(Utils){co,h5,ine,H5} = Hur5M[jac] yields fifth-order Hurwitz analysis";
H4::usage="(Utils)H4[co] gives the 4th Hurwitz determinant, needed in Routh-Hurwitz theory";
H6::usage="(Utils)H6[co] computes the 6th Hurwitz determinant for stability analysis";
makeLPM::usage="(Utils)makeLPM[mat] yields the leading principal minors";
ACM::usage="(Utils)ACM[A,k] yields additive compound matrix";
perR::usage="(Utils)perR[M,i,j] swaps rows i and j in matrix M";
perC::usage="(Utils)perC[matrix,cycle] performs a cyclic permutation on the rows (or columns) of the input matrix";
Stab::usage="(Utils)Stab[mod,cfp,cn] analyzes stability at fixed point";
Sta::usage="(Utils)Sta[jac,X,Xv] numeric stability analysis";
JTD::usage="(Utils)JTD[mod,cn] performs JTD analysis";
JTDP::usage="(Utils)JTDP[mod,zeta,cn] performs JTDP analysis with parameter zeta";
Grobpol::usage="(Utils)Grobpol[RHS,var,par,ind,cn] computes a reduced polynomial system by eliminating variables using Gr\[ODoubleDot]bner basis methods";
RUR::usage="(Utils)RUR[mod,ind,cn] attempts to reduce the fixed point system to one with variables specified by the list ind";
GBH::usage="(Utils)GBH[pol,var,sc,cn] performs Gr\[ODoubleDot]bner basis analysis";
L1Planar::usage="(Utils)L1Planar[fg,eq] performs L1 planar analysis";
DerEq::usage="(Utils)DerEq[fg,eq] computes derivative equations";
Stodola::usage="(Utils)Stodola[pol,var] implements Stodola method for polynomial analysis";
red::usage="(Utils)red[re,cond:{}] erases from the output of a Reduce all the conditions in cond";
reCL::usage="(Utils)reCL[re] erases from the output of a Reduce all the conditions in cond";
onePR::usage="(Utils)onePR[cof,cp] outputs conditions that the first and last coefs of a list have different signs";
expon::usage="(Utils)expon computes the maximum power of an expanded form";
posM::usage="(Utils) posM[matrix] keeps all syntactically positive terms.";
whenP::usage="(Utils) whenP[list] reduces inequalities ensuring positivity of all components of the list.";
FposEx::usage="(Utils)FposEx[matrix] extracts first syntactically positive term in a nonnegative matrix";
GetVec::usage="(Utils)GetVec[A,om] extracts vectors, used in L13,L23";
L13::usage="(Utils)L13[A,B,CC,cp] computes L13 coefficient";
L23::usage="(Utils)L23[A,B,CC,DD,EE] computes L23 coefficient";
Res1F::usage="(Utils)Res1F[mod,csr,pol,in,cn] computes first residue form";
Deg::usage="(Utils)Deg[poly,var] computes degree of polynomial";
Hirono::usage="(Utils)Hirono[S,intRows,intCols] implements Hirono method";
DerL::usage="(Utils)DerL[expr] computes derivatives in L format";

(* Visualization functions *)
EucFHJ::usage="(Visualization)EucFHJ[reactions] draws the Euclidean Feinberg-Horn-Jackson graph for a two-species reaction network, showing complexes as red points, reactions as blue arrows with red heads, and the Newton polytope in gray. Example: EucFHJ[{{0,\"A\"},{\"A\"+\"B\",2*\"B\"}}]";
NewtPol::usage="(Visualization)NewtPol[complexCoords] generates the Newton polytope (convex hull) graphics primitive from an association of complex names to 2D coordinates. Returns {LightGray, EdgeForm[{Black,Thin}], ConvexHullMesh[points]} or {} if insufficient points";
FHJ::usage="(Visualization)FHJ[comp,edges,rates,ver,groups] generates the Feinberg-Horn-Jackson graph with vertex coloring and labeling options";
endo::usage="(Visualization)endo[reactions] analyzes a reaction network for endotactic properties. Options: \"ShowPlot\" -> True/False/Automatic (default: Automatic), \"Verbose\" -> True/False (default: False)";
isEndotactic::usage="(Visualization)isEndotactic[reactions, speciesList] checks if network is endotactic";
isStronglyEndotactic::usage="(Visualization)isStronglyEndotactic[reactions, speciesList] checks if network is strongly endotactic";
phasePl2::usage="(Visualization)phasePl2[RHS,var,par,p0val,opts] creates 2D phase plane plots for dynamical systems";

Begin["`Private`"];

root=If[StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[$InputFileName],
    NotebookDirectory[]
];

(* Load in proper order - simple to complex *)
Get[FileNameJoin[{root, "Core.wl"}]];
Get[FileNameJoin[{root, "Siphons.wl"}]];
Get[FileNameJoin[{root, "Boundary.wl"}]];
Get[FileNameJoin[{root, "Bifurcation.wl"}]];
Get[FileNameJoin[{root, "Conv.wl"}]];
Get[FileNameJoin[{root, "Utils.wl"}]];
Get[FileNameJoin[{root, "CRNT.wl"}]];
Get[FileNameJoin[{root, "Visualize.wl"}]];
Get[FileNameJoin[{root, "Extra.wl"}]];

End[];
EndPackage[];





