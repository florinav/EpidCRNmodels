(* ::Package:: *)

(* EpidCRN Main Package - Loader for all subpackages *)
BeginPackage["EpidCRN`"];

(* Global variables *)
Global`ome;

(* Core functions *)
extMat::usage = "{spe, al, be, gamma, Rv, RHS, def} = extMat[reactions] 
 Returns species list, 
alpha matrix (reactants), beta matrix (products), gamma matrix (net stoichiometric), 
reaction rate vector, RHS of mass action ODEs, 
and deficiency as {formula, terms, result}.";

compToAsso::usage = "compToAsso[side] parses a reaction side (left or right) and returns an association of species names to stoichiometric coefficients. Example: compToAsso[k*\"i\" + 2*\"s\"] returns <|\"i\"->k, \"s\"->2|>";

extSpe::usage = "extSpe[reactions] extracts all species names from a reaction network. Returns a list of unique species strings.";

asoRea::usage = "asoRea[RN] transforms classic reaction network format into association format with \"Substrates\" and \"Products\" keys.";

convertReactionFormat::usage = "convertReactionFormat[reactions] converts arrow format reactions {A->B, C->D} to pair format {{A,B}, {C,D}}. Returns input unchanged if already in pair format.";

strToSymb::usage = "strToSymb[reactions] converts string format reactions to symbolic format for internal processing.";

symbToStr::usage = "symbToStr[complex] converts symbolic expressions to string format.";

stoichiometricMatrices::usage = "{alpha, beta, gamma, species} = stoichiometricMatrices[reactions] creates stoichiometric matrices for a reaction network.";

reaToRHS::usage = "{RHS, species, Rv} = reaToRHS[reactions] generates the right-hand side of the ODE system for a reaction network using mass action kinetics.";

expM::usage = "expM[var,expo] gives the vector var at power in matrix expo using Inner[OperatorApplied[Power],#2,#1,Times]&";


(* CRNT functions *)
IaFHJ::usage="IaFHJ[vert,edg] analyzes Feinberg-Horn-Jackson graph structure, returns {incidenceMatrix, tableForm}";
IkFHJ::usage="IkFHJ[vert,edg,tk] computes Ik matrix for FHJ analysis";
SpeComInc::usage="SpeComInc[comp,spec] creates species-complex incidence matrix";
getComE::usage="getComE[RN] extracts complexes and edges from reaction network, returns {complexes, edges}";
lapK::usage="lapK[RN, rates] computes Laplacian matrix of reaction network";


(* Boundary analysis functions *)
NGM::usage="{Jx,F,V,K,Jy,Jxy,Jyx,chp,Kd} = NGM[mod,inf] yields NGM analysis components";
DFE::usage="DFE[mod,inf] yields the disease-free equilibrium of the model";
mRts::usage="mRts[RN,ks] creates mass action rates with names ks";
JR0::usage="{R0,co} = JR0[pol,u] computes basic reproduction number from polynomial";
extHD::usage="{hd,li} = extHD[pol,u] factors polynomial, extracts high degree and linear factors";

bd1::usage="bd1[RN,rts] boundary analysis for single strain: expects DFE and one endemic equilibrium. Returns {RHS,var,par,cp,mSi,Jx,Jy,E0,ngm,R0A,EA,E1}";

bd2::usage="bd2[RN,rts] boundary analysis for two strains: expects DFE, two boundary points, and one coexistence equilibrium. Returns {RHS,var,par,cp,mSi,Jx,Jy,E0,ngm,R0A,EA,E1t,E2t}";


(* Bifurcation analysis functions *)
fpHopf::usage="fpHopf[RHS,var,par,p0val] finds positive fixed points and analyzes for Hopf bifurcations. Returns {posSols,complexEigs,angle,eigs} where angle is ArcTan[Re/Im]*180/Pi of complex eigenvalues";

simpleOptHopf::usage="simpleOptHopf[RHS,var,par,coP,optInd,numTries] simple optimization to maximize angle (find Hopf bifurcations) by varying parameters at optInd positions. Returns {bestAngle,bestValues,finalP0Val}";

optHopf::usage="optHopf[RHS,var,par,coP,optInd,timeLimit,method,accGoal,precGoal,maxIter] sophisticated optimization using NMaximize to find Hopf bifurcations. Returns {bestAngle,bestValues,finalP0Val}";

cont::usage="cont[RHS,var,par,p0val,stepSize,plotInd,bifInd,analyticalPlot] continuation analysis tracking equilibria along parameter bifInd, plotting in parameter space plotInd. Returns curve data with stability analysis";

hopfD::usage="hopfD[curve] analyzes continuation curve for Hopf bifurcation points by detecting stability changes. Returns list of Hopf points with parameters";

scanPar::usage="scanPar[RHS,var,par,p0val,plotInd,gridRes,plot,hTol,delta,wRan,hRan,R01,R02,R21,R12] comprehensive parameter space scanning with equilibrium classification. Returns {finalPlot,errors,res}";

scanPR::usage="scanPR[RHS,var,par,coP,plotInd,R01,R02,R21,R12,rangeB1,rangeB2,stepSize] simple parameter range scanning for debugging";

pertIC::usage="pertIC[equilibrium,var,factor,minq,n] generates n perturbed initial conditions around equilibrium with perturbation factor";

TS::usage="TS[RHS,var,par,p0val,tmax] time series simulation using NDSolve with BDF method. Returns solution functions";

intEq::usage="intEq[RHS,var,par,p0val,att] finds equilibria using perturbed initial conditions around attractor att. Returns list of positive solutions";


(* Siphon and persistence functions *)
minSiph::usage="minSiph[species,reactions] finds minimal siphons in a reaction network";
isSiph::usage="isSiph[species,reactions,siphon] checks if a given set forms a siphon in the reaction network";

isDrainable::usage="isDrainable[reactions, speciesSet] checks if speciesSet is drainable for the given reaction network. A set is drainable if there exists a reaction pathway that decreases all species in the set";

isSelfReplicable::usage="isSelfReplicable[reactions, speciesSet] checks if speciesSet is self-replicable for the given reaction network. A set is self-replicable if there exists a reaction pathway that increases all species in the set";

isCritical::usage="isCritical[reactions, speciesSet] checks if speciesSet is critical (no positive conservation law has support contained in the set)";

siphonAnalysis::usage="siphonAnalysis[reactions] provides comprehensive siphon classification for a reaction network. Returns an association with each siphon classified as drainable, self-replicable, critical, and categorized by type";

persistenceAnalysis::usage="persistenceAnalysis[reactions] analyzes network persistence based on drainable siphons. Returns detailed information about persistence, drainable siphons, and extinction threats. Networks without drainable siphons are persistent";

catalysisAnalysis::usage="catalysisAnalysis[reactions] identifies autocatalytic behavior in reaction networks by finding self-replicable critical siphons. Returns information about catalytic sets and autocatalytic potential";

autocatalysisReport::usage="autocatalysisReport[reactions] provides comprehensive analysis of autocatalytic behavior and persistence properties. Returns detailed association with network info, persistence analysis, catalysis analysis, and theoretical insights";

checkPersistence::usage="checkPersistence[RN] determines persistence status of reaction network RN. Returns {status, analysis} where status is 'Persistent', 'Unknown', or 'Non-persistent'";

persistenceReport::usage="persistenceReport[RN] provides comprehensive persistence analysis of reaction network RN with detailed output";

isCatalytic::usage="isCatalytic[RN] checks if the reaction network contains catalytic sets";

findCatalyticSets::usage="findCatalyticSets[RN] identifies all catalytic sets in the reaction network";

findMinimalCriticalSiphons::usage="findMinimalCriticalSiphons[reactions] finds all minimal critical siphons in a reaction network. These are the siphons that pose potential threats to persistence according to the theory";

cons::usage="cons[mat,cp] parametrizes positively the left kernel of mat, using also conditions cp; cp is not necessary if mat is numeric";

invFacet::usage="invF = invFacet[reactions,maxCodim] finds invariant facets of the positive orthant for reaction networks up to specified codimension";

isInvariantFacet::usage="isInvariantFacet[facetSet,reactions] checks if a given set of species forms an invariant facet";

CofRH::usage="CofRH[A] yields coefficients of CharacteristicPolynomial, as required by Routh-Hurwitz theory, normalized so the free coefficient is 1";
CofP::usage="CofP[list] yields coefficients of a polynomial as required by Routh-Hurwitz theory, normalized so the free coefficient is 1";
Hur2::usage="Hur2[co] yields stability conditions for second-order system";
Hur3M::usage="{co,h3,ine} = Hur3M[A] yields third-order Hurwitz analysis";
Hur4M::usage="{co,h4,ine} = Hur4M[A] yields fourth-order Hurwitz analysis";
Hur5M::usage="{co,h5,ine,H5} = Hur5M[jac] yields fifth-order Hurwitz analysis";
H4::usage="H4[co] gives the 4th Hurwitz determinant, needed in Routh-Hurwitz theory";
H6::usage="H6[co] computes the 6th Hurwitz determinant for stability analysis";

(* =========================================================================*)
(*MATRIX AND LINEAR ALGEBRA UTILITIES*)
(* =========================================================================*)

makeLPM::usage="makeLPM[mat] yields the leading principal minors";
ACM::usage="ACM[A,k] yields additive compound matrix";
perR::usage="perR[M,i,j] swaps rows i and j in matrix M";
perC::usage="perC[matrix,cycle] performs a cyclic permutation on the rows (or columns) of the input matrix";

(* =========================================================================*)
(*STABILITY ANALYSIS*)
(* =========================================================================*)

Stab::usage="Stab[mod,cfp,cn] analyzes stability at fixed point";
Sta::usage="Sta[jac,X,Xv] numeric stability analysis";
JTD::usage="JTD[mod,cn] performs JTD analysis";
JTDP::usage="JTDP[mod,zeta,cn] performs JTDP analysis with parameter zeta";

(* =========================================================================*)
(*POLYNOMIAL AND ALGEBRAIC TOOLS*)
(* =========================================================================*)

Grobpol::usage="Grobpol[RHS,var,par,ind,cn] computes a reduced polynomial system by eliminating variables using Gr\[ODoubleDot]bner basis methods";
RUR::usage="RUR[mod,ind,cn] attempts to reduce the fixed point system to one with variables specified by the list ind";
GBH::usage="GBH[pol,var,sc,cn] performs Gr\[ODoubleDot]bner basis analysis";
L1Planar::usage="L1Planar[fg,eq] performs L1 planar analysis";
DerEq::usage="DerEq[fg,eq] computes derivative equations";
Stodola::usage="Stodola[pol,var] implements Stodola method for polynomial analysis";

(* =========================================================================*)
(*CONVERSION AND FORMAT UTILITIES*)
(* =========================================================================*)

mat2Matl::usage="mat2Matl[matrix] converts Mathematica matrix to MATLAB format string";
matl2Mat::usage="matl2Mat[string] converts MATLAB format string to Mathematica matrix";
matlr2Mat::usage="matlr2Mat[string] converts MATLAB row vector string to Mathematica list";
l2L::usage="l2L[list] converts lowercase list format to uppercase format";
m2toM::usage="m2toM[expr] converts m2 expressions to M format";
toSum::usage="toSum[expr] converts expression to sum format";
toProd::usage="toProd[expr] converts expression to product format";

(* =========================================================================*)
(*GENERAL UTILITIES*)
(* =========================================================================*)

remZ::usage="remZ[li] removes zeroes from list";
rtS::usage="rtS[RHS] extracts reaction terms from RHS";
albe::usage="albe[RHS,var] extracts stoichiometric matrices al,be,ga from RHS";
RHS2RN::usage="RHS2RN[RHS,var] extracts reactions representation of an ODE";
selZR::usage="selZR[con] selects zero rules from conditions";
red::usage="red[re,cond] erases from the output of a Reduce all the conditions in cond";
reCL::usage="reCL[re,cond] erases from the output of a Reduce all the conditions in cond";
seZF::usage="seZF[so] removes in a list of lists those with a 0";
onePR::usage="onePR[cof,cp] outputs conditions that the first and last coefs of a list have different signs";
expon::usage="expon computes the maximum power of an expanded form";
posM::usage="posM[matrix] keeps all syntactically positive terms";
FposEx::usage="FposEx[matrix] extracts first syntactically positive term in a nonnegative matrix";
onlyP::usage="onlyP[m] checks whether all the coefficients of the numerator of a rational expression m are nonnegative";
CreateMassActionRate::usage="CreateMassActionRate[reactants,kParam] creates mass action rate expression from reactant association and rate parameter";
Par::usage="Par[RHS,var] extracts parameters from dynamics";
strEdg::usage="strEdg[edges] processes edge structures for graph operations";
countMS::usage="countMS[matrix] counts negative coefficients in a matrix";
convNum::usage="convNum[expr] converts numerical expressions";

(* =========================================================================*)
(*SPECIALIZED ANALYSIS FUNCTIONS*)
(* =========================================================================*)

GetVec::usage="GetVec[A,om] extracts vectors, used in L13,L23";
L13::usage="L13[A,B,CC,cp] computes L13 coefficient";
L23::usage="L23[A,B,CC,DD,EE] computes L23 coefficient";
Res1F::usage="Res1F[mod,csr,pol,in,cn] computes first residue form";
Deg::usage="Deg[poly,var] computes degree of polynomial";
Hirono::usage="Hirono[S,intRows,intCols] implements Hirono method";
DerL::usage="DerL[expr] computes derivatives in L format";


Begin["`Private`"];
root=Which[StringQ[$InputFileName]&&
$InputFileName=!="",DirectoryName[$InputFileName],
StringQ[NotebookFileName[]]&&NotebookFileName[]=!="",
NotebookDirectory[],True,Directory[]];Get/@FileNameJoin[{root,#<>".wl"}]&/@
{"Core","CRNT","Boundary","Bifurcation","Siphons","Utils"};



(*Get["EpidCRN`Visualization`"];*)

End[];
EndPackage[];
