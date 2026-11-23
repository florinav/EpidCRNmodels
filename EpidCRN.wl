(* ::Package:: *)

(* EpidCRN1 Main Package - Loader for all subpackages *)
BeginPackage["EpidCRN`"];
(* Global variables *)
Global`ome;

(*Core funtions*)
reS::usage = "reS[ine,cp]:=FullSimplify[Reduce[Append[cp,ine]]//FullSimplify,cp]";
cSol::usage = "summary = cSol[RHS,var] or summary = cSol[RHS,var,prF]  solves RHS==0, classifies solutions by algebraic degree. summary = {{\"Total\",n},{\"Rational\",n,idx},{\"SquareRoots\",n,idx},{\"Degree3\",n,idx},{\"Degree>3\",n,idx},{\"AllNonzero\",n,idx},{\"Solutions\",so}}. Optional prF applies pretty-printing to the solutions.";
(*chkSol;*)
minSiph::usage = "(Core){mSi,cDFE,cE0,nonm}=minSiph[var,RN,RHS] or {mSi,nonm}=minSiph[var,RN] finds minimal siphons. var=variable list, RN=reaction network, RHS=ODE right-hand side (optional). mSi=list of minimal siphons (each a list of species strings), cDFE=DFE substitution rules (all siphon species->0), cE0=DFE equilibrium solution (requires RHS), nonm=species not in any minimal siphon.";
extMat::usage = "{spe, al, be, gamma, Rv, RHS, def} = extMat[reactions]
Returns species list, alpha matrix (reactants), beta matrix (products),
gamma matrix (net stoichiometric), reaction rate vector, RHS of mass action ODEs,
and deficiency as {formula, terms, result}.";
NGM::usage = "NGM[mod, infVars, F] computes Next Generation Matrix.
mod={RHS,var,par}. If F provided, validates Inverse[V] all positive under par>0.
If not or invalid, uses posM. Returns {Jx,F,V,K,Jy,Jxy,fx,Kd} where fx=RHS[[inf]]/.{infVars->0}.
Admissible opaque functions satisfy f[...,0,...]=0 (ODE positivity). zeroOpaque zeroes
f-values when any siphon arg is 0, and derivatives Derivative[ords][f][args] when
some j has ords[[j]]=0 and args[[j]]=0 (differentiating the admissibility identity
in non-j directions). Derivatives in direction j at args[[j]]=0 (e.g. df/dB1 at B1=0)
are left symbolic as they are the non-trivial invasion rates.";
staPP::usage = "{lSta,qSta,hDeg,lFac,qFac}=staPP[pol,par]
analyzes polynomial stability by factoring. Automatically assumes par>0.
Linear factors (au+b): requires b>0. Quadratic factors (au^2+bu+c): requires b>0 && c>0 (Routh-Hurwitz).
Returns lSta, qSta as logical expressions (simplified using par>0), hDeg as list of higher-degree factors (deg>=3),
lFac as list of linear factors, qFac as list of quadratic factors.";
Hur3::usage="ine=Hur3[pol,var] Routh-Hurwitz for degree-3 polynomial. Assumes parameters positive (cp=Complement[Variables[pol],{var}]). Prints Reduce result if it evaluates within 5s. Returns ine (list of 3 inequalities).";
Hur4::usage="ine=Hur4[pol,var] Routh-Hurwitz for degree-4 polynomial. Assumes parameters positive (cp=Complement[Variables[pol],{var}]). Prints Reduce result if it evaluates within 5s. Returns ine (list of 4 inequalities).";
resHopf::usage="{co,res}=resHopf[M] monic char poly coefficients (c0=1,..,cn) and Hopf resultant (Resultant of Re,Im parts of p(I*w) w.r.t. w).";
Hur4M::usage="{co,h4,ine}=Hur4M[mat] computes Hurwitz determinant for 4x4 matrix. Returns coefficients co, determinant h4, inequalities ine.";
cons::usage="cons[gam] finds positive conservation laws from gamma matrix. Returns list of positive left nullspace vectors w such that w.gam=0.";
ODE2RN::usage="{RN,rts,spe,alp,bet,gam}=ODE2RN[RHS,var] converts ODE system to reaction network. RHS is list of expressions, var is list of variables. Returns reaction network RN, rates vector rts, species list spe, reactant matrix alp (alpha), product matrix bet (beta), and net stoichiometric matrix gam (gamma). Algorithm: identifies monomials, uses negative terms to define rates, constructs stoichiometric matrices where sources use alpha=0, then beta=gamma+alpha. Inverse operation to extMat";
ODE2RNp::usage="{RN,rts,spe,alp,bet,gam,rnRed}=ODE2RNp[RHS,var,prF] 
converts ODE system with rational rates to reaction network RN. 
Separates numerators/denominators, analyzes polynomial model, returns 
true rates multiplied by 1/denominators. Prints also rnRed, a version of RN where common factors 
appearing
both in alp and bet are printed above the transition arrow";

actRea::usage = "(Organizations)actRea[RN,spe] returns the sublist of RN whose reactant species are all in spe (React(rho) subseteq spe). Inflow reactions (lhs=0) are always activated.";
speClo::usage = "(Organizations)speClo[RN,spe] returns the species closure of spe: fixed-point iteration adding products of every reaction activated by the current set.";
isClo::usage = "(Organizations)isClo[RN,spe] returns True iff spe equals its own species closure (speClo[RN,spe]===Union[spe]).";
ercClo::usage = "(Organizations)ercClo[RN,rxn] returns the elementary reaction closure of a single reaction rxn: speClo[RN,React(rxn)].";
allErc::usage = "(Organizations)allErc[RN] returns {{rxn,closure},...} for every reaction in RN, deduplicated by closure (distinct ERCs, first generator kept).";
isReac::usage = "(Organizations)isReac[RN,spe] returns True iff spe is closed and actRea[RN,spe] is nonempty (reactive closed set).";
isSM::usage = "(Organizations)isSM[RN,spe] returns True iff spe is self-maintaining: exists v>=1 on every reaction activated by spe with gamma.v>=0 on the species of spe (LP feasibility via FindInstance).";
isSSM::usage = "(Organizations)isSSM[RN,spe] returns True iff spe is semi-self-maintaining: exists v>=0, not all zero, on the reactions activated by spe with gamma.v>=0 on the species of spe.";
isOrg::usage = "(Organizations)isOrg[RN,spe] returns True iff spe is a chemical organization: isClo[RN,spe] && isSM[RN,spe].";
allOrg::usage = "(Organizations)allOrg[RN,Method->\"Exhaustive\"|\"ERCGenerated\"] returns the list of chemical organizations of RN. Exhaustive checks all species subsets; ERCGenerated checks closures of joins of distinct ERCs (allErc).";
synQ::usage = "(Organizations)synQ[RN,ercList] returns {synergyQ,witnessRxns,joinClosure}: True iff some reaction is activated by the closure of the union of ercList but by none of the individual ERCs; witnessRxns are those reactions.";
cmpERC::usage = "(Organizations)cmpERC[RN,e1,e2] returns {complementaryQ,supply12,supply21}: supply12 = products of reactions activated by e1 intersected with the reactant species of e2's own generating reaction(s) (and symmetrically supply21).";
isOutC::usage = "(Buffering)isOutC[RN,m,r] returns True iff (m,r) is output-complete: every reaction of RN whose reactant intersects species set m lies in r.";
locIdx::usage = "(Buffering)locIdx[RN,m,r] returns rank(gamma[[All,r]])-Length[m], using ALL species rows of the net stoichiometric matrix gamma (from extMat), never only the rows in m.";
isBuff::usage = "(Buffering)isBuff[RN,m,r] returns True iff (m,r) is an Okada-Mochizuki buffering structure: isOutC[RN,m,r] && locIdx[RN,m,r]==0.";
allBuff::usage = "(Buffering)allBuff[RN,Method->\"Exhaustive\"|\"ERCGenerated\"] returns the list of {m,r} buffering structures of RN, using the forced-minimal r for each candidate m. General form of buffS (buffS additionally minimizes m by inclusion).";
compToAsso::usage = "(Core1)compToAsso[side] parses a reaction side (left or right) and returns an association of species names to stoichiometric coefficients. Example: compToAsso[k*\"i\" + 2*\"s\"] returns <|\"i\"->k, \"s\"->2|>";
extSpe::usage = "(Core1)extSpe[reactions] extracts all species names from a reaction network. Returns a list of unique species strings.";
comp2Asso::usage = "(Core)comp2Asso[complex] converts a reaction complex to an association species->stoichiometry. Declared here rather than left private because Buffering.wl calls it, and Buffering.wl is loaded on demand from user files rather than by this loader.";
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
chk::usage="(Core1)chk[RHS] checks if RHS variables conflict with package names (V,F,K,mSi,RN,alp,bet,gam,etc). Returns list of conflicts and prints warning.";
allT::usage="(Core1)allT[expr] extracts all terms from expanded polynomial as list";
isN::usage="(Core1)isN[term] returns True if term has negative sign, False otherwise";
symRM::usage="(Core)symRM[alp,n,m] builds m x n symbolic reactivity matrix of r[j,k] or 0";
rSet::usage="(Core)rSet[term] extracts list of r[j,k] symbols from a monomial";
csInfo::usage="(Core)csInfo[k,term,gam] returns {k,term,sIdx,rIdx,Sk,det} for a child selection, or Nothing if non-square";
upfQ::usage="(Core)upfQ[{k,...,detSk}] returns True if child selection has unstable positive feedback";
hurQ::usage="(Core)hurQ[Sk] returns True if matrix Sk is Hurwitz-stable (all eigenvalues have negative real part)";
negCS::usage="(Core)negCS[co,gam,n] extracts all negative child selections from char poly coefficients co";
hasseCS::usage="(Core){gH,edges}=hasseCS[nCS] builds Hasse diagram of r-set inclusions among child selections";
unfQ::usage="(Core)unfQ[{k,...,detSk}] returns True if child selection has unstable negative feedback";
upfCS::usage="(Core)upfCS[nCS,gH,spe,rts] finds UPF feedbacks with Hurwitz-stable supersets, prints species and rates";
oscRI::usage="(Core)oscRI[nCS,gH,spe,rts] Recipe I: finds minimal stable supersets for each UPF (Oscillatory Cores of Class I), prints whether Recipe I holds and the networks that achieve it";
oscRII::usage="(Core)oscRII[nCS,gH,spe,rts] Recipe II: finds minimal UNFs with a Hurwitz-stable principal submatrix (Oscillatory Cores of Class II)";
fixSgn::usage="(Core)fixSgn[nCS,n] fixed-sign test [Cor 5.7/5.8]: for each order k=1..n, checks if all k-CS determinants share the same sign. Returns {{k,isFixed,signs},...}";
csExamples::usage="(Core)csExamples[nCS,spe,rts] prints one illustrative UPF and one NF child selection (preferring k=2): species, reactions, stoichiometric matrix Sk, det, sign class (-1)^k*det, and eigenvalues";
ODE2RNs::usage="(Core){RN,rts,spe,alp,bet,gam,rnR}=ODE2RNs[RHS,var] like ODE2RN but handles fully symbolic functional rates f[S,i], T[i]. Detects opaque function calls (excludes NumericFunction heads), extracts stoichiometry via Coefficient, delegates polynomial remainder to ODE2RN.";

ODE2WY::usage="(Core){W,Y}=ODE2WY[RHS,var] converts ODE system to W,Y matrix format for wr0RealODE. W is coefficient matrix, Y is exponent matrix. Satisfies RHS==W.x^Y.";
ab2RN::usage="(Core)RN=ab2RN[alp,bet,spe] converts alpha/beta matrices to reaction network format. alp is reactant matrix, bet is product matrix, spe is species list.";
ODE2WR0::usage="(Core){success,RN,rts}=ODE2WR0[RHS,var] finds WR0 (weakly reversible deficiency zero) realization of ODE system if it exists. Returns {True,RN,rts} on success, {False,reason,{}} on failure.";
defi::usage="(Core){Nc,l,s,delta}=defi[RN] computes deficiency of reaction network. Nc=complexes, l=linkage classes, s=rank(gamma), delta=Nc-l-s.";
wrQ::usage="(Core)wrQ[RN] returns True if RN is weakly reversible (each linkage class is strongly connected).";
subsCon::usage="(Core)subsCon[str] converts variable names to subscripted form. Converts \"i12\" to Subscript[i,12], \"i1\" to Subscript[i,1], etc. Works with strings or symbols. Returns original if no digits found.";

FLaws::usage="{y0,mR,R0,Dw,ye,k,xe} = FLaws[A,vbe,ba,AS,La]  Four Laws (Fall/Bonzi) endemic equilibrium for ABa models. A: n x n Metzler-Hurwitz infected-block. vbe: m x n WAIFW matrix. ba: n-entry new-infection vector. AS: m x m Metzler-Hurwitz susceptible-block. La: m-recruitment vector. Returns y0=DFE susceptibles, mR=replacement vector (Law1: R0=y0.mR), Dw=dwell vector (Law2: xe=k*Dw), ye=endemic susceptibles (Law3), k=epidemic intensity (Law4: ye.mR=1), xe=endemic infected.";

extABa::usage="{La,AS,vbe,ba,A} = extABa[RHS,var,yvar,xvar]  Extract ABa model matrices from a bilinear ODE. RHS,var: full ODE system. yvar: susceptible variables (m-vector). xvar: infected variables (n-vector). La: constant recruitment vector (m). AS: susceptible transition matrix (m x m). vbe: WAIFW matrix (m x n), extracted via Coefficient[yRHS[[i]],yvar[[i]]*xvar[[k]]]. ba: entry-probability vector (n). A: infected transition matrix (n x n). Emits extABa::notABa if any y-equation has non-zero ABa residual; extABa::novbe if vbe is all zeros.";
extABa::notABa="extABa: susceptible equation `` (variable ``) has non-zero ABa residual: ``. The model is NOT pure ABa with the given y-variables.";
extABa::novbe="extABa: vbe matrix is all zeros \[LongDash] no susceptible class has bilinear incidence. Cannot determine entry vector ba.";
extABa::notrank1="extABa: new-infection matrix F (at y=1) is not rank 1 \[LongDash] model has multiple independent infection pathways. F = ``.";

extABP::usage="{La,AS,B,A,P,yvar} = extABP[RHS,var,xvar]  General DSI bilinear decomposition. var: all variables. xvar: infected/x-class variables (yvars = complement). Returns La (m-recruitment), AS (m x m susceptible matrix), B (m x n WAIFW matrix), A (n x n infected matrix), P (n x m insertion matrix), yvar (susceptible variables in var order). Algorithm: B[i,j]=-Coeff(yRHS_i, y_i*x_j); A[k,l]=Coeff(xRHS_k,x_l); P[k,i]=Coeff(rem_k,y_i*x_j)/B[i,j] for nonzero B[i,j]. Prints 'not bilinear' and returns $Failed if model is not DSI. Use MatrixRank[P]==1 or MatrixRank[B]==1 to detect the two rank-one simplifying cases (see rem.tex).";
extABP::notbil="extABP: equation `` is not purely bilinear (residual non-zero). Model is not DSI.";
extABP::badP="extABP: inconsistent P at k=``, i=``: ratios=`` differ. Model is not DSI bilinear.";

(* CRNT functions *)
IaFHJ::usage="(CRNT)IaFHJ[vert,edg] analyzes Feinberg-Horn-Jackson graph structure, returns {incidenceMatrix, tableForm}";
IkFHJ::usage="(CRNT)IkFHJ[vert,edg,tk] computes Ik matrix for FHJ analysis";
SpeComInc::usage="(CRNT)SpeComInc[comp,spec] creates species-complex incidence matrix";
getComE::usage="(CRNT)getComE[RN] extracts complexes and edges from reaction network, returns {complexes, edges}";
lapK::usage="(CRNT)lapK[RN, rates] computes Laplacian matrix of reaction network";
whichInv::usage="(CRNT)whichInv[R0A,kv] returns list of R0 values for each invasion block. R0A from NGMRN may be scalar, list (one per block), or matrix K (returns eigenvalues sorted descending). R0>=1 means that block can invade. NOTE: invasion condition is heuristic for Kdt membership; use vbCheck for the vertex-balancing question.";
vbCheck::usage="(CRNT){isVB,defect}=vbCheck[rts,RN,kxv] checks vertex-balancing of flux at a specific (k,x*) point. rts=rate vector from ODE2RN, RN=reaction network, kxv=rules for both k and x*. Returns True if max|defect|<1e-8. For Kdt membership, there must exist x*>0 making this True.";
kdtPoint::usage="(CRNT){isKdt,xStar,defect}=kdtPoint[RHS,var,kv] checks Kdt membership at a numerical parameter point kv. Calls ODE2RN[RHS,var] internally. Finds a positive equilibrium x* via NSolve[RHS/.kv==0,var,Reals], then calls vbCheck. Returns {True,xStar,defect~0} if vertex-balanced (confirmed Kdt member at this point), {False,xStar,defect} if not, {False,{},{}} if no positive equilibrium exists. isKdt=True certifies that the Horn-Jackson function is a global Lyapunov function at this kv.";
GetMassAction::usage="(CRNT){f,rates,G,monomials,vars,flux}=GetMassAction[reactions] builds mass-action ODE from E-graph. reactions=directed edges {src,tgt} between stoichiometry vectors (e.g. {2,0}->{0,1}). Returns ODE f, rate constants rates={\[Kappa]_i}, graph G, monomials \[Kappa]_i x^\[Alpha]_i, species vars, flux symbols {\[Beta]_i}.";
GetEQ::usage="(CRNT)GetEQ[reactions] vertex-balance equilibrium equations: substitutes \[Kappa]->\[Beta] and all vars->1 in the mass-action ODE.";
GetFlux::usage="(CRNT)GetFlux[reactions] returns the flux symbol vector {\[Beta]_1,...,\[Beta]_m} for reactions.";
GetDEandVB::usage="(CRNT){DE,VB}=GetDEandVB[rxnG,rxnH] dynamical equivalence + vertex-balance conditions for E-graphs G and H. DE: coefficient equations forcing f_G(\[Beta])=f_H(\[Gamma]). VB: IncidenceMatrix[H].\[Gamma]=0.";
Get\[Beta]2\[Kappa]::usage="(CRNT){\[Beta]2\[Kappa],\[Kappa]s,vars}=Get\[Beta]2\[Kappa][reactions] returns substitution rules mapping flux \[Beta]_i to mass-action monomials \[Kappa]_i x^\[Alpha]_i, plus rate constants and species.";
DrawEGraph::usage="(CRNT)DrawEGraph[rxn] draws E-graph with stoichiometry vectors as vertex coordinates.";
maxEg::usage="(CRNT)Eg=maxEg[Eg] computes G^max: edges in G^comp for which G admits a realization w.r.t. some weakly reversible H \[SubsetEqual] G^comp. Enumerates strongly connected subgraphs containing each edge (largest first), checks DE feasibility with kappa>0, lambda>0.";
Eg2RN::usage="(CRNT){RN,rts}=Eg2RN[Eg] converts BCH E-graph (directed edges {src,tgt} between stoichiometry vectors) to EpidCRN format. RN = list of {srcMonomial,tgtMonomial}, rts = Subscript[\[Kappa],i]*x^src_i. Compatible with vbCheck, getComE, kdtPoint. Rate constants use BCH convention Subscript[\[Kappa],i]; pass {Subscript[\[Kappa],1]->v1,...,x->x0,...} as kv.";
RN2Eg::usage="(CRNT)Eg=RN2Eg[alp,bet] inverse of Eg2RN. alp,bet = stoichiometry matrices from ODE2RN. Returns directed edges alp[[i]]->bet[[i]].";

(* Boundary1 analysis functions *)
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
invN::usage = "invN[E1t, E2t, R0A, E0, par, cp, in1, in2, fval, ins] computes invasion numbers and persistence conditions for multi-strain Coredemic models. E1t and E2t are lists of boundary fixed points for strains 1 and 2. R0A contains basic reproduction numbers. E0 is the disease-free equilibrium. par and cp are model parameters and constraints. in1 and in2 are indices selecting equilibria from E1t and E2t (use -1 for irrational equilibria). Optional: fval gives values for fixed parameters, ins specifies which parameters to fix. Returns {E1, E2, R12, R21, coP} where R12, R21 are invasion numbers and coP contains coexistence parameter values. For irrational equilibria, returns 'nonRat' 
for the equilibrium and invasion number, 'unknown' for coP.";bdCo::usage="(Boundary){RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ElTRat}=bdCo[RN,rts,var] boundary analysis for two strains with nondisjoint minSiph. Pass var from ODE2RN to preserve variable order. Returns RHS, var, par, cp, mSi, Jx, Jy, cDFE (DFE conditions), E0 (DFE), K (NGM), R0A (reproduction numbers), infVars (infection variable symbols labeling K rows/columns), ElTRat (boundary equilibrium solutions for each siphon).";
bdCom::usage="(Boundary){RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, Esys, EA}=bdCom[RN,rts,var] boundary analysis using bdFp. Uses complementarity principle: for siphon constraints, non-siphon variables are positive so equations can be divided by them, reducing polynomial degree (e.g., degree 2 instead of 4). Esys contains unsolved systems (each {equations, variables}), EA contains rational boundary equilibrium solutions.";
(*NGMRN::usage="(Core){RHS, var, par, cp, mSi, Jx, Jy, Jxy, fx, cDFE, cE0, F, V, K, R0A, infV}=
NGMRN[RN,rts,var]  Pass var and gam from ODE2RN to preserve variable order.
{RN,rts,spe,alp,bet,gam}=ODE2RN[RHS1,var]; fx=RHS[[inf]]/.{infVars->0} (inf flow on siphon).
 infV are infection variable symbols corresponding to rows/columns of K).";*)
comp2Asso1::usage="(Cl)comp2Asso1[expr] parses complex expression with parametric coefficients. Handles arbitrary Times like be*et*\"Vd\".";
extMat1::usage="(Cl)extMat1[RN,var] like extMat but uses comp2Asso1 for parametric stoichiometry.";
NGMRN1::usage="(Cl)NGMRN1[RN,rts,var] like NGMRN but uses extMat1 for parametric stoichiometry.";
bdAnC::usage="(Boundary){RHS,var,par,cp,mSi,Jx,Jy,cDFE,E0,K,R0A,infVars,ngm}=bdAnC[RN,rts,var] 
boundary analysis for closed systems (N conserved). Checks if Sum[RHS]==0. If yes, assumes N=1 and eliminates first variable using conservation law, then applies bdAn to reduced system. Returns same outputs as bdAn but for reduced system.";
bdFp::usage = "bdFp[RHS, var, mSi] computes boundary equilibria on siphon facets. RHS is the ODE right-hand side, var is the list of variables, mSi is the list of minimal siphons (strings). For each siphon, sets siphon variables to 0, solves the system, and filters using onlyNN to remove DFE and negative solutions. Returns list where bdFp[[j]] is either positive solutions (list of rule lists) or {equations, variables} if no positive solutions exist. Typically bdFp[[j,1]] is the biologically relevant endemic equilibrium.";

(* Equilibria classification functions *)
sipLat::usage="(Equilibria)sipLat[mSi,vars] constructs lattice of minimal siphons ordered by decreasing length. mSi is list of minimal siphons from bdAn, vars is variable list. Builds all unions of minimal siphons, removes cDFE (union of all), converts to rules (siphon variables -> 0). Returns list of rules ordered longest first.";
slvLvl::usage="(Equilibria)slvLvl[eqs,vars,rules,solved,dfe] solves equilibrium equations at one siphon length level. eqs are equilibrium equations (RHS==0), vars is variable list, rules are siphon rules of same length, solved is list of already-found equilibria, dfe is disease-free equilibrium. Attempts Solve with 300s timeout for each rule, filters out DFE and duplicates, keeps only rational solutions. Returns {ratSols,timeOuts}.";
slvAll::usage="(Equilibria)slvAll[eqs,vars,sipRules,dfe] iteratively solves all siphon length levels from longest to shortest. Maintains accumulator of solved equilibria to avoid duplicates. For timeouts, attempts RUR (Real Univariate Representation) computation. Returns {ratSols,rurSols,timeOuts}.";
clsEqb::usage="(Equilibria)clsEqb[RHS,var,mSi,E0] classifies equilibria per siphon facet. Returns list of {siphon, nRat, nRur, {ratSols, rurSols}} for each siphon (excluding cDFE). nRat includes DFE, nRur counts equations not solutions. If Sp/Time limit reached, returns \"limit\". DFE filtered from stored solutions.";
clsEq::usage="(Equilibria)clsEq[RHS,var,mSi,E0,prF,tim,keepVar] first runs clsEqb to print all siphon results, then finds endemic equilibrium using unrestricted Solve with timeout tim. If timeout, attempts RUR approach. keepVar (optional) is variable to keep for RUR extraction (default: first var not in any siphon). Prints nSol with Short solutions if successful, nRUR with RUR polynomials if RUR used, or TIME-OUT if both timeout. Filters out DFE and siphon-zero equilibria. Returns list of solutions.";
clsEq1::usage="(Equilibria)clsEq1[RHS,var,mSi,E0,prF,tim,keepVar] same as clsEq. keepVar (optional) allows user to suggest which variable to keep for RUR polynomial extraction.";
clsEqR::usage="(Equilibria)clsEqR[RHS,var,mSi,E0,prF,tim,S1,S2] restricted clsEq excluding siphons overlapping S1 or S2. S1,S2 are siphons from mSi (e.g. {v1,y1}). S2 optional.";

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
R01,R02,R12,R21] performs  parameter space scanning for multi-strain Coredemic models 
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
Uses specialized tolerance controls for multi-scale Coredemic dynamics.";

scanPar::usage = "(Bifurcation)scanPar[RHS,var,par,p0val,plotInd,gridRes,plot,hTol,delta,
wRan,hRan,R01,R02,R21,R12] performs comprehensive parameter space scanning with 
flexible grid resolution and range controls. Takes ODE system RHS, variables var, 
parameters par, reference parameter values p0val (as list), plot parameter indices plotInd,
 grid resolution gridRes (or Automatic for range mode), base plot to overlay on, 
Hopf tolerance hTol for bifurcation detection, step size delta for range mode, width range wRan and height range hRan as fractions of reference values. Reproduction numbers R01,R02 for individual strains and invasion numbers R21,R12 (note reversed order from scan). Returns {finalPlot,errors,res} where finalPlot is publication-ready plot with color-coded equilibrium types and counts, errors contains diagnostic information, res is scanning results with equilibrium stability classification. Supports both grid mode (fixed NxN resolution) and range mode (adaptive stepping). Optimized for general dynamical systems analysis with automatic equilibrium detection and stability analysis via eigenvalue computation. Does not require siphon analysis, 
making it suitable for broader classes of models beyond Coredemic systems.";

(* Siphon and persistence functions *)
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
IGMS::usage="(Siphons){edges,elabs,g,cycles}=IGMS[RN,mSi] builds the Invasion Graph of Minimal Siphons. Nodes are minimal siphons T_1,...,T_n (labeled as subscripts). A directed edge T_i->T_j is drawn when the species in siphon T_j can invade the boundary face where T_i species are absent (i.e. some reaction in RN has a reactant in T_j but a product touching T_i). Prints siphon index table and edge list. Returns edge list, edge label rules, the Graph object, and detected cycles.";
findCores::usage="(Siphons)minimal=findCores[RN,maxSize] finds all minimal autocatalytic cores in RN up to size maxSize. A core is a nonempty set of species W such that some subset of reactions is autocatalytic for W (each species in W is produced at a higher net rate when all species in W are present than when absent). Returns the minimal such sets (no proper autocatalytic subset exists). Uses extSpe and stoichiometric matrices alpha, beta internally.";

(* Conv functions - Simple one-liner conversions *)
toSum::usage="(Conv)toSum[expr] converts expression to sum format";
toProd::usage="(Conv)toProd[expr] converts expression to product format";
l2L::usage="(Conv)l2L[list] converts lowercase list format to uppercase format";
m2toM::usage="(Conv)m2toM[expr] converts m2 expressions to M format";
remZ::usage="(Conv)remZ[li]=Select[li,#=!=0&] removes zeroes from list. See also selZF for removing sublists containing 0.";
remN::usage="(Utils)remN[expr] removes negative terms from expression. Applies DeleteCases to remove all negative numbers and terms with negative coefficients at all levels.";
redTim::usage="(Utils)redTout[con,var,time] applies Reduce[con,var] with timeout. Returns 'timeOut' if computation exceeds time limit (default 300 seconds = 5 minutes). var defaults to empty list, time defaults to 300.";
selZR::usage="(Conv)selZR[con] selects zero rules from conditions";
selZF::usage="(Conv)selZF[so] removes in a list of lists those with a 0";
strEdg::usage="(Conv)strEdg[edges] processes edge structures for graph operations";
countMS::usage="(Conv)countMS[matrix] counts negative coefficients in a matrix";
onlyP::usage="(Conv)onlyP[m] checks whether all the coefficients of the numerator of a rational expression m are nonnegative";
onlyN::usage="(Conv)onlyN[m] checks whether all the coefficients of the numerator of a rational expression m are negative";
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
Hur3M::usage="(Boundary){co,h3,ine} = Hur3M[A] yields third-order Hurwitz analysis";
Hur4M::usage="(Boundary){co,h4,ine} = Hur4M[A] yields fourth-order Hurwitz analysis";
Hur5M::usage="(Boundary){co,h5,ine,H5} = Hur5M[jac] yields fifth-order Hurwitz analysis";
H4::usage="(Utils)H4[co] gives the 4th Hurwitz determinant, needed in Routh-Hurwitz theory";
H6::usage="(Utils)H6[co] computes the 6th Hurwitz determinant for stability analysis";
makeLPM::usage="(Utils)makeLPM[mat] yields the leading principal minors";
ACM::usage="(Utils)ACM[A,k] yields additive compound matrix";
dulP::usage="(Utils){TJ,trm}=dulP[RHS,var,{a,b}] planar Dulac test for the monomial metric p=x^a y^b. In two dimensions J^[2]=tr J, so TJ=pdot/p+tr J=div(p f)/p. TJ<0 on the open quadrant excludes periodic orbits and homoclinic loops; trm is the expanded summand list, giving a term-by-term sufficient sign test, while Reduce[ForAll[var,pos,TJ<0]] gives the exact region.";
Pdot::usage="(Utils)Pdot[P,f,x] state-dependent derivative of a matrix metric along the flow, P_f=Sum[D[P,x[[i]]]*f[[i]],{i,Length[x]}]. D acts elementwise, so P may be a full matrix, not only diagonal.";
TACM::usage="(Utils)B=TACM[J,k,P,f,x] the k-th transformed additive compound matrix, B_k=P_f.P^-1+P.J^[k].P^-1. Generalizes dulP's planar TJ (a 1x1, k=2 special case) to any order k and any invertible state-dependent metric P. Reuses ACM for J^[k] and Pdot for P_f.";
matMu::usage="(Utils)matMu[B,p] matrix measure (logarithmic norm) of B for p in {1,2,Infinity}: mu_1 is the max column sum with off-diagonal entries in absolute value, mu_infty the max row sum, mu_2 the largest eigenvalue of (B+B^T)/2. Compose with TACM to test sup_x matMu[TACM[J,k,P,f,x],p]<0, the k-contraction certificate.";
qLMI::usage="(Utils){feas,P}=qLMI[Jlist,eta] fixed-rate feasibility test for a common quadratic contraction metric on the matrix polytope Jlist={J1,...,JN} (single matrix: Jlist={J}). Finds symmetric P with P>=I and Jl^T.P+P.Jl<=-2*eta*P-I for every vertex via SemidefiniteOptimization. Fixing eta avoids the bilinear eta*P term. Returns {True,P} or {False,Null}.";
qRate::usage="(Utils){eta,P}=qRate[Jlist,etaMax,tol] largest common quadratic-contraction rate for the matrix polytope Jlist, found by bisection on eta using qLMI as the feasibility oracle. etaMax defaults to 10, tol to 10^-4.";
jPoly::usage="(Utils)verts=jPoly[A0,As,hb] matrix polytope vertices for J(x)=A0+Sum hj(x)*Aj, hb={{lb1,ub1},...} the bounds on each hj. Returns the 2^m corner matrices (Tuples[hb]); feed to qRate for a polytope-robust quadratic certificate. Does not derive the affine decomposition or bounds from J,x -- that step is the caller's.";
sosCon::usage="(Utils){cons,vars,Q}=sosCon[M,x,d] builds (without solving) the SOS-matrix feasibility problem for symmetric polynomial matrix M(x): v.M.v is SOS in (x,v) iff it equals z.Q.z for a PSD Gram matrix Q on the monomial basis z=v_i*(monomials of x up to degree d). Separated from solving so several SOS-matrix conditions sharing outer unknowns can be joined into one SDP, as polyMet does.";
sosMat::usage="(Utils){feas,Q}=sosMat[M,x,d] standalone SOS-matrix certificate for a symmetric polynomial matrix M(x) (a bare polynomial p(x) is the 1x1 case {{p}}), via sosCon+SemidefiniteOptimization. Returns {True,Q} or {False,Null}. Following Parrilo's SOS relaxation of polynomial matrix inequalities.";
polyP::usage="(Utils){P,vars}=polyP[x,n,d] symmetric n x n matrix ansatz whose entries are independent degree<=d polynomials in x with fresh unknown coefficients. vars is the free coefficient list.";
polyMet::usage="(Utils){feas,P}=polyMet[J,f,x,d,eta] searches for a degree<=d polynomial contraction metric P(x)>0 certifying -(P_f+J^T.P+P.J+2*eta*P)>0 via SOS (compu.tex's PolynomialContractionMetric), both conditions joined into one SDP so they share the same metric coefficients. At d=0 (constant P, P_f=0) reduces exactly to qLMI[{J},eta]. Fixing eta avoids the bilinear eta*P term, as in qLMI.";
optMon::usage="(Utils){found,{a,b},iters}=optMon[RHS,var,box,pts0,maxIter] exponent search for a monomial contraction metric p=x^a y^b via dulP, box={{aLo,aHi},{bLo,bHi}}. dulP's TJ is affine in (a,b) at fixed (x,y), so this is a semi-infinite linear program, not a nonconvex bi-level search: cutting-plane loop over LinearOptimization (maximize margin over a growing sample set) alternated with exact Reduce[ForAll[...]] verification and FindInstance-extracted counterexamples. Returns {False,Null,iters} if no certificate exists in the box.";
kTau::usage="(Utils)val=kTau[A,k,p,T] the k-shifted log norm tau_{p,k}(A)=tr(A)+(n-k)*mu_q(T.(-A).T^-1) (1/p+1/q=1, T defaults to identity) of Dalin-Ofir-BarShalom-Ovseevich-Bullo-Margaliot, IEEE TAC 2024 (arXiv:2209.01046), Theorem 17. Satisfies mu_{p,T^(k)}(A^[k])<=tau_{p,k}(A) without ever forming the k-additive compound A^[k]; tau_{p,k}(A)<=-eta<0 certifies k-contraction directly. Convex in A, so Max[kTau[#,k,p,T]&/@vertices] tests a whole polytope/box of Jacobians. This is compu.tex's CompoundFreeKContraction, Method->\"DalinEtAl\"; see cecKC for Method->\"GeneralizedLyapunov\" (Cecilia et al. 2026), constant-matrix case only.";
exKC::usage="(Utils)exKC[val,tol] numeric-to-exact reconstruction (elementwise Rationalize to within tol) of an SDP-solver candidate; touches only numeric coefficients, symbols are left alone. compu.tex's \"P_hat -> P_rational\" step.";
posDefQ::usage="(Utils)posDefQ[M] exact positive-definiteness test for a numeric/rational constant matrix via Sylvester's criterion (all leading principal minors positive). No Reduce, no floating point once M is exact.";
verQ::usage="(Utils)verQ[Pex,Jlist,eta] exact verification of a constant quadratic k-contraction certificate: Pex>0 and -(Jl^T.Pex+Pex.Jl+2*eta*Pex)>0 for every vertex Jl in Jlist, via posDefQ. Companion to qLMI/qRate, which return the same certificate shape but numeric.";
exVerQ::usage="(Utils){label,Pex}=exVerQ[P,Jlist,eta,tols] tries exKC at a sequence of tolerances, keeps those verQ confirms exactly, returns the simplest (smallest max denominator). {\"Exact\",Pex} if some tolerance verified, {\"Sampled\",P} (original numeric, unverified) otherwise. compu.tex's ExactifyKCertificate+VerifyKCertificate for the constant case; \"RigorousNumerical\" (interval arithmetic) is not implemented.";
verPoly::usage="(Utils)verPoly[M,x,cond] exact positive-definiteness of a symmetric polynomial matrix M(x) for all x satisfying cond (True=unrestricted), via Sylvester's leading-principal-minor test quantified over x with Reduce[ForAll[...]].";
exVerPoly::usage="(Utils){label,Pex}=exVerPoly[J,f,x,eta,Pnum,cond,tols] the polyMet analogue of exVerQ: rebuilds P-I and the differential contraction-inequality matrix from each exactified candidate and checks both via verPoly, returning {\"Exact\",Pex} for the first tolerance that verifies globally, {\"Sampled\",Pnum} otherwise. Reuses Pdot unchanged.";
cecKC::usage="(Utils){ok,P0,P1,mu0,mu1}=cecKC[A,k] generalized Lyapunov k-contraction certificate for constant xdot=A.x, Andreu/Cecilia-Zoboli-Astolfi-Serres-Andrieu, IEEE TAC 2026, Theorem k_contraction_LMI_cLTI at l=2: finds symmetric P0 (In=(0,0,n)), P1 (In=(k-1,0,n-k+1)), and mu0,mu1 with A^T.Pi+Pi.A<2*mui*Pi (i=0,1) and mu1+(k-1)*mu0<0 -- necessary AND sufficient for k-contraction of a constant A, so sharper than kTau's compound-free but sufficient-only test whenever the identity metric is a bad fit for a non-normal A. Built via LyapunovSolve on the shifted Sylvester equation (A-mu I)^T.P+P.(A-mu I)=-I (their Remark[algorithm]+Lemma[inertia_shifting] item 4, the Ostrowski-Schneider inertia theorem), never an inertia-constrained SDP. Requires 2<=k<=n. State-dependent/polytope case (their Theorem k_contraction_LMI proper) not attempted here.";
inCnt::usage="(Utils){neg,zero,pos}=inCnt[M] exact inertia triple of a symmetric matrix M's eigenvalue signs, matching Cecilia et al.'s In(P)=(pi_-,pi_0,pi_+) convention (negative count first). Exact when M has exact/rational entries.";
exCecKC::usage="(Utils){label,P0,P1,mu0,mu1}=exCecKC[A,k,P0,P1,mu0,mu1,tols] numeric-to-exact reconstruction and independent verification of a cecKC certificate: exactifies elementwise at each tolerance and keeps the first confirming inCnt[P0]===(0,0,n), inCnt[P1]===(k-1,0,n-k+1), both LMI residuals via posDefQ, and the scalar sum condition, all independently of cecKC's own reasoning. Returns {\"Exact\",...} or {\"Sampled\",P0,P1,mu0,mu1}.";
jntCec::usage="(Utils){ok,P0,P1,mu0,mu1,i}=jntCec[Jlist,k] the state-dependent/polytope extension of cecKC: a single constant (P0,P1,mu0,mu1) valid at EVERY vertex of a matrix polytope Jlist (e.g. from jPoly), matching what Cecilia et al.'s Theorem k_contraction_LMI actually asks for on a forward-invariant set, not just one frozen point. Tries cecKC's own exact single-matrix certificate at each vertex in turn against every OTHER vertex's LMI residual, valid by continuity whenever the polytope is not too large; returns the first that survives everywhere, or {False,Null,Null,Null,Null,Null} if none does. No claim of completeness -- a polytope with no single vertex whose certificate extends everywhere may still admit some other joint certificate this strategy does not find.";
perR::usage="(Utils)perR[M,i,j] swaps rows i and j in matrix M";
perC::usage="(Utils)perC[matrix,cycle] performs a cyclic permutation on the rows (or columns) of the input matrix";
JTD::usage="(Utils)JTD[mod,cn] performs JTD analysis";
JTDP::usage="(Utils)JTDP[mod,zeta,cn] performs JTDP analysis with parameter zeta";
Grobpol::usage="(Utils)Grobpol[RHS,var,par,ind,cn] computes a reduced polynomial system by eliminating variables using Gr\[ODoubleDot]bner basis methods";
RUR::usage="(Utils)RUR[mod,ind,cn] attempts to reduce the fixed point system to one with variables specified by the list ind";
GBH::usage="(Utils)GBH[pol,var,sc,cn] performs Gr\[ODoubleDot]bner basis analysis";
L1Planar::usage="(Utils)L1Planar[fg,eq] performs L1 planar analysis";
DerEq::usage="(Utils)DerEq[fg,eq] computes derivative equations";
Stodola::usage="(Utils)Stodola[pol,var] implements Stodola method for polynomial analysis";
red::usage="(Utils)red[re,cond:{}] erases from the output of a Reduce all the conditions in cond";
(*reCL::usage="(Utils)reCL[re] erases from the output of a Reduce all the conditions in cond";*)
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
phase2::usage="(Visualization){Xs,eig,stab,sad,pFP,pStr,rhN}=phase2[RHS,var,cN,tMax,nTraj] creates 2D phase plane with nullclines and streamplot. RHS is list of 2 ODEs, var is list of 2 variables, cN is parameter substitutions, tMax is max time, nTraj is number of trajectories. Returns sorted positive fixed points Xs, eigenvalues eig, stabilities stab, saddles sad, fixed point plot pFP, streamplot pStr, numeric RHS rhN";
phN::usage="(Visualization){Xs,eig,stab,sad,pFP,pStr}=phN[rhN,tMax,nTraj] creates 2D phase portrait for fully numeric RHS. rhN is numeric RHS (all parameters substituted), variables extracted automatically. Optional: tMax (default 50), nTraj (default 15). Returns sorted positive fixed points Xs, eigenvalues eig, stabilities stab, saddles sad, fixed point plot pFP, streamplot pStr";
proj2::usage="(Visualization){fpProj,eig,stab,pFP,pStream}=proj2[RHS,var,projVars,cN] projects n-dimensional system onto 2D plane. RHS is n-dimensional ODE system, var is list of all n variables, projVars is list of 2 variables to project onto, cN is parameter substitutions. Solves nullclines of other (n-2) variables to create projected 2D phase portrait. Returns projected fixed points fpProj, eigenvalues eig (from full nD system), stabilities stab, fixed point plot pFP, streamplot pStream";
phase3::usage="(Visualization){Xs,eig,stab,pFP,pStream}=phase3[rhN,var,tMax,pts,rad,opts] creates 3D phase portrait for fully numeric RHS. rhN is numeric 3D RHS (all parameters substituted), var is list of 3 variables IN THE SAME ORDER as the RHS equations (e.g., {s,i2,i12} if equation 1 is ds/dt, equation 2 is di2/dt, equation 3 is di12/dt). Optional: tMax (default 50), pts is StreamPoints density (default 10), rad is sphere radius for fixed points (default 0.015), opts are StreamPlot3D options (e.g., StreamScale->0.05). Default uses Arrow3D markers with rainbow colors. Returns positive fixed points Xs, eigenvalues eig, stabilities stab, fixed point plot pFP (Graphics3D with colored spheres: Green=stable, Orange=saddle, Red=unstable), streamplot pStream (StreamPlot3D). Use Show[pFP,pStream,ViewPoint->{1.5,-2,1.5}] to combine";
invGr::usage="(Visualization){gra,ver,edg,col}=invGr[mSi,bdFps,R0A,RHS,var,E0] constructs invasion graph for multi-strain models. mSi are minimal siphons (list of species name lists), bdFps are boundary fixed points, R0A are basic reproduction numbers, RHS is ODE right-hand side, var is variable list, E0 is disease-free equilibrium. Returns Graph object gra with invasion edges, vertex list ver (communities as sorted index lists), edge list edg showing invasion pathways, vertex colors col (Green=successful invasion, Yellow=failed invasion). Nodes represent communities (siphon complements), edges show strain invasions when Rj(ySc)>1";
hsGr::usage = "(Boundary){vtx,eqL,rts,edg,acy,inv}=hsGr[RHS,var,pop,par] computes the Hofbauer-Schreiber invasion graph at numerical parameter rules par, for populations pop written as x_i'=x_i f_i. vtx lists the supports K of the boundary and interior equilibria (one entry per equilibrium), eqL the equilibria, rts the invasion rates f_i of all populations at each, edg the edges {a,b} with vtx[[a]] strictly inside vtx[[b]] and every population of vtx[[b]] not in vtx[[a]] invading, acy whether the graph is acyclic, inv whether each vertex is invadable (or full).";
invSen::usage = "(Boundary)sen=invSen[RHS,var,pop,par,prs,hs] lists {K,i,d} for every boundary equilibrium of hs=hsGr[...] and missing population i, where d is the vector of derivatives of the invasion rate f_i(E_K) with respect to the parameters prs, by implicit differentiation on the face, at the numbers par. A zero vector means the rate is protected, as in Theorem t:bufInv.";
hsPlot::usage = "(Visualize)hsPlot[hs] draws the invasion graph from hsGr: supports layered by size, green if invadable or full, red if not, arrows K->L.";
senPlot::usage = "(Visualize)senPlot[sen,prs] draws invSen as a table: rows (resident support : missing population), columns parameters prs, white where the invasion rate does not depend on the parameter, red/blue for positive/negative derivative.";
DSR::usage="(Visualization)DSR[RN] builds the directed species-reaction (SR) bipartite graph for reaction network RN. Returns an association with keys: \"Graph\" -> the Graph object with signed, stoichiometric edge labels, \"EdgeData\" -> list of associations {<|\"edge\"->_, \"sign\"->_, \"stoich\"->_|>, ...}, \"ReactionNodes\" -> list of reaction node labels {\"R1\", \"R2\", ...}, \"SpeciesNodes\" -> list of species node labels";

(*Cl - experimental versions*)
minSiph1::usage="(Cl) Copy of minSiph0. {mSi,cDFE,cE0,nonm}=minSiph1[var,RN,RHS] or {mSi,cDFE,nonm}=minSiph1[var,RN]. Outputs cE0 when RHS present.";
minSiph0::usage="(Cl) Exact copy of minSiph. {mSi,cDFE,cE0,nonm}=minSiph0[var,RN,RHS] or {mSi,cDFE,nonm}=minSiph0[var,RN]. Always outputs cDFE even without RHS.";

NGMRN::usage = "{RHS,par,cp,mSi,Jx,Jy,Jxy,fx,cDFE,cE0,F,V,K,R0A,infV}=NGMRN[RN,rts,var] computes NGM from a reaction network. R0A is an array of per-siphon reproduction numbers: R0A[[k]] = spectral radius of the NGM restricted to minimal siphon mSi[[k]], as a symbolic function of the resident variables. Length[R0A]=Length[mSi].";

(* NGM for Frequency-Dependent incidence (NGMFD.wl) -- for models where a rate is
   normalized by a population total (e.g. beta*S*I/N) whose support spans both the
   siphon species and the resident species; NGMRN's all-zero two-stage linearization
   can hit a literal 0 denominator on such models (Power::infy), since NGMRN itself
   is not to be modified for this without explicit sign-off. NGMFD instead evaluates
   the transversal Jacobian at the true (nonzero-resident) boundary equilibrium. *)
NGMFD::usage = "{JxE,F,V,K,evals,R0}=NGMFD[{RHS,var,par},sigVars,faceEq] computes the NGM on the siphon face sigVars==0, evaluating the transversal Jacobian at the true boundary equilibrium faceEq (resident variables at their solved nonzero values, not 0), correctly handling frequency-dependent incidence where the population-total denominator spans both sides of the siphon partition. Reuses posM and minSiph unchanged.";
NGMFDsub::usage = "faceEq=NGMFDsub[{RHS,var,par},sigVars] solves RHS==0 on the boundary face sigVars==0 for the remaining (resident) variables, like minSiph's own cE0 but for an arbitrary siphon face, not only the total DFE siphon -- use to get faceEq for NGMFD at a one-strain or other boundary equilibrium.";
NGMFDRN::usage = "{RHS,par,mSi,cDFE,cE0,R0A}=NGMFDRN[RN,rts,var] is the frequency-dependent-incidence analog of NGMRN: reuses minSiph unchanged for the siphon lattice and DFE, then computes a per-minimal-siphon reproduction number via NGMFD using minSiph's own cE0 as the evaluation point.";

(* Okada-Mochizuki buffering structures (Buffering.wl) *)
bufLat::usage = "(Buffering)lat=bufLat[RN,spe] returns the COMPLETE list of buffering structures {m,r} of RN (species strings m, reaction indices r): for every species subset m, r=Out(m) (reactions with a reactant in m) when rank(gam[[All,Out(m)]])==|m|, together with every enlargement of Out(m) by reactions whose gam column lies in the span of gam[[All,Out(m)]]. Includes the bottom {{},{}}. Returns {} when rank(gam)<number of species (conservation laws need a reduced formulation). spe (optional) fixes the species order, as in extMat.";
latHas::usage = "(Buffering)cov=latHas[lat] returns the Hasse cover edges {i,j} (lat[[i]] strictly below lat[[j]], nothing in between) of the inclusion order of a list lat of {m,r} pairs, e.g. from bufLat.";
latBlk::usage = "(Buffering){blk,inf}=latBlk[lat,cov] reconstructs the canonical (CCF) blocks from the buffering lattice by Birkhoff duality: blk lists {node,newSpecies,newReactions} for every join-irreducible node (exactly one lower cover); inf lists {a,b} (positions in blk) meaning that perturbing a reaction of block a can move the species of block b, transitively reduced.";
rcMat::usage = "(Buffering){L,dL,fL}=rcMat[RN,spe] builds the symbolic rate-cycle matrix L=[T|C] of RN (T: formal entries FormalT[j,i] where reaction j has reactant i; C: basis of ker gam) and returns its factored determinant dL and the list fL of nonconstant irreducible factors (one per species-bearing canonical block). dL is $Aborted after 120 s.";
ccfJl::usage = "(Buffering){blk,ord,cov,C0,Cinf,R0T,RinfT}=ccfJl[Rpat,gam] computes the combinatorial canonical form of [Rpat|C], C a basis of ker gam, with LayeredMixedMatrices.jl (Kuhrs-Telek-Vassena) via ccfRun.jl; polynomial time, no 2^n loop. Rpat: reactions x species 0/1 rate dependence (any kinetics); gam: species x reactions net stoichiometry. blk={{reactions},{species},{cycleRows}} per block; ord: pairs {k,l} = every buffering structure containing block k contains block l; cov: its Hasse edges.";
ccfRN::usage = "(Buffering)ccfRN[RN,spe] = ccfJl with the mass-action rate pattern of RN (rate depends on reactants).";
ccfIdl::usage = "(Buffering)idl=ccfIdl[nb,ord,cap] lists the order ideals (block index lists) of the CCF block poset of nb blocks, i.e. the buffering structures; stops after cap (default 100000).";
latPlot::usage = "(Buffering)latPlot[lat,cov] draws the Hasse diagram of a buffering lattice lat (from bufLat) with covers cov (from latHas), bottom element at the bottom; each node shows its species set and |R|.";
blkPlot::usage = "(Buffering)blkPlot[blk,inf] draws the canonical blocks blk (from latBlk) with an arrow a->b when perturbing a reaction of block a can move the species of block b (inf from latBlk).";
rnPlot::usage = "(Buffering)rnPlot[RN,spe,{m,r}] draws the species-reaction graph of RN (species disks, reaction squares numbered as in RN) and highlights one buffering structure: species m in orange, reactions r in blue, e.g. rnPlot[RN,spe,lat[[2]]].";
locChk::usage ="(Buffering)locChk[RN,rts,var,m] returns True iff the law of localization holds for the species set m: at a random exact rational point, the rows of J^-1 gam[[All,Out(m)]] for the species outside m vanish, J=gam.D[rts,var]. A False refutes localization; a True is a probabilistic certificate.";
buffS::usage = "{mBuffS,allBuffS}=buffS[RN]: finds Okada-Mochizuki buffering structures of a reaction network RN (ODE2RN/extMat format). For each nonempty species subset m, sets r to the reaction set with a reactant in m (the forced minimal output-complete r), computes the localization index lam=rank(gamma[[All,r]])-Length[m], and keeps (m,r) when lam==0. Returns mBuffS (minimal buffering structures, i.e. no other found structure's species set is a proper subset) and allBuffS (all found), each entry {species,reactions}. Reuses extMat unchanged; does not touch minSiph or NGM.";

(* CEP decision algorithm (CEP.wl) *)
cepLat::usage = "(CEP)cepLat[mSi]: cover skeleton of the siphon lattice. mSi = list of minimal siphons (sorted string lists). Returns {{parentSiphon, childSiphon, invaderBlock},...} for all cover relations Sigma'=Sigma\\sigma with sigma a minimal siphon.";
nStbQ::usage = "(CEP)nStbQ[nd, cv, ps] or nStbQ[nd, cv, ps, tol]: True iff every outgoing T_k = R_k-1 < -tol at node nd. Vacuously True at bottom nodes. cv = {{parent,child,block,R-expr},...}.";
hotCh::usage = "(CEP)hotCh[nd, cv, ps] or hotCh[nd, cv, ps, tol]: outgoing covers at nd with T_k = R_k-1 > tol.";
cepDec::usage = "(CEP)cepDec[cv, topNode, ps, opts]: strict chain walk (alg:CEP numerical mode). Options: tol->10^-10, strict->True. Returns {\"stable\",nd}, {\"threshold\",nd}, or {\"cepViolation\",nd,blocks}.";
cepDecX::usage = "(CEP)cepDecX[cv, topNode, ps, opts]: BFS extension for non-strict/coexistence case. Options: tol->10^-10. Returns {\"stable\",nd} or {\"none\",seen}.";
cepDecSym::usage = "(CEP)cepDecSym[cv, topNode]: symbolic mode (alg:CEP symbolic mode). Enumerates all 2^n sign patterns for cover thresholds T_k=R_k-1. Returns list of {\"region\", conditions, node} for stable patterns. Patterns with cepViolation or threshold are omitted.";
cepCovers::usage = "(CEP)cepCovers[mSi, RHS, var]: builds the full covers list from model data alone. Calls cepLat[mSi] for the skeleton, then for each cover computes the face equilibrium (Solve), transversal Jacobian (Outer[D,...]), regular splitting (posM), and face reproduction number R_sigma = rho(F_sigma.Inverse[V_sigma]). Returns {{parentSiphon, childSiphon, invaderStr, R_sigma},...}. Node identifiers are sorted siphon lists (absent-species sets).";

Sh::usage="{shMat,redRHS,nSh}=Sh[RHS,var] classifies each monomial of Expand/@RHS by sign pattern across equations. Returns: shMat = n x n out-degree Laplacian of the species graph; redRHS = RHS with inflow/outflow monomials removed; nSh = list of NonShuai monomials (empty iff model is a pure Shuai model).";

(* BoolNet - Boolean network analysis and EpidCRN interface *)
mkNet::usage = "net=mkNet[genes,rules] builds Boolean network Association from string gene names and Mathematica Boolean rule expressions (using &&,||,!). Regulators inferred automatically. Use readBnet or readBoolODE to load from file.";
readBnet::usage = "net=readBnet[file] reads Bnet-format file (lines: GeneName, BoolRule). Converts & -> &&, | -> ||. Returns network Association with keys genes, rules, regs.";
readBoolODE::usage = "net=readBoolODE[file] reads BoolODE tab-separated format (header: Gene<tab>Rule). Returns network Association with keys genes, rules, regs.";
evalNet::usage = "state=evalNet[net,state] synchronous Boolean update: evaluates all rules at 0/1 state list (indexed by net[\"genes\"]). Returns next 0/1 state list.";
steadyBool::usage = "fps=steadyBool[net] returns list of Boolean fixed points (states s where synchronous update maps s->s). Each fixed point is a 0/1 list.";
boolSTG::usage = "g=boolSTG[net] builds synchronous state transition graph as Mathematica Graph. Vertices are binary-integer state labels (FromDigits[state,2]).";
boolMultiRHS::usage = "{var,RHS}=boolMultiRHS[net] or boolMultiRHS[net,gam] converts Boolean network to multilinear (BoolCube) polynomial ODE: dxi/dt=gam*(fi^multi(x)-xi). Output plugs directly into ODE2RN[RHS,var] for CRN/NGM analysis. gam defaults to 1.";
boolHill2LRHS::usage = "{var2L,RHS2L}=boolHill2LRHS[net] or boolHill2LRHS[net,pars] converts Boolean network to two-layer Hill ODE (BoolODE style): mRNA xG and protein pG per gene. pars: Association overriding defaults m->1,lx->0.5,r->2,dp->0.8,k->0.5,n->2. Variables ordered {xA,pA,xB,pB,...}. Output plugs into ODE2RN[RHS2L,var2L].";

(*old version*)
Sta::usage="(old)Stab[mod,cfp,cn] analyzes stability at fixed point";
Stab::usage="(old)Sta[jac,X,Xv] numeric stability analysis";

(* Plot functions *)
crnGr::usage="(Plot)g=crnGr[RN] bipartite species-reaction graph. Species nodes (blue), reaction nodes R1,R2,... (orange). Directed edges: reactant->Ri->product.";
intGr::usage="(Plot)g=intGr[RN] interaction graph (species only). Edges: self-loop on each reactant, reactant->product, bidirectional between co-reactants. Opposing directed edges merged to undirected.";
qGr::usage="(Plot)g=qGr[RN] or g=qGr[RN,groups] quotient graph. groups={\"G1\"->{\"s1\",\"s2\"},...} collapses species into group nodes. Without groups identical to intGr.";

(* BCH flux/rate symbols -- declared public so they display without context prefix
   and remain distinct from Global`\[Beta] used in multi-strain epidemic models *)
\[Kappa]; \[Beta]; \[Gamma];

(* Per-siphon NGM functions (ClaudeFuckups.wl) *)
NGMsiph2::usage = "{mSi,JsigList,Klist,R0list} = NGMsiph2[RN,rts,var]  Like NGMsiph but extracts R_k via the factored characteristic polynomial of K_k. Trivial lam^k factors (zero eigenvalues) are excluded from Max. Nonzero-eigenvalue factors: degree <= 2 yield closed-form roots (quadratic formula, Simplify not FullSimplify to avoid Root rewriting); degree >= 3 yield R_k = {\"fac\", {poly,...}} — the irreducible polynomials, no Root objects. Helpers: eigsFac, rhoFac.";
NGMsiph::usage  = "{mSi,JsigList,Klist,ChpList,R0list} = NGMsiph[RN,rts,var]  Per-siphon Next Generation Matrices. Computes one K_k per minimal siphon Sigma_k via D(RHS[Sigma_k],Sigma_k)|_{Sigma_k=0}, split posM[J]/posM[-J]. Only Sigma_k is zeroed; all other variables remain free. Returns: symbolic K_k matrices, Factor[CharacteristicPolynomial[K_k,Lambda]] for each siphon, and spectral abscissae R_k as functions of free variables. Use R0siph to evaluate at a specific equilibrium.";
R0siph::usage   = "invM = R0siph[mSi,JsigList,var,EA]  Invasion numbers at boundary equilibria. EA = output of bdFp. Entry invM[[k,j]] = R_k(E_j) = spectral abscissa of K_k evaluated at the Sigma_j=0 face equilibrium EA[[j]]. Biologically: strain k invades the resident equilibrium on face j.";
invTable::usage = "tab = invTable[mSi,R0list,JsigList,var,EA]  Combined invasion table (ns x ns matrix). Diagonal = R_k(E0) from NGMsiph; off-diagonal [k,j] = R_k(E_j) from R0siph.";

allRatSub::usage="(Core)allRatSub[RHS,var,prF] enumerates all proper non-empty subsets Z of var as candidate fast variables for quasi-steady-state (QSS) elimination. For each Z, solves RHS_Z==0 for the variables in Z, keeps only rational solutions (no Sqrt/Root/complex), then computes the reduced drift b(y)=RHS_Y evaluated at z=z*(y) and keeps only those where b is also rational. Returns a list of Associations <|\"fast\"->Z,\"slow\"->Y,\"zStar\"->z*(y),\"b\"->b(y)|>, one per valid rational QSS. prF is optional (default Identity). Used in Bin.wl and MM.wl to list all Class A (algebraic) reductions with rational rates.";
pPrint::usage="pPrint[mat,mSi] displays mat as a Grid with border and with vertical+horizontal dividers after each siphon block. mSi is a list of variable-name lists (strings or symbols) as returned by NGMRN; var and prF must be defined in the calling context. Block boundaries are determined by the maximum position in var of each siphon's members.";
zQ::usage="zQ[v,s] returns True iff variable v evaluates to 0 under substitution s (list of rules). Used to identify which compartments are zero at a candidate equilibrium solution.";
stabPrint::usage="stabPrint[c] prints stability analysis at equilibrium substitution c. The label is derived automatically from the symbol name of c via SymbolName (HoldFirst). Requires jac (Jacobian matrix), par (parameter list), prF (pretty-printer) defined in calling context. Factors the characteristic polynomial via staPP; applies Hur3/Hur4 to any remaining higher-degree factor.";
LyapPV::usage="LyapPV[p, svar, zBlocks, blkData, cEp, RHS, var] computes the Perron-Volterra Lyapunov function L_p and its time derivative for the rank-one block competition model at the block-p boundary equilibrium. p: resident block index (1-based, matching zBlocks order). svar: susceptible variable symbol. zBlocks: list of compartment variable lists per block (e.g. {{i},{z1,z2}}). blkData: list of {wInc, wEnt, Vmat} per block, where phi_k = wInc_k . z_k (scalar incidence) and new infections enter as wEnt_k * phi_k * svar. cEp: equilibrium substitution rule list (from Solve). RHS: ODE right-hand side. var: variable list. Returns {Lp, dLpdt, mRList, cList, zpstar}: symbolic Lyapunov function, its derivative, per-block reproduction coefficients mR_k = wInc_k^T V_k^{-1} wEnt_k, canonical weight vectors c_k = wInc_k^T V_k^{-1}/mR_k, and resident equilibrium coordinates zpstar.";

Begin["`Private`"];

root=If[StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[$InputFileName],
    NotebookDirectory[]
];

(* Loaded later may override doubles loaded former *)
Get[FileNameJoin[{root, "Siphons.wl"}]];
Get[FileNameJoin[{root, "Boundary.wl"}]];
Get[FileNameJoin[{root, "Equilibria.wl"}]];
Get[FileNameJoin[{root, "Bifurcation.wl"}]];
Get[FileNameJoin[{root, "Conv.wl"}]];
Get[FileNameJoin[{root, "Utils.wl"}]];
Get[FileNameJoin[{root, "Core.wl"}]];
Get[FileNameJoin[{root, "CRNT.wl"}]];
Get[FileNameJoin[{root, "Visualize.wl"}]];
Get[FileNameJoin[{root, "Extra.wl"}]];
(*Get[FileNameJoin[{root, "bdCom.wl"}]];*)
Get[FileNameJoin[{root, "Cl.wl"}]];  (* Experimental versions for testing *)
Get[FileNameJoin[{root, "Core.wl"}]];  (* ODE2RN, extMat with rnRed,... *)
Get[FileNameJoin[{root, "NGMFD.wl"}]];  (* NGM for Frequency-Dependent incidence: NGMFD, NGMFDsub, NGMFDRN *)
Get[FileNameJoin[{root, "Organizations.wl"}]];  (* Chemical Organization Theory, ERCs: actRea,speClo,isClo,ercClo,allErc,isReac,isSSM,isSM,isOrg,allOrg,synQ,cmpERC *)
(* Buffering.wl is NOT loaded here.  Nothing in EpidCRN calls it, and buffering
   is studied in a handful of files; loading it made its self-test run in every
   session.  The files that need it load it themselves, right after <<EpidCRN`:
     Get[FileNameJoin[{$HomeDirectory, "EpidCRNmodels", "applications", "Buffering.wl"}]];
   The usage statements for buffS,isOutC,locIdx,isBuff,allBuff stay above, so the
   definitions still attach to the EpidCRN` symbols when that Get runs. *)
Get[FileNameJoin[{root, "CEP.wl"}]];  (* CEP lattice-walk decision algorithm *)
Get[FileNameJoin[{root, "ClaudeFuckups.wl"}]];  (* Per-siphon NGM: NGMsiph, R0siph, invTable *)
Get[FileNameJoin[{root, "Plot.wl"}]];
Get[FileNameJoin[{root, "BoolNet.wl"}]];
(*

*)

End[];
EndPackage[];




