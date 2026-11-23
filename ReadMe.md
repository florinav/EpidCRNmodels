# EpidCRN Package - ReadMe

## Important: rts = full rates (monomials * constants), cts (constants alone) basically never used

## Standard Notebook Startup (for RN-input notebooks)

```mathematica
ClearAll["Global`*"]
SetDirectory[NotebookDirectory[]]; SetOptions[$FrontEndSession, NotebookAutoSave -> True];
NotebookSave[];
AppendTo[$Path, FileNameJoin[{$HomeDirectory, "Dropbox", "EpidCRNmodels"}]];
<<EpidCRN`;
(* 0: Latex dictionary *)
prF[ex_] := ex /. {
  x1->Subscript[x,1], x2->Subscript[x,2], x3->Subscript[x,3],
  S1->Subscript[S,1], S2->Subscript[S,2],
  B1->Subscript[B,1], B2->Subscript[B,2],
  et -> η, de->δ, be->β, la -> λ, La->Λ, mu -> μ,
  mu1 -> Subscript[μ, 1], mu2 -> Subscript[μ, 2],
  ga1 -> Subscript[γ, 1], ga2 -> Subscript[γ, 2],
  a1 -> Subscript[α, 1], a2 -> Subscript[α, 2], al -> α,
  be1 -> Subscript[β, 1], be2 -> Subscript[β, 2], be3 -> Subscript[β, 3],
  ep1 -> Subscript[ε, 1], ep2 -> Subscript[ε, 2],
  om -> ω};

(* Then define RN and var *)
RN = {...};
var = {...};
```

**Key points:**
- Always reload package with `<<EpidCRN\`` (no caching issues)
- prF is model-specific LaTeX dictionary for pretty printing
- Modify prF rules for each model's variables

---

## CRITICAL RULES FOR AI AGENTS

1. **⚠️ NEVER MODIFY ANY PACKAGE WITHOUT EXPLAINING WHAT YOU DO (EXCEPT ClaudeFuckups.wl) ⚠️**
   - All subpackages (epi.wl, Core.wl, Equilibria.wl, etc.) contain trusted functions
   - Before modifying any package file, you must explain what you plan to change and why
   - ClaudeFuckups.wl is the ONLY exception - you can freely modify it for experimental versions
   - If you need to create a new version of a function, create it in ClaudeFuckups.wl first

2. **⚠️ ALWAYS TEST ClaudeFuckups.wl COMPILES BEFORE PRESENTING ⚠️**
   - Before offering any code in ClaudeFuckups.wl to the user, verify the package loads without errors
   - Test by loading: `<<EpidCRN\`` and checking for compilation errors
   - Only modify ClaudeFuckups.wl (not other subpackages) until user approves promotion

## Coding Conventions

**Variable Naming:**
- Use global names consistently across notebooks: RN, var, alp, bet, gam, jac, eig, par
- Never use numbered variants (RN1, RN2, etc.) - not necessary since we don't run multiple networks simultaneously
- Always use 3-letter abbreviations when printing: "RN=", "gam=", "eig="
- Always print variable name before value: `Print["RN=",RN," gam=",gam," eig=",eig]`
- Never use Print statements without variable names

**Function Outputs:**
- All important functions must output lists (not associations)
- Usage statements appear only in loader (EpidCRN.wl) and only once per function

**Testing:**
- Each new function/script must have uncommented test at end of subpackage file (before End[])
- Tests run automatically when package loads
- Follow pattern in Core.wl for ODE2RN test
- Provide self-contained tests for each new submodule

## Package Structure

**Subpackages:**
- **epi.wl** - Foundational functions: minSiph, ODE2RN, NGM, etc. **NEVER modify without permission.**
- **Siphons.wl** - Siphon analysis functions
- **CRNT.wl** - Chemical Reaction Network Theory
- **Utils.wl** - Utility functions
- **Boundary.wl** - Boundary equilibria analysis: Hur3M, Hur4M, Hur5M
- **Equilibria.wl** - Equilibria classification: clsEq, clsEqb, sipLat
- **ClaudeFuckups.wl** - Experimental/untested versions of functions. Use for testing new versions before promoting to main subpackages.
- **Core.wl** - (Legacy, separate from epi.wl)

**Known Issue - minSiph Usage Display:**
- When clicking on `minSiph` in the ?EpidCRN`* list, the full code definition is displayed instead of just the usage statement
- This occurs even though the usage statement is correctly defined in EpidCRN.wl (line 9)
- All other functions display only their usage statement as expected
- This is a package structure issue that needs investigation
- **Workaround:** Use `?minSiph` directly to see the usage statement
- **TODO:** Fix package structure to ensure minSiph displays usage like other functions

## Key Functions

**Boundary Analysis:**
- `bdFp[RHS, var, mSi]` - Computes boundary equilibria on siphon facets. For each minimal siphon in mSi: sets siphon variables to 0, solves the resulting system, filters out DFE and negative solutions using onlyNN. Returns list where `bdFp[[j]]` is either positive solutions (list of rule lists) or `{equations, variables}` for siphon j when no positive solutions exist.

**Bifurcation Analysis:**
- `scan[RHS, var, par, persRule, plotInd, mSi, gridRes, steadyTol, stabTol, chopTol, R01, R02, R12, R21, rangeExtension]` - Scans parameter space and classifies equilibria.
  - **Reproduction numbers**: R01 = R0 for strain 1 at DFE, R02 = R0 for strain 2 at DFE, R12 = invasion number for strain 1 into endemic E2, R21 = invasion number for strain 2 into endemic E1
  - **Classification logic**:
    - DFE: R01 < 1 and R02 < 1
    - E1: R01 > 1 and R02 < 1 (strain 1 only)
    - E2: R02 > 1 and R01 < 1 (strain 2 only)
    - Coexistence region: R01 > 1 and R02 > 1
      - If R12 < 1 and R21 < 1: competitive exclusion (winner determined by R01 vs R02)
      - If R12 > 1 and R21 < 1: strain 1 wins (E1)
      - If R12 < 1 and R21 > 1: strain 2 wins (E2)
      - If R12 > 1 and R21 > 1: potential coexistence - solve for equilibrium and check eigenvalues
  - **Performance optimization**: Stability (eigenvalue) analysis performed only in coexistence region (R01>1, R02>1, R12>1, R21>1) to reduce computation time. Other regions classified by reproduction numbers only.
  - **rangeExtension**: Parameter controlling plot range extension beyond intersection point (default 2.5 = 250%)

## Recent Work: Two-Reduction Regime Analysis (December 2024)

### Motivation
Integrated Feliu-Lax-Walcher-Wiuf (FLWW) asymptotic reduction theory with Saez-Feliu-Wiuf (SFW) algebraic elimination theory to validate QSS reductions via singular perturbation analysis.

### Key Theoretical Contributions

**FLWW Paper ("Quasi-steady state and singular perturbation reduction for reaction networks"):**
- Provides necessary/sufficient conditions for Tikhonov-Fenichel reduction validity
- **Blanket condition (i):** B₀(x) invertible ⟺ spanning trees exist in G_U (Lemma 11)
- **Blanket condition (ii):** Critical manifold stationary ⟺ cycle conditions (Lemma 12)
- **Corollary 13:** Setting inflow rates (∗→U reactions) to zero guarantees TFPV
- **Key insight:** Same eliminable set can be valid in different parameter regimes

**Comparison: SFW vs FLWW**
| Aspect | SFW (2017) | FLWW |
|--------|------------|------|
| Purpose | Identify eliminable sets | Validate QSS asymptotically |
| Conditions | Noninteracting, U-linear, spanning tree | + Blanket conditions (i), (ii) |
| Output | Reduced network structure | Parameter regimes where valid |
| Validation | Algebraic (steady state) | Asymptotic (Tikhonov-Fenichel) |

### Implementation: John.nb

**File:** `C:\Users\flori\Dropbox\EpidCRNmodels\multiStrain\John.nb`

Comprehensive single-file analysis of John's two-strain SIR model demonstrating:

**Model:** 5D system {s, i1, i2, r1, r2} with:
- Cross-immunity (σ₁, σ₂ parameters)
- Demographic turnover (Λ, μ rates)
- Strain-specific infection/recovery (β₁, β₂, γ₁, γ₂)
- Waning immunity (θ₁, θ₂)

**Maximal Eliminable Sets Found:**
1. {i1, i2} → reduces to 3D system in {s, r1, r2}
2. {r2, r1, s} → reduces to 2D system in {i1, i2}
3. {i1, r1} → reduces to 3D system

**Two Reduction Regimes Illustrated:**

*Regime 1: Fast Infection Dynamics*
- **Parameters:** β₁=2.0, β₂=1.5, γ₁=1.0, γ₂=0.8 (fast) vs Λ=0.05, μ=0.002, θ₁=0.05, θ₂=0.04 (slow)
- **Elimination:** {i1, i2} via QSS assumption (i₁' = 0, i₂' = 0)
- **Reduced system:** 3D in {s, r1, r2}
- **Reconstruction:** i₁ = φ₁(s, r1, r2), i₂ = φ₂(s, r1, r2)
- **Biological interpretation:** Epidemic on fast timescale → track slow demographic changes

*Regime 2: Fast Demographic Turnover*
- **Parameters:** Λ=2.0, μ=0.8, θ₁=1.0, θ₂=0.8 (fast) vs β₁=0.2, β₂=0.15, γ₁=0.1, γ₂=0.08 (slow)
- **Elimination:** {s, r1, r2} via QSS assumption (s' = 0, r₁' = 0, r₂' = 0)
- **Reduced system:** 2D in {i1, i2}
- **Reconstruction:** s = ψ₁(i1, i2), r₁ = ψ₂(i1, i2), r₂ = ψ₃(i1, i2)
- **Biological interpretation:** Demography on fast timescale → track slow epidemic spread

**File Structure:**
```mathematica
(* Section 1: MODEL DEFINITION *)
RHS = {...};
var = {s, i1, i2, r1, r2};
{RN, rts, spe, alp, bet, gam} = ODE2RN[RHS, var];

(* Section 2: BOUNDARY ANALYSIS *)
{RHS1, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, ng} = bdAn[RN, rts, var];
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, K, R0A, infVars, EA} = bdCo[RN, rts, var];

(* Section 3: REDUCTION ANALYSIS *)
allMaxElim = maxElim[RN, rts];

(* Section 4: GRAPHICAL ILLUSTRATION - automatic if ≥2 sets *)
If[Length[allMaxElim] >= 2,
  U1 = allMaxElim[[1]];
  U2 = allMaxElim[[2]];
  (* Regime 1 comparison plots *)
  (* Regime 2 comparison plots *)
];
```

**Design Philosophy:**
- ONE file per model (named after the model)
- User changes RHS/var at top → everything updates automatically
- Automatically finds all maximal eliminable sets
- If ≥2 sets exist, automatically illustrates both regimes
- Verifies FSW conditions for each set
- Creates comparison plots: Blue = Full system, Red = Reduced system

**Graphical Output:**
- Slow variables: Full vs Reduced trajectories overlay
- Fast variables: Reconstructed via QSS vs Full system
- Validates asymptotic agreement in appropriate parameter regimes

### Related Files

**Analysis summaries:**
- `TwoReductionsSummary.md` - Detailed comparison of both regimes
- `illustrateReductions.nb` - Generic tool (works for any CRN)
- `illustrateTwoReductions.nb` - Earlier working prototype
- `compareReductionsFinal.nb` - Development version

**Recommendation for New Models:**
1. Copy John.nb to NewModel.nb
2. Change RHS and var definitions
3. Run → automatic two-regime analysis if multiple eliminable sets exist

---

## Current Plan

Integration of FLWW asymptotic theory with SFW elimination theory complete for John's model. Framework ready for application to other multi-strain/multi-pathway epidemic models.

---

\subsection{Equilibria Classification Module (January 2025)}

\subsubsection{Motivation}
After \texttt{bdAn} computes minimal siphons, DFE, reproduction functions and $R_0$, we need systematic classification of equilibria into: rational, non-rational with explicit RUR, and time-outs.

\subsubsection{Algorithm Overview}
\textbf{Input:} Output from \texttt{bdAn}

\begin{verbatim}
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars}
\end{verbatim}

\textbf{Steps:}
\begin{enumerate}
\item Construct lattice of minimal siphons, order by decreasing length
\item Remove cDFE rule (DFE already found by \texttt{bdAn}, pruned from future results)
\item Process rules by length level (longest first)
\item For each siphon: attempt \texttt{Solve} with 300s timeout and 1GB memory limit
\item Classify solutions: rational vs non-rational (detecting \texttt{Sqrt}, \texttt{Root}, \texttt{AlgebraicNumber})
\item Filter DFE from stored solutions, but include in counts
\item Count RUR as equations (not solutions): even number of solutions counted as half (likely $\pm$ pairs)
\end{enumerate}

\textbf{Output:} Per-siphon list: \texttt{\{siphon, nRat, nRur, \{ratSols, rurSols\}\}} or \texttt{"limit"} if Sp/Time limit exceeded

\subsubsection{Implementation: Equilibria.wl}

\textbf{Location:} \texttt{C:\textbackslash Users\textbackslash flori\textbackslash Dropbox\textbackslash EpidCRNmodels\textbackslash Equilibria.wl}

\textbf{Function Hierarchy:}

\textbf{Level 1:} \texttt{sipLat[mSi, vars]} -- Builds siphon lattice
\begin{itemize}
\item Input: minimal siphons from \texttt{bdAn}, variable list
\item Output: list of siphon rules ordered by decreasing length
\item Removes cDFE rule (union of all minimal siphons)
\end{itemize}

\textbf{Level 2:} \texttt{slvLvl[eqs, vars, rules, solved, dfe]} -- Solves one length level
\begin{itemize}
\item For each rule: \texttt{Solve[eqs /. rule, vars]} with 300s timeout
\item Filters: DFE, duplicates from solved list
\item Separates rational (polynomial ratios) from non-rational (contains \texttt{Sqrt}, \texttt{Root})
\item Returns: \texttt{\{ratSols, nonRatSols, timeOuts\}}
\end{itemize}

\textbf{Level 3:} \texttt{slvAll[eqs, vars, sipRules, dfe]} -- Iterates all levels
\begin{itemize}
\item Processes longest siphons first
\item Maintains solved accumulator to avoid duplicates
\item Returns: \texttt{\{ratSols, rurSols, timeOuts\}}
\end{itemize}

\textbf{Level 4:} \texttt{clsEqb[RHS, var, mSi, E0, prF]} -- Main interface
\begin{itemize}
\item Creates equilibrium equations: \texttt{Thread[RHS == 0]}
\item Solves for each siphon separately
\item \texttt{prF}: Optional print flag (default \texttt{False}). If \texttt{True}, prints results per siphon with rational solutions tracked to avoid reprinting
\item Returns per-siphon classification: list of \texttt{\{siphon, nRat, nRur, \{ratSols, rurSols\}\}}
\item \textbf{Counting logic:}
  \begin{itemize}
  \item nRat: Total rational solutions (including DFE if present)
  \item nRur: Number of distinct RUR equations (not solutions). Heuristic: if even number of solutions, likely $\pm$ pairs, count as half
  \item Stored solutions: Exclude DFE from storage
  \end{itemize}
\item \textbf{RUR handling:} If solutions are extremely long ($>$5000 chars), treats as RUR and solves reduced system keeping $x_1$ (substrate) as parameter
\end{itemize}

\subsubsection{Usage Example}

\begin{verbatim}
<<EpidCRN`;

(* Define model *)
RHS = {La - be1*i1*s - be2*i2*s - mu*s,
       be1*i1*s - (ga1 + mu)*i1,
       be2*i2*s - (ga2 + mu)*i2};
var = {s, i1, i2};

(* Convert to reaction network *)
{RN, rts, spe, alp, bet, gam, rnRed} = ODE2RN[RHS, var];

(* Boundary analysis *)
{RHS, var, par, cp, mSi, Jx, Jy, cDFE, E0, F, V, K, R0A, infVars} =
  bdAn[RN, rts, var];

(* Classify equilibria per siphon *)
perSiphon = clsEqb[RHS, var, mSi, E0];

(* With printing enabled *)
perSiphon = clsEqb[RHS, var, mSi, E0, True];

(* Each entry: {siphon, nRat, nRur, {ratSols, rurSols}} *)
(* If prF=False, manually print results *)
Do[
  {sip, nRat, nRur, sols} = perSiphon[[i]];
  Print["Siphon ", i, ": ", sip, " rat=", nRat, " rur=", nRur];
, {i, Length[perSiphon]}];
\end{verbatim}

\subsubsection{Test Results}

\textbf{Two-strain SIS model:}
\begin{itemize}
\item Variables: \texttt{\{s, i1, i2\}}
\item Minimal siphons: \texttt{\{\{i1\}, \{i2\}\}}
\item Per-siphon results:
  \begin{itemize}
  \item Siphon 1 \texttt{\{i2\}}: \texttt{rat=2, rur=0} (1 DFE + 1 endemic with strain 1 only)
  \item Siphon 2 \texttt{\{i1\}}: \texttt{rat=2, rur=0} (1 DFE + 1 endemic with strain 2 only)
  \end{itemize}
\end{itemize}

\textbf{OSN model (Online Social Network with botnets):}
\begin{itemize}
\item Variables: \texttt{\{x1, x2, x3, B1, S1, B2, S2, R\}} (8 variables)
\item Minimal siphons: \texttt{\{\{x2\}, \{B1, S1\}, \{B2, S2\}\}}
\item Siphon lattice: 6 siphons (excluding cDFE)
\item Per-siphon results:
  \begin{itemize}
  \item Siphon 1 \texttt{\{B1,S1,B2,S2\}}: \texttt{rat=3, rur=0} (1 DFE + 2 endemic, stores 2)
  \item Siphons 2,3 (len=3): \texttt{rat=1, rur=0} (DFE only)
  \item Siphons 4,5 \texttt{\{B2,S2\}}, \texttt{\{B1,S1\}}: \texttt{rat=3, rur=1} (1 DFE + 2 endemic + 1 RUR equation with 2 roots)
  \item Siphon 6 \texttt{\{x2\}}: \texttt{rat=1, rur=0} (DFE only)
  \end{itemize}
\end{itemize}

\textbf{Key Implementation Details:}

\textbf{1. Detection of non-rational solutions:}
\begin{verbatim}
isRational[s_] := FreeQ[s, Root | RootOf | AlgebraicNumber] &&
  FreeQ[s, _Power?((!IntegerQ[#[[2]]] && Head[#[[2]]] =!= Symbol)&)]
\end{verbatim}
This correctly identifies \texttt{Sqrt[x]} as \texttt{Power[x, 1/2]} (non-integer exponent) and classifies it as non-rational.

\textbf{2. RUR equation counting:}
RUR count represents number of distinct algebraic equations, not individual roots. When \texttt{Solve} returns multiple solutions differing by $\pm$ signs (e.g., $x = a \pm \sqrt{b}$), these are roots of a single equation. Heuristic: if even number of RUR solutions, count as half (likely $\pm$ pairs). For epidemic models, only positive roots are biologically valid.

\textbf{3. RUR safeguard - keeping substrate as parameter:}
CRITICAL: If a solution passes \texttt{isRational} but has extremely long expressions ($>$5000 characters), it is likely irrational. For RUR cases, the proper approach is:
\begin{itemize}
\item Keep $x_1$ (substrate, first variable in \texttt{var}) as a free parameter
\item Delete the first equation (equation for $x_1$)
\item Solve the reduced system for remaining variables
\item Result: RUR solutions expressed in terms of $x_1$
\end{itemize}
This prevents spurious "rational" solutions with extremely long denominators that are actually irrational. The code automatically tries this reduced approach when RUR is detected and uses it if it produces cleaner output.

\subsubsection{Bug Fix: Complementarity Filtering (January 2025)}

\textbf{Problem:} Siphon 1 was not printing its boundary equilibria, and Siphons 4-5 were reprinting solutions already shown in Siphon 1.

\textbf{Root Cause:} The \texttt{hasOtherSiphonZero} helper function builds \texttt{antiSiph} -- the union of OTHER minimal siphons. A solution on one siphon (e.g., $\{B_2, S_2\} = 0$) should have variables from other minimal siphons non-zero (e.g., $\{B_1, S_1, x_2\} \neq 0$). The bug occurred because:
\begin{enumerate}
\item Initially used \texttt{Complement[Rest[vars], siphonVars]} which included ALL non-siphon variables (like R which isn't in any minimal siphon)
\item This incorrectly filtered out valid boundary equilibria with $R=0$
\item Later fix used \texttt{SubsetQ[siphonVars, mSi\_element]} to find other siphons, but \texttt{SubsetQ} failed because mSi elements might be strings while siphonVars contains symbols
\end{enumerate}

\textbf{Solution:} Convert mSi elements to symbols before SubsetQ test:
\begin{verbatim}
hasOtherSiphonZero[s_] := Module[{otherSiphons, antiSiph, mSiSymbols},
  (* Convert mSi elements to symbols first *)
  mSiSymbols = Map[
    Map[If[StringQ[#], ToExpression[#], #]&, #]&,
    mSi
  ];
  (* Find minimal siphons NOT contained in current siphon *)
  otherSiphons = Select[mSiSymbols, !SubsetQ[siphonVars, #]&];
  (* Union of other minimal siphons *)
  antiSiph = Flatten[otherSiphons, 1];
  antiSiph = DeleteDuplicates[antiSiph];
  (* Check if any are zero *)
  !nonZat[s, antiSiph]
];
\end{verbatim}

\textbf{Example (OSN model):}
\begin{itemize}
\item \texttt{mSi = \{\{x2\}, \{B1, S1\}, \{B2, S2\}\}}
\item For Siphon 1 with \texttt{siphonVars = \{B1, S1, B2, S2\}}:
  \begin{itemize}
  \item Before fix: \texttt{antiSiph = \{x2, x3, R\}} (WRONG - includes R which isn't in any minimal siphon)
  \item After fix: \texttt{antiSiph = \{x2\}} (CORRECT - only the other minimal siphon)
  \item Solutions with $x_2 > 0$ now correctly pass filter and print
  \end{itemize}
\end{itemize}

\textbf{Related fixes:}
\begin{itemize}
\item \texttt{prF} parameter added to clsEqb signature: formatting function applied with \texttt{//prF}
\item Return structure changed from \texttt{\{siphonVars, nRat, nRur, sols\}} to \texttt{\{nRat, nRur, sols\}}
\item \texttt{printedRat} list tracks already-printed solutions to prevent duplication
\item RUR filtering uses \texttt{Join[solRed[[i]], rule]} to complete solution before checking complementarity
\end{itemize}

\subsubsection{Integration}
\begin{enumerate}
\item File added: \texttt{Equilibria.wl} in root directory
\item Loader updated: \texttt{Get[FileNameJoin[\{root, "Equilibria.wl"\}]]} in \texttt{EpidCRN.wl} (line 219)
\item Usage statements added to \texttt{EpidCRN.wl} (lines 69-72)
\item Test files: \texttt{testEquilibria2.wl}, \texttt{testOSNfull.wl}, \texttt{testSqrt.wl}
\end{enumerate}

\subsection{clsEq: Endemic Equilibrium Solver with Timeout (2026-01-14)}

\subsubsection{Purpose}
New function \texttt{clsEq[RHS, var, mSi, E0, prF, tim]} specifically tackles endemic equilibrium (no siphon variables zero) with configurable timeout and fallback to RUR.

\subsubsection{Algorithm}
\begin{enumerate}
\item First run \texttt{clsEqb[RHS, var, mSi, E0, prF]} and print all siphon results
\item Then attempt unrestricted \texttt{Solve[eqs, vars]} with timeout \texttt{tim} for endemic equilibrium
\item If timeout: Attempt RUR approach (drop first equation, solve reduced system, compute polynomial in x1)
\item If RUR also times out: Print "TIME-OUT" and return empty list
\item Filter: Remove DFE, solutions with any siphon zero, all-zero solutions
\end{enumerate}

\subsubsection{Output}
\begin{itemize}
\item Success (Solve): Prints \texttt{nSol=N} followed by \texttt{E1=Short[...]} through \texttt{EN=Short[...]}
\item Success (RUR): Prints \texttt{nRUR=N} followed by \texttt{ratS=Short[...]} and \texttt{RUR=Short[...]} for each polynomial
\item Failure: Prints "TIME-OUT"
\item Returns: List of solution rules (empty if timeout or no solutions)
\end{itemize}

\subsubsection{Usage}
\texttt{endemic = clsEq[RHS, var, mSi, E0, prF, 60]}
Test file: \texttt{testClsEq.wl} with 60-second timeout.

## Reserved Variable Names

**IMPORTANT**: The following variable names are reserved by the package and should NOT be used in model RHS definitions:

- **E0**: Equilibrium values (numeric/symbolic values of variables at equilibrium)
- **cE0**: Equilibrium rule (list of rules var -> value at DFE, returned by minSiph)
- **cDFE**: Disease-Free Equilibrium rule (siphon variables -> 0)
- **V, F, K**: NGM matrices (Transition, New infections, Next Generation Matrix)
- **Jx, Jy, Jxy, Jyx, Kd**: Jacobian matrices and blocks
- **RN, rts, var, par**: Reaction network, rates, variables, parameters
- **alp, bet, gam, spe**: Stoichiometric matrices (alpha, beta, gamma, species)
- **mSi, inf, infc, mod, cp, sol**: Minimal siphons, infection vars, complements, model, constraints, solutions
- **I, E**: Reserved to avoid conflict with built-in Mathematica symbols (imaginary unit, exponential)

Use the \texttt{chk[RHS]} function immediately after defining RHS to detect conflicts:
```mathematica
RHS = {...};
chk[RHS];  (* Warns if any variables conflict with reserved names *)
```

**Function Signatures**: All package functions MUST use reserved names in their input/output signatures (e.g., return `V` not `VV`). Internal implementations may use alternative names (e.g., `VV` internally) but must map back to reserved names before returning.

---

## Future Plans

**Python Port**: Core.wl (containing the fundamental CRN functions: ODE2RN, extMat, minSiph, NGM, etc.) is planned to be implemented in Python. The Mathematica version serves as the reference implementation.

---

## Four Laws Plan (2026-03-07, updated 2026-03-10)

**Reference:** DSI.tex (Mathematical_Epidemiology), Theorem t:end; Fall/Bonzi papers.

**Files:** `C:\Users\flori\Dropbox\Cholera\Bonzi.wl`, `C:\Users\flori\Dropbox\Cholera\SEI2RB.nb`

### What was done
`Bonzi.wl` converted to a proper package `Bonzi`` exporting one function:
```mathematica
{y0, mR, R0, Dw, ye, k, xe} = fourLaws[A, vbe, ba, AS, La]
```

### The four laws implemented in fourLaws
| Output | Formula | Law |
|--------|---------|-----|
| `y0` | `(-AS)^{-1} La` | DFE |
| `mR` | `vbe (-A)^{-1} ba` | replacement vector |
| `R0` | `y0 . mR` | Law 1 |
| `Dw` | `(-A)^{-1} ba` | dwell times |
| `ye` | `(k*Diag[mR] - AS)^{-1} La` | Law 3 |
| `k`  | solve `ye . mR = 1` | Law 4 normalization |
| `xe` | `k * Dw` | Law 2 |

**Inputs:** `A` (n×n MH, infected), `vbe` (m×n WAIFW), `ba` (n-vec entry probs), `AS` (m×m MH, susceptible), `La` (m-vec recruitment).

### Testing plan
1. **Test 1 (done in Bonzi.wl):** n=3 staged progression, single S. Run `<<Bonzi.wl` and verify symbolic + numeric output.
2. **Test 2 (next):** SEI2RB model in `SEI2RB.nb`. After NGMRN runs, extract A, vbe, ba, AS, La and call `fourLaws`. Verify R0 matches `R0A` from NGMRN. Add cells to SEI2RB.nb for this.

### Implementation (2026-03-10)
`FLaws` and `extABa` implemented in **`Cl.wl`** (loaded as experimental subpackage by EpidCRN.wl).
Usage statements already in `EpidCRN.wl` lines 83-88.

#### extABa extraction algorithm
Given `RHS, var, yvar (susceptible), xvar (infected)`:

**From y-equations** (for each susceptible class i):
- `vbe[i,j] = -Coefficient(yRHS_i, y_i * x_j)` — WAIFW from bilinear incidence loss
- `phi_i = yRHS_i + y_i*(vbe_i . xvar)` — remove bilinear term
- `La_i = phi_i |_{y=x=0}` — constant recruitment
- `AS[i,j] = Coefficient(phi_i - La_i, y_j)` — transition matrix

**From x-equations** (for each infected class k):
- `A[k,l] = Coefficient(xRHS_k, x_l)` — infected transition matrix
- `rem_k = xRHS_k - A_k . xvar` — bilinear remainder = `ba[k] * y^t vbe x`
- `ba[k] = Coefficient(rem_k, y_{i0} * x_{j0}) / vbe[i0,j0]` — for any nonzero vbe entry

**Checks:** notABa message if y-equation residual nonzero; novbe if vbe all zero.

#### General P case: extABP (implemented 2026-03-10)
`extABP[RHS, var, xvar]` → `{La, AS, B, A, P, yvar}` in Cl.wl, usage in EpidCRN.wl.

Algorithm: `bil[k,i,j] = Coeff(rem_k, y_i*x_j) = P[k,i]*B[i,j]`.
Recover `P[k,i]` for each (k,i): `ratios = bil[k,i,j]/B[i,j]` for all nonzero `B[i,j]`.
If ratios are not all equal → `extABP::badP` message (model not DSI bilinear).

From rem.tex: `rank(F(y))=1 ∀y` iff `rank(P)=1` or `rank(B)=1` (two independent cases).
- `rank(P)=1` → P-rank-one model (P = ba*ones^t)
- `rank(B)=1` → B-rank-one = ABa model (B = ba*vbe^t)
Check: `MatrixRank[P]`, `MatrixRank[B]` after calling `extABP`.

### Testing plan
1. **Test 1 (coded in Bonzi.wl):** n=3 staged progression, single S. Run `Bonzi.wl`, verify extABa recovers A, vbe, ba, AS, La matching hardcoded values, and FLaws gives same R0.
2. **Test 2 (next):** SEI2RB (`SEI2RB.nb`). Call `extABa[RHS, var, {S}, {e,iM,iA,B}]` and verify R0 matches NGMRN output.

---

## NCB Oscillatory Core Analysis (February 2025)

**File:** `SIRWS.wl` -- NCB (Negative Coefficient Bundle) approach for detecting Hopf/oscillatory cores in CRNs, applied to the SIRWS (Susceptible-Infected-Recovered-Waning) model with 4 species and 12 reactions.

**Method (following Bi.R.Ne/ChatGPT formulation):**
1. Build symbolic reactivity matrix RM where RM[j,k] = r[j,k] if species k is reactant of reaction j
2. Form symbolic Jacobian G = gam . RM (stoichiometric matrix times reactivity)
3. Expand characteristic polynomial coefficients into summands; each monomial = a k-child-selection (bijection from species subset to reactions subset)
4. Negative summands identify destabilizing feedbacks. A k-CS is an unstable positive feedback (UPF) if sign(det(gam[kappa])) = (-1)^(k-1)
5. Build Hasse diagram by r-set inclusion, apply transitive reduction
6. Oscillatory core = UPF whose Hasse descendants include a Hurwitz-stable CS-matrix (checked via Hur3/Hur4)

**Reuses (read-only):** `allT`, `isN`, `Hur3`, `Hur4` from Core.wl; `ODE2RN`, `NGMRN` from epi.wl

**Expected result:** Oscillatory cores involving the infection feedback loop {S, i} and the immunity waning cycle {R, W}.

---

## Per-Siphon Reproduction Functions — Plan (2026-05-08)

### Problem

`NGMRN` (Core.wl:733) computes a **single global NGM** using
`infVarsStr = Union[Flatten[mSi]]` — the union of all minimal siphon variables.
For models with overlapping siphons (e.g., BG two-strain vector-borne model,
minimal siphons Σ1={I1,I21,V1,V12}, Σ2={I2,I12,V2,V12} sharing V12),
this produces a single spectral abscissa that does not correspond to any
biologically meaningful reproduction number for either strain.

The correct reproduction function for siphon Σ_k is:

    R_k = ρ( K_{Σ_k} ),   K_{Σ_k} = F_{Σ_k} · V_{Σ_k}^{-1}

where F_{Σ_k}, V_{Σ_k} are the positive/non-positive parts of Jx restricted to
the infection variables of Σ_k, evaluated at the DFE.

Existing correct per-siphon logic: `cepCovers` (CEP.wl:56–82) — reuse its
approach, do NOT copy or modify it.

### Rule

**Do NOT modify** `NGMRN`, `NGM`, `bdFp`, `cepCovers`, or any existing function.
Prototype new functions in `ClaudeFuckups.wl` first. Promote to `Core.wl` only
after user approval.

---

### New Function 1: `NGMsiph`

**Signature:** `NGMsiph[RN_, rts_, var_]`

Parallel to `NGMRN`. Returns per-siphon NGM matrices and reproduction functions
evaluated at the DFE.

**Algorithm:**

```
1. Run standard pipeline to obtain: mSi, E0 (DFE), Jx (at E0)
2. For each siphon Σ_k in mSi:
   a. infVars_k = infection variables of Σ_k
      (same detection as cepCovers: variables whose Jx diagonal can be positive at E0)
   b. Extract Jx_k = submatrix of Jx restricted to rows/cols indexed by infVars_k
   c. Split Jx_k = F_k - V_k  via posM  (same as Core.wl:NGM)
   d. K_k = F_k . Inverse[V_k]
   e. R_k = Max[Re[Eigenvalues[K_k]]]
3. Return {mSi, Klist, R0list}
```

**Output** (tuple):
- `mSi`     — list of minimal siphons (strings)
- `Klist`   — list of per-siphon NGM matrices {K_1, K_2, ...} (symbolic)
- `R0list`  — list of per-siphon spectral abscissae {R_1, R_2, ...} (symbolic)

---

### New Function 2: `R0siph`

**Signature:** `R0siph[R0list_, var_, BElist_]`

Given pre-computed `R0list` (from `NGMsiph`) and boundary equilibria `BElist`
(from `bdFp`), evaluate R_k at E0 and at each boundary equilibrium.

**Algorithm:**

```
1. R_k(E0)  = R0list[[k]] at DFE                   (basic reproduction number)
2. For each boundary eq E_j in BElist (face Σ_j = 0):
   R_k(E_j) = R0list[[k]] /. Thread[var -> E_j],  k ≠ j  (invasion number)
3. Return Table[{k, j, R_k(E_j)}]  (* invasion matrix *)
```

**Output:** matrix/table, one entry per (invader k, resident face j) pair.

---

### New Function 3: `invTable`

**Signature:** `invTable[R0list_, var_, BElist_]`

Thin wrapper: calls `R0siph`, formats output as a Mathematica Table ready for
`Grid` display. No computation, pure formatting.

---

### Implementation Notes

- **posM split**: reuse `posM` from Core.wl (internal helper).
- **infVars detection**: same logic as `cepCovers` — diagonal of Jx restricted to
  Σ_k variables, check sign pattern at E0.
- **Symbolic safety**: all operations symbolic; `FullSimplify` only on final R_k.
- **No side effects**: functions are pure — no `Print`, no global assignments.
- **Usage declarations**: add `NGMsiph::usage`, `R0siph::usage`, `invTable::usage`
  in `EpidCRN.wl` public section (before `Begin["\`Private\`"]`). HARD RULE.

### Files to Touch (after approval only)

| File | Action |
|------|---------|
| `EpidCRNmodels/ClaudeFuckups.wl` | Prototype all three functions here first |
| `EpidCRNmodels/Core.wl` | Promote after approval, add after NGMRN block (~line 782) |
| `EpidCRNmodels/EpidCRN.wl` | Add `::usage` declarations for all three |

### Validation Target

BG model (`Burcu/BG.wl`), 14 variables, two overlapping siphons:
- Σ1 = {I1, I21, V1, V12}, Σ2 = {I2, I12, V2, V12}
- Expected: R_1 = R01 from mod1.tex, R_2 = R02
- Expected: R_2(E1) = invasion number of strain 2 into strain-1 boundary eq
- Cross-check: `NGMRN` output should differ (it conflates Σ1∪Σ2)

---

## New module: Organizations.wl (Chemical Organization Theory, ERCs)

`EpidCRNmodels/Organizations.wl` — Chemical Organization Theory (closed,
self-maintaining species sets) and Elementary Reaction Closures (ERCs),
plus synergy/complementarity between ERCs, for SubFF's Section 1 rewrite
(`CocktailArXiv-updated/COT.tex`).

**Functions (all list output, no Associations):**
`actRea[RN,spe]`, `speClo[RN,spe]`, `isClo[RN,spe]`, `ercClo[RN,rxn]`,
`allErc[RN]`, `isReac[RN,spe]`, `isSM[RN,spe]`, `isSSM[RN,spe]`,
`isOrg[RN,spe]`, `allOrg[RN,Method->...]`, `synQ[RN,ercList]`,
`cmpERC[RN,e1,e2]`.

**`Buffering.wl` additions** (do not touch `buffS`): `isOutC[RN,m,r]`,
`locIdx[RN,m,r]`, `isBuff[RN,m,r]`, `allBuff[RN,Method->...]` — general
forms exposing what `buffS` kept internal (arbitrary output-complete `r`,
not just the forced-minimal one).

Reuses `comp2Asso`/`asoRea`/`extMat` from Core.wl; RN representation
unchanged (list of `lhs->rhs` rules over string species). Loaded via
`EpidCRN.wl` between `NGMFD.wl` and `Buffering.wl`.

**Verified examples:** `EpidCRNmodels/OrganizationsExamples.wl` — bare
two-strain network (organizations vs. siphons vs. buffering all differ),
`Okada.tex`'s S,I1,I2,P1,P2 network (ERC/organization angle added),
coinfection synergy, sequential-infection complementarity, vector-host
cycle (with/without demographic inflow), and the dengue model reusing
`TwoStrainDengueBuffering.wl`'s exact 38-reaction network (finds ERC
synergy between the two strains witnessing the co-infected vector class
`V12`). All numbers transcribed into `COT.tex` come from this file's
output, not hand-derivation. `SubFF.tex` compiles cleanly with the new
`COT.tex` (verified via `pdflatex`, two passes, zero new undefined
references).

`EpidCRNmodels/OrganizationsGrid.wl` extends this to every method
(`allErc`,`allOrg`,`minSiph`,`synQ`,`cmpERC`,`allBuff`/`buffS`) on every
example uniformly, backing the summary table in the new preceding section
`CocktailArXiv-updated/Summary.tex` (`\input` before `COT` in `SubFF.tex`).
Notable verified findings from the full grid, not previously in `COT.tex`:
`cmpERC` returns a trivial witness (`{S}`, the shared demographic
susceptible) between almost any two ERCs in an open network — the
interesting cases are the ones with a *non-trivial* witness (Example D's
`{S1}`, the closed vector-host model's `{IH}`/`{IV}`); the closed
(no-inflow) vector-host network has **zero** chemical organizations
(demographic recruitment is structurally necessary for persistence, not
just typical); and `Method->"ERCGenerated"` on the dengue network misses
the known `{R}` buffering structure that `Method->"Exhaustive"` finds
(concrete evidence the heuristic mode is not complete).

Added two more rows to the same grid/table: the OSN rumor-propagation
model of Fakih, Halanay and Avram (`CholeraShuai13/mo3.tex`), both
`omega=0` and `omega>0` variants (built from its reaction table directly —
the methods don't depend on its non-mass-action fractional rates).
`minSiph` independently reproduces that source's own claimed siphons
exactly (`{U},{S1,B1},{S2,B2}`, plus `{W}` only when `omega=0`). Notable:
`cmpERC` returns `False` (no complementarity at all, not even the
trivial-`S` kind) for both variants, because this model's `U` class is
never freely supplied by inflow (unlike the epidemic examples' shared
`S`); and `omega=0` roughly doubles organizations (6→12) and more than
triples nontrivial buffering structures (7→24), a quantitative
confirmation of that source's own informal claim that `R` "does not feed
back" when `omega=0`.

---

## k-Contraction Package — Plan (2026-08-19)

`GLVbook/Rap/compu.tex` proposes a `KContraction`` package implementing the
chain `f -> J -> J^[k] -> transformed compound metric -> matrix measure`,
plus SOS/LMI metric search and compound-free criteria. Rather than a
separate package, its purely symbolic layer belongs in `Utils.wl`, next to
`ACM` and `dulP`, which already implement two pieces of it.

**What already existed, unrecognized by compu.tex:** `ACM[A,k]` (additive
compound `J^[k]`, via the standard `D[Minors[I+t*A,k],t]/.t->0` identity) is
exactly compu.tex's proposed `AdditiveCompound[J,k]`. `dulP[RHS,var,{a,b}]`
is exactly the planar (`n=2,k=2`, monomial `P`) special case of the
"Transformed compound metrics" section.

**Phase 1 (done, this entry) — the symbolic layer, no external solver:**
`Pdot[P,f,x]` (`P_f`, elementwise `D` so `P` need not be diagonal), and
`TACM[J,k,P,f,x]` (`B_k=P_f.P^{-1}+P.J^{[k]}.P^{-1}`, generalizing `dulP`'s
`TJ` to any order `k` and any invertible state-dependent `P`), and
`matMu[B,p]` (matrix measure/logarithmic norm, `p` in `{1,2,Infinity}`).
The `k`-contraction certificate itself is the two-line composition
`matMu[TACM[J,k,P,f,x],p]<0`, matching the division of labor `dulP` already
uses (it returns `TJ`; the caller applies `Reduce[ForAll[...]]`). Reuses
`ACM` unchanged; does not touch `posM`, `NGM`, or `minSiph`.

**Verified:** `TACM` at `k=2` with the 1x1 metric `P={{dulx^-1}}` on `dulP`'s
own Bogdanov–Takens test system reproduces `dulP`'s `TJ` term-for-term
(checked independently in Python/SymPy before writing the Wolfram code, then
as the uncommented in-file test). `matMu` reduces to the max-column-sum
formula on a Metzler matrix, matching the note already in `Utils.wl`'s
`Hur4M`-adjacent documentation.

**Phase 2a (done, 2026-08-19) — the constant/polytope quadratic-LMI layer,**
using Mathematica's native `SemidefiniteOptimization` directly, no YALMIP:
`qLMI[Jlist,eta]` (fixed-rate feasibility test on a matrix polytope, symmetric
`P` built from fresh `Unique[]` variables, `P>=I` normalized, `-I` margin;
fixing `eta` avoids the bilinear `eta*P` term, exactly as compu.tex's own
"Quadratic contraction metrics as LMIs" subsection recommends), `qRate`
(bisection on `eta` using `qLMI` as the oracle, returns `{eta,P}`), and
`jPoly[A0,As,hb]` (matrix polytope vertex construction for
`J(x)=A0+Sum hj(x) Aj` given known bounds `hb` on each `hj` — only the
mechanical `2^m`-vertex enumeration is automated, matching how Phase 1 also
deferred symbolic auto-derivation; feed `jPoly`'s output straight into
`qRate` for the polytope-robust certificate compu.tex describes).

**Verified:** the `SemidefiniteOptimization[0,cons,vars]` feasibility-SDP
pattern (`VectorGreaterEqual[{M,0},"SemidefiniteCone"]` constraints) was
tried interactively against a known-feasible and a known-infeasible toy LMI
before being written into the package, so the infeasibility signal
(`Indeterminate` values + `SemidefiniteOptimization::nsolc`, not an
exception) is handled deliberately, not guessed. `qRate` on `A=diag(-1,-2)`
recovers the exact known best rate `eta=1` (the slower eigenvalue) to
bisection tolerance, both alone and paired with a second polytope vertex
that only flips the off-diagonal sign. `jPoly` reproduces hand-computed
vertices exactly at `m=1` and `m=2`. All four new functions' in-file tests
were run through a real kernel (`wolframscript`, loading the actual
`EpidCRN\`` package the normal way) before and after being embedded, with
no regression in the pre-existing `dulP`/`ACM`/`TACM`/`matMu` tests.

**Phase 2b (done, 2026-08-19) — the SOS polynomial-metric layer,** following
Parrilo's Gram-matrix relaxation ([AylwardParriloSlotine2008], cited by
compu.tex): `sosCon[M,x,d]` (builds, without solving, the linear-equality-
plus-PSD-cone system certifying that `v.M.v` is SOS in `(x,v)` for a
symmetric polynomial matrix `M(x)`, kept separate from solving so several
such conditions can share one SDP), `sosMat[M,x,d]` (standalone solver
wrapping `sosCon`), `polyP[x,n,d]` (symmetric polynomial matrix ansatz with
fresh unknown coefficients), and `polyMet[J,f,x,d,eta]` — compu.tex's
`PolynomialContractionMetric`: joins the SOS conditions for `P(x)>0` and the
differential contraction inequality (reusing Phase 1's `Pdot` for `P_f`) into
one SDP over the shared metric coefficients, at a fixed `eta` (avoiding the
same bilinear `eta*P` term `qLMI` avoids).

**Verified:** the equality-constraint-plus-cone `SemidefiniteOptimization`
pattern (needed for Gram coefficient-matching, new API surface beyond Phase
2a) was tried on `p(x)=x^4+2x^2+1=(x^2+1)^2` before being written into the
package. `sosMat` correctly certifies a hand-verified SOS matrix
(`[[1,x],[x,x^2+1]]`, whose `v.M.v=(v1+xv2)^2+v2^2`) and correctly rejects a
matrix that provably cannot be one (`[[x,0],[0,x]]`, odd in `x`). `polyMet`
at `d=0` (constant `P`, so `P_f` vanishes identically) was checked to agree
exactly with `qLMI[{J},eta]` on `A=diag(-1,-2)` at both a feasible
(`eta=0.9`) and an infeasible (`eta=1.5`) rate — ties Phase 2b back to the
already-validated Phase 2a rather than testing it in isolation. On a
genuinely nonlinear system (`f=-x-x^3`, so `Pdot` does real work) `polyMet`
found `P(x)~35.4+20.1x^2` (an even polynomial — matching the problem's own
even/odd symmetry, a sign the solver found something structurally sensible,
not a numerical accident), independently re-verified by direct numeric
substitution of both `P(x)-I` and the contraction-inequality matrix at seven
sample points across `[-2,2]`, not by trusting the SDP's own internal Gram
bookkeeping. Full package reload after embedding: all nine self-tests pass
(the six from Phases 1-2a plus these three), load time unchanged
(~2.5-2.8s), no regressions.

One environment note for whoever runs this next: `wolframscript -file`
sometimes silently fails to `Get` the package on a reused kernel session
(`Get::noopen`) even when the exact same script worked moments before: prefer
an explicit absolute path over `<<EpidCRN\`` + `$Path`, and prefer forward
slashes over `\\`-escaped Windows paths in any shell-heredoc-embedded `.wl`
script — backslash characters have repeatedly been silently eaten when
written this way.

**`optMon` (done, 2026-08-19) — compu.tex's `OptimizeKMetric`, but not the
feared bi-level search.** For a *monomial* ansatz `p=x^a y^b`, `dulP`'s `TJ`
turns out to be **affine in `(a,b)` at any fixed `(x,y)`** (no `a*b` cross
term — visible directly in `LVA-monomial-formula`, `Rap/LVA3.tex`), so
"search over exponents" is a semi-infinite *linear* program, not a nonconvex
one: `optMon[RHS,var,box,pts0]` runs a cutting-plane loop
(`LinearOptimization` maximizing the margin over a growing sample set),
verifies the winning `(a,b)` *exactly* via `Reduce[ForAll[...]]`, and on
failure extracts a genuine counterexample via `FindInstance` to add to the
sample set for the next round. Reuses `dulP` unchanged.

**Tested against the real system named in this task** — `Rap/LVA3.tex`'s own
`LVA-counterexample`, at `kappa^*=(1,1,4,1,1/4000,1/500)`, where the theorem
*proves* the `b=-1` family has no certifying `a in [-3,0]` (witness
`TJ(1/4,1/20)>=9/400` for every such `a`) but `p=x^-2 y^-2` (outside that
family) works. `optMon` reproduces both facts: restricted to `b=-1` it
correctly reports infeasible, seeded at the paper's own witness point — and
at `(a,b)=(0,-1)`, dulP's `TJ(1/4,1/20)` came out to exactly `9/400`,
matching the paper's hand computation bit for bit, an independent check that
`dulP`+`optMon` are reading the same system the paper is. The unrestricted
2D search correctly finds a certificate — `(a,b)=(0,-5)`, not the paper's own
`(-2,-2)` but a genuine different one, confirmed independently by hand
(`sup_{x,y>0} TJ` is attained only in the limit `x->0+,y=1/10`, value
`-0.199<0` there, via elementary calculus, not by trusting the code).

The refinement loop specifically (not just the answer) was also exercised
and checked: every attempt on the real system landed on a box corner in one
shot (`iters=0`), so a hand-built adversarial single-exponent case was added
(`RHS1/x=-x^2+6x-8`, non-monotone: `x=3` forces `a<-1`, `x=1/2` forces
`a>-0.524`, contradictory), seeded with only the misleading point `x=3` —
`optMon` found the locally-good-but-globally-wrong `a` first, `Reduce`
caught it, `FindInstance` extracted the real counterexample, and the
resolved LP correctly went infeasible (`iters=1`), matching the
hand-derived contradiction. All three checks embedded as the in-file test;
full package reload after embedding: ten self-tests, all `True`, no
regressions, load time ~3.3s (up from ~2.8s, the cost of two `Reduce` calls).

## Phase 3 — compound-free criteria and orchestration (started 2026-08-19)

`kTau` (done) — compu.tex's `CompoundFreeKContraction`, `Method->"DalinEtAl"`
only, faithfully implementing Theorem 17 of Dalin, Ofir, Bar-Shalom,
Ovseevich, Bullo & Margaliot, *Verifying k-Contraction without Computing
k-Compounds*, IEEE TAC 69(3):1492-1506, 2024 — full text pulled from
arXiv:2209.01046 (`wolframscript`/PDF rendering can't display pages on this
machine, so the PDF was fetched and its text extracted with Python `pypdf`
instead, then read directly) specifically **because** compu.tex insists on
"the hypotheses of each theorem literally, rather than replacing them by one
generic 'LMI test'" — a generic stand-in was never on the table for this
piece. `kTau[A,k,p,T]` computes the k-shifted log norm
`tau_{p,k}(A)=tr(A)+(n-k)*mu_q(T.(-A).T^-1)` (their Definition 4) directly
from `A`, reusing `matMu` unchanged; their Theorem 17 gives
`mu_{p,T^(k)}(A^[k])<=tau_{p,k}(A)`, so `tau_{p,k}(A)<=-eta<0` certifies
k-contraction **without ever forming the k-additive compound** `A^[k]`.
`tau_{p,k}` is convex in `A` (their Prop. 19(7)), so exactly as `qLMI`/
`qRate`/`jPoly` already exploit convexity via vertex enumeration, a Jacobian
ranging over a polytope or box is handled by `Max[kTau[#,k,p,T]&/@vertices]`.

No `Method->` dispatcher was built: with only two modes implemented (see
`cecKC` below), and the second covering only the constant-matrix case, a
dispatcher would still be premature.

**Verified against the paper's own stated exact results, not paraphrased
ones:** `A=I_n` gives `tau_{p,k}(I_n)=k` exactly for any `n,k,p` (their
Example 10, stated as the case where the bound is *tight*, i.e. equality,
not merely an inequality); the general diagonal-LTI closed form (their
Eq. 41) matched by direct substitution; and their own **published worked
numerical example** (Example 12): an `n=3` Hopfield network with all weights
`w_ij=1`, `phi=tanh`, `r_i=r` for all `i` — provably *not* contractive in any
norm (three equilibria: `0` and `+-1.2447*(1,1,1)`, stated exactly in the
paper) — yet satisfying their sufficient condition for 2-contraction exactly
on their own stated threshold `r<1/2`: `kTau` at `r=0.49` gives the exact
rational `-4/49<0` (certified) and at `r=0.51` gives `+4/51>0` (not
certified), the sign flipping precisely at their threshold, both values
matching hand substitution into their closed form. `J(x)` is affine in the
activation derivatives `phi_i'(x_i) in [0,1]`, so by convexity the worst
case is at a corner of `[0,1]^3`, checked by brute force (`n=3` makes this
trivial — `Tuples[{0,1},3]`, 8 corners). Full reload: 11 self-tests, all
`True`, no regressions, load time unchanged (`kTau` is pure algebra, no
solver call).

**`exKC`/`posDefQ`/`verQ`/`exVerQ`/`verPoly`/`exVerPoly` (done, 2026-08-19) —
compu.tex's `ExactifyKCertificate`+`VerifyKCertificate`,** in two parallel
tracks matching the two shapes of certificate Phase 2 actually produces:

- *Constant certificates* (`qLMI`/`qRate`): `exKC[val,tol]` rationalizes a
  numeric candidate (elementwise, so it works equally on a pure numeric
  matrix or one with symbolic entries and numeric coefficients — only the
  numbers move); `posDefQ[M]` is exact positive-definiteness via Sylvester's
  criterion (leading principal minors, no `Reduce`, no floating point once
  `M` is exact); `verQ` applies `posDefQ` to both `Pex>0` and the LMI;
  `exVerQ` tries a sequence of tolerances (coarsest first), keeps the ones
  that verify exactly, returns the simplest.
- *Polynomial certificates* (`polyMet`): `verPoly[M,x,cond]` is the same
  Sylvester idea but quantified over `x` via `Reduce[ForAll[...]]` — the
  same exact-verification idiom `optMon` already uses, generalized from a
  scalar to a matrix; `exVerPoly` is `exVerQ`'s counterpart, rebuilding both
  `P-I` and the differential contraction-inequality matrix (reusing `Pdot`)
  from each exactified candidate.

Both return `{"Exact",Pex}` or `{"Sampled",P}` (the unverified numeric
original) — **`"RigorousNumerical"` (interval-arithmetic fallback) is
explicitly not implemented**, named as such rather than silently absorbed
into `"Sampled"`, so nothing claims a guarantee it doesn't have.

**A real syntax pitfall hit and fixed while building `verPoly`:**
`ForAll[{x},cond,expr]` (a *list* containing one variable) is not valid —
`ForAll` needs a bare symbol for one variable, or an explicit list of two-or-
more; passing `{x}` produces a `Reduce::naqs` failure that silently looks
like a legitimate "not verified" result if you're not watching for it (it
isn't — it's a syntax error dressed as one). A second, distinct pitfall:
`Reals` is not a `ForAll` *condition* (its 2nd argument is a boolean, and
`Reals` isn't one) — the domain restriction belongs on the *outer* `Reduce`
call, not stuffed into the quantifier itself. Both are noted directly in
`verPoly`'s comment so they aren't rediscovered.

**Verified two ways, tying this new layer back to already-validated
certificates rather than testing it in isolation:** on `qLMI`'s own
`A=diag(-1,-2)`, `eta=0.9` result, `exVerQ` returns the clean
`Pex={{339/8,0},{0,229/12}}`, `label="Exact"` — confirmed a second time by
direct `Eigenvalues` of the exact `Pex` substituted into the LMI, bypassing
`verQ`'s own code entirely. On `polyMet`'s own nonlinear test system
(`f=-x-x^3`, degree 2, `eta=0.5`), `exVerPoly` returns
`Pex=248/7+(181/9)x^2`, `label="Exact"` — confirmed a second time by direct
numeric evaluation of both matrix inequalities at seven sample points across
`[-5,5]`, bypassing `verPoly`'s own `Reduce`-based code entirely. Full
reload: 13 self-tests, all `True`, no regressions, load time ~3.7s (up from
~3.25s, entirely accounted for by the added `Reduce` calls).

**`cecKC`/`inCnt`/`exCecKC` (done, 2026-08-19) — `Method->"GeneralizedLyapunov"`,
constant-matrix case, faithfully implementing Andreu/Cecilia, Zoboli, Astolfi,
Serres & Andrieu, *Generalized Lyapunov Conditions for k-Contraction: Analysis
and Feedback Design*, IEEE TAC 71(2):1022-1037, 2026 — full source (`2ndVersion.tex`,
supplied directly by the user, not fetched) read specifically for the exact
theorem statement, not a paraphrase.** For `xdot=A.x`, their
Theorem `k_contraction_LMI_cLTI` (linear case, specialized to `l=2`) asks for
constant symmetric `P0` (`In=(0,0,n)`, positive definite), `P1`
(`In=(k-1,0,n-k+1)`, exactly `k-1` negative eigenvalues), and scalars
`mu0,mu1` with `A^T.Pi+Pi.A<2*mui*Pi` (`i=0,1`) and `mu1+(k-1)*mu0<0` — and,
unlike `kTau`, this is stated as **necessary and sufficient** for
k-contraction of a constant `A`, not merely sufficient.

The construction is direct, not an SDP search. Their own Remark `algorithm`
says the prescribed-inertia requirement "cannot be represented as
semidefinite constraints... [but] can be dropped without significant impact,"
because Lemma `inertia_shifting` item (4) is exactly the classical
Ostrowski-Schneider Lyapunov inertia theorem: for `Q=I` and any `mu` not
equal to any `Re(eigenvalue of A)`, the Sylvester equation
`(A-mu I)^T.P+P.(A-mu I)=-I` has a *unique* symmetric solution `P`
(`LyapunovSolve` in Wolfram Language), and expanding the equation gives
`A^T.P+P.A=2*mu*P-I` — exactly the required strict LMI with margin `I` — while
`P`'s inertia is *forced* to `(p,0,n-p)` where `p=#{eigenvalues of A with
Re>mu}`. So choosing `mu0` above the largest `Re(eigenvalue)` (`p=0`) and
`mu1` strictly between the `(k-1)`th and `k`th largest `Re(eigenvalue)`
(`p=k-1`) automatically yields `P0,P1` of the *required* inertia, with no
inertia-constrained (nonconvex) search at all — the crux insight that makes
this tractable. `cecKC[A,k]` tries several `mu0`/`mu1` candidates spanning the
valid range (`lyaTry`, internal helper) rather than the exact midpoint, after
hitting a genuine degeneracy while building this: the arithmetic midpoint of
two real eigenvalues makes `(A-mu I)`'s eigenvalues sum to exactly zero,
which is precisely `LyapunovSolve`'s non-uniqueness condition
(`LyapunovSolve::nuniq`) — not a bug, a real boundary case now documented in
`cecKC`'s own comment. `exCecKC` mirrors `exVerQ`'s exactify-then-verify
pattern, independently confirming all four of the theorem's own conditions
(`inCnt` — new, exact inertia triple via eigenvalue sign counts, matching the
paper's own `In(P)=(pi_-,pi_0,pi_+)` convention, negative count first,
confirmed directly from the paper's own `def.tex` source — plus `posDefQ` on
both LMI residuals, plus the scalar sum condition) rather than trusting
`cecKC`'s own reasoning.

**Verified three ways.** (1) `A=diag(-1,-2)`, `k=2`, `n=2`: since `k=n`,
`In(P1)=(1,0,1)` is honestly indefinite, not a degenerate corner case;
`exCecKC` returns exact `P0` with eigenvalues `{50,1/2}` (both positive) and
`P1` with eigenvalues `{5/3,-5/7}` (one of each sign), confirmed independently
via `inCnt`. (2) A genuinely discriminating example demonstrating the
theorem's own advertised advantage over `kTau`: a triangular (hence exactly
known eigenvalues `-3/10,-9/10,-3/2`, top-2 sum `-6/5<0`, genuinely
2-contractive) but strongly non-normal matrix (off-diagonal shear `200`) —
`kTau` at its default identity metric gives `tau=12.35>0`, **failing to
certify** a fact that is true, exactly the known weakness of a fixed-metric
sufficient-only test on a non-normal matrix — while `cecKC` still certifies
it exactly, `exCecKC` reconstructing exact rational `mu0=-29/100`,
`mu1=-18/25` and algebraic (`Root`-object) `P0,P1` eigenvalues with exactly
the required signs. (3) `diag(1,-1/10,-1/10)`, top-2 eigenvalue sum `9/10>0`,
genuinely **not** 2-contractive — `cecKC` correctly and cleanly returns
`False` rather than a spurious certificate. Full reload: 15 self-tests, all
`True`, no regressions.

**`jntCec` (done, 2026-08-19) — the state-dependent/polytope case of
Theorem `k_contraction_LMI` proper (`J(x)` varying over a forward-invariant
set `S`, "for all `x` in `S`"), previously deferred.** The obvious approach
— a fresh feasibility-only SDP search for a genuinely indefinite `P1`
shared across several vertices (dropping the inertia constraint per
Remark `algorithm`, exactly as `cecKC`'s own single-matrix construction
does) — was tried at length and abandoned: the solver, minimizing nothing
in particular, consistently returns a **positive definite** point whenever
anything is feasible at all, never an indefinite one, even with an added
linear constraint biasing toward a specific vertex's own unstable
eigenvector (tried on the SEIR2 example from the six-example sweep below;
see `SEIR2cert.tex`'s own open remark for that history). A cleaner
diagnosis followed: for the same SEIR2 example, one vertex has a tied
complex-conjugate top-2 eigenvalue pair, which can *only* ever contribute
an even eigenvalue-count crossing to any real threshold — a hard
structural reason no vertex-by-vertex eigenvalue-matching recipe could
ever produce a common `P1` there, independent of search method.

What worked instead needed no new search at all: `jntCec[Jlist,k]` tries
`cecKC`'s own *exact* single-matrix certificate at each vertex in turn
against every *other* vertex's LMI residual, returning the first that
survives everywhere. Validated two ways. (1) A small perturbation of
`cecKC`'s own `A=diag(-1,-2)` test system — the exact certificate found at
the unperturbed matrix survives at the perturbed one by continuity,
confirmed independently by direct eigenvalue evaluation of the same
residuals, bypassing `jntCec`'s own code. (2) The four-dimensional network
of `4Pcert.tex`, restricted to the box `x1,x2,x3∈[2/5,8/5]` (eight
vertices) — `jntCec` finds a certificate at one vertex that, exactified via
`exCecKC`, is confirmed **exact and independently verified via `posDefQ`
at all eight vertices simultaneously**: a genuine forward-invariant-region
certificate, not a single point, closing exactly the gap `4Pcert.tex`
itself left open (a constant certificate cannot extend to the *whole*
unbounded orthant, but a substantial *bounded* region was open, and now
isn't). The same network's full box `[0,2]^3` correctly returns `False`:
two of its own eight vertices individually fail `cecKC` (exact eigenvalue
ties from the network's permutation symmetry, already documented in
`4Pcert.tex`), so no single vertex-certificate could possibly extend to
the whole box, and none does. Full reload: 16 self-tests, all `True`, no
regressions. No claim of completeness: a polytope where no single vertex's
certificate extends everywhere may still admit some other joint
certificate this strategy does not find — the harder, genuinely
open-ended construction (search jointly rather than reuse a single
vertex) remains unexplored.

**`jntCec` sweep across all six Phase-3 examples (done, 2026-08-19+).**
Once `jntCec` existed, every one of the six examples in `Compound.tex`
(the four-dimensional network, split-exposed SEIR2, reversible LVA,
Bogdanov–Takens, classical mass-action SEIR, and the structured two-strain
quotient) was retested on a shrunk box around a real-eigenvalue point,
each written up as its own `*cert.tex`/`*cert.wl` pair with an exact,
independently `posDefQ`-verified certificate. Two of the six (LVA, BT)
have a genuinely quadratic Jacobian entry, handled by treating the square
as a third independent box parameter — a valid but conservative
relaxation (some resulting vertices mismatch the linear and quadratic
parameters and correspond to no real point on the curve); the other four
are affine in the relevant state variables directly, no relaxation
needed. Beyond the spiral obstruction already known from the pointwise
sweep, this pass surfaced two further, genuinely distinct reasons a joint
box can be much smaller than expected even with no complex eigenvalues
anywhere: (i) `cecKC`'s own default `mu0` sits barely one percent above
its vertex's own top eigenvalue, essentially no margin, so whenever the
top eigenvalue is sensitive to the box parameters (BT, classical SEIR,
two-strain all hit this to varying degrees) the surviving box is far
smaller than the box `cecKC` succeeds on pointwise at every corner; (ii)
in the two-strain example, a box straddling the invariant face `I2=0`
failed badly even at tiny width, because the block-triangular structure
that makes the *single-point* certificate exact and clean there
generalizes poorly to nearby points — centering away from the face
instead fixed it immediately. Every one of the six exact certificates
found is recorded in its own file: `4Pcert.tex` (already covered above),
`SEIR2cert.tex`, `LVA3cert.tex`, `BTcert.tex`, `SEIRcert.tex`,
`exaI2cert.tex`, all in the `GLVbook/Rap` project folder, each paired
with the `.wl` script that reproduces its own claimed numbers exactly.

**Deferred, not yet implemented:**
- *`KContractionSpectrum`* — should be nearly free now: for k=1..n, compare
  `qRate`/`polyMet` (needs compound `ACM[J,k]`, already in Phase 1) against
  `kTau` (compound-free) at each order, tabulating margins. Straightforward
  orchestration of already-validated pieces; not yet written only for lack
  of a real multi-order example lined up to validate it against.
- *`KContractionRegion`* — boundary tracing over a parameter, combining
  symbolic elimination with numerical continuation; the most open-ended
  remaining piece, deserves its own session.
- *`KContractionAnalyze`* — final orchestration layer; depends on everything
  above existing first.

`compu.tex` itself is untouched; this section is the authoritative plan.

---
